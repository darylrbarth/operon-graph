import pandas as pd
import networkx as nx
import re
from itertools import combinations
from typing import List, Dict, Set, Tuple

class GraphBuilderProst:
    """
    Builds a graph from PROST cluster results and operons.
    Cleans raw PROST DataFrame on initialization.
    """

    def __init__(self, prost_tsv, header_file, operonic_distance=40):
        """
        prost_tsv: Path to PROST TSV file
        header_file: Path to gene headers FASTA file
        operonic_distance: Max distance to group genes into an operon
        """
        self.prost_df = self._load_and_clean_prost(prost_tsv)
        self.header_file = header_file
        self.operonic_distance = operonic_distance

        self.genes = []
        self.operons = []

        self.member_to_best_head = {}
        self.head_to_members = {}
        self.member_to_cluster = {}

        self.nodes = set()
        self.gene_pairs = []
        self.edges = []
        self.edge_properties = {}
        self.weighted_edges = []

        self.G = nx.Graph()

    # ---------------- Load & clean PROST results ----------------

    def _load_and_clean_prost(self, prost_tsv):
        """
        Cleans the raw PROST DataFrame:
        - Strips whitespace
        - Cleans IDs
        - Converts score and e_value to numeric
        - Keeps columns: ['head','member','score','e_value']
        """
        df = pd.read_csv(prost_tsv, sep="\t", header=None)
        
        df = df.dropna(axis=1, how='all')  # drop empty columns

        # Assign default column names if none exist
        if df.shape[1] >= 4:
            df = df.iloc[:, :4]
            df.columns = ['head', 'member', 'score', 'e_value']
        else:
            raise ValueError("PROST DataFrame must have at least 4 columns: head, member, score, e_value")

        # Strip whitespace
        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

        # Convert score/e_value to numeric
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df['e_value'] = pd.to_numeric(df['e_value'], errors='coerce')

        # Clean IDs (remove spaces and # annotations)
        df['head'] = df['head'].apply(lambda x: str(x).split()[0].split('#')[0])
        df['member'] = df['member'].apply(lambda x: str(x).split()[0].split('#')[0])

        return df
    
     # ---------------- Gene header parsing ----------------

    def parse_headers(self):
        """
        Parse headers file to extract ORF information: orf_id, contig, start, end, strand
        """
        genes = []
        with open(self.header_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    continue
                parts = [p.strip() for p in line.split('#')]
                orf_id = parts[0].split()[0][1:]
                contig = re.sub(r'_\d+$', '', orf_id)
                start = int(parts[1])
                end = int(parts[2])
                strand = '+' if parts[3].strip() == '1' else '-'

                genes.append({
                    "orf_id": orf_id,
                    "contig": contig,
                    "start": min(start, end),
                    "end": max(start, end),
                    "strand": strand
                })
        self.genes = genes
        return genes

    # ---------------- Operon grouping ----------------

    @staticmethod
    def gene_distance(g1: Dict, g2: Dict) -> int:
        if g1['contig'] != g2['contig'] or g1['strand'] != g2['strand']:
            return None
        g1, g2 = sorted([g1, g2], key=lambda x: x['start'])
        return g2['start'] - g1['end']

    def group_operons(self):
        """
        Group genes into operons based on operonic_distance
        """
        genes_sorted = sorted(self.genes, key=lambda g: (g["contig"], g["strand"], g["start"]))
        current_operon = [genes_sorted[0]]

        for g in genes_sorted[1:]:
            last_gene = current_operon[-1]
            dist = self.gene_distance(last_gene, g)
            if dist is not None and dist <= self.operonic_distance:
                current_operon.append(g)
            else:
                self.operons.append(current_operon)
                current_operon = [g]

        self.operons.append(current_operon)
        return

    # ---------------- Load operons from headers ----------------
    def load_operons(self):
        """
        Parse headers file and automatically group genes into operons.
        Stores genes in self.genes and operons in self.operons
        """
        self.parse_headers()
        self.group_operons()
        return
    
    # ---------------- PROST cluster processing ----------------

    def resolve_clusters(self):
        """
        Assign each member to the head with the highest score to remove overlaps.
        """
        for idx, row in self.prost_df.iterrows():
            member = row['member']
            head = row['head']
            score = row['score']

            if member not in self.member_to_best_head or score > self.member_to_best_head[member][1]:
                self.member_to_best_head[member] = (head, score)

        # Build head_to_members
        for member, (head, _) in self.member_to_best_head.items():
            if head not in self.head_to_members:
                self.head_to_members[head] = set()
            self.head_to_members[head].add(member)

        # Build member -> head map (for fast lookup)
        for head, members in self.head_to_members.items():
            for m in members:
                self.member_to_cluster[m] = head

    # ---------------- Nodes ----------------
    
    def create_nodes(self):
        """Create nodes from heads of clusters and prepare basic node attributes."""
        self.nodes = set(self.head_to_members.keys())
        # Basic attributes you can use in Cytoscape:
        # - cluster_size: number of member proteins in this PROST cluster
        self.node_attributes = {h: {'cluster_size': len(members)} for h, members in self.head_to_members.items()}
        return
    
    # ---------------- Operon-based edges ----------------

    def create_gene_pairs(self):
        """Create gene pairs only between adjacent genes within each operon (same contig + strand by construction)."""
        self.gene_pairs = []
        for operon in self.operons:
            if len(operon) < 2:
                continue
            # operon list is already sorted by genomic coordinate (see group_operons)
            for i in range(len(operon) - 1):
                g1 = operon[i]
                g2 = operon[i + 1]
                self.gene_pairs.append((g1['orf_id'], g2['orf_id']))
        return
        
    def create_edges(self):
        """Map gene pairs to cluster head edges, skipping genes not present in PROST clustering."""
        self.edges = []
        missing = 0
        for g1, g2 in self.gene_pairs:
            c1 = self.member_to_cluster.get(g1)
            c2 = self.member_to_cluster.get(g2)
            if c1 is None or c2 is None:
                missing += 1
                continue
            if c1 == c2:
                continue
            self.edges.append((c1, c2))
        if missing:
            print(f"[create_edges] Skipped {missing} gene-pairs due to missing cluster assignments.")
        return

    # ---------------- Edge weighting ----------------

    def calculate_edge_properties(self):
        """
        Add absolute and normalized weights for edges.
        """
        for edge in self.edges:
            n1 = len(self.head_to_members[edge[0]])
            n2 = len(self.head_to_members[edge[1]])
            possible_connections = n1 * n2

            if edge in self.edge_properties:
                props = self.edge_properties[edge]
                props['abs_weight'] += 1
                props['norm_weight'] += 1 / possible_connections
            else:
                self.edge_properties[edge] = {
                    'abs_weight': 1,
                    'norm_weight': 1 / possible_connections
                }

    def make_weighted_edges(self):
        """
        Convert edge properties into weighted edge tuples for networkx.
        """
        self.weighted_edges = [
            (e[0], e[1], self.edge_properties[e])
            for e in self.edge_properties
        ]

    # ---------------- Build Graph ----------------

    def run(self):
        """
        Runs the full graph construction pipeline and returns nx.Graph
        """
        self.load_operons()
        
        self.resolve_clusters()
        # print(f'Number of cluster heads: {len(self.head_to_members)}, example: {list(self.head_to_members.keys())[:5]}')
        
        self.create_nodes()
        # print(f'Number of nodes: {len(self.nodes)}, example: {list(self.nodes)[:5]}')
        
        self.create_gene_pairs()
        self.create_edges()
        self.calculate_edge_properties()
        self.make_weighted_edges()

        self.G.add_nodes_from(self.nodes)
        # attach basic node attributes
        if hasattr(self, 'node_attributes'):
            nx.set_node_attributes(self.G, self.node_attributes)

        self.G.add_edges_from(self.weighted_edges)
        return self.G


    # ---------------- Export helpers ----------------

    def export_gml(self, out_gml: str):
        """Write the built graph to a .gml file for Cytoscape Desktop."""
        if self.G is None or self.G.number_of_nodes() == 0:
            raise ValueError("Graph is empty. Run builder.run() before exporting.")
        nx.write_gml(self.G, out_gml)
        return out_gml

    def export_tsv(self, nodes_tsv: str, edges_tsv: str):
        """Export node and edge tables as TSVs (good for Cytoscape Web or import as tables).

        Produces:
          - nodes_tsv: id, cluster_size
          - edges_tsv: source, target, abs_weight, norm_weight
        """
        if self.G is None or self.G.number_of_nodes() == 0:
            raise ValueError("Graph is empty. Run builder.run() before exporting.")

        # Nodes
        node_rows = []
        for n in self.G.nodes():
            attrs = self.G.nodes[n]
            node_rows.append({
                'id': n,
                'cluster_size': attrs.get('cluster_size', None),
            })
        pd.DataFrame(node_rows).to_csv(nodes_tsv, sep='\t', index=False)

        # Edges
        edge_rows = []
        for u, v, attrs in self.G.edges(data=True):
            edge_rows.append({
                'source': u,
                'target': v,
                'abs_weight': attrs.get('abs_weight', None),
                'norm_weight': attrs.get('norm_weight', None),
            })
        pd.DataFrame(edge_rows).to_csv(edges_tsv, sep='\t', index=False)
        return nodes_tsv, edges_tsv



 # ---------------- USAGE EXAMPLE ----------------
if __name__ == "__main__":
    builder = GraphBuilderProst(
        prost_tsv="./data/prost_results/ecoli.tsv",
        header_file="./data/genomes/ecoli_genome.fna.faa",
        operonic_distance=40
    )

    G = builder.run()
    print(f"Graph built with {len(G.nodes())} nodes and {len(G.edges())} edges")
    print(f"Node examples:\n {list(G.nodes)[:5]}")
    print(f"Edge examples:\n {list(G.edges(data=True))[:5]}")

    # Exports for Cytoscape (Desktop + Web)
    builder.export_gml("operon_network.gml")
    builder.export_tsv("operon_nodes.tsv", "operon_edges.tsv")
    print("Wrote: operon_network.gml, operon_nodes.tsv, operon_edges.tsv")
