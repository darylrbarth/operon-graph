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

    def __init__(self, prost_tsv, header_file, operonic_distance=200):
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
        """Create nodes from heads of clusters."""
        self.nodes = set(self.head_to_members.keys())
        return
    
    # ---------------- Operon-based edges ----------------

    def create_gene_pairs(self):
        """Map gene pairs to cluster head edges."""
        for operon in self.operons:
            # connect all pairs of genes within an operon
            for gene1, gene2 in combinations(operon, 2):
                self.gene_pairs.append((gene1['orf_id'], gene2['orf_id']))
        return
        
    def create_edges(self):
        """Map gene pairs to cluster head edges."""
        for pair in self.gene_pairs:
            # match the gene pair to their cluster heads
            self.edges.append((self.member_to_cluster[pair[0]], self.member_to_cluster[pair[1]]))
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
        self.G.add_edges_from(self.weighted_edges)
        return self.G



 # ---------------- USAGE EXAMPLE ----------------
if __name__ == "__main__":
    builder = GraphBuilderProst(
        prost_tsv="./data/prost_results/ecoli.tsv",
        header_file="./data/genomes/ecoli_genome.fna.faa",
        operonic_distance=200
    )

    G = builder.run()
    print(f"Graph built with {len(G.nodes())} nodes and {len(G.edges())} edges")
    print(f"Node examples:\n {list(G.nodes)[:5]}")
    print(f"Edge examples:\n {list(G.edges(data=True))[:5]}")