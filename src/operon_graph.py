
import networkx as nx 

from src.utils.prot_clust_utils import mmseq_cluster
from src.graph_builder import GraphBuilder
from src.config import OperonGraphConfig


class OperonGraph():

    """
    Operon Graph. Top level class that owns the graph object (nx Graph), uses the clustering tools and GraphBuilder to help.
    """

    def __init__(self, config: OperonGraphConfig):
        self.threshold = config.threshold
        self.evalue = config.evalue
        self.input_format = config.input_format
        self.operonic_distance = config.operonic_distance
        self.graph_output_path = config.graph_output_path

        self.G = None
        self.prot_node_tsv = None


    def run(self, input_fasta: str):
        #Step 1: Cluster Proteins ~ make ze nodes!
        self.prot_node_tsv = mmseq_cluster(input_fasta, self.threshold, self.evalue)
        return self.run_nocluster(self.prot_node_tsv)

    def run_nocluster(self, input_tsv: str):
        G = nx.DiGraph()

        #Step 2: Create Protein Chains ~ make ze edges!
        graph_builder = GraphBuilder(input_tsv, self.input_format, self.operonic_distance)
        graph_builder.run()
        G.add_nodes_from(graph_builder.nodes)
        G.add_edges_from(graph_builder.weighted_edges)

        self.G = G
        return G
    
    def save_graph(self):
        print(f"Saving graph to {self.graph_output_path}")
        nx.write_gml(self.G, self.graph_output_path)
    
    def draw(self, G):
        nx.draw(G, with_labels=True, font_weight='bold') 

    #Function that takes in G, desired node, and degree and returns the subgraph
    def get_subgraph(self, G, node, degree):
        return nx.ego_graph(G, node, radius=degree)
