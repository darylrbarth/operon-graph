#Dis be 
import src
from src.prot_clust import mmseq_cluster
from src.prot_chain_sequence import EdgeConstructor
import networkx as nx 

class GraphConstructor():
    def __init__(self, threshold, evalue, input_format):
        self.threshold = threshold
        self.evalue = evalue
        self.input_format = input_format

    def run(self, input_fasta):
        G = nx.DiGraph()

        #Step 1: Cluster Proteins ~ make ze nodes!
        self.prot_node_tsv = mmseq_cluster(input_fasta, self.threshold, self.evalue)

        #Step 2: Create Protein Chains (list of )
        edge_constructor = EdgeConstructor(self.prot_node_tsv, self.input_format)
        edge_constructor.run()

        #Add nodes and edges to graph
        G.add_nodes_from(edge_constructor.nodes)
        G.add_edges_from(edge_constructor.weighted_edges)
        #Save graph to class instance
        self.G = G
        return G
    
    def run_nocluster(self, input_tsv):
        G = nx.DiGraph()
        #Step 2: Create Protein Chains (list of )
        edge_constructor = EdgeConstructor(input_tsv, self.input_format)
        edge_constructor.run()
        #Add nodes and edges to graph
        G.add_nodes_from(edge_constructor.nodes)
        G.add_edges_from(edge_constructor.weighted_edges)
        #Save graph to class instance
        self.G = G
        return G

    def draw(self, G):
        nx.draw(G, with_labels=True, font_weight='bold') 

    #Function that takes in G, desired node, and degree and returns the subgraph
    def get_subgraph(self, G, node, degree):
        return nx.ego_graph(G, node, radius=degree)