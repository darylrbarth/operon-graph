#Dis be 
import src
from src.prot_clust.py import mmseq_cluster, get_prot_names
from src.prot_chain_sequence.py import EdgeConstructor
import networkx as nx 

class GraphConstructor():
    def __init__(self, threshold, evalue):
        self.threshold = threshold
        self.evalue = evalue

    def run(self, input_fasta):
        G = nx.DiGraph()

        #Step 1: Cluster Proteins ~ make ze nodes!
        prot_node_tsv = mmseq_cluster(input_fasta, self.threshold, self.evalue)
        prot_node_names = get_prot_names(prot_node_tsv)
        #Add nodes to graph
        G.add_nodes_from(prot_node_names)

        #Step 2: Create Protein Chains (list of )
        edge_constructor = EdgeConstructor(prot_node_tsv)
        edge_constructor.run()

        #Add edges to graph
        G.add_edges_from(edge_constructor.weighted_edges)
        #Save graph to class instance
        self.G = G
        return G

    def draw(self, G):
        nx.draw(G, with_labels=True, font_weight='bold') 