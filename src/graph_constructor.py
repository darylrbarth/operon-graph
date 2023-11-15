#Dis be 
import src
from src.prot_clust.py import mmseq_cluster, get_prot_names
import networkx as nx 

class GraphConstructor():
    def __init__(self):
        pass

    def run():
        G = nx.DiGraph()

        #Step 1: Cluster Proteins ~ make ze nodes!
        prot_node_tsv = mmseq_cluster()
        prot_node_names = get_prot_names(prot_node_tsv)
        #Add nodes to graph
        G.add_nodes_from(prot_node_names)

        #Step 2: Create Protein Chains
        fasta_parse()
        #Step 3: Determine Edge Weights 

        return 