#Python script to set parameters and run graph_constructor.py
import src
from src.graph_constructor.py import GraphConstructor

#Set parameters
threshold = 0.5
evalue = 0.001
fasta = 'data/interim/prot_clust_test.fasta'

#Run graph constructor
graph_constructor = GraphConstructor(threshold, evalue)
G = graph_constructor.run(fasta)

#Save graph
nx.write_gml(G, 'data/interim/prot_clust_test.gml')
