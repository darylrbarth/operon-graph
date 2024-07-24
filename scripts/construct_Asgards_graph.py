#Python script to set parameters and run graph_constructor.py
import src
from src.graph_constructor import GraphConstructor
import networkx as nx
import matplotlib.pyplot as plt

#Set parameters
threshold = 0.6
evalue = 1.000E-05
fasta = 'data/All_Asgards_042423.faa'
graph_output = 'data/prot_clust_test.fasta.gml'

#Already run mmseqs2? Put in tsv here: 
#tsv = 'data/All_Asgards_wTMs_Oct212023_mmseqsDB_clu.tsv'

#Run graph constructor
graph_constructor = GraphConstructor(threshold, evalue)
G = graph_constructor.run(fasta)
#G = graph_constructor.run_nocluster(tsv)


#Save graph
nx.write_gml(G, graph_output)

#Save graph as png
#pos = nx.spring_layout(G)
#nx.draw_networkx_nodes(G, pos)
#nx.draw_networkx_edges(G, pos)
#nx.draw_networkx_labels(G, pos)
#plt.show()


#Draw graph
#graph_constructor.draw(G)

print('Done!')