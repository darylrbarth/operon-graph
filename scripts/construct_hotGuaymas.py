#Python script to set parameters and run graph_constructor.py
import src
from src.graph_constructor import GraphConstructor
import networkx as nx
import matplotlib.pyplot as plt
import time
start_time = time.time()

#Set parameters
threshold = 0.5
evalue = 1.000E-05
fasta = '/stor/work/Marcotte/project/drbarth/plastics/reference/Guaymas_mining/Guaymas2020/05.reference/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.faa'
graph_output = '../data/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.gml'

#Already run mmseqs2? Put in tsv here: 
tsv = '../data/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C_mmseqsDB_clu_50.tsv'

#Run graph constructor
graph_constructor = GraphConstructor(threshold, evalue)
#G = graph_constructor.run(fasta)
#If mmseqs clustering is already done: uncomment the following line
G = graph_constructor.run_nocluster(tsv)


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
elapsed_time = time.time() - start_time
print(f"Execution time: {elapsed_time} seconds")
print('Done!')