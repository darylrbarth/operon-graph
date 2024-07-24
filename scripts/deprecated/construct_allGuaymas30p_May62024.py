#Python script to set parameters and run graph_constructor.py
import src
from scripts.src.operon_graph import GraphConstructor
import networkx as nx
import matplotlib.pyplot as plt
import time
start_time = time.time()

#Set parameters
threshold = 0.5
evalue = 1.000E-05
#input string format
input_format = 'Guaymas_fasta_temp_scaffold' #Options are 'Guaymas_fasta' that looks like this: D4944_C32_H1_scaffold_1932_9 or 'Guaymas_fasta_temp_scaffold' that looks like this: D4944_C32_H1_scaffold_1932_9_100 and D4944_C32_H1-scaffold_1932_9_100
operonic_distance = 10 
#fasta = '/stor/work/Marcotte/project/drbarth/plastics/reference/Guaymas_mining/Guaymas2020/05.reference/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.faa'
graph_output = '../data/GRAPH_ALL_Guaymas2020_hottest_clu30_May62024.gml'


#Already run mmseqs2? Put in tsv here: 
tsv = '../data/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_mmseqs_clu30.tsv'

#Run graph constructor
graph_constructor = GraphConstructor(threshold, evalue, input_format, operonic_distance)
#G = graph_constructor.run(fasta)
#If mmseqs clustering is already done: uncomment the following line
G = graph_constructor.run_nocluster(tsv)


#Save graph
nx.write_gml(G, graph_output)

elapsed_time = time.time() - start_time
print(f"Time to make graph: {elapsed_time} seconds")
print(f'Graph output here: {graph_output}')

#Save graph as tsv files

import networkx as nx

# Load the graph from a GML file
G = nx.read_gml(graph_output)

# Open a file to write the node data
output_nodes_tsv = '../data/nodes_GRAPH_ALL_Guaymas2020_hottest_clu30_May62024.tsv'
output_edges_tsv = '../data/edges_GRAPH_ALL_Guaymas2020_hottest_clu30_May62024.tsv'


with open(output_nodes_tsv, 'w') as node_file:
    node_file.write("NodeId\tlabel\n")  # Modify based on your node attributes
    for node, data in G.nodes(data=True):
        # Write node and attributes to the file, ensure attributes match what's in your graph
        node_file.write(f"{node}\t{data.get('label', '')}\n")
print(f'saved nodes to {output_nodes_tsv}')
# Open a file to write the edge data
with open(output_edges_tsv, 'w') as edge_file:
    edge_file.write("Source\tTarget\tWeight\n")  # Modify if you have different or additional attributes
    for source, target, data in G.edges(data=True):
        # Write edge and attributes to the file, ensure attributes match what's in your graph
        edge_file.write(f"{source}\t{target}\t{data.get('weight', '')}\n")
print(f'saved edges to {output_edges_tsv}')

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