
import networkx as nx
import matplotlib.pyplot as plt
import time

from src.operon_graph import OperonGraph
from src.config import OperonGraphConfig

def test_operon_graph():
    start_time = time.time()

    config = OperonGraphConfig(
        threshold=0.5,
        evalue=1.000E-05,
        input_format='Guaymas_fasta_temp_scaffold',
        graph_output='data/test_data/test_GRAPH_Guaymas2020_hottest_clu30.gml',
    )

    input_fasta = '/stor/work/Marcotte/project/drbarth/plastics/reference/Guaymas_mining/Guaymas2020/05.reference/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.faa' # TODO: Choose smaller test file, possibly E. coli
    operon_graph = OperonGraph(config)
    G = operon_graph.run(input_fasta)

    #Save graph
    operon_graph.save_graph()

    elapsed_time = time.time() - start_time
    print(f"Execution time: {elapsed_time} seconds")
    print('Done!')


if __name__ == '__main__':
    test_operon_graph()
