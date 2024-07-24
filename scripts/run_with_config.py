
import time
import sys

from src.operon_graph import OperonGraph
from src.config import OperonGraphConfig

"""
Simple script to
- Instantiate OperonGraph, mandatorily taking in a json config file
- Build the graph on Guaymas data
- Save the graph to the specified path
- TODO what else needed in the standard workflow?

Usage example (from repo root dir): python scripts/run_with_config.py config/sample.json
"""

if __name__ == "__main__":
    start_time = time.time()

    config_path = sys.argv[1]
    config = OperonGraphConfig.from_json(config_path)

    input_fasta = '/stor/work/Marcotte/project/drbarth/plastics/reference/Guaymas_mining/Guaymas2020/05.reference/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.faa' # TODO: Choose smaller test file, possibly E. coli
    operon_graph = OperonGraph(config)
    G = operon_graph.run(input_fasta)

    #Save graph
    operon_graph.save_graph()

    elapsed_time = time.time() - start_time
    print(f"Execution time: {elapsed_time} seconds")
    print('Done!')