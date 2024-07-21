
from dataclasses import dataclass

"""
Parameters needed for constructing the operon graph

input_format options:
  - 'Guaymas_fasta': D4944_C32_H1_scaffold_1932_9
  - 'Guaymas_fasta_temp_scaffold': D4944_C32_H1_scaffold_1932_9_100 and D4944_C32_H1-scaffold_1932_9_100
"""

@dataclass
class OperonGraphConfig:
    threshold: float = 0.9
    evalue: float = 1.000E-05
    operonic_distance: int = 10
    input_format: str = 'Guaymas_fasta_temp_scaffold'
    graph_output_path: str = '../data/graphs/tmp_operon_graph.gml'

    def __post_init__(self):
        assert self.threshold > 0, 'Threshold must be greater than 0'
        assert self.evalue > 0, 'Evalue must be greater than 0'
        assert self.operonic_distance > 0, 'Operonic distance must be greater than 0'