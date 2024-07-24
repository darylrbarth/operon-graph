
from src.config import OperonGraphConfig

asgard_config1 = OperonGraphConfig(
    threshold=0.5,
    evalue=1.000E-05,
    operonic_distance=10,
    input_format='Guaymas_fasta_temp_scaffold',
    graph_output_path='data/graphs/All_Asgards_042423.gml',
)
