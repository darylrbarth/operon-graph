
from src.config import OperonGraphConfig

asgard_config1 = OperonGraphConfig(
    threshold=0.5,
    evalue=1.000E-05,
    operonic_distance=20,
    input_format='Guaymas_fasta_temp_scaffold',
    graph_output_path='data/graphs/20_operonic_updated_undirected_All_Asgards_042423.gml',
)
