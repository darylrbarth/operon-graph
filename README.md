# operon-graph

Given a dataset of proteins, we build an "Operon Graph" with nodes as clusters of similar proteins, and weighted edges indicating how frequently proteins from the given clusters are found colocated on operons within the dataset. This can help find operationally compatible proteins, and map our the more evolutionarily conserved operons.

###
Data extracted from each protein
- Scaffold
- Position on scaffold
- Whether forward or reverse

### Todos

Functionality
PR 1
- Get tests working
- Assure edge weights working, merge

PR 2
- Protein dataclass: scaffold, location, forward and reverse strand, etc.
- Data adapters for Ecoli, prodigal, truncated prodigal, maybe others (support multiple formats)

Later
- Viz tools?
- Notebook where we walk through E. coli

### Example settings (scratch)
threshold = 0.5
evalue = 1.000E-05
#input string format
input_format = 'Guaymas_fasta_temp_scaffold' #Options are 'Guaymas_fasta' that looks like this: D4944_C32_H1_scaffold_1932_9 or 'Guaymas_fasta_temp_scaffold' that looks like this: D4944_C32_H1_scaffold_1932_9_100 and D4944_C32_H1-scaffold_1932_9_100
#fasta = '/stor/work/Marcotte/project/drbarth/plastics/reference/Guaymas_mining/Guaymas2020/05.reference/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_greaterthan80C.faa'
graph_output = 'data/GRAPH_ALL_Guaymas2020_hottest_clu30_May62024.gml'

#Already run mmseqs2? Put in tsv here: 
tsv = 'data/Guaymas2020_Scaffolds_Bins_deduplicated_hottest_mmseqs_clu30.tsv'

Run with python -m call_file_name
