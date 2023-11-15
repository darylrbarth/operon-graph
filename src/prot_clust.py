import subprocess

def mmseq_cluster(fasta, thresh, eval):
    filename = fasta.split(".")[0]
    db_name = f'data/interim/{filename}_mmseqsDB' 
    subprocess.run(['mmseqs', 'createdb', fasta, db_name])
    clu_name = f'{db_name}_clu'
    #Make a tmp directory 
    subprocess.run(['mkdir', 'data/interim/tmp'])
    #Cluster the input protein sequences 
    subprocess.run(['mmseqs', 'cluster', db_name, clu_name, 'data/interim/tmp', '--min-seq-id', thresh, '-e', eval, '--cluster-mode', '2', '--threads', '20'])
    #Make the tsv
    tsv_name=f'data/{filename}_mmseqsDB_clu.tsv'
    subprocess.run(['mmseqs', 'createtsv', db_name, db_name, clu_name, tsv_name])
    #Remove the now unnecessary files and folders
    subprocess.run(['rm', 'data/interim/{db_name}*'])
    subprocess.run(['rm', '-r', 'data/interim/tmp'])
    return tsv_name 

def get_prot_names(tsv):
    nodes = set()
    with open(tsv, 'r') as f:
        for line in f:
            nodes.add(line.split('\t')[0])
    return nodes