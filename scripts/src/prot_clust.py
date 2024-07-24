import subprocess

def mmseq_cluster(fasta, thresh, eval):
    #print working directory
    print(subprocess.run(['pwd']))
    filename = fasta.split("/")[-1].split(".")[0]
    subprocess.run(['mkdir', '../data/interim/'])
    db_name = f'../data/interim/{filename}_mmseqsDB'
    subprocess.run(['mmseqs', 'createdb', fasta, db_name])
    subprocess.run(['mkdir', '../data/interim/clusters'])
    clu_name = f'../data/interim/clusters/{filename}_clu'
    #Make a tmp directory 
    subprocess.run(['mkdir', '../data/interim/tmp'])
    #Cluster the input protein sequences 
    subprocess.run(['mmseqs', 'cluster', db_name, clu_name, '../data/interim/tmp', '--min-seq-id', str(thresh), '-e', str(eval), '--cluster-mode', '2', '--threads', '20'])
    #Make the tsv
    tsv_name=f'../data/{filename}_mmseqsDB_clu.tsv'
    subprocess.run(['mmseqs', 'createtsv', db_name, db_name, clu_name, tsv_name])
    #Remove the now unnecessary files and folders
    #subprocess.run(['rm', f'{db_name}*'])
    #subprocess.run(['rm', '-r', 'data/interim'])
    print(f'Clustered {fasta} with threshold {thresh} and evalue {eval} and saved as {tsv_name}')
    return tsv_name 