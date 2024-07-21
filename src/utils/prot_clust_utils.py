
import subprocess

def mmseq_cluster(fasta, thresh, eval, skip_if_done=True, cleanup=False):
    """
    Runs MMSeqs via calls to subprocess, returns protein clusters as tsv
    File paths assume running from repo root directory.

    If skip_if_done is True, will check if the tsv file already exists and skip if it does.
    """
    filename = fasta.split("/")[1].split(".")[0]
    db_name = f'data/interim/{filename}_mmseqsDB'
    clu_name = f'data/interim/clusters/{filename}_clu'
    tsv_name=f'data/{filename}_mmseqsDB_clu.tsv'
    print(subprocess.run(['pwd']))

    #Check if the tsv file already exists
    if skip_if_done:
        try:
            with open(tsv_name, 'r') as f:
                print(f'{tsv_name} already exists, skipping clustering')
                return tsv_name
        except FileNotFoundError:
            pass

    #Make DB and temporary directories
    subprocess.run(['mkdir', 'data/interim/'])
    subprocess.run(['mmseqs', 'createdb', fasta, db_name])
    subprocess.run(['mkdir', 'data/interim/clusters'])
    subprocess.run(['mkdir', 'data/interim/tmp'])

    #Cluster the input protein sequences 
    subprocess.run(['mmseqs', 'cluster', db_name, clu_name, 'data/interim/tmp', '--min-seq-id', str(thresh), '-e', str(eval), '--cluster-mode', '2', '--threads', '20'])
    
    #Make the tsv
    subprocess.run(['mmseqs', 'createtsv', db_name, db_name, clu_name, tsv_name])

    #Cleanup
    if cleanup:
        subprocess.run(['rm', f'{db_name}*'])
        subprocess.run(['rm', '-r', 'data/interim'])

    print(f'Clustered {fasta} with threshold {thresh} and evalue {eval} and saved as {tsv_name}')
    return tsv_name 