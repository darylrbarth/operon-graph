#Read in tsv and make hash table connecting cluster heads to cluster members pairwise
# {prot76:prot1, prot985:prot1, prot1:prot1}
#Also make hash table connecting cluster heads to a list of cluster members
# {prot1:[prot76, prot985, prot1]}

def prot_hash(tsv):
    member_to_cluster = {}
    head_to_members = {}
    with open(tsv, 'r') as f:
        for line in f:
            cluster_head, cluster_member = line.split('\t')
            member_to_cluster[cluster_member] = cluster_head
            if cluster_head in head_to_members:
                head_to_members[cluster_head].add(cluster_member)
            else:
                head_to_members[cluster_head] = set(cluster_member)
    return member_to_cluster, head_to_members

#Parse through cluster members and make a list of proteins that share the same scaffold
def scaffold(tsv):
    members = prot_hash(tsv)[0].keys()
    #scaffolds = set()
    scaffolds_to_members = {}
    for member in members: 
        #Remove the last '_' and everything after it
        scaffold = member.split('_')[:-1]
        #scaffolds.add(scaffold)
        if scaffold in scaffolds_to_members:
            scaffolds_to_members[scaffold].add(member)
        else:
            scaffolds_to_members[scaffold] = set(member)
    #return scaffolds, scaffolds_to_members
    return scaffolds_to_members

def scaffold_pairs(tsv):
    scaffolds_to_members = scaffold(tsv)
    scaffold_pairs = set()
    for scaffold in scaffolds_to_members.keys():
        members = scaffolds_to_members[scaffold]
        for member in members:
            for other_member in members:
                if member != other_member:
                    scaffold_pairs.add((member, other_member))
    return scaffold_pairs


        

    
    










#Record pairs as a list of tuples indicated by the indices of the nodes in the tsv, ie [(1,2), (7,84)]

def get_prot_names(tsv):
    nodes = set()
    with open(tsv, 'r') as f:
        for line in f:
            nodes.add(line.split('\t')[0])
    return nodes