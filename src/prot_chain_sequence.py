#Read in tsv and make hash table connecting cluster heads to cluster members pairwise
# {prot76:prot1, prot985:prot1, prot1:prot1}
#Also make hash table connecting cluster heads to a list of cluster members
# {prot1:[prot76, prot985, prot1]}

class EdgeConstructor():
    def __init__(self, tsv):
        self.tsv = tsv
        self.member_to_cluster = {}
        self.head_to_members = {}
        self.scaffolds = set()
        self.scaffolds_to_members = {}
        self.positions = {}
        self.ordered_scaffolds = {}
        self.scaffold_pairs = []
        self.edges = []
        self.edge_weights = {}
        self.weighted_edges = []

    def run(self):
        self.prot_hash()
        self.scaffold()
        self.ordered_scaffold_members()
        self.scaffold_pairs()
        self.create_edges()
        self.calculate_edge_weights()
        self.weighted_edges()

    #Make the dictionary of cluster heads to cluster members
    def prot_hash(self):
        #member_to_cluster = {}
        #head_to_members = {}
        with open(self.tsv, 'r') as f:
            for line in f:
                cluster_head, cluster_member = line.split('\t')
                self.member_to_cluster[cluster_member] = cluster_head
                if cluster_head in self.head_to_members:
                    self.head_to_members[cluster_head].add(cluster_member)
                else:
                    self.head_to_members[cluster_head] = set(cluster_member)
        return 
        #return member_to_cluster, head_to_members

    #Parse through cluster members and make a list of proteins that share the same scaffold
    def scaffold(self):
        members = self.member_to_cluster.keys()
        for member in members: 
            #Remove the last '_' and everything after it
            scaffold = member.split('_')[:-1]
            member_position = member.split('_')[-1]
            #scaffolds.add(scaffold)
            if scaffold in self.scaffolds_to_members:
                self.scaffolds_to_members[scaffold].add(member)
                self.positions[scaffold].append(member_position)
            else:
                self.scaffolds_to_members[scaffold] = set(member)
                self.positions[scaffold] = [member_position]
        return

    #Sort the members of each scaffold by position
    def ordered_scaffold_members(self):
        #scaffolds_to_members, positions = scaffold(tsv)
        #ordered_scaffolds = {}
        for scaffold in self.scaffolds_to_members.keys():
            scaffold_positions = self.positions[scaffold].sort()
            for member in self.scaffolds_to_members[scaffold]:
                position = member.split('_')[-1]
                index = scaffold_positions.index(position)
                scaffold_positions[index] = member
            self.ordered_scaffolds[scaffold] = scaffold_positions
        return 


    def scaffold_pairs(self):   
        #list of tuples of scaffold pairs, ie [(prot1, prot2), (prot2, prot3), (prot3, prot4)]
        #scaffold_pairs = []
        for scaffold in self.ordered_scaffolds.keys():
            member_list = self.ordered_scaffolds[scaffold]
            #Make a tuple of each member and the member after it in the list
            for i in range(len(member_list)-1):
                self.scaffold_pairs.append((member_list[i], member_list[i+1]))
        return


    #Make a list of edges between cluster heads by replacing the scaffold pairs with the cluster heads
    def create_edges(self):
        #edges = []
        #scaffold_pairs = scaffold_pairs(tsv)
        #member_to_cluster = prot_hash(tsv)[0]
        for pair in self.scaffold_pairs:
            #Replace the scaffold pair with the cluster head pair
            self.edges.append((self.member_to_cluster[pair[0]], self.member_to_cluster[pair[1]]))
        #match the scaffold pairs to the cluster heads
        return

    #Add edge weights to the edges
    #Syntax looks like this: (2, 3, {'weight': 3.1415})
    def calculate_edge_weights(self):
        #edges = create_edges(tsv)
        #hash table of edges and their weights
        #edge_weights = {}
        #head_to_members = prot_hash(tsv)[1]
        for edge in self.edges:
            num_members_n1 = len(self.head_to_members[edge[0]])
            num_members_n2 = len(self.head_to_members[edge[1]])
            num_possible_connections = num_members_n1 * num_members_n2
            if edge in self.edge_weights:
                self.edge_weights[edge] += 1/num_possible_connections
            else:
                self.edge_weights[edge] = 1/num_possible_connections
        return 


    def weighted_edges(self):
        #edge_weights = calculate_edge_weights(tsv)
       # weighted_edges = []
        for edge in self.edge_weights.keys():
            self.weighted_edges.append((edge[0], edge[1], {'weight': self.edge_weights[edge]}))
        return 










#Record pairs as a list of tuples indicated by the indices of the nodes in the tsv, ie [(1,2), (7,84)]