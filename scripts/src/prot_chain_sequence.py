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
        self.nodes = set()

    def run(self):
        self.get_nodes()
        print(f'Number of nodes: {len(self.nodes)}')
        self.prot_hash()
        print(f'Number of cluster heads: {len(self.head_to_members)}')
        self.scaffold()
        print(f'Number of scaffolds: {len(self.scaffolds_to_members)}, example: {list(self.scaffolds_to_members.keys())[0]}, {list(self.scaffolds_to_members.values())[0][:5]}')
        self.ordered_scaffold_members()
        print(f'Number of ordered scaffolds: {len(self.ordered_scaffolds)}, example: {list(self.ordered_scaffolds.keys())[0]}, {list(self.ordered_scaffolds.values())[0][:5]}')
        self.make_scaffold_pairs()
        print(f'Number of scaffold pairs: {len(self.scaffold_pairs)}, example: {self.scaffold_pairs[:5]}')
        self.create_edges()
        print(f'Number of edges: {len(self.edges)}, example: {self.edges[:5]}')
        self.calculate_edge_weights()
        print(f'Number of edge weights: {len(self.edge_weights)}, example: {list(self.edge_weights.keys())[:5]}')
        self.make_weighted_edges()
        print(f'Number of weighted edges: {len(self.weighted_edges)}, example: {self.weighted_edges[:5]}')

    def get_nodes(self):
        with open(self.tsv, 'r') as f:
            for line in f:
                self.nodes.add(line.split('\t')[0])
        return
    
    #Make the dictionary of cluster heads to cluster members
    def prot_hash(self):
        with open(self.tsv, 'r') as f:
            for line in f:
                cluster_head, cluster_member = line.split('\t')
                #Remove the \n
                cluster_member = cluster_member[:-1]
                self.member_to_cluster[cluster_member] = cluster_head
                if cluster_head in self.head_to_members:
                    self.head_to_members[cluster_head].add(cluster_member)
                else:
                    self.head_to_members[cluster_head] = set(cluster_member)
        return 

    #Parse through cluster members and make a list of proteins that share the same scaffold
    def scaffold(self):
        members = self.member_to_cluster.keys()
        for member in members: 
            #Remove the last '_' and everything after it
            scaffold = member.split('_')[:-1]
            scaffold = '_'.join(scaffold)
            member_position = member.split('_')[-1]
            #Drop the \n
            member_position = member_position.split('\n')[0]
            #scaffolds.add(scaffold)
            if scaffold in self.scaffolds_to_members.keys():
                self.scaffolds_to_members[scaffold].append(member)
                self.positions[scaffold].append(member_position)
            else:
                self.scaffolds_to_members[scaffold] = [member]
                self.positions[scaffold] = [member_position]
        return

    #Sort the members of each scaffold by position
    def ordered_scaffold_members(self):
        for scaffold in self.scaffolds_to_members.keys():
            self.positions[scaffold].sort()
            scaffold_positions = [None] * len(self.positions[scaffold])
            for member in self.scaffolds_to_members[scaffold]:
                #print(member)
                position = member.split('_')[-1]
                #print(position)
                index = self.positions[scaffold].index(position)
                scaffold_positions[index] = member
            self.ordered_scaffolds[scaffold] = scaffold_positions
        return 


    def make_scaffold_pairs(self):   
        #list of tuples of scaffold pairs, ie [(prot1, prot2), (prot2, prot3), (prot3, prot4)]
        for scaffold in self.ordered_scaffolds.keys():
            member_list = self.ordered_scaffolds[scaffold]
            #Make a tuple of each member and the member after it in the list
            for i in range(len(member_list)-1):
                self.scaffold_pairs.append((member_list[i], member_list[i+1]))
        return


    #Make a list of edges between cluster heads by replacing the scaffold pairs with the cluster heads
    def create_edges(self):
        for pair in self.scaffold_pairs:
            #Replace the scaffold pair with the cluster head pair
            self.edges.append((self.member_to_cluster[pair[0]], self.member_to_cluster[pair[1]]))
        #match the scaffold pairs to the cluster heads
        return

    #Add edge weights to the edges
    #Syntax looks like this: (2, 3, {'weight': 3.1415})
    def calculate_edge_weights(self):
        for edge in self.edges:
            num_members_n1 = len(self.head_to_members[edge[0]])
            num_members_n2 = len(self.head_to_members[edge[1]])
            num_possible_connections = num_members_n1 * num_members_n2
            if edge in self.edge_weights:
                self.edge_weights[edge] += 1/num_possible_connections
            else:
                self.edge_weights[edge] = 1/num_possible_connections
        return 


    def make_weighted_edges(self):
        for edge in self.edge_weights.keys():
            self.weighted_edges.append((edge[0], edge[1], {'weight': self.edge_weights[edge]}))
        return 
