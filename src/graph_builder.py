#Read in tsv and make hash table connecting cluster heads to cluster members pairwise
# {prot76:prot1, prot985:prot1, prot1:prot1}
#Also make hash table connecting cluster heads to a list of cluster members
# {prot1:[prot76, prot985, prot1]}

from typing import List, Dict, Set

"""
Core logic for constructing the nodes, edges, and edge weights for an OperonGraph
"""

class GraphBuilder():
    def __init__(self, tsv, input_format, operonic_distance):
        self.tsv = tsv
        self.input_format = input_format
        self.operonic_distance = operonic_distance
        self.member_to_cluster = {}
        self.head_to_members = {}
        self.scaffolds = set()
        self.scaffolds_to_members = {}
        self.positions : Dict[str, List[int]] = {}
        self.ordered_scaffolds = {}
        self.scaffold_pairs = []
        self.edges = []
        self.edge_properties = {}
        self.weighted_edges = []
        self.nodes = set()

    def run(self):
        self.get_nodes()
        print(f'Number of nodes: {len(self.nodes)}, example: {list(self.nodes)[:5]}')
        self.prot_hash()
        print(f'Number of cluster heads: {len(self.head_to_members)}, example: {list(self.head_to_members.keys())[:5]}')
        self.scaffold()
        print(f'Number of scaffolds: {len(self.scaffolds_to_members)}, example: {list(self.scaffolds_to_members.keys())[0]}, {list(self.scaffolds_to_members.values())[0][:5]}')
        self.ordered_scaffold_members()
        print(f'Number of ordered scaffolds: {len(self.ordered_scaffolds)}, example: {list(self.ordered_scaffolds.keys())[0]}, {list(self.ordered_scaffolds.values())[0][:5]}')
        self.make_scaffold_pairs()
        print(f'Number of scaffold pairs: {len(self.scaffold_pairs)}, example: {self.scaffold_pairs[:5]}')

        # Edge construction
        
        self.create_edges()
        print(f'Number of edges: {len(self.edges)}, example: {self.edges[:5]}')
        self.calculate_edge_properties()
        print(f'Number of edge weights: {len(self.edge_properties)}, example: {list(self.edge_properties.keys())[:5]}')
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
                if self.input_format == 'Guaymas_fasta_temp_scaffold':
                    cluster_head = cluster_head.split('_')[:-1]
                    cluster_head = '_'.join(cluster_head)
                    #print(cluster_head)
                    cluster_member = cluster_member.split('_')[:-1]
                    cluster_member = '_'.join(cluster_member)
                    #print(cluster_member)
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
            #'Guaymas_fasta' that looks like this: D4944_C32_H1_scaffold_1932_9 or 'Guaymas_fasta_temp_scaffold' that looks like this: D4944_C32_H1_scaffold_1932_9_100 and D4944_C32_H1-scaffold_1932_9_100
            #if self.input_format == 'Guaymas_fasta':
            scaffold = member.split('_')[:-1]
            scaffold = '_'.join(scaffold)
            member_position = member.split('_')[-1]
            #print(f'Scaffold: {scaffold}, member_position: {member_position}')
            #elif self.input_format == 'Guaymas_fasta_temp_scaffold':
            #    #print(member)
            #    scaffold = member.split('_')[:-2]
            #    #print(scaffold)
            #    scaffold = '_'.join(scaffold)
            #    member_position = member.split('_')[-2]
                #print(f'Scaffold: {scaffold}, member_position: {member_position}')
            #Drop the \n
            member_position = member_position.split('\n')[0]
            #Add scaffolds to members list! 
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
        """
        Connect all pairs of proteins within self.operonic distance (Warning: n^2 operation)
        Operonic distance controls the density of the graph
        """
        #list of tuples of scaffold pairs, ie [(prot1, prot2), (prot2, prot3), (prot3, prot4)]
        for scaffold in self.ordered_scaffolds.keys():
            member_list = self.ordered_scaffolds[scaffold]
            #Make a tuple of each member and the member after it, if the member after has a position within operonic_distance of each other
            for i in range(len(member_list)):
                member_pos = member_list[i].split('_')[-1]
                for j in range(i+1, len(member_list)):
                    next_member_pos = member_list[j].split('_')[-1]
                    if int(next_member_pos) - int(member_pos) > self.operonic_distance:
                        break 
                    self.scaffold_pairs.append((member_list[i], member_list[j]))
        return
    
    # def make_adjacent_scaffold_pairs(self):
    #     """
    #     Old function for generating scaffold pairs purely by adjacency in operon
    #     """
    #     #list of tuples of scaffold pairs, ie [(prot1, prot2), (prot2, prot3), (prot3, prot4)]
    #     for scaffold in self.ordered_scaffolds.keys():
    #         print(scaffold)
    #         member_list = self.ordered_scaffolds[scaffold]
    #         print(member_list)
    #         #Make a tuple of each member and the member after it, if the member after has a position within operonic_distance of each other
    #         for i in range(len(member_list)-1):
    #             member_pos = member_list[i].split('_')[-1]
    #             next_member_pos = member_list[i+1].split('_')[-1]
    #             if int(next_member_pos) - int(member_pos) <= self.operonic_distance:
    #             #print(member_list[i], member_list[i+1])
    #                 self.scaffold_pairs.append((member_list[i], member_list[i+1]))
    #     return


    # --------------- Edge construction ----------------

    #Make a list of edges between cluster heads by replacing the scaffold pairs with the cluster heads
    def create_edges(self):
        for pair in self.scaffold_pairs:
            #Replace the scaffold pair with the cluster head pair
            self.edges.append((self.member_to_cluster[pair[0]], self.member_to_cluster[pair[1]]))
        #match the scaffold pairs to the cluster heads
        return
    

    #Add edge weights to the edges
    #Syntax looks like this: (2, 3, {'abs_weight': 3.1415, norm_weight: 0.18235})
    def calculate_edge_properties(self):
        for edge in self.edges:
            num_members_n1 = len(self.head_to_members[edge[0]])
            num_members_n2 = len(self.head_to_members[edge[1]])
            num_possible_connections = num_members_n1 * num_members_n2
            if edge in self.edge_properties:
                edge_properties = self.edge_properties[edge]
                edge_properties['abs_weight'] += 1
                edge_properties['norm_weight'] += 1/num_possible_connections
            else:
                edge_properties = {
                    'abs_weight': 1,
                    'norm_weight': 1/num_possible_connections,
                }
                self.edge_properties[edge] = edge_properties
        return 


    def make_weighted_edges(self):
        for edge in self.edge_properties.keys():
            self.weighted_edges.append((edge[0], edge[1], self.edge_properties[edge]))
        return 
