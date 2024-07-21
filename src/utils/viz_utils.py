

def convert_graph_to_tsv(G):
    """
    Currently slow, there's likely a faster way to use the GML directly
    """
    # Open a file to write the node data
    output_nodes_tsv = 'data/test_GRAPH_nodes.tsv'
    output_edges_tsv = 'data/test_GRAPH_edges.tsv'

    with open(output_nodes_tsv, 'w') as node_file:
        node_file.write("NodeId\tlabel\n")  # Modify based on your node attributes
        for node, data in G.nodes(data=True):
            # Write node and attributes to the file, ensure attributes match what's in your graph
            node_file.write(f"{node}\t{data.get('label', '')}\n")
    print(f'saved nodes to {output_nodes_tsv}')
    # Open a file to write the edge data
    with open(output_edges_tsv, 'w') as edge_file:
        edge_file.write("Source\tTarget\tWeight\n")  # Modify if you have different or additional attributes
        for source, target, data in G.edges(data=True):
            # Write edge and attributes to the file, ensure attributes match what's in your graph
            edge_file.write(f"{source}\t{target}\t{data.get('weight', '')}\n")
    print(f'saved edges to {output_edges_tsv}')
