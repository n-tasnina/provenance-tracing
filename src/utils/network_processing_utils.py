def get_edges_and_rank_from_list_of_paths(paths):
    '''
    This will return a list of edges givena list of paths.
    Input: paths: list of lists. Each inner list contains the nodes along a path.
    output: a dict where edge is a key and list of rank(s) of the path(s) that edge belong to is the value

    CAUTION: Thi rank is meaningful only when we give all top k paths in an ordered way as input.
    '''
    all_edges = {}
    cur_rank=1
    for path_nodes in paths:
        # if cur_rank==648:
        #     print(path_nodes)
        new_edges = path_nodes_to_set_of_edges(path_nodes)
        for edge in new_edges:
            if edge not in all_edges:
                all_edges[edge]=[cur_rank]
            else:
                all_edges[edge].append(cur_rank)
        cur_rank+=1
    return all_edges


def get_edges_from_list_of_paths(paths):
    '''
    This will return a list of edges givena list of paths.
    Input: paths: list of lists. Each inner list contains the nodes along a path.
    output: a dict where edge is a key and list of rank(s) of the path(s) that edge belong to is the value
    '''
    all_edges = []
    for path_nodes in paths:
        new_edges = path_nodes_to_set_of_edges(path_nodes)
        for edge in new_edges:
            all_edges.append(edge)
    return set(all_edges)


def path_nodes_to_set_of_edges(path_nodes):
    '''
    convert a list of nodes along a path into a set of edges.
    Return a set of tuples where each tuple in an edge.
    '''
    edges = set()
    for i in range(0, len(path_nodes)-1):
        edges.add(tuple([path_nodes[i],path_nodes[i+1]]))
    return edges

def path_nodes_to_ordered_edges(path_nodes):
    '''
    convert a list of nodes along a path into a list of edges.
    Return a set of tuples where each tuple in an edge.
    '''
    edges = []
    for i in range(0, len(path_nodes)-1):
        edges.append(tuple([path_nodes[i],path_nodes[i+1]]))
    return edges

