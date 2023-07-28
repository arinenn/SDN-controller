import sys
import json
import argparse
import xml.etree.ElementTree as ElementTree
from math import acos, sin, cos, inf



# --- CONSTANTS ---

SOURCE_PATH = './input'
OUTPUT_PATH = './output'
namespaces = {
    'ns': "http://graphml.graphdrawing.org/xmlns",
    'xsi': "http://www.w3.org/2001/XMLSchema-instance",
    'schema': "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"
}  # header info from TopologyZOO



# --- SELECTORS, GETTERS, HELPERS ---

# return translation of key_name in attributes
def translate_key(key_name, attributes):
    # example "d32" translates into "Longitude"
    return attributes[key_name]['name']

# returns attributes of given topology which are necessary for computing
def is_needed_attributes(node):
    needed = {
        'id',
        'label',
        'Latitude',
        'Longitude'
    }
    for key in needed:
        if key not in node:
            return False
    return True

# returns array with ids of edges of node
def get_adjacent_edges_of_node(node_id, edges):
    # accepts id of node
    result = []
    for edge in edges:
        if edge['source'] == node_id:
            result.append(edge)
    return result

# search node by id in nodes
def get_node_by_id(node_id, nodes):
    for node in nodes:
        if node['id'] == node_id:
            return node
    # never reached if parser works correctly
    return None

# gets index in array of node with node_id
def get_node_index_by_id(node_id, nodes):
    for i in range(0, len(nodes)):
        if nodes[i]['id'] == node_id:
            return i

# returns adjacency matrix for given graph
def get_adjacency_matrix(root_id, nodes, edges):
    adjacency_matrix = [[inf]*len(nodes) for i in range(0,len(nodes))] # two-dim array
    vertex_ids = []  # for bijection between vertexes and ids of nodes
    # get position of root_id in nodes
    start = 0
    for i in range(0, len(nodes)):
        if nodes[i]['id'] == root_id:
            start = i
            break
    # fill matrix
    for i in range(start, len(nodes) + start):
        index = i % len(nodes)
        temp = i - start
        adjacency_matrix[temp][temp] = 0  # no delay in node
        # fill adjacency_matrix[index] row and memorize id
        vertex_ids.append(nodes[index]['id'])
        adjacent_edges = get_adjacent_edges_of_node(nodes[index]['id'], edges)
        for edge in adjacent_edges:
            adjacency_matrix[temp][get_node_index_by_id(edge['target'], nodes)] = edge['Delay']
    return adjacency_matrix, vertex_ids

# returns paths from root node to other nodes by spanning tree created by dijkstras algorithm
def get_main_paths(delay_table, optimum_connections, vertex_ids):
    # returns two-dim array with path, delay, type for every other node starting from root node
    main_paths = []
    size = len(optimum_connections)
    vertex_start = 0
    for vertex_end_ in range(vertex_start + 1, size):
        if optimum_connections[vertex_end_] == inf:
            # can't return
            continue
        path = [vertex_end_]
        vertex_end = vertex_end_
        while vertex_end != vertex_start:
            # while haven't return
            vertex_end = optimum_connections[path[-1]]
            path.append(vertex_end)  # memorize id of node to go through
        path = list(map(lambda p: vertex_ids[p], path[::-1])) # transform to ids and invert path
        # get path delay
        main_paths.append({
            'Type': 'main',
            'Path': path,
            'Delay': delay_table[vertex_end_]
        })
    return main_paths



# --- TOPOLOGY and GRAPH PARSERS ---

# returns array with fields as ids of parsed XML nodes
def parse_graph_nodes(xml_nodes, attributes):
    # accepts array of XML elements
    # requires id in xml nodes
    nodes = []
    banned_nodes_ids = []  # nodes without needed information
    for node in xml_nodes:
        node_id = node.attrib['id']
        # get data of node
        new_node = {}
        node_data = node.findall('ns:data', namespaces)
        for data in node_data:
            new_node[translate_key(data.attrib['key'], attributes)] = data.text
        # check if all needed information is present
        if is_needed_attributes(new_node):
            nodes.append(new_node)
        else:
            # ban node, save it's id
            banned_nodes_ids.append(node_id)
    return [nodes, banned_nodes_ids]

# returns array with fields as ids of parsed XML edges
def parse_graph_edges(xml_edges, attributes, banned_nodes_ids):
    # accepts array of XML elements
    # requires id in data nodes of given xml nodes
    edges = []
    for edge in xml_edges:
        # check if source or target is banned
        if (edge.attrib['source'] in banned_nodes_ids) or (edge.attrib['target'] in banned_nodes_ids):
            # skip current edge
            continue
        # add source and target
        new_edge = {
            'source': edge.attrib['source'],
            'target': edge.attrib['target']
        }
        edge_data = edge.findall('ns:data', namespaces)
        for data in edge_data:
            element = translate_key(data.attrib['key'], attributes)
            new_edge[element] = data.text
        edges.append(new_edge)  # save data of edge
    return edges

# returns parsed topology
def parse_topology(filename):
    # get links from file
    try:
        tree = ElementTree.parse(filename)
    except OSError as e:
        print(f"Unable to open {filename}", file=sys.stderr)
        exit(1)
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    root_keys = root.findall('ns:key', namespaces)
    root_graph = root.find('ns:graph', namespaces)
    graph_data = root_graph.findall('ns:data', namespaces)
    graph_nodes = root_graph.findall('ns:node', namespaces)
    graph_edges = root_graph.findall('ns:edge', namespaces)
    # get information for tags of nodes, fill attributes by root_keys, graph_data
    attributes = {}
    for key in root_keys:
        d = dict(key.attrib)
        attributes[d['id']] = {
            'name': d['attr.name']
        }
    for data in graph_data:
        # get attribute by data_id and append it to element with data_id
        attributes[data.attrib['key']]['value'] = data.text
    # attributes are finished, now parse graph
    [nodes, banned_nodes_ids] = parse_graph_nodes(graph_nodes, attributes)  # get nodes
    edges = parse_graph_edges(graph_edges, attributes, banned_nodes_ids)  # get edges
    return nodes, edges



# --- CALCULATORS ---

# calculates distance between two nodes by ids
def calculate_distance(first_node_id, second_node_id, nodes):
    earth_radius = 6371
    first_node = get_node_by_id(first_node_id, nodes)
    second_node = get_node_by_id(second_node_id, nodes)
    lat1 = float(first_node['Latitude'])
    lon1 = float(first_node['Longitude'])
    lat2 = float(second_node['Latitude'])
    lon2 = float(second_node['Longitude'])
    # celestial mechanics formulae
    return earth_radius * acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1))

# modifies weights of graph edges
def calculate_distance_and_delay(nodes, edges):
    # modifies edges
    delay_by_km = 4.8  # calculating in microseconds
    for node in nodes:
        adjacent_edges = get_adjacent_edges_of_node(node['id'], edges)
        for edge in adjacent_edges:
            edge['Distance'] = calculate_distance(edge['source'], edge['target'], nodes)
            edge['Delay'] = edge['Distance'] * delay_by_km



# --- ALGORITHMS ---

# returns tree and maximum Delay
def kruskals_algorithm(edges):
    # zero step of algorithm: sort edges by Delay
    sorted_edges = sorted(edges, key=lambda e: e['Delay'])
    unity = set()  # set of connected nodes
    groups = {}  # dictionary of isolated node (ids) groups of
    new_edges = []  # minimum spanning tree edges
    # first step of algorithm
    for edge in sorted_edges:
        first_node, second_node = edge['source'], edge['target'] # ids of nodes
        if first_node not in unity or second_node not in unity:
            if first_node not in unity and second_node not in unity:
                # connect source and target nodes of edge
                groups[first_node] = [first_node, second_node]
                groups[second_node] = groups[first_node]
            else:
                if not groups.get(first_node):
                    groups[second_node].append(first_node)
                    groups[first_node] = groups[second_node]
                else:
                    groups[first_node].append(second_node)
                    groups[second_node] = groups[first_node]
            new_edges.append(edge)
            unity.add(first_node)
            unity.add(second_node)
    # second step of algorithm
    for edge in sorted_edges:
        first_node, second_node = edge['source'], edge['target']  # ids of nodes
        if second_node not in groups[first_node]:
            new_edges.append(edge)
            temp = groups[first_node].copy()
            groups[first_node] += groups[second_node]
            groups[second_node] += temp
    return new_edges

# realisation of dijkstras algorithm for graph given in matrix form
def dijkstras_algorithm(nodes, adjacency_matrix):
    size = len(nodes)
    delay_table = [inf] * size  # minimum delay of path from root_id to vertex
    optimum_connections = [inf] * size  # for each node there is best neighbour to visit
    optimum_connections[0] = 0  # inf indicates that we can't visit this node in new spanning tree
    vertex = 0  # starting from node with 0 index
    visited_nodes = {vertex}
    delay_table[0] = 0  # zero delay for starting node
    while vertex != -1:
        for j, delay in enumerate(adjacency_matrix[vertex]):
            if j not in visited_nodes:
                new_delay = delay_table[vertex] + delay
                if new_delay < delay_table[j]:
                    delay_table[j] = new_delay
                    optimum_connections[j] = vertex
        # find new vertex to visit
        temp = -1
        maximum_delay = max(delay_table)
        for i, delay in enumerate(delay_table):
            if delay < maximum_delay and i not in visited_nodes:
                maximum_delay = delay
                temp = i
        vertex = temp
        # memorize
        if vertex >= 0:
            visited_nodes.add(vertex)
    return delay_table, optimum_connections

# the recursive part of depth_first_search
def dfs_recursive(visited, nodes, edges, curr_node, end_node, path, container, delay):
    # nodes are represented by ids
    if curr_node in visited:
        return
    if curr_node not in visited:
        path.append(curr_node)
        visited.add(curr_node)
        if curr_node == end_node:
            container.append({
                'Type': 'reserve',
                'Path': path,
                'Delay': delay
            })
            return
        adjacent_edges = get_adjacent_edges_of_node(curr_node, edges)
        for edge in adjacent_edges:
            neighbour = edge['target']
            dfs_recursive(visited, nodes, edges, neighbour, end_node, path, container, delay + edge['Delay'])

# returns all the paths in given topology
def depth_first_search(root_id, nodes, edges):
    paths = []  # return value
    start = get_node_index_by_id(root_id, nodes)
    # run for each node and return every path from root_id to nodes[i]
    # nodes[i] is end node
    for i in range(start + 1, len(nodes) + start):
        index_i = i % len(nodes)
        visited = set()
        current_path = []
        delay_on_path = 0
        dfs_recursive(visited, nodes, edges, root_id, nodes[index_i]['id'], current_path, paths, delay_on_path)
    # filter paths from duplicates
    duplicates = []
    for i in range(len(paths)):
        for j in range(i+1, len(paths)):
            if paths[i]['Path'] == paths[j]['Path'] and paths[i]['Path'] not in list(map(lambda p: p['Path'], duplicates)):
                duplicates.append(paths[i])
    paths = list(filter(lambda p: p['Path'] not in list(map(lambda p_d: p_d['Path'], duplicates)), paths))
    paths += duplicates
    return paths


# --- FILE GENERATORS ---

# creates file with data of given topology
def generate_topology_description(filename, nodes, edges):
    try:
        csv_file = open(f'{filename}_topo.csv', 'w+')
    except OSError as e:
        print(f"Unable to open {filename}", file=sys.stderr)
        exit(1)
    labels = [
        'Node 1 (id)',
        'Node 1 (label)',
        'Node 1 (longitude)',
        'Node 1 (latitude)',
        'Node 2 (id)',
        'Node 2 (label)',
        'Node 2 (longitude)',
        'Node 2 (latitude)',
        'Distance (km)',
        'Delay (mks)',
    ]
    print(','.join(labels), file=csv_file)
    for node in nodes:
        # get every edge of current node
        adjacent_edges = get_adjacent_edges_of_node(node['id'], edges)
        # walk through edges
        for edge in adjacent_edges:
            source_node = get_node_by_id(edge['source'], nodes)
            target_node = get_node_by_id(edge['target'], nodes)
            first_node_info = [
                source_node['id'],
                source_node['label'],
                source_node['Longitude'],
                source_node['Latitude'],
            ]
            second_node_info = [
                target_node['id'],
                target_node['label'],
                target_node['Longitude'],
                target_node['Latitude'],
            ]
            edge_info = [
                f'{edge["Distance"]:.2f}',
                f'{edge["Delay"]:.2f}',
            ]
            print(','.join(first_node_info + second_node_info + edge_info), file=csv_file)
    csv_file.close()

# creates file with description of given paths
def generate_spanning_tree_description(filename, controller_id, main_paths, reserve_paths):
    # main_paths include information about every other node included in tree
    try:
        csv_file = open(f'{filename}_topo.csv', 'w+')
    except OSError as e:
        print(f"Unable to open {filename}", file=sys.stderr)
        exit(1)
    labels = [
        'Node 1 (id)',
        'Node 2 (id)',
        'Path type',
        'Path',
        'Delay (mks)',
    ]
    print(','.join(labels), file=csv_file)
    for path in main_paths:
        # Node 1 (id) is always equal to controller_id
        # find main paths
        print(','.join([
            controller_id,
            path['Path'][-1],
            path['Type'],
            json.dumps(path['Path']),
            f"{path['Delay']:.2f}"
        ]), file=csv_file)
        # find reserve paths with same endings
        for r_path in reserve_paths:
            if r_path['Path'][-1] == path['Path'][-1]:
                print(','.join([
                    controller_id,
                    r_path['Path'][-1],
                    r_path['Type'],
                    json.dumps(r_path['Path']),
                    f"{r_path['Delay']:.2f}"
                ]), file=csv_file)
    csv_file.close()



# --- PROGRAM IMPLEMENTATION ---

# implements error handling for user input, returns given parameters
def parse_program_arguments():
    parser = argparse.ArgumentParser(
        description='Analyse topology from TopologyZOO.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-k', help='enter criteria for file')
    parser.add_argument('-t', help='enter filename of topology from TopologyZOO')
    argv = vars(parser.parse_args())
    # error handling
    if argv['k'] is None:
        print('Criteria is not given.', file=sys.stderr)
        exit(1)
    if argv['t'] is None:
        print("Filename of topology not given.", file=sys.stderr)
        exit(1)
    K = int(argv['k'])
    if K != 1 and K != 2:
        print("Invalid value of parameter 'k' given.", file=sys.stderr)
        exit(1)
    filename = argv['t']
    f_l = filename.lower()
    if ('.graphml' not in f_l) or ('.graphml' in f_l and '.graphml' != f_l[len(f_l) - 8 : len(f_l)]):
        print("Invalid file extension.", file=sys.stderr)  # take only .graphml files
        exit(1)
    return K, filename

# runs first criteria (type) of program
def run_first_criteria(NODES, EDGES):
    min_max_delay = 0
    maximum_nodes = -inf
    controller_id = NODES[0]['id']
    spanning_tree = {
        'table': [],
        'connections': [],
        'ids': []
    }
    for node in NODES:
        # search place for controller
        root_id = node['id']
        adjacency_matrix, vertex_ids = get_adjacency_matrix(root_id, NODES, EDGES)
        delay_table, optimum_connections = dijkstras_algorithm(NODES, adjacency_matrix)
        filtered_delays = list(filter(lambda i: i != inf and i != 0, delay_table))
        if len(filtered_delays) == 0 or len(filtered_delays) <= maximum_nodes:
            # second condition is for searching the best spanning tree in means of connection count
            continue
        maximum_nodes = len(filtered_delays)  # then len(filtered_delays) > maximum_nodes
        if min_max_delay == 0 or max(filtered_delays) <= min_max_delay:
            # memorize this spanning tree
            min_max_delay = max(delay_table)
            controller_id = root_id
            spanning_tree['table'] = delay_table
            spanning_tree['connections'] = optimum_connections
            spanning_tree['ids'] = vertex_ids
    # delay_table[i] indicates what delay accumulates between nodes with ids of root_id and vertex_ids[i]
    # optimum_connections[i] indicates if vertex_ids[i] node connected to spanning tree, generated by dijkstras algo
    # therefore, this condition is equivalent to optimum_connections[i] != inf
    main_paths = get_main_paths(
        spanning_tree['table'],
        spanning_tree['connections'],
        spanning_tree['ids']
    )  # get main paths
    return controller_id, main_paths


if __name__ == '__main__':
    K, filename = parse_program_arguments()  # get program variables
    NODES, EDGES = parse_topology(f'{SOURCE_PATH}/{filename}')  # parses given topology
    calculate_distance_and_delay(NODES, EDGES)  # makes graph weighted
    new_filename = filename[:-8]  # remove .graphml
    csv_path = f'{OUTPUT_PATH}/{new_filename}'  # where output files are generated
    if K == 1:
        # first criteria (K1)
        # generate first .csv file
        generate_topology_description(csv_path, NODES, EDGES)
        # run K1
        controller_id, main_paths = run_first_criteria(NODES, EDGES)
        reserve_paths = depth_first_search(controller_id, NODES, EDGES)
        # generate second .csv file
        generate_spanning_tree_description(csv_path, controller_id, main_paths, reserve_paths)
    elif K == 2:
        # second criteria + first criteria (K2 + K1)
        # get minimum spanning tree
        MIN_SPANNING_TREE_EDGES = kruskals_algorithm(EDGES)
        # generate first .csv file
        generate_topology_description(csv_path, NODES, MIN_SPANNING_TREE_EDGES)
        # run K1 for spanning tree
        controller_id, main_paths = run_first_criteria(NODES, MIN_SPANNING_TREE_EDGES)
        reserve_paths = depth_first_search(controller_id, NODES, MIN_SPANNING_TREE_EDGES)
        # generate second .csv file
        generate_spanning_tree_description(csv_path, controller_id, main_paths, reserve_paths)