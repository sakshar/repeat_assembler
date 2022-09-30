from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv


def get_assembly_maps(naive, cc, kmeans):
    naive_map = dict()
    cc_map = dict()
    kmeans_map = dict()
    for record in naive:
        naive_map[str(record.id)] = record.seq
    for record in cc:
        cc_map[str(record.id)] = record.seq
    for record in kmeans:
        kmeans_map[str(record.id)] = record.seq
    return naive_map, cc_map, kmeans_map


def get_final_assembly(paths, edges, contig_map, path_to_output):
    for path in paths:
        output_file = "".join(path)
        assembly_seq = ""
        flag = 0  # only for resolving the first overlap
        for i in range(1, len(path)):
            u, v = path[i-1], path[i]
            if flag == 0:
                assembly_seq = str(contig_map[u][:-(edges[(u, v)][0][1] - edges[(u, v)][0][0])]) + str(contig_map[v])
                flag = 1
            else:
                assembly_seq = str(assembly_seq[:-(edges[(u, v)][0][1] - edges[(u, v)][0][0])]) + str(contig_map[v])
        print(output_file, ":", len(assembly_seq))
        seqs = []
        record = SeqRecord(Seq(assembly_seq), id=output_file, description="path_"+output_file)
        seqs.append(record)
        SeqIO.write(seqs, path_to_output + output_file + ".fasta", "fasta")


def paf_reader(infile, slack):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        removed_nodes = set()
        for row in reader:
            if row[0] != row[5] and row[0] not in removed_nodes and row[5] not in removed_nodes:
                new = [str(row[0])] + [int(i) for i in row[1:4]] + [str(row[5])] + [int(i) for i in row[6:9]]
                if new[5] > new[1] >= new[3] - new[2] >= new[1] - slack:
                    removed_nodes.add(new[0])
                elif new[1] > new[5] >= new[7] - new[6] >= new[5] - slack:
                    removed_nodes.add(new[4])
                elif row[4] == "-":
                    continue
                else:
                    data.append(new)
    filterer_data = []
    for row in data:
        if row[0] not in removed_nodes and row[4] not in removed_nodes:
            filterer_data.append(row)
    return filterer_data


def get_overlap_graph(overlap_data, repeat_length, tolerance):
    edges = dict()
    nodes = dict()
    for row in overlap_data:
        length1, start1, end1, length2, start2, end2 = row[1], row[2], row[3], row[5], row[6], row[7]
        delta1, delta2 = end1 - start1, end2 - start2
        if delta1 < repeat_length and delta2 < repeat_length:
            if 0 <= start1 <= tolerance and length2 - tolerance <= end2 <= length2:
                edges[(row[4], row[0])] = ((start2, end2), (start1, end1))
                nodes[row[0]] = length1
                nodes[row[4]] = length2
            elif 0 <= start2 <= tolerance and length1 - tolerance <= end1 <= length1:
                edges[(row[0], row[4])] = ((start1, end1), (start2, end2))
                nodes[row[0]] = length1
                nodes[row[4]] = length2
    return edges, nodes


def get_adjacency_list(edges):
    adj_list = dict()
    for key in edges.keys():
        v1, v2 = key[0], key[1]
        if v1 not in adj_list:
            adj_list[v1] = [v2]
        else:
            adj_list[v1].append(v2)
    return adj_list


def get_adjacency_list_with_all_nodes(edges, nodes):
    adj_list = dict()
    for node in nodes.keys():
        adj_list[node] = []
    for key in edges.keys():
        v1, v2 = key[0], key[1]
        adj_list[v1].append(v2)
    return adj_list


# functions to check for cycles
def isCyclicUtil(graph, v, visited, recStack):

    # Mark current node as visited and
    # adds to recursion stack
    visited[v] = True
    recStack[v] = True

    # Recur for all neighbours
    # if any neighbour is visited and in
    # recStack then graph is cyclic
    for neighbour in graph[v]:
        if not visited[neighbour]:
            if isCyclicUtil(graph, neighbour, visited, recStack):
                return True
        elif recStack[neighbour]:
            return True

    # The node needs to be popped from
    # recursion stack before function ends
    recStack[v] = False
    return False


# Returns true if graph is cyclic else false
def isCyclic(graph, nodes):
    visited, recStack = dict(), dict()
    for node in nodes.keys():
        visited[node] = False
        recStack[node] = False
    for node in nodes.keys():
        if not visited[node]:
            if isCyclicUtil(graph, node, visited, recStack):
                return True
    return False


# code for getting all possible paths in a DAG
def dfs(data, path, paths):
    datum = path[-1]
    if datum in data:
        for val in data[datum]:
            new_path = path + [val]
            paths = dfs(data, new_path, paths)
    else:
        paths += [path]
    return paths


def enumerate_paths(graph):
    nodes = list(graph.keys())
    all_paths = []
    for node in nodes:
        node_paths = dfs(graph, [node], [])
        all_paths += node_paths
    return all_paths


paf = "simulated/k21/output/sim5_depth50_max20000_hifi/naive2/overlaps.paf"
path_to_output = "simulated/k21/output/sim5_depth50_max20000_hifi/naive2/final_assembly/"
naive = SeqIO.parse('simulated/k21/output/sim5_depth50_max20000_hifi/naive2/naive2.asm.fasta', 'fasta')
cc = SeqIO.parse('simulated/k21/output/sim5_depth50_max20000_hifi/cc/cc.asm.fasta', 'fasta')
kmeans = SeqIO.parse('simulated/k21/output/sim5_depth50_max20000_hifi/kmeans/kmeans.asm.fasta', 'fasta')
slack = 100
overlap_data = paf_reader(paf, slack)
#naive_map, cc_map, kmeans_map = get_assembly_maps(naive, cc, kmeans)
#get_final_assembly(overlap_data, naive_map, path_to_output)
#print(len(overlap_data))
#for row in overlap_data:
#    print(row)
repeat_length = 20000
tolerance = 1000
edges, nodes = get_overlap_graph(overlap_data, repeat_length, tolerance)
print("----- contigs with length -----")
print(nodes)
print("----- overlaps with suffix and prefix windows -----")
for edge in edges.keys():
    print(edge, ":", edges[edge])
graph_all_nodes = get_adjacency_list_with_all_nodes(edges, nodes)
graph = get_adjacency_list(edges)
print("----- overlap graph as adjacency list -----")
print(graph)
print("----- check for DAG -----")
if isCyclic(graph_all_nodes, nodes):
    print("The overlap graph contains cycle :'(")
else:
    print("The overlap graph is a DAG :)")
    paths = enumerate_paths(graph)
    print("----- all possible paths in the overlap graph -----")
    for path in paths:
        print(path)
    naive_map, cc_map, kmeans_map = get_assembly_maps(naive, cc, kmeans)
    print("----- all possible assemblies with length -----")
    get_final_assembly(paths, edges, naive_map, path_to_output)
"""
    naive_map, cc_map, kmeans_map = get_assembly_maps(naive, cc, kmeans)
    print("----- all possible assemblies with length -----")
    get_final_assembly(paths, edges, naive_map, path_to_output)
print("------- naive assembly -------")
print(naive_map)
print("------- cc assembly -------")
print(cc_map)
print("------- kmeans assembly -------")
print(kmeans_map)
"""