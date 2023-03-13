from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sys

def get_assembly_map(assembly):
    assembly_map = dict()
    for record in assembly:
        assembly_map[str(record.id)] = record.seq
    return assembly_map


def get_final_assembly_without_paths(contig_map, path_to_output):
    max_output_file = ""
    max_assembly_seq = ""
    for id in contig_map.keys():
        if len(contig_map[id]) > len(max_assembly_seq):
            max_output_file = id
            max_assembly_seq = contig_map[id]
    print(max_output_file + "_final", ":", len(max_assembly_seq))
    seqs = []
    record = SeqRecord(Seq(max_assembly_seq), id=max_output_file + "_final",
                       description="size_" + str(len(max_assembly_seq)))
    seqs.append(record)
    SeqIO.write(seqs, path_to_output + max_output_file + "_final.fasta", "fasta")


def get_final_assembly(paths, edges, contig_map, path_to_output):
    max_output_file = ""
    max_assembly_seq = ""
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
        if len(assembly_seq) > len(max_assembly_seq):
            max_output_file = output_file
            max_assembly_seq = assembly_seq
        seqs = []
        record = SeqRecord(Seq(assembly_seq), id=output_file, description="size_"+str(len(assembly_seq)))
        seqs.append(record)
        SeqIO.write(seqs, path_to_output + output_file + ".fasta", "fasta")
    print(max_output_file+"_final", ":", len(max_assembly_seq))
    seqs = []
    record = SeqRecord(Seq(max_assembly_seq), id=max_output_file+"_final", description="size_" + str(len(max_assembly_seq)))
    seqs.append(record)
    SeqIO.write(seqs, path_to_output + max_output_file + "_final.fasta", "fasta")



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
        #if delta1 < repeat_length and delta2 < repeat_length:
        if 0 <= start1 <= tolerance and length2 - tolerance <= end2 <= length2:
            if (row[4], row[0]) not in edges:
                edges[(row[4], row[0])] = ((start2, end2), (start1, end1))
                nodes[row[0]] = length1
                nodes[row[4]] = length2
            else:
                u, v = edges[(row[4], row[0])][0], edges[(row[4], row[0])][1]
                prev_delta1, prev_delta2 = v[1] - v[0], u[1] - u[0]
                if delta1 < prev_delta1 and delta2 < prev_delta2:
                    edges[(row[4], row[0])] = ((start2, end2), (start1, end1))
        elif 0 <= start2 <= tolerance and length1 - tolerance <= end1 <= length1:
            if (row[0], row[4]) not in edges:
                edges[(row[0], row[4])] = ((start1, end1), (start2, end2))
                nodes[row[0]] = length1
                nodes[row[4]] = length2
            else:
                u, v = edges[(row[0], row[4])][0], edges[(row[0], row[4])][1]
                prev_delta1, prev_delta2 = u[1] - u[0], v[1] - v[0]
                if delta1 < prev_delta1 and delta2 < prev_delta2:
                    edges[(row[0], row[4])] = ((start1, end1), (start2, end2))
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


k = 21
slack = 100
tolerance = 1000
methods = ["naive", "cc"]
repeat_size, copy, snp, depth = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
ref_size = 100000 + (copy * repeat_size)
parent_dir = str(repeat_size) + "_" + str(copy) + "_" + str(snp)
paf_file = "../output/"+parent_dir+"/"+str(depth)+"/"+methods[0]+"/overlaps.paf"
path_to_output = "../output/"+parent_dir+"/"+str(depth)+"/"+methods[0]+"/final_assembly/"
contigs_file = SeqIO.parse("../output/"+parent_dir+"/"+str(depth)+"/"+methods[0]+"/asm.fasta", 'fasta')
assembly_map = get_assembly_map(contigs_file)
overlap_data = paf_reader(paf_file, slack)
edges, nodes = get_overlap_graph(overlap_data, repeat_size, tolerance)
if len(list(nodes.keys())) == 0:
    print("----- the best possible assembly with length -----")
    get_final_assembly_without_paths(assembly_map, path_to_output)
else:
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
        print("----- the best possible assembly with length -----")
        get_final_assembly_without_paths(assembly_map, path_to_output)
    else:
        print("The overlap graph is a DAG :)")
        paths = enumerate_paths(graph)
        print("----- all possible paths in the overlap graph -----")
        for path in paths:
            print(path)
        if len(paths) == 0:
            print("----- the best possible assembly with length -----")
            get_final_assembly_without_paths(assembly_map, path_to_output)
        else:
            print("----- all possible assemblies and the best one with length -----")
            get_final_assembly(paths, edges, assembly_map, path_to_output)