from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sys

cycle_start, cycle_end = -1, -1


def dfs(v, color, parent, adj):
    global cycle_start, cycle_end
    color[v] = 1
    for u in adj[v]:
        if color[u] == 0:
            parent[u] = v
            if dfs(u, color, parent, adj):
                return True
        elif color[u] == 1:
            cycle_end = v
            cycle_start = u
            return True
    color[v] = 2
    return False


def find_cycle(adj, nodes):
    node_list = list(nodes.keys())
    color = dict()
    parent = dict()
    for node in node_list:
        color[node] = 0
        parent[node] = -1
    global cycle_start, cycle_end
    cycle_start = -1

    for v in node_list:
        if color[v] == 0 and dfs(v, color, parent, adj):
            break
    cycle = list()
    if cycle_start == -1:
        # print("Acyclic")
        return cycle
    else:
        cycle.append(cycle_start)
        v = cycle_end
        while v != cycle_start:
            cycle = [v] + cycle
            v = parent[v]
        cycle = [cycle_start] + cycle

        # print("Cycle found: ", cycle)
        return cycle


def get_assembly_map(assembly):
    assembly_map = dict()
    for record in assembly:
        seq = record.seq
        assembly_map[str(record.id)] = (seq, Seq(seq).reverse_complement())
    return assembly_map


def get_final_assebmly_without_overlap_graph(contig_map, path_to_output, ref_size):
    current_assembly = dict()
    current_dist = ref_size
    sorted_contigs = dict(sorted(contig_map.items(), key=lambda x: len(x[1][0]), reverse=True))
    # for id in sorted_contigs:
    #    print(id, len(sorted_contigs[id][0]))
    for id in sorted_contigs:
        current_contig = sorted_contigs[id][0]
        if abs(current_dist - len(current_contig)) < abs(current_dist):
            current_assembly[id] = current_contig
            current_dist -= len(current_contig)
        # else:
        #    break
    seqs = []
    # print("Best assembly:", list(current_assembly.keys()), current_dist, ref_size)
    # for id in current_assembly:
    #    print(id, len(current_assembly[id]))
    for id in current_assembly:
        record = SeqRecord(Seq(current_assembly[id]), id=id, description=str(len(current_assembly[id])))
        seqs.append(record)
    SeqIO.write(seqs, path_to_output + "rambler.fasta", "fasta")


def get_final_assembly_for_DAGs(contig_map, path_to_output, ref_size, edges, nodes):
    print(list(nodes.keys()))
    paths = enumerate_paths(get_adjacency_list(edges))
    candidate_assemblies = list()
    for path in paths:
        current_nodes = list(nodes.keys())
        current_contigs = dict()
        current_contig = ""
        current_contig_id = path[0]
        used_nodes = [path[0]]
        current_dist = ref_size
        current_orientations = ""
        last_edge = ""
        flag = False
        for i in range(1, len(path)):
            u, v = path[i - 1], path[i]
            current_orientations += edges[(u, v)][2]
            # when merging the first edge on a path
            if not flag:
                if edges[(u, v)][2] == "+":
                    current_contig = str(contig_map[u][0]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                    last_edge = "+"
                elif edges[(u, v)][2] == "-":
                    current_contig = str(contig_map[u][1]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                    last_edge = "-"
                elif edges[(u, v)][2] == "*":
                    current_contig = str(contig_map[u][0]) + str(contig_map[v][1][edges[(u, v)][1][1]:])
                    last_edge = "*"
                current_contig_id += v
                flag = True
            else:
                if edges[(u, v)][2] == "+":
                    # handling "++", "-+": continue with the current contig
                    if last_edge in ["+", "-"]:
                        current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    # handling "*+": break the current one and start a new contig
                    elif last_edge == "*":
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    last_edge = "+"
                elif edges[(u, v)][2] == "-":
                    # handling "+-", "--": break the current one and start a new contig
                    if last_edge in ["+", "-"]:
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    # handling "*-": continue with the current contig
                    elif last_edge == "*":
                        current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    last_edge = "-"
                elif edges[(u, v)][2] == "*":
                    # handling "+*", "-*": continue with the current contig
                    if last_edge in ["+", "-"]:
                        current_contig += str(contig_map[v][1][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    # handling "**": break the current one and start a new contig
                    elif last_edge == "*":
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][1][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    last_edge = "*"
            used_nodes.append(v)
        if current_contig_id not in current_contigs.keys():
            current_contigs[current_contig_id] = current_contig
            current_dist -= len(current_contig)
        remaining_nodes = [x for x in current_nodes if x not in used_nodes]
        for node in remaining_nodes:
            current_contigs[node] = contig_map[node][0]
            current_dist -= nodes[node]
        candidate_assemblies.append((current_contigs, current_dist, current_orientations))
    sorted_candidate_assemblies = sorted(candidate_assemblies, key=lambda x: abs(x[1]))
    seqs = []
    best_assembly, best_dist, best_orientations = sorted_candidate_assemblies[0][0], sorted_candidate_assemblies[0][1], \
                                                  sorted_candidate_assemblies[0][2]
    print("Best assembly:", list(best_assembly.keys()), best_orientations, best_dist, ref_size)
    print("All assemblies:")
    for assembly in sorted_candidate_assemblies:
        print(list(assembly[0].keys()), assembly[2], assembly[1])
    for id in best_assembly.keys():
        record = SeqRecord(Seq(best_assembly[id]), id=id, description=str(len(best_assembly[id])))
        seqs.append(record)
    SeqIO.write(seqs, path_to_output + "rambler.fasta", "fasta")


def get_assembly_for_all_rotations_of_a_cycle(contig_map, ref_size, cycle, edges):
    all_rotations = [cycle[:-1]]
    for i in range(1, len(cycle) - 1):
        all_rotations.append(cycle[:-1][i:] + cycle[:-1][:i])
    cyclic_assemblies = []
    cyclic_contigs_list = []
    for current_rotation in all_rotations:
        current_contigs = dict()
        current_orientations = ""
        current_contig = ""
        current_contig_id = current_rotation[0]
        current_dist = ref_size
        last_edge = ""
        flag = False
        for i in range(1, len(current_rotation)):
            u, v = current_rotation[i - 1], current_rotation[i]
            current_orientations += edges[(u, v)][2]
            # when merging the first edge on a path
            if not flag:
                if edges[(u, v)][2] == "+":
                    current_contig = str(contig_map[u][0]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                    last_edge = "+"
                elif edges[(u, v)][2] == "-":
                    current_contig = str(contig_map[u][1]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                    last_edge = "-"
                elif edges[(u, v)][2] == "*":
                    current_contig = str(contig_map[u][0]) + str(contig_map[v][1][edges[(u, v)][1][1]:])
                    last_edge = "*"
                current_contig_id += v
                flag = True
            else:
                if edges[(u, v)][2] == "+":
                    # handling "++", "-+": continue with the current contig
                    if last_edge in ["+", "-"]:
                        current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    # handling "*+": break the current one and start a new contig
                    elif last_edge == "*":
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    last_edge = "+"
                elif edges[(u, v)][2] == "-":
                    # handling "+-", "--": break the current one and start a new contig
                    if last_edge in ["+", "-"]:
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    # handling "*-": continue with the current contig
                    elif last_edge == "*":
                        current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    last_edge = "-"
                elif edges[(u, v)][2] == "*":
                    # handling "+*", "-*": continue with the current contig
                    if last_edge in ["+", "-"]:
                        current_contig += str(contig_map[v][1][edges[(u, v)][1][1]:])
                        current_contig_id += v
                    # handling "**": break the current one and start a new contig
                    elif last_edge == "*":
                        current_contigs[current_contig_id] = current_contig
                        current_dist -= len(current_contig)
                        current_contig = str(contig_map[v][1][edges[(u, v)][1][1]:])
                        current_contig_id = v
                    last_edge = "*"
        if current_contig_id not in current_contigs.keys():
            current_contigs[current_contig_id] = current_contig
            current_dist -= len(current_contig)
        cyclic_assemblies.append((current_contigs, current_dist, current_orientations))
        cyclic_contigs_list += [list(current_contigs.keys())]
    return cyclic_assemblies, all_rotations, cyclic_contigs_list


def get_final_assembly_for_cycles(contig_map, path_to_output, ref_size, edges, nodes):
    cycle = find_cycle(get_adjacency_list_with_all_nodes(edges, nodes), nodes)
    # j = 0
    # print(j, ":", cycle)
    true_edges = edges.copy()
    while True:
        cyclic_assemblies, all_rotations, cyclic_contigs_list = get_assembly_for_all_rotations_of_a_cycle(contig_map,
                                                                                                          ref_size,
                                                                                                          cycle, edges)
        cycle_nodes = all_rotations[0]
        non_cycle_nodes = [x for x in list(nodes.keys()) if x not in cycle_nodes]
        cycle_edges = dict()
        non_cycle_edges = edges.copy()
        for i in range(1, len(cycle)):
            cycle_edges[(cycle[i - 1], cycle[i])] = non_cycle_edges.pop((cycle[i - 1], cycle[i]))
        # print(non_cycle_nodes)
        # print(list(non_cycle_edges.keys()))
        new_cycle = find_cycle(get_adjacency_list_with_all_nodes(non_cycle_edges, nodes), nodes)
        # j += 1
        # print(j, ":", new_cycle)
        if len(new_cycle) == 0:
            # print("DAG")
            break
        elif len(new_cycle) > len(cycle):
            cycle = new_cycle.copy()
            # print(j, ":", new_cycle)
        else:
            # non_cycle_edges.pop((new_cycle[-2], new_cycle[-1]))
            non_cycle_edges.pop((new_cycle[1], new_cycle[2]))
            temp_cycle = find_cycle(get_adjacency_list_with_all_nodes(non_cycle_edges, nodes), nodes)
            if len(temp_cycle) == 0:
                # print("no more cycles")
                break
            else:
                # print("still cycles, think!!!", temp_cycle)
                # edges.pop((new_cycle[1], new_cycle[2]))
                for k in range(1, len(temp_cycle)):
                    edges.pop((temp_cycle[k - 1], temp_cycle[k]))

    # print("cycle:", cycle)
    cyclic_assemblies, all_rotations, cyclic_contigs_list = get_assembly_for_all_rotations_of_a_cycle(contig_map,
                                                                                                      ref_size, cycle,
                                                                                                      edges)
    cycle_nodes = all_rotations[0]
    non_cycle_nodes = [x for x in list(nodes.keys()) if x not in cycle_nodes]
    cycle_edges = dict()
    # non_cycle_edges = edges.copy()
    for i in range(1, len(cycle)):
        cycle_edges[(cycle[i - 1], cycle[i])] = edges[(cycle[i - 1], cycle[i])]
    print("cycle_nodes:", cycle_nodes)
    print("non-cycle_nodes:", non_cycle_nodes)
    paths = enumerate_paths(get_adjacency_list(non_cycle_edges))
    # print("all paths:")
    # for path in paths:
    #    print(path)
    adjusted_paths = []
    # fragmented_paths = []
    for path in paths:
        if path[0] in cycle_nodes:
            new_path = [path[0]]
            for i in range(1, len(path)):
                if path[i] in non_cycle_nodes:
                    new_path += [path[i]]
                else:
                    # if len(path[i:]) > 1 and path[i:] not in fragmented_paths:
                    #    fragmented_paths.append(path[i:])
                    break
            if len(new_path) > 1 and new_path not in adjusted_paths:
                adjusted_paths.append(new_path)
        elif path[0] in non_cycle_nodes:
            new_path = [path[0]]
            for i in range(1, len(path)):
                if path[i] in non_cycle_nodes:
                    new_path += [path[i]]
                else:
                    new_path += [path[i]]
                    # if len(path[i:]) > 1 and path[i:] not in fragmented_paths:
                    #    fragmented_paths.append(path[i:])
                    break
            if len(new_path) > 1 and new_path not in adjusted_paths:
                adjusted_paths.append(new_path)
    print("adjusted paths:")
    for path in adjusted_paths:
        print(path)
    # print("fragmented paths:")
    # for path in fragmented_paths:
    #    print(path)
    candidate_assemblies = []
    if len(non_cycle_nodes) == 0:
        sorted_candidate_assemblies = sorted(cyclic_assemblies, key=lambda x: abs(x[1]))
        seqs = []
        best_assembly, best_dist, best_orientations = sorted_candidate_assemblies[0][0], sorted_candidate_assemblies[0][
            1], sorted_candidate_assemblies[0][2]
        print("Best assembly:", list(best_assembly.keys()), best_orientations, best_dist, ref_size)
        print("All assemblies:")
        for assembly in sorted_candidate_assemblies:
            print(list(assembly[0].keys()), assembly[2], assembly[1])
        for id in best_assembly.keys():
            record = SeqRecord(Seq(best_assembly[id]), id=id, description=str(len(best_assembly[id])))
            seqs.append(record)
        SeqIO.write(seqs, path_to_output + "rambler.fasta", "fasta")
        return
    if len(adjusted_paths) == 0:
        for assembly in cyclic_assemblies:
            current_contigs, current_dist, current_orientations = assembly[0], assembly[1], "||" + assembly[2] + "||"
            for node in non_cycle_nodes:
                current_contigs[node] = contig_map[node][0]
                current_dist -= nodes[node]
            candidate_assemblies.append((current_contigs, current_dist, current_orientations))
        sorted_candidate_assemblies = sorted(candidate_assemblies, key=lambda x: abs(x[1]))
        seqs = []
        best_assembly, best_dist, best_orientations = sorted_candidate_assemblies[0][0], sorted_candidate_assemblies[0][
            1], sorted_candidate_assemblies[0][2]
        print("Best assembly:", list(best_assembly.keys()), best_orientations, best_dist, ref_size)
        print("All assemblies:")
        for assembly in sorted_candidate_assemblies:
            print(list(assembly[0].keys()), assembly[2], assembly[1])
        for id in best_assembly.keys():
            record = SeqRecord(Seq(best_assembly[id]), id=id, description=str(len(best_assembly[id])))
            seqs.append(record)
        SeqIO.write(seqs, path_to_output + "rambler.fasta", "fasta")
        return
    l = 0
    for assembly in cyclic_assemblies:
        true_current_contigs, true_current_dist, true_current_orientations = assembly[0], assembly[1], assembly[2]
        current_rotation = all_rotations[l]
        true_current_contigs_ids = cyclic_contigs_list[l]
        # print(current_orientations)
        # print(current_rotation)
        # print(current_contigs_ids)
        for path in adjusted_paths:
            current_contigs = true_current_contigs.copy()
            current_contigs_ids = true_current_contigs_ids.copy()
            current_dist = true_current_dist
            current_orientations = true_current_orientations[:]
            # extending the cycle from the end
            if path[0] == current_rotation[-1]:
                used_nodes = []
                current_contig_id = current_contigs_ids[-1]
                current_contig = current_contigs.pop(current_contig_id)
                current_dist += len(current_contig)
                current_orientation = current_orientations[:] + "||"
                last_edge = edges[(current_rotation[-2], current_rotation[-1])][2]
                for p in range(1, len(path)):
                    u, v = path[p - 1], path[p]
                    current_orientation += edges[(u, v)][2]
                    if edges[(u, v)][2] == "+":
                        # handling "++", "-+": continue with the current contig
                        if last_edge in ["+", "-"]:
                            current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                            current_contig_id += v
                        # handling "*+": break the current one and start a new contig
                        elif last_edge == "*":
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                            current_contig_id = v
                        last_edge = "+"
                    elif edges[(u, v)][2] == "-":
                        # handling "+-", "--": break the current one and start a new contig
                        if last_edge in ["+", "-"]:
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                            current_contig_id = v
                        # handling "*-": continue with the current contig
                        elif last_edge == "*":
                            current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                            current_contig_id += v
                        last_edge = "-"
                    elif edges[(u, v)][2] == "*":
                        # handling "+*", "-*": continue with the current contig
                        if last_edge in ["+", "-"]:
                            current_contig += str(contig_map[v][1][edges[(u, v)][1][1]:])
                            current_contig_id += v
                        # handling "**": break the current one and start a new contig
                        elif last_edge == "*":
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[v][1][edges[(u, v)][1][1]:])
                            current_contig_id = v
                        last_edge = "*"
                    used_nodes += [v]
                if current_contig_id not in current_contigs.keys():
                    current_contigs[current_contig_id] = current_contig
                    current_dist -= len(current_contig)
                # extend from the front here, think about it later
                """
                current_remaining_nodes = [x for x in non_cycle_nodes if x not in used_nodes]
                for opposite_path in adjusted_paths:
                    if opposite_path[-1] == current_rotation[0]:
                        skip = False
                        for n in opposite_path[:-1]:
                            if n not in current_remaining_nodes:
                                skip = True
                        if not skip:
                            current_contig_id = current_contigs_ids[0]
                            current_contig = current_contigs.pop(current_contig_id)
                            current_dist += len(current_contig)
                            current_orientation = "||" + current_orientations[:]
                            last_edge = edges[(current_rotation[0], current_rotation[1])][2]
                """
                remaining_nodes = [x for x in non_cycle_nodes if x not in used_nodes]
                for node in remaining_nodes:
                    current_contigs[node] = contig_map[node][0]
                    current_dist -= nodes[node]
                candidate_assemblies.append((current_contigs, current_dist, current_orientation))
            # extending the cycle from the front
            elif path[-1] == current_rotation[0]:
                used_nodes = []
                current_contig_id = current_contigs_ids[0]
                current_contig = current_contigs.pop(current_contig_id)
                current_dist += len(current_contig)
                current_orientation = "||" + current_orientations[:]
                last_edge = edges[(current_rotation[0], current_rotation[1])][2]
                for p in range(len(path) - 2, -1, -1):
                    u, v = path[p], path[p + 1]
                    current_orientation = edges[(u, v)][2] + current_orientation
                    if edges[(u, v)][2] == "+":
                        # handling "++", "-+": continue with the current contig
                        if last_edge in ["+", "*"]:
                            current_contig = str(contig_map[u][0][:edges[(u, v)][0][0]]) + current_contig
                            current_contig_id = u + current_contig_id
                        # handling "*+": break the current one and start a new contig
                        elif last_edge == "-":
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[u][0][:edges[(u, v)][0][0]])
                            current_contig_id = u
                        last_edge = "+"
                    elif edges[(u, v)][2] == "*":
                        # handling "+-", "--": break the current one and start a new contig
                        if last_edge in ["+", "*"]:
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[u][0][:edges[(u, v)][0][0]])
                            current_contig_id = u
                        # handling "*-": continue with the current contig
                        elif last_edge == "-":
                            current_contig = str(contig_map[u][0][:edges[(u, v)][0][0]]) + current_contig
                            current_contig_id = u + current_contig_id
                        last_edge = "*"
                    elif edges[(u, v)][2] == "-":
                        # handling "+*", "-*": continue with the current contig
                        if last_edge in ["+", "*"]:
                            current_contig = str(contig_map[u][1][:edges[(u, v)][0][0]]) + current_contig
                            current_contig_id = u + current_contig_id
                        # handling "**": break the current one and start a new contig
                        elif last_edge == "-":
                            current_contigs[current_contig_id] = current_contig
                            current_dist -= len(current_contig)
                            current_contig = str(contig_map[u][1][:edges[(u, v)][0][0]])
                            current_contig_id = u
                        last_edge = "-"
                    used_nodes += [u]
                if current_contig_id not in current_contigs.keys():
                    current_contigs[current_contig_id] = current_contig
                    current_dist -= len(current_contig)
                # extend from the end here,  think about it later
                remaining_nodes = [x for x in non_cycle_nodes if x not in used_nodes]
                for node in remaining_nodes:
                    current_contigs[node] = contig_map[node][0]
                    current_dist -= nodes[node]
                candidate_assemblies.append((current_contigs, current_dist, current_orientation))
            # the path has no cycle nodes
            elif path[0] in non_cycle_nodes and path[-1] in non_cycle_nodes:
                used_nodes = [path[0]]
                current_contig_id = path[0]
                current_contig = ""
                # current_dist += len(current_contig)
                current_orientation = "||" + current_orientations[:] + "||"
                last_edge = ""
                flag = False
                for p in range(1, len(path)):
                    u, v = path[p - 1], path[p]
                    current_orientation += edges[(u, v)][2]
                    # when merging the first edge on a path
                    if not flag:
                        if edges[(u, v)][2] == "+":
                            current_contig = str(contig_map[u][0]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                            last_edge = "+"
                        elif edges[(u, v)][2] == "-":
                            current_contig = str(contig_map[u][1]) + str(contig_map[v][0][edges[(u, v)][1][1]:])
                            last_edge = "-"
                        elif edges[(u, v)][2] == "*":
                            current_contig = str(contig_map[u][0]) + str(contig_map[v][1][edges[(u, v)][1][1]:])
                            last_edge = "*"
                        current_contig_id += v
                        flag = True
                    else:
                        if edges[(u, v)][2] == "+":
                            # handling "++", "-+": continue with the current contig
                            if last_edge in ["+", "-"]:
                                current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                                current_contig_id += v
                            # handling "*+": break the current one and start a new contig
                            elif last_edge == "*":
                                current_contigs[current_contig_id] = current_contig
                                current_dist -= len(current_contig)
                                current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                                current_contig_id = v
                            last_edge = "+"
                        elif edges[(u, v)][2] == "-":
                            # handling "+-", "--": break the current one and start a new contig
                            if last_edge in ["+", "-"]:
                                current_contigs[current_contig_id] = current_contig
                                current_dist -= len(current_contig)
                                current_contig = str(contig_map[v][0][edges[(u, v)][1][1]:])
                                current_contig_id = v
                            # handling "*-": continue with the current contig
                            elif last_edge == "*":
                                current_contig += str(contig_map[v][0][edges[(u, v)][1][1]:])
                                current_contig_id += v
                            last_edge = "-"
                        elif edges[(u, v)][2] == "*":
                            # handling "+*", "-*": continue with the current contig
                            if last_edge in ["+", "-"]:
                                current_contig += str(contig_map[v][1][edges[(u, v)][1][1]:])
                                current_contig_id += v
                            # handling "**": break the current one and start a new contig
                            elif last_edge == "*":
                                current_contigs[current_contig_id] = current_contig
                                current_dist -= len(current_contig)
                                current_contig = str(contig_map[v][1][edges[(u, v)][1][1]:])
                                current_contig_id = v
                            last_edge = "*"
                    used_nodes.append(v)
                if current_contig_id not in current_contigs.keys():
                    current_contigs[current_contig_id] = current_contig
                    current_dist -= len(current_contig)
                remaining_nodes = [x for x in non_cycle_nodes if x not in used_nodes]
                for node in remaining_nodes:
                    current_contigs[node] = contig_map[node][0]
                    current_dist -= nodes[node]
                candidate_assemblies.append((current_contigs, current_dist, current_orientation))
        l += 1
    sorted_candidate_assemblies = sorted(candidate_assemblies, key=lambda x: abs(x[1]))
    seqs = []
    best_assembly, best_dist, best_orientations = sorted_candidate_assemblies[0][0], sorted_candidate_assemblies[0][1], \
                                                  sorted_candidate_assemblies[0][2]
    print("Best assembly:", list(best_assembly.keys()), best_orientations, best_dist, ref_size)
    print("All assemblies:")
    for assembly in sorted_candidate_assemblies:
        print(list(assembly[0].keys()), assembly[2], assembly[1])
    for id in best_assembly.keys():
        record = SeqRecord(Seq(best_assembly[id]), id=id, description=str(len(best_assembly[id])))
        seqs.append(record)
    SeqIO.write(seqs, path_to_output + "rambler.fasta", "fasta")


def paf_reader(infile, slack):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        removed_nodes = set()
        for row in reader:
            if row[0] != row[5] and row[0] not in removed_nodes and row[5] not in removed_nodes:
                new = [str(row[0])] + [int(i) for i in row[1:4]] + [str(row[5])] + [int(i) for i in row[6:9]] + [
                    str(row[4])] + [int(i) for i in row[9:11]]
                if new[5] > new[1] >= new[3] - new[2] >= new[1] - slack:
                    removed_nodes.add(new[0])
                elif new[1] > new[5] >= new[7] - new[6] >= new[5] - slack:
                    removed_nodes.add(new[4])
                else:
                    data.append(new)
    filterer_data = []
    for row in data:
        if row[0] not in removed_nodes and row[4] not in removed_nodes:
            filterer_data.append(row)
    return filterer_data


def get_overlap_graph(overlap_data, tolerance):
    edges = dict()
    nodes = dict()
    for row in overlap_data:
        length1, start1, end1, length2, start2, end2 = row[1], row[2], row[3], row[5], row[6], row[7]
        strand = row[8]
        exact_match, matched_block = row[9], row[10]
        percent_identity = 1.0 * exact_match / matched_block
        if strand == "+":
            if 0 <= start1 <= tolerance and length2 - tolerance <= end2 <= length2:
                if (row[4], row[0]) not in edges:
                    edges[(row[4], row[0])] = ((start2, end2), (start1, end1), strand, exact_match, matched_block)
                    nodes[row[0]] = length1
                    nodes[row[4]] = length2
                else:
                    prev_exact_match, prev_matched_block = edges[(row[4], row[0])][3], edges[(row[4], row[0])][4]
                    prev_percent_identity = 1.0 * prev_exact_match / prev_matched_block
                    if percent_identity > prev_percent_identity or (
                            percent_identity == prev_percent_identity and matched_block < prev_matched_block):
                        edges[(row[4], row[0])] = ((start2, end2), (start1, end1), strand, exact_match, matched_block)
            elif 0 <= start2 <= tolerance and length1 - tolerance <= end1 <= length1:
                if (row[0], row[4]) not in edges:
                    edges[(row[0], row[4])] = ((start1, end1), (start2, end2), strand, exact_match, matched_block)
                    nodes[row[0]] = length1
                    nodes[row[4]] = length2
                else:
                    prev_exact_match, prev_matched_block = edges[(row[0], row[4])][3], edges[(row[0], row[4])][4]
                    prev_percent_identity = 1.0 * prev_exact_match / prev_matched_block
                    if percent_identity > prev_percent_identity or (
                            percent_identity == prev_percent_identity and matched_block < prev_matched_block):
                        edges[(row[0], row[4])] = ((start1, end1), (start2, end2), strand, exact_match, matched_block)
        elif strand == "-":
            if (0 <= start1 <= tolerance and 0 <= start2 <= tolerance) or (
                    length2 - tolerance <= end2 <= length2 and length1 - tolerance <= end1 <= length1):
                if (row[0], row[4]) not in edges:
                    edges[(row[0], row[4])] = ((start1, end1), (start2, end2), strand, exact_match, matched_block)
                    nodes[row[0]] = length1
                    nodes[row[4]] = length2
                else:
                    prev_exact_match, prev_matched_block = edges[(row[0], row[4])][3], edges[(row[0], row[4])][4]
                    prev_percent_identity = 1.0 * prev_exact_match / prev_matched_block
                    if percent_identity > prev_percent_identity or (
                            percent_identity == prev_percent_identity and matched_block < prev_matched_block):
                        edges[(row[0], row[4])] = ((start1, end1), (start2, end2), strand, exact_match, matched_block)
    edges = get_strand_adjusted_edges(edges, nodes)
    self_loop_removed_edges = dict()
    edge_list = list(edges.keys())
    for edge in edges:
        u, v = edge[0], edge[1]
        percent_identity = 1.0 * edges[edge][3] / edges[edge][4]
        if (v, u) in edge_list:
            reverse_percent_identity = 1.0 * edges[(v, u)][3] / edges[(v, u)][4]
            if percent_identity > reverse_percent_identity or (
                    percent_identity == reverse_percent_identity and edges[edge][4] < edges[(v, u)][4]):
                self_loop_removed_edges[edge] = edges[edge]
            else:
                self_loop_removed_edges[(v, u)] = edges[(v, u)]
        else:
            self_loop_removed_edges[edge] = edges[edge]
    # return edges, nodes
    # adjusted_edges = get_strand_adjusted_edges(self_loop_removed_edges, nodes)
    return self_loop_removed_edges, nodes
    # return adjusted_edges, nodes


def get_strand_adjusted_edges(edges, nodes):
    adjusted_edges = dict()
    for edge in edges:
        if edges[edge][2] == "+":
            adjusted_edges[edge] = edges[edge]
        elif edges[edge][2] == "-":
            # code here for adjusting orientation
            start1, end1, start2, end2 = edges[edge][0][0], edges[edge][0][1], edges[edge][1][0], edges[edge][1][1]
            length1, length2 = nodes[edge[0]], nodes[edge[1]]
            new_start1, new_end1 = length1 - end1, length1 - start1
            dist_to_end1, dist_to_end2 = length1 - new_end1, length2 - end2
            if dist_to_end1 < dist_to_end2 and start2 < start1:
                adjusted_edges[edge] = (
                (new_start1, new_end1), (start2, end2), edges[edge][2], edges[edge][3], edges[edge][4])
            else:
                adjusted_edges[(edge[1], edge[0])] = (
                (start2, end2), (new_start1, new_end1), "*", edges[edge][3], edges[edge][4])
    return adjusted_edges


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


# code for getting all possible paths in a DAG
def dfs_for_DAG(data, path, paths):
    datum = path[-1]
    if datum in data:
        for val in data[datum]:
            new_path = path + [val]
            paths = dfs_for_DAG(data, new_path, paths)
    else:
        paths += [path]
    return paths


def enumerate_paths(graph):
    nodes = list(graph.keys())
    all_paths = []
    for node in nodes:
        node_paths = dfs_for_DAG(graph, [node], [])
        all_paths += node_paths
    return all_paths


def cycle_info_writer():
    k = 21
    slack = 100
    tolerance = 1000
    methods = ["naive", "cc"]
    rows = [["repeat_size", "copy", "snp", "depth", "label", "nodes", "edges", "cycle"]]
    for repeat_size in ["5000", "10000", "15000", "20000"]:
        for copy in ["2", "5", "10"]:
            for snp in ["100", "250", "500", "1000", "2000"]:
                for depth in ["10", "20", "30", "40"]:
                    ref_size = 100000 + (int(copy) * int(repeat_size))
                    parent_dir = repeat_size + "_" + copy + "_" + snp
                    paf_file = "../output/" + parent_dir + "/" + depth + "/" + methods[0] + "/overlaps.paf"
                    path_to_output = "../output_reproduced/default/" + parent_dir + "/" + depth + "/"
                    contigs_file = SeqIO.parse(
                        "../output/" + parent_dir + "/" + depth + "/" + methods[0] + "/asm.fasta", 'fasta')
                    assembly_map = get_assembly_map(contigs_file)
                    overlap_data = paf_reader(paf_file, slack)
                    edges, nodes = get_overlap_graph(overlap_data, tolerance)
                    if (len(list(nodes.keys()))) == 0:
                        # print("--------------------")
                        # print(repeat_size, copy, snp, depth)
                        rows.append([repeat_size, copy, snp, depth, "None", "0", "0", "None"])
                        # get_final_assebmly_without_overlap_graph(assembly_map, path_to_output, ref_size)
                        # print("--------------------")
                    else:
                        graph_all_nodes = get_adjacency_list_with_all_nodes(edges, nodes)
                        cycle = find_cycle(graph_all_nodes, nodes)
                        if len(cycle) == 0:
                            # print("--------------------")
                            # print(repeat_size, copy, snp, depth)
                            rows.append([repeat_size, copy, snp, depth, "DAG", nodes, edges, cycle])
                            # get_final_assembly_for_DAGs(assembly_map, path_to_output, ref_size, edges, nodes)
                            # print("--------------------")
                        else:
                            print("--------------------")
                            print(repeat_size, copy, snp, depth)
                            rows.append([repeat_size, copy, snp, depth, "Cycle", nodes, edges, cycle])
                            get_final_assembly_for_cycles(assembly_map, path_to_output, ref_size, edges, nodes)
                            print("--------------------")
    # cycle_info_file = "../output/cycle_info_reproduced_withAdjustedEdges_v2.csv"
    # with open(cycle_info_file, "w") as csvfile:
    #    writer = csv.writer(csvfile, delimiter="\t")
    #    [writer.writerow(r) for r in rows]


cycle_info_writer()