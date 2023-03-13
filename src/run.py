from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import csv
import numpy as np


def get_sunk_map(sunk_file):
    sunk_map = dict()
    counter = 0
    while True:
        line = sunk_file.readline()
        if not line:
            break
        sunk_map[line.strip()] = counter
        counter += 1
    return sunk_map


def write_sunk_map(out_file, sunk_map):
    output = open(out_file, "w")
    for sunk in sunk_map.keys():
        output.write(str(sunk) + "," + str(sunk_map[sunk]) + "\n")
    output.close()


def get_read_to_sunk_and_sunk_to_read_maps(sunk_map, read_file):
    read_to_sunk_map = dict()
    sunk_to_read_map = dict()
    sunks = sunk_map.keys()
    for record in read_file:
        read = record.seq
        id = record.id
#        print("HiFi:", str(id))
        for i in range(len(read) - k + 1):
            # for sunk in hifi_sunk_map.keys():
            kmer = read[i:i + k]
            str_kmer = str(kmer)
            str_rev_kmer = str(kmer.reverse_complement())
            if str_kmer in sunks:
                if id not in read_to_sunk_map.keys():
                    read_to_sunk_map[id] = []
                read_to_sunk_map[id] += [(sunk_map[str_kmer], i)]  # (sunk_id, pos)
                if sunk_map[str_kmer] not in sunk_to_read_map.keys():
                    sunk_to_read_map[sunk_map[str_kmer]] = []
                sunk_to_read_map[sunk_map[str_kmer]] += [(id, i)]  # (read_id, pos)
            elif str_rev_kmer in sunks:
                if id not in read_to_sunk_map.keys():
                    read_to_sunk_map[id] = []
                read_to_sunk_map[id] += [(sunk_map[str_rev_kmer], i)]  # (sunk_id, pos)
                if sunk_map[str_rev_kmer] not in sunk_to_read_map.keys():
                    sunk_to_read_map[sunk_map[str_rev_kmer]] = []
                sunk_to_read_map[sunk_map[str_rev_kmer]] += [(id, i)]  # (read_id, pos)
    return read_to_sunk_map, sunk_to_read_map


def write_read_to_sunk_and_sunk_to_read_maps(out_r2s, out_s2r, read_to_sunk_map, sunk_to_read_map):
    out_r2s_file = open(out_r2s, "w")
    out_s2r_file = open(out_s2r, "w")
    for read in read_to_sunk_map.keys():
        out_r2s_file.write(str(read) + ">")
        for tup in read_to_sunk_map[read]:
            out_r2s_file.write(str(tup) + ";")
        out_r2s_file.write("\n")
    out_r2s_file.close()
    for sunk in sunk_to_read_map.keys():
        out_s2r_file.write(str(sunk) + ">")
        for tup in sunk_to_read_map[sunk]:
            out_s2r_file.write(str(tup) + ";")
        out_s2r_file.write("\n")
    out_s2r_file.close()


def get_read_to_sunk_map(infile):
    read_to_sunk_map = dict()
    read_to_sunk_list = dict()
    read_to_sunkpos_map = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        words = line.strip().split(">")
        read = int(words[0].split("_")[1])
        read_to_sunk_map[read] = dict()
        read_to_sunk_list[read] = []
        read_to_sunkpos_map[read] = dict()
        pairs = words[1].split(";")
        pos = 0
        for pair in pairs[:-1]:
            values = pair[1:-1].split(",")
            read_to_sunk_map[read][int(values[0])] = int(values[1])
            read_to_sunk_list[read] += [(int(values[0]), int(values[1]))]
            read_to_sunkpos_map[read][int(values[1])] = pos
            pos += 1
    return read_to_sunk_map, read_to_sunk_list, read_to_sunkpos_map


def get_sunk_to_read_map(infile):
    sunk_to_read_map = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        words = line.strip().split(">")
        sunk_to_read_map[int(words[0])] = dict()
        pairs = words[1].split(";")
        for pair in pairs[:-1]:
            values = pair[1:-1].split(",")
            sunk_to_read_map[int(words[0])][int(values[0].split("_")[1][:-1])] = int(values[1])
    return sunk_to_read_map


def get_pairwise_profiles(sunk_to_read_map, read_to_sunk_map, read_to_sunk_list, read_to_sunkpos_map, tolerance=5):
    pairwise_profiles = dict()
    sunks = sunk_to_read_map.keys()
    count = 0
    for sunk in sunks:
#        print(count, ":", sunk)
        reads = list(sunk_to_read_map[sunk].keys())
        for i in range(len(reads)):
            r1 = reads[i]
            for j in range(i+1, len(reads)):
                r2 = reads[j]
                pos1, pos2 = read_to_sunk_map[r1][sunk], read_to_sunk_map[r2][sunk]
                pairwise_profiles[(r1, r2)] = [(sunk, pos1, pos2)]
                sunk_index1, sunk_index2 = read_to_sunkpos_map[r1][pos1], read_to_sunkpos_map[r2][pos2]
                r1_list, r2_list = read_to_sunk_list[r1], read_to_sunk_list[r2]
                while sunk_index1 + 1 < len(r1_list) and sunk_index2 + 1 < len(r2_list):
                    s1, s2 = r1_list[sunk_index1+1][0], r2_list[sunk_index2+1][0]
                    delta1, delta2 = r1_list[sunk_index1+1][1] - pos1, r2_list[sunk_index2+1][1] - pos2
                    if s1 == s2 and delta1 - tolerance <= delta2 <= delta1 + tolerance:
                        sunk_index1 += 1
                        sunk_index2 += 1
                        pos1, pos2 = r1_list[sunk_index1][1], r2_list[sunk_index2][1]
                        pairwise_profiles[(r1, r2)] += [(s1, pos1, pos2)]
                    else:
                        break
        count += 1
    return pairwise_profiles


def write_pairwise_profiles(pairwise_profiles, file_path):
    output_pairwise_profiles = open(file_path, "w")
    for pairs in pairwise_profiles.keys():
        output_pairwise_profiles.write("(" + str(pairs[0]) + "," + str(pairs[1]) + ")>")
        for profile in pairwise_profiles[pairs]:
            output_pairwise_profiles.write(
                "(" + str(profile[0]) + "," + str(profile[1]) + "," + str(profile[2]) + ");")
        output_pairwise_profiles.write("\n")
    output_pairwise_profiles.close()


def maf_parser(infile):
    maf = AlignIO.parse(infile, "maf")
    maf_profile = [["id", "start", "refStrand", "refSize", "readStrand", "readSize"]]
    read = 0
    for multiple_alignment in maf:
        read += 1
        start, refStrand, refSize, readStrand, readSize = -1, -1, -1, -1, -1
        i = 0
        for seqrec in multiple_alignment:
            if i == 0:
                start, refStrand, refSize = seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["size"]
                i += 1
            else:
                readStrand, readSize = seqrec.annotations["strand"], seqrec.annotations["size"]
        maf_profile.append([read, start, refStrand, refSize, readStrand, readSize])
    return maf_profile


def maf_writer(outfile, maf_profile):
    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in maf_profile]


def read_pairwise_profiles(infile):
    pairwise_profiles = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        key_values = line.strip().split(">")
        read_pair = key_values[0][1:-1].split(",")
        pairwise_profiles[(int(read_pair[0]), int(read_pair[1]))] = []
        triplets = key_values[1].split(";")
        for triplet in triplets[:-1]:
            values = triplet[1:-1].split(",")
            pairwise_profiles[(int(read_pair[0]), int(read_pair[1]))] += [(int(values[0]), int(values[1]), int(values[2]))]
    return pairwise_profiles


# sunk_delta = sunkPos_r1 - sunkPos_r2, true_delta = start_r1 - start_r2
def get_pairwise_location_data(maf_data, pw_profiles, start, end, threshold):
    overlap_locations = [["r1", "start_r1", "r2", "start_r2", "true_delta", "sunk_delta_mean", "sunk_delta_sd", "#_of_SUNKs"]]
    read_set = set()
    for key in pw_profiles.keys():
        r1, r2 = key[0], key[1]
        shared_sunks = pw_profiles[(r1, r2)]
        sunk_delta_mean = 0
        no_shared_sunk = len(shared_sunks)
        sunk_deltas = np.zeros(no_shared_sunk)
        for i in range(no_shared_sunk):
            sunk_delta = shared_sunks[i][1] - shared_sunks[i][2]
            sunk_delta_mean += sunk_delta/no_shared_sunk
            sunk_deltas[i] = sunk_delta
        sunk_delta_sd = np.sqrt(np.sum((sunk_deltas - sunk_delta_mean)**2)/no_shared_sunk)
        start_r1, start_r2 = maf_data[r1][1], maf_data[r2][1]
        if start <= start_r1 <= end and start <= start_r2 <= end and len(shared_sunks) > threshold:
            read_set.add(r1)
            read_set.add(r2)
            overlap_locations.append([r1, start_r1, r2, start_r2, start_r1 - start_r2, np.round(sunk_delta_mean, 2), np.round(sunk_delta_sd, 6), no_shared_sunk])
#            print(sunk_delta_mean, sunk_delta_sd)
#    print(len(read_set))
    return overlap_locations


def pairwise_location_data_writer(outfile, overlap_locations):
    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in overlap_locations]


class DisjointSet:
    def __init__(self, vertices, parent):
        self.vertices = vertices
        self.parent = parent

    def find(self, item):
        if self.parent[item] == item:
            return item
        else:
            res = self.find(self.parent[item])
            self.parent[item] = res
            return res

    def union(self, set1, set2):
        root1 = self.find(set1)
        root2 = self.find(set2)
        self.parent[root1] = root2


def csv_reader(infile):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    return data


def get_distance_array(data, size, read_dict):
    distance = np.zeros((size, size))
    for i in range(size):
        distance[i, i] = 1000.0
    for i in range(1, len(data)):
        r1, r2, dist = int(data[i][0]), int(data[i][2]), float(data[i][-1])
        distance[read_dict[r1], read_dict[r2]], distance[read_dict[r2], read_dict[r1]] = dist, dist
    X = []
    for i in range(size):
        for j in range(i+1, size):
            X.append([distance[i, j]])
    return distance, X


def get_adjacency_list(data, size, read_dict):
    adj_list = dict()
    for i in range(1, size):
        r1, r2, dist = read_dict[int(data[i][0])], read_dict[int(data[i][2])], int(data[i][-1])
        if r1 not in adj_list:
            adj_list[r1] = dict()
        if r2 not in adj_list:
            adj_list[r2] = dict()
        adj_list[r1][r2], adj_list[r2][r1] = dist, dist
    return adj_list


def dfs(graph, temp, v, visited):
    # Mark the current vertex as visited
    visited[v] = True
    # Store the vertex to list
    temp.append(v)
    # Repeat for all vertices adjacent
    # to this vertex v
    for i in graph[v].keys():
        if not visited[i]:
            # Update the list
            temp = dfs(graph, temp, i, visited)
    return temp


def get_connected_components(graph):
    visited = dict()
    cc = []
    vertices = sorted(list(graph.keys()))
    for i in vertices:
        visited[i] = False
    for v in vertices:
        if not visited[v]:
            temp = []
            cc.append(dfs(graph, temp, v, visited))
    return cc


def get_clusters_cc(ccs, reverse_read_dict):
    clusters = dict()
    count = 0
    for cc in ccs:
        cluster = []
        for r in cc:
            cluster.append(reverse_read_dict[r])
        clusters[count] = cluster
        count += 1
    return clusters


def get_read_set(data):
    reads = set()
    for i in range(1, len(data)):
        r1, r2 = int(data[i][0]), int(data[i][2])
        reads.add(r1)
        reads.add(r2)
    return reads


def convert_data_to_float(data, sorted):
    new_data = []
    size = len(data)
    for i in range(1, size):
        new_row = [float(r) for r in data[i]]
        new_data.append(new_row)
    if sorted:
        new_data.sort(key=lambda x: x[-1], reverse=True)
    return new_data


def get_disjoint_set(reads, data, greedy):
    data = convert_data_to_float(data, greedy)
    parent = {}
    for r in reads:
        parent[r] = r
    ds = DisjointSet(reads, parent)
    size = len(data)
    for i in range(size):
        r1, r2 = int(data[i][0]), int(data[i][2])
        ds.union(r1, r2)
    return ds


def get_clusters_naive(reads, ds):
    clusters = dict()
    for r in reads:
        parent = ds.find(r)
        if parent not in clusters:
            clusters[parent] = [r]
        else:
            clusters[parent] += [r]
    return clusters


def index_reads(reads):
    count = 0
    read_dict = dict()
    reverse_read_dict = dict()
    for r in reads:
        read_dict[r] = count
        reverse_read_dict[count] = r
        count += 1
    return read_dict, reverse_read_dict


def write_clusters(outfile, clusters):
    outfile.write(str(len(clusters.keys()))+"\n")
    for key in clusters.keys():
        outfile.write(','.join(str(i) for i in clusters[key]))
        outfile.write("\n")


def read_clusters(cluster_file):
    no_of_clusters = int(cluster_file.readline().strip())
    clusters = []
    for i in range(no_of_clusters):
        reads = cluster_file.readline().strip().split(",")
        clusters.append([int(j) for j in reads])
    return clusters, no_of_clusters


def read_alignment(alignment_file):
    data = []
    with open(alignment_file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    read_position_dict = dict()
    for i in range(1, len(data)):
        read_position_dict[int(data[i][0])] = int(data[i][1])
    return read_position_dict


def write_read_clusters_to_fasta(cluster_file, read_file, path_to_output):
    clusters, no_of_clusters = read_clusters(cluster_file)
    read_dict = dict()
    for record in read_file:
        id = int(str(record.id).split("_")[1])
        read_dict[id] = record.seq
    for i in range(no_of_clusters):
        seqs = []
        for j in range(len(clusters[i])):
            record = SeqRecord(Seq(read_dict[clusters[i][j]]), id=str(clusters[i][j]),
                               description="cluster_"+str(i)+"_read_"+str(clusters[i][j]))
            seqs.append(record)
        SeqIO.write(seqs, path_to_output + str(i) + ".fasta", "fasta")


def validate_clusters(cluster_file, alignment_file, path_to_output):
    clusters, no_of_clusters = read_clusters(cluster_file)
    read_position_dict = read_alignment(alignment_file)
    cluster_profiles = []
    cluster_ranges = dict()
    for i in range(no_of_clusters):
        cluster_profiles.append([["id", "pos"]])
        max_, min_ = 0, 200000
        for read in clusters[i]:
            pos = read_position_dict[int(read)]
            if pos > max_:
                max_ = pos
            if pos < min_:
                min_ = pos
            cluster_profiles[-1].append([read, str(pos)])
        cluster_ranges[i] = (min_, max_)
    for i in range(no_of_clusters):
        with open(path_to_output + str(i) + ".csv", "w") as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            [writer.writerow(r) for r in cluster_profiles[i]]
    return cluster_ranges


def get_clusters_from_sunks_and_reads(sunk_file, read_file1, read_file2, ref_size, depth, parent_dir, method):
    sunk_map = get_sunk_map(sunk_file)
    out_sunk_map = inters_path+parent_dir+"/"+str(depth)+"/sunk_map.txt"
    write_sunk_map(out_sunk_map, sunk_map)
    r2s, s2r = get_read_to_sunk_and_sunk_to_read_maps(sunk_map, read_file1)
    out_r2s = inters_path+parent_dir+"/"+str(depth)+"/read_to_sunk_map.txt"
    out_s2r = inters_path+parent_dir+"/"+str(depth)+"/sunk_to_read_map.txt"
    write_read_to_sunk_and_sunk_to_read_maps(out_r2s, out_s2r, r2s, s2r)
    print("--- k-mer barcoding done ---")
    in_r2s = open(out_r2s, "r")
    in_s2r = open(out_s2r, "r")
    r2s_map, r2s_list, r2sp_map = get_read_to_sunk_map(in_r2s)
    s2r_map = get_sunk_to_read_map(in_s2r)
    pw_profiles = get_pairwise_profiles(s2r_map, r2s_map, r2s_list, r2sp_map)
    out_pw_profiles = inters_path+parent_dir+"/"+str(depth)+"/pairwise_profiles.txt"
    write_pairwise_profiles(pw_profiles, out_pw_profiles)
    print("--- pairwise profiling done ---")
    in_maf = reads_path + parent_dir + "/" + str(depth) + ".maf"
    maf_data = maf_parser(in_maf)
    out_alignment = inters_path+parent_dir+"/"+str(depth)+"/alignment.csv"
    maf_writer(out_alignment, maf_data)
    start, end, threshold = 0, ref_size, 20
    overlap_data = get_pairwise_location_data(maf_data, pw_profiles, start, end, threshold)
    out_overlap = inters_path+parent_dir+"/"+str(depth)+"/overlap.csv"
    pairwise_location_data_writer(out_overlap, overlap_data)
    print("--- alignment data parsing for validation done ---")
    reads = list(get_read_set(overlap_data))
#    print("# of reads:", len(reads))
    out_path = clusters_path + parent_dir + "/" + str(depth) + "/" + method + "/"
    if method == "naive":
        greedy = False
        ds = get_disjoint_set(reads, overlap_data, greedy)
        clusters = get_clusters_naive(reads, ds)
        if greedy:
            out_path = out_path + "greedy.txt"
            outfile = open(out_path, "w")
        else:
            out_path = out_path+"basic.txt"
            outfile = open(out_path, "w")
    else:
        sz = len(overlap_data)
        read_dict, reverse_read_dict = index_reads(reads)
        graph = get_adjacency_list(overlap_data, sz, read_dict)
        cc = get_connected_components(graph)
        clusters = get_clusters_cc(cc, reverse_read_dict)
        out_path = out_path + "basic.txt"
        outfile = open(out_path, "w")
    write_clusters(outfile, clusters)
    outfile.close()
    print("--- clustering based on overlaps done ---")
    #out_path = clusters_path + parent_dir + "/" + str(depth) + "/" + method + "/basic.txt"
    cluster_file = open(out_path, "r")
    out_clusters = clusters_path + parent_dir + "/" + str(depth) + "/" + method + "/"
    write_read_clusters_to_fasta(cluster_file, read_file2, out_clusters)
    print("--- writing clusters to fasta files done ---")


k = 21
repeat_sizes = [5000, 10000, 15000, 20000]
copies = [2, 5, 10, 20]
snps = [100, 250, 500, 1000, 2000]
depths = [10, 20, 30, 40, 50]
methods = ["naive", "cc"]
reads_path = "../data/HiFi/"
sunks_path = "../data/kmers/"
inters_path = "../intermediates/"
clusters_path = "../clusters/"
for repeat_size in [5000]:
    for copy in copies:
        for snp in snps:
            parent_dir = str(repeat_size) + "_" + str(copy) + "_" + str(snp)
            for depth in depths:
                print("Start:", parent_dir+"_"+str(depth))
                sunk_file = open(sunks_path + parent_dir + "/" + str(depth) + ".kmers")
                read_file1 = SeqIO.parse(reads_path + parent_dir + "/" + str(depth) + ".fasta", "fasta")
                read_file2 = SeqIO.parse(reads_path + parent_dir + "/" + str(depth) + ".fasta", "fasta")
                ref_size = 100000 + (copy * repeat_size)
                get_clusters_from_sunks_and_reads(sunk_file, read_file1, read_file2, ref_size, depth, parent_dir, methods[0])
                print("Done:", parent_dir+"_"+str(depth))