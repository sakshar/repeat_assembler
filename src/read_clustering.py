# cluster reads based on the similarity of pairwise k-mer profiling
# ML-based: using one-hot encoding for sunks per read
# Distance-based: using pairwise profiling
import csv
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


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


def get_pairwise_profiles(infile):
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


def get_clusters_kmeans(data, k, iter, sz, reverse_read_dict):
    centroids, labels = kmeans2(data, k, iter=iter)
    clusters = dict()
    for i in range(k):
        clusters[i] = []
    for i in range(sz):
        clusters[labels[i]] += [reverse_read_dict[i]]
    return centroids, labels, clusters


def write_clusters(outfile, clusters):
    outfile.write(str(len(clusters.keys()))+"\n")
    for key in clusters.keys():
        outfile.write(','.join(str(i) for i in clusters[key]))
        outfile.write("\n")


#infile = open("simulated/k21/output/sim5_depth50_nano_pairwise_profiles_optimized.txt", "r")
#pw_profiles = get_pairwise_profiles(infile)
#print(pw_profiles[(21, 557)])
"""
size = 739 # number of reads
infile = "simulated/HiFi/hifi.sim5_depth50_max20000_overlap.csv"
data = csv_reader(infile)
distance, X = get_distance_array(data, size)
"""
#Z = linkage(X, 'single')
#fig = plt.figure(figsize=(25, 10))
#dn = dendrogram(Z)
#plt.show()
infile = "simulated/HiFi/hifi.sim5_depth50_max20000_overlap2.csv"
data = csv_reader(infile)
reads = list(get_read_set(data))
print("# of reads:", len(reads))
greedy = False
ds = get_disjoint_set(reads, data, greedy)
clusters_naive = get_clusters_naive(reads, ds)
print("--------- naive clustering statistics ---------")
print("# of clusters:", len(clusters_naive.keys()))
for key in clusters_naive:
    print("parent:", key, ">")
    print("size:", len(clusters_naive[key]))
    print(clusters_naive[key])
if greedy:
    outfile_naive = open("simulated/k21/output/sim5_depth50_max20000_hifi_clusters_greedy_naive2.txt", "w")
else:
    outfile_naive = open("simulated/k21/output/sim5_depth50_max20000_hifi_clusters_naive2.txt", "w")
write_clusters(outfile_naive, clusters_naive)
outfile_naive.close()

"""
print("--------- k-means clustering statistics ---------")
sz = len(reads)
read_dict, reverse_read_dict = index_reads(reads)
dist, X = get_distance_array(data, sz, read_dict)
k, iter = 7, 100
centroids, labels, clusters_kmeans = get_clusters_kmeans(dist, k, iter, sz, reverse_read_dict)
print("# of clusters:", k)
for i in range(k):
    print("cluster ID:", i)
    print("size: ", len(clusters_kmeans[i]))
    print(clusters_kmeans[i])
outfile_kmeans = open("simulated/k21/output/sim5_depth50_max20000_hifi_clusters_kmeans.txt", "w")
write_clusters(outfile_kmeans, clusters_kmeans)
outfile_kmeans.close()
"""
print("--------- connected-components based clustering statistics ---------")
sz = len(data)
read_dict, reverse_read_dict = index_reads(reads)
graph = get_adjacency_list(data, sz, read_dict)
print(len(list(graph.keys())))
cc = get_connected_components(graph)
print("# of clusters:", len(cc))
clusters_cc = get_clusters_cc(cc, reverse_read_dict)
for i in range(len(cc)):
    print("cluster ID:", i)
    print("size: ", len(clusters_cc[i]))
    print(clusters_cc[i])
outfile_cc = open("simulated/k21/output/sim5_depth50_max20000_hifi_clusters_cc2.txt", "w")
write_clusters(outfile_cc, clusters_cc)
outfile_cc.close()