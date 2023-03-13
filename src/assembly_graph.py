# codes for building overlap graph from shared sunk profile and then traversing a hamiltonian path through it
# to obtain a reconstructed assembly with improved accuracy at repeat regions
# add a function to write clusters of reads to separate fasta files, then use HiCanu & HiFiasm on them

import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def copy_contigs(contig_file, coverage, path_to_output):
    contig_dict = dict()
    for record in contig_file:
        id = str(record.id)
        contig_dict[id] = record.seq
    seqs = []
    for i in contig_dict.keys():
        for j in range(coverage):
            record = SeqRecord(Seq(contig_dict[i]), id=i+"_"+str(j),
                               description="contig_" + i + "_copy_" + str(j))
            seqs.append(record)
    SeqIO.write(seqs, path_to_output + "_" + str(coverage) + ".fasta", "fasta")


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


clustering_methods = ["naive", "kmeans", "cc"]
cluster_file = open("simulated/k21/output/sim5_depth50_max20000_hifi_clusters_" + clustering_methods[0] + ".txt", "r")
read_file = SeqIO.parse("simulated/HiFi/hifi.sim5_depth50_max20000.fasta", "fasta")
alignment_file = "simulated/HiFi/hifi.sim5_depth50_max20000_alignment.csv"
path_to_output = "simulated/k21/output/sim5_depth50_max20000_hifi/" + clustering_methods[0]
#write_read_clusters_to_fasta(cluster_file, read_file, path_to_output)
ranges = validate_clusters(cluster_file, alignment_file, path_to_output)
for i in range(len(ranges.keys())):
    print(i, ":", ranges[i])
"""
coverage = 30
contig_file = SeqIO.parse("simulated/k21/output/sim5_depth50_max20000_hifi/naive/naive.asm.fasta", "fasta")
path_to_output = "simulated/k21/output/sim5_depth50_max20000_hifi/naive/naive"
copy_contigs(contig_file, coverage, path_to_output)
"""