import numpy as np
from Bio import AlignIO
import csv


def maf_parser(infile):
    maf = AlignIO.parse(infile+".maf", "maf")
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


def maf_writer(infile, maf_profile):
    with open(infile+"_alignment.csv", "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in maf_profile]


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
            print(sunk_delta_mean, sunk_delta_sd)
    print(len(read_set))
    return overlap_locations


def pairwise_location_data_writer(infile, overlap_locations):
    with open(infile+"_overlap.csv", "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in overlap_locations]


infile_maf = "simulated/HiFi/hifi.sim5_depth50_max20000"
maf_data = maf_parser(infile_maf)
infile_pw = open("simulated/k21/output/sim5_depth50_max20000_hifi_pairwise_profiles_optimized.txt", "r")
pw_profiles = get_pairwise_profiles(infile_pw)
#maf_writer(infile, maf_data)
start, end, threshold = 0, 200000, 20
overlap_data = get_pairwise_location_data(maf_data, pw_profiles, start, end, threshold)
pairwise_location_data_writer(infile_maf, overlap_data)