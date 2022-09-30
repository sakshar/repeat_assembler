# outputs the pairwise kmer profiles of reads (SUNK_id, position and conserved relative distance)
# improved, efficient and hopefully no bugs
k = 21


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
        print(count, ":", sunk)
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


infile_read_to_sunk = open("simulated/k21/output/sim5_depth50_max20000_hifi_read_to_sunk_map.txt", "r")
infile_sunk_to_read = open("simulated/k21/output/sim5_depth50_max20000_hifi_sunk_to_read_map.txt", "r")
r2s_map, r2s_list, r2sp_map = get_read_to_sunk_map(infile_read_to_sunk)
s2r_map = get_sunk_to_read_map(infile_sunk_to_read)
pw_profiles = get_pairwise_profiles(s2r_map, r2s_map, r2s_list, r2sp_map)
outfile = "simulated/k21/output/sim5_depth50_max20000_hifi_pairwise_profiles_optimized.txt"
write_pairwise_profiles(pw_profiles, outfile)