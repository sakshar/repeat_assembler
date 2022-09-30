import random

base_dic = {0:'A', 1:'C', 2:'G', 3:'T'}
reverse_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}


def reverse_complement(dna_seq):
    seq_len = len(dna_seq)
    rev = ""
    #print(seq_len)
    for i in range(seq_len):
        #print(dna_seq[-(i+1)], reverse_map[dna_seq[-(i+1)]])
        rev += reverse_map[dna_seq[-(i+1)]]
    #print(rev)
    return rev


# Returns index of x in arr if present, else -1
def binary_search(arr, low, high, x):
    # Check base case
    if high >= low:

        mid = (high + low) // 2

        # If element is present at the middle itself
        if arr[mid] == x:
            return mid

        # If element is smaller than mid, then it can only
        # be present in left subarray
        elif arr[mid] > x:
            return binary_search(arr, low, mid - 1, x)

        # Else the element can only be present in right subarray
        else:
            return binary_search(arr, mid + 1, high, x)

    else:
        # Element is not present in the array
        return -1


"""
repeat = 12
rd = "A"*repeat + 'C' + "G"*repeat + 'T' + "A"*repeat

#print(rd)

genome = ""
f = open("./data", 'r')
while True:
    line = f.readline()
    if not line:
        break
    if line[0] != '>':
        line = line.strip()
        genome += line
genome = genome.upper()
f.close()
#print(genome[:5])
print("Genome size:", len(genome))

kmers = {}

k = 20

sz = len(genome)

for i in range(sz - k + 1):
    kmer = genome[i:i+k]
    if kmer in kmers:
        kmers[kmer] += 1
        continue
    kmers[kmer] = 1

#print(len(kmers))

unique_kmers = []
for key in kmers:
    if kmers[key] == 1:
        unique_kmers.append(key)

print(len(unique_kmers))

f = open("rice_MH63_unique_kmers_"+str(k)+".txt", "w")
for ukmer in unique_kmers:
    f.write(ukmer+'\n')
f.close()
"""
parts = 12
unique_kmers = []
for i in range(1, parts+1):
    part_kmers = []
    f1 = open("./encoded/drosophila/k25/genome/"+str(i), 'r')
    while True:
        line = f1.readline()
        if not line:
            break
        part_kmers.append(int(line.strip()))
    total_unique_kmers = len(part_kmers)
    print(total_unique_kmers)
    unique_kmers.append(part_kmers)

kmer_files = ["178_178", "163_177", "179_193", "147_162", "194_209", "132_146", "210_224", "116_131", "225_240"]
tp_fp_pairs = []
last_tp, last_fp = 0, 0
for file in kmer_files:
    tp, fp = 0, 0
    f = file.split("_")
    for m in range(int(f[0]), int(f[1])+1):
        f2 = open("./encoded/drosophila/k25/pacbio/"+str(m), 'r')
        while True:
            line = f2.readline()
            if not line:
                break
            line = int(line.strip())
            res = -1
            for i in range(parts):
                res = binary_search(unique_kmers[i], 0, len(unique_kmers[i]) - 1, line)
                if res != -1:
                    tp += 1
                    break
            if res == -1:
                fp += 1
    print("count range: " + file)
    print('true positive:', tp)
    print('false positive:', fp)
    last_tp += tp
    last_fp += fp
    tp_fp_pairs.append((last_tp, last_fp))

print(tp_fp_pairs)
