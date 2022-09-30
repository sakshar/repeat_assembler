# iterates through every read and every SUNK and outputs SUNK_id and position in each read
from Bio import SeqIO
from Bio.Seq import Seq


k = 21
input_sunks_hifi = open("simulated/k21/hifi.sim5_depth50_max20000.20_75.kmers")
#input_sunks_ont = open("simulated/k21/nano.sim5_depth50.24_66.kmers")
input_reads_hifi = SeqIO.parse("simulated/hifi.sim5_depth50_max20000.fasta", "fasta")
#input_reads_ont = SeqIO.parse("simulated/nano.sim5_depth50.fasta", "fasta")

hifi_sunk_map = dict()
#ont_sunk_map = dict()
counter = 0
while True:
    line = input_sunks_hifi.readline()
    if not line:
        break
    hifi_sunk_map[line.strip()] = counter
    counter += 1
"""
counter = 0
while True:
    line = input_sunks_ont.readline()
    if not line:
        break
    ont_sunk_map[line.strip()] = counter
    counter += 1
print("ONT:", len(ont_sunk_map.keys()))
"""
print("HiFi:", len(hifi_sunk_map.keys()))

#output_ont_sunk_map = open("simulated/k21/output/sim5_depth50_nano_sunk_map.txt", "w")
output_hifi_sunk_map = open("simulated/k21/output/sim5_depth50_max20000_hifi_sunk_map.txt", "w")

"""
for sunk in ont_sunk_map.keys():
    output_ont_sunk_map.write(str(sunk)+","+str(ont_sunk_map[sunk])+"\n")
output_ont_sunk_map.close()
"""
for sunk in hifi_sunk_map.keys():
    output_hifi_sunk_map.write(str(sunk)+","+str(hifi_sunk_map[sunk])+"\n")
output_hifi_sunk_map.close()

print("HiFi SUNKs mapping done.")

#ont_read_to_sunk_map = dict()
hifi_read_to_sunk_map = dict()
#ont_sunk_to_read_map = dict()
hifi_sunk_to_read_map = dict()

"""
ont_sunks = ont_sunk_map.keys()
for record in input_reads_ont:
    read = record.seq
    id = record.id
    print("ONT:", str(id))
    for i in range(len(read) - k + 1):
        #for sunk in ont_sunk_map.keys():
        kmer = read[i:i + k]
        str_kmer = str(kmer)
        str_rev_kmer = str(kmer.reverse_complement())
        if str_kmer in ont_sunks:
            if id not in ont_read_to_sunk_map.keys():
                ont_read_to_sunk_map[id] = []
            ont_read_to_sunk_map[id] += [(ont_sunk_map[str_kmer], i)]  #(sunk_id, pos)
            if ont_sunk_map[str_kmer] not in ont_sunk_to_read_map.keys():
                ont_sunk_to_read_map[ont_sunk_map[str_kmer]] = []
            ont_sunk_to_read_map[ont_sunk_map[str_kmer]] += [(id, i)]  #(read_id, pos)
        elif str_rev_kmer in ont_sunks:
            if id not in ont_read_to_sunk_map.keys():
                ont_read_to_sunk_map[id] = []
            ont_read_to_sunk_map[id] += [(ont_sunk_map[str_rev_kmer], i)]  #(sunk_id, pos)
            if ont_sunk_map[str_rev_kmer] not in ont_sunk_to_read_map.keys():
                ont_sunk_to_read_map[ont_sunk_map[str_rev_kmer]] = []
            ont_sunk_to_read_map[ont_sunk_map[str_rev_kmer]] += [(id, i)]  #(read_id, pos)
"""
hifi_sunks = hifi_sunk_map.keys()
for record in input_reads_hifi:
    read = record.seq
    id = record.id
    print("HiFi:", str(id))
    for i in range(len(read) - k + 1):
        #for sunk in hifi_sunk_map.keys():
        kmer = read[i:i + k]
        str_kmer = str(kmer)
        str_rev_kmer = str(kmer.reverse_complement())
        if str_kmer in hifi_sunks:
            if id not in hifi_read_to_sunk_map.keys():
                hifi_read_to_sunk_map[id] = []
            hifi_read_to_sunk_map[id] += [(hifi_sunk_map[str_kmer], i)]  # (sunk_id, pos)
            if hifi_sunk_map[str_kmer] not in hifi_sunk_to_read_map.keys():
                hifi_sunk_to_read_map[hifi_sunk_map[str_kmer]] = []
            hifi_sunk_to_read_map[hifi_sunk_map[str_kmer]] += [(id, i)]  # (read_id, pos)
        elif str_rev_kmer in hifi_sunks:
            if id not in hifi_read_to_sunk_map.keys():
                hifi_read_to_sunk_map[id] = []
            hifi_read_to_sunk_map[id] += [(hifi_sunk_map[str_rev_kmer], i)]  # (sunk_id, pos)
            if hifi_sunk_map[str_rev_kmer] not in hifi_sunk_to_read_map.keys():
                hifi_sunk_to_read_map[hifi_sunk_map[str_rev_kmer]] = []
            hifi_sunk_to_read_map[hifi_sunk_map[str_rev_kmer]] += [(id, i)]  # (read_id, pos)

output_hifi_read_to_sunk_map = open("simulated/k21/output/sim5_depth50_max20000_hifi_read_to_sunk_map.txt", "w")
#output_ont_read_to_sunk_map = open("simulated/k21/output/sim5_depth50_nano_read_to_sunk_map.txt", "w")
output_hifi_sunk_to_read_map = open("simulated/k21/output/sim5_depth50_max20000_hifi_sunk_to_read_map.txt", "w")
#output_ont_sunk_to_read_map = open("simulated/k21/output/sim5_depth50_nano_sunk_to_read_map.txt", "w")
"""
for read in ont_read_to_sunk_map.keys():
    output_ont_read_to_sunk_map.write(str(read)+">")
    for tup in ont_read_to_sunk_map[read]:
        output_ont_read_to_sunk_map.write(str(tup)+";")
    output_ont_read_to_sunk_map.write("\n")
output_ont_read_to_sunk_map.close()
"""
for read in hifi_read_to_sunk_map.keys():
    output_hifi_read_to_sunk_map.write(str(read)+">")
    for tup in hifi_read_to_sunk_map[read]:
        output_hifi_read_to_sunk_map.write(str(tup) + ";")
    output_hifi_read_to_sunk_map.write("\n")
output_hifi_read_to_sunk_map.close()
"""
for sunk in ont_sunk_to_read_map.keys():
    output_ont_sunk_to_read_map.write(str(sunk)+">")
    for tup in ont_sunk_to_read_map[sunk]:
        output_ont_sunk_to_read_map.write(str(tup) + ";")
    output_ont_sunk_to_read_map.write("\n")
output_ont_sunk_to_read_map.close()
"""
for sunk in hifi_sunk_to_read_map.keys():
    output_hifi_sunk_to_read_map.write(str(sunk)+">")
    for tup in hifi_sunk_to_read_map[sunk]:
        output_hifi_sunk_to_read_map.write(str(tup) + ";")
    output_hifi_sunk_to_read_map.write("\n")
output_hifi_sunk_to_read_map.close()

print("All output files writing done.")