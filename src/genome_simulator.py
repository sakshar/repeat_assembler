from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint
import math

# code for genome simulation
base_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def genome_simulator(repeat_size, no_of_copies, snp_rate):
    sim_seq = ""
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    repeat = ""
    for i in range(repeat_size):
        repeat += base_map[randint(0, 3)]
    offsets = []
    for j in range(repeat_size//snp_rate):
        offset = randint(0, snp_rate - 1)
        offsets.append(offset)
    if snp_rate == 2000 and repeat_size in [5000, 15000, 25000]:
        offsets.append(randint(0, 999))
    sim_seq += repeat
    for i in range(no_of_copies - 1):
        new_repeat = ""
        new_repeat += repeat
        for j in range(int(math.ceil(repeat_size//snp_rate))):
            char = randint(0, 3)
            new_repeat = new_repeat[:j * snp_rate + offsets[j]] + base_map[char] + new_repeat[j * snp_rate + offsets[j] + 1:]
        sim_seq += new_repeat
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    filename = str(repeat_size)+"_"+str(no_of_copies)+"_"+str(snp_rate)
    seqs = []
    print("Genome size:", len(sim_seq))
    for i in range(1):
        record = SeqRecord(Seq(sim_seq), id=filename,
                           description="simulated with " + str(no_of_copies) + " repeats of length " + str(repeat_size) + " with SNP/" + str(snp_rate) + " bp")
        seqs.append(record)
    SeqIO.write(seqs, "../data/ref/"+filename+".fasta", "fasta")
    offset_file = open("../data/ref/"+filename+"_offset.txt", "w")
    offset_file.write(','.join([str(i) for i in offsets]))
    offset_file.close()


def single_copy_repeat_extractor(infile):
    genome = SeqIO.parse("../data/ref/"+infile, "fasta")
    genome_seq = ""
    for record in genome:
        genome_seq = str(record.seq)
    filename = infile.split(".")[0]
    tokens = filename.split("_")
    repeat_size, no_of_copies, snp_rate = int(tokens[0]), int(tokens[1]), int(tokens[2])
    genome_seq = genome_seq[:unique_size+repeat_size] + genome_seq[unique_size + (no_of_copies * repeat_size):]
    print("Extracted genome size:", len(genome_seq))
    seqs = []
    for i in range(1):
        record = SeqRecord(Seq(genome_seq), id=filename+"_extracted",
                           description="extracted first repeat from " + str(no_of_copies) + " repeats of length " + str(repeat_size) + " with SNP/" + str(snp_rate) + " bp")
        seqs.append(record)
    SeqIO.write(seqs, "../data/ref/"+filename+"_extracted.fasta", "fasta")


def run():
    repeat_sizes = [20000]
    no_of_copies = [5]
    snp_rates = [100, 250, 500, 2000]
    for i in repeat_sizes:
        for j in no_of_copies:
            for k in snp_rates:
                print("repeat size:", i, "copies:", j, "snp_rate:", k)
                genome_simulator(i, j, k)
                filename = str(i) + "_" + str(j) + "_" + str(k)
                single_copy_repeat_extractor(filename + ".fasta")


file = "../arabidopsis/ref/chr5_11300k-11350k.fasta"
genome = SeqIO.parse(file, "fasta")
for record in genome:
    id = str(record.id)
    seq = record.seq
new_seq = seq[:31500] + seq[20500:31500]*3 + seq[31500:]
print(len(new_seq))
seqs = []
for i in range(1):
    record = SeqRecord(Seq(new_seq), id=id+"_replicated")
    seqs.append(record)
SeqIO.write(seqs, "../arabidopsis/ref/chr5_11300k-11350k_replicated.fasta", "fasta")

unique_size = 50000
"""
repeat_size, no_of_copies, snp_rate = 10000, 5, 1000
filename = str(repeat_size)+"_"+str(no_of_copies)+"_"+str(snp_rate)
genome_simulator(repeat_size, no_of_copies, snp_rate)
single_copy_repeat_extractor(filename+".fasta")
"""
#run()
