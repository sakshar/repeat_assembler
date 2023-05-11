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


def run(indel=0.0):
    print(indel)
    repeat_sizes = [15000, 20000, 25000]
    no_of_copies = [5, 10]
    snp_rates = [250, 500, 1000]
    for i in repeat_sizes:
        for j in no_of_copies:
            for k in snp_rates:
                print("repeat size:", i, "copies:", j, "snp_rate:", k)
                genome_simulator_indel_with_random_snp(i, j, k, indel)
                #filename = str(i) + "_" + str(j) + "_" + str(k)
                #single_copy_repeat_extractor(filename + ".fasta")


def genome_simulator_indel(repeat_size, no_of_copies, snp_rate):
    sim_seq = ""
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    repeat = ""
    for i in range(repeat_size):
        repeat += base_map[randint(0, 3)]
    offsets = []
    repeat_sizes = []
    for j in range(repeat_size//snp_rate):
        offset = randint(0, snp_rate - 1)
        offsets.append(offset)
    if snp_rate == 2000 and repeat_size in [5000, 15000, 25000]:
        offsets.append(randint(0, 999))
    sim_seq += repeat
    repeat_sizes.append(len(repeat))
    insertion = ""
    insert_offsets = []
    for i in range(repeat_size//10):
        insertion += base_map[randint(0, 3)]
    for j in range((repeat_size//10)//snp_rate):
        insert_offset = randint(0, snp_rate - 1)
        insert_offsets.append(insert_offset)
    if snp_rate == 1000 and repeat_size//10 in [1500, 2500]:
        insert_offsets.append(randint(0, 499))
    for i in range(no_of_copies - 1):
        new_repeat = ""
        new_repeat += repeat
        for j in range(int(math.ceil(repeat_size//snp_rate))):
            char = randint(0, 3)
            new_repeat = new_repeat[:j * snp_rate + offsets[j]] + base_map[char] + new_repeat[j * snp_rate + offsets[j] + 1:]
        in_del = randint(0, 1)
        in_del_size = int(randint(0, 10)/100*repeat_size)
        if in_del == 0: #deletion
            #del_pos = randint(0, repeat_size - in_del_size - 1)
            #new_repeat = new_repeat[:del_pos] + new_repeat[del_pos+in_del_size:]
            new_repeat = new_repeat[:repeat_size-in_del_size]
        elif in_del == 1: #insertion
            in_pos = randint(0, repeat_size - 1)
            new_insertion = ""
            new_insertion += insertion
            for k in range(int(math.ceil((repeat_size // 10) // snp_rate))):
                char = randint(0, 3)
                new_insertion = new_insertion[:k * snp_rate + insert_offsets[k]] + base_map[char] + new_insertion[
                                                                                       k * snp_rate + insert_offsets[k] + 1:]
            #new_repeat = new_repeat[:in_pos] + insertion + new_repeat[in_pos:]
            new_repeat = new_repeat + new_insertion[:in_del_size]
        sim_seq += new_repeat
        repeat_sizes.append(len(new_repeat))
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    filename = str(repeat_size)+"_"+str(no_of_copies)+"_"+str(snp_rate)
    seqs = []
    print("Genome size:", len(sim_seq))
    for i in range(1):
        record = SeqRecord(Seq(sim_seq), id=filename,
                           description="simulated with " + str(no_of_copies) + " repeats of length " + str(repeat_size) + " with SNP/" + str(snp_rate) + " bp")
        seqs.append(record)
    SeqIO.write(seqs, "../data/ref/indel2/"+filename+".fasta", "fasta")
    offset_file = open("../data/ref/indel2/"+filename+"_offset.txt", "w")
    offset_file.write(','.join([str(i) for i in offsets])+'\n')
    offset_file.write(','.join([str(i) for i in insert_offsets]) + '\n')
    offset_file.write(','.join([str(i) for i in repeat_sizes]))
    offset_file.close()


def genome_simulator_indel_with_random_snp(repeat_size, no_of_copies, snp_rate, indel):
    sim_seq = ""
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    repeat = ""
    for i in range(repeat_size + int(repeat_size * indel)):
        repeat += base_map[randint(0, 3)]
    snp_sites = []
    repeat_sizes = []
    parent = ""
    snps = []
    for j in range(repeat_size + int(repeat_size * indel)):
        snp_site = randint(0, snp_rate - 1)
        if snp_site == 0:
            snp_sites.append((j, repeat[j]))
            parent += repeat[j]
    print("total snps:", len(snp_sites))
    print(snp_sites)
    #if snp_rate == 2000 and repeat_size in [5000, 15000, 25000]:
    #    offsets.append(randint(0, 999))
    sim_seq += repeat[:repeat_size]
    repeat_sizes.append(repeat_size)
    #insertion = ""
    #insert_offsets = []
    #for i in range(repeat_size//10):
    #    insertion += base_map[randint(0, 3)]
    #for j in range((repeat_size//10)//snp_rate):
    #    insert_offset = randint(0, snp_rate - 1)
    #    insert_offsets.append(insert_offset)
    #if snp_rate == 1000 and repeat_size//10 in [1500, 2500]:
    #    insert_offsets.append(randint(0, 499))
    for i in range(no_of_copies - 1):
        new_repeat = ""
        new_repeat_size = randint(repeat_size - int(repeat_size * indel), repeat_size + int(repeat_size * indel))
        print("repeat_copy:", i+2, "size:", new_repeat_size)
        new_repeat += repeat[:new_repeat_size]
        child = ""
        j = 0
        while j < len(snp_sites) and snp_sites[j][0] < new_repeat_size:
            char = randint(0, 3)
            while base_map[char] == snp_sites[j][1]:
                char = randint(0, 3)
            new_repeat = new_repeat[:snp_sites[j][0]] + base_map[char] + new_repeat[snp_sites[j][0] + 1:]
            child += base_map[char]
            j += 1
        sim_seq += new_repeat
        repeat_sizes.append(new_repeat_size)
        snps.append(child)
    for i in range(unique_size):
        sim_seq += base_map[randint(0, 3)]
    filename = str(repeat_size)+"_"+str(no_of_copies)+"_"+str(snp_rate)
    seqs = []
    print("Genome size:", len(sim_seq))
    for i in range(1):
        record = SeqRecord(Seq(sim_seq), id=filename, description="")
        seqs.append(record)
    SeqIO.write(seqs, "../data/ref/indel5/"+filename+".fasta", "fasta")
    snpsite_file = open("../data/ref/indel5/"+filename+"_snpsites.txt", "w")
    snpsite_file.write('total snps:'+str(len(snp_sites))+'\n')
    snpsite_file.write(','.join(["("+str(i[0])+","+str(i[1])+")" for i in snp_sites])+'\n')
    snpsite_file.write(','.join([str(i) for i in repeat_sizes])+'\n')
    snpsite_file.write(parent+'\n')
    snpsite_file.write('\n'.join([i for i in snps]))
    snpsite_file.close()


"""
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
"""

unique_size = 50000
repeat_size, no_of_copies, snp_rate = 15000, 10, 1000
#filename = "indel_"+str(repeat_size)+"_"+str(no_of_copies)+"_"+str(snp_rate)
#genome_simulator_indel(repeat_size, no_of_copies, snp_rate)
#single_copy_repeat_extractor(filename+".fasta")
indel = 0.05
run(indel)
