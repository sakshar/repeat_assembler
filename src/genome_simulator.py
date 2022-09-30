from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint

# code for genome simulation
base_map = {0: 'A', 1: 'C', 2:'G', 3: 'T'}


def genome_simulator(unique_seq_size, repeat_size, no_of_copies, snp_rate):
    sim_seq = ""
    for i in range(unique_seq_size):
        sim_seq += base_map[randint(0, 3)]
    repeat = ""
    for i in range(repeat_size):
        repeat += base_map[randint(0, 3)]
    offsets = []
    for j in range(repeat_size//snp_rate):
        offset = randint(0, snp_rate - 1)
        offsets.append(offset)
    sim_seq += repeat
    for i in range(no_of_copies - 1):
        new_repeat = ""
        new_repeat += repeat
        for j in range(repeat_size//snp_rate):
            char = randint(0, 3)
            new_repeat = new_repeat[:j * snp_rate + offsets[j]] + base_map[char] + new_repeat[j * snp_rate + offsets[j] + 1:]
        sim_seq += new_repeat
    for i in range(unique_seq_size):
        sim_seq += base_map[randint(0, 3)]
    filename = "unique_"+str(unique_seq_size)+"_repeat_"+str(repeat_size)+"_copies_"+str(no_of_copies)+"_snp_"+str(snp_rate)
    seqs = []
    print("Genome size:", len(sim_seq))
    for i in range(1):
        record = SeqRecord(Seq(sim_seq), id=filename,
                           description="simulated with " + str(no_of_copies) + " repeats of length " + str(repeat_size) + " with SNP/" + str(snp_rate) + " bp")
        seqs.append(record)
    SeqIO.write(seqs, "simulated/genomes/"+filename+".fasta", "fasta")
    offset_file = open("simulated/genomes/"+filename+"_offset.txt", "w")
    offset_file.write(','.join([str(i) for i in offsets]))
    offset_file.close()


def single_copy_repeat_extractor(infile):
    genome = SeqIO.parse("simulated/genomes/"+infile, "fasta")
    genome_seq = ""
    for record in genome:
        genome_seq = str(record.seq)
    filename = infile.split(".")[0]
    tokens = filename.split("_")
    unique_size, repeat_size, no_of_copies, snp_rate = int(tokens[1]), int(tokens[3]), int(tokens[5]), int(tokens[7])
    genome_seq = genome_seq[:unique_size+repeat_size] + genome_seq[unique_size + (no_of_copies * repeat_size):]
    print("Extracted genome size:", len(genome_seq))
    seqs = []
    for i in range(1):
        record = SeqRecord(Seq(genome_seq), id=filename+"_extracted",
                           description="extracted first repeat from " + str(no_of_copies) + " repeats of length " + str(repeat_size) + " with SNP/" + str(snp_rate) + " bp")
        seqs.append(record)
    SeqIO.write(seqs, "simulated/genomes/"+filename+"_extracted.fasta", "fasta")


unique_size, repeat_size, no_of_copies, snp_rate = 50000, 10000, 5, 500
filename = "unique_"+str(unique_size)+"_repeat_"+str(repeat_size)+"_copies_"+str(no_of_copies)+"_snp_"+str(snp_rate)
genome_simulator(unique_size, repeat_size, no_of_copies, snp_rate)
single_copy_repeat_extractor(filename+".fasta")