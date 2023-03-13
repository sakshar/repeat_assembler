# This is a sample Python script_illu.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import xlsxwriter


def kmer_histo_generator(filename):
    kmer_map = dict()
    kmer_histo = dict()
    file = open(filename, 'r')
    while True:
        row = file.readline()
        if not row:
            break
        row = row.strip().split()
        freq = int(row[1])
        if freq not in kmer_histo:
            kmer_map[freq] = [row[0]]
            kmer_histo[freq] = 1
        else:
            kmer_map[freq].append(row[0])
            kmer_histo[freq] = kmer_histo[freq] + 1
    return kmer_histo, kmer_map


def write_histo_to_xcel(input):
    workbook = xlsxwriter.Workbook(input + '.xlsx')
    worksheet = workbook.add_worksheet()

    histo = open(input + '.histo', 'r')
    i = 1
    for line in histo:
        args = line.strip().split()
        worksheet.write("A" + str(i), args[0])
        worksheet.write("B" + str(i), args[1])
        i += 1
    workbook.close()


def write_read_stat_to_xcel(input):
    workbook = xlsxwriter.Workbook(input + '.xlsx')
    worksheet = workbook.add_worksheet()
    file = open(input + ".txt", "r")
    for i in range(1, 126):
        num = file.readline().strip().split(":")[1]
        depth = file.readline().strip().split(":")[1]
        mean_sd = file.readline().strip().split(":")[1]
        min_val = file.readline().strip().split(":")[1]
        max_val = file.readline().strip().split(":")[1]
        worksheet.write("A" + str(i), num)
        worksheet.write("B" + str(i), mean_sd)
        worksheet.write("C" + str(i), max_val)
        worksheet.write("D" + str(i), min_val)
        i += 1
    workbook.close()


repeat_sizes = [5000, 10000, 15000, 20000, 25000]
no_of_copies = [2, 5, 10, 20, 50]
snp_rates = [100, 250, 500, 1000, 2000]
depths = [10, 20, 30, 40, 50]
sizes = ["10k", "15k", "20k", "25k"]
path = '../read_stats/read_stat_'
#for size in sizes:
#    input_file = path + size
#    write_read_stat_to_xcel(input_file)
"""
for repeat_size in repeat_sizes:
    for copy in no_of_copies:
        for snp in snp_rates:
            for depth in depths:
                input = path + str(repeat_size) + "_" + str(copy) + "_" + str(snp) + "/" + str(depth)
                write_histo_to_xcel(input)
"""
write_histo_to_xcel("../arabidopsis/ERR6210723")
"""
histo_, map_ = kmer_histo_generator('./yeast_complex/k21/hifiasm.histo')
workbook = xlsxwriter.Workbook('./yeast_complex/k21/hifiasm.xlsx')
worksheet = workbook.add_worksheet()
kmer_file = open('./yeast_complex/k21/hifiasm.histo.txt', 'w')
i = 1
for key in sorted(histo_):
    worksheet.write("A"+str(i), str(key))
    worksheet.write("B"+str(i), str(histo_[key]))
    i += 1
    kmer_file.write(str(key)+","+str(histo_[key])+"\n")
    for val in map_[key]:
        kmer_file.write(val+"\n")
workbook.close()
kmer_file.close()
"""
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
