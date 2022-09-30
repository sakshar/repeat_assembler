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


path = './simulated/k21/hifi.sim5_depth50_max20000.'
workbook = xlsxwriter.Workbook(path+'xlsx')
worksheet = workbook.add_worksheet()

histo = open(path+'histo', 'r')
i = 1
for line in histo:
    args = line.strip().split()
    worksheet.write("A"+str(i), args[0])
    worksheet.write("B"+str(i), args[1])
    i += 1
workbook.close()

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
