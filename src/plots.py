from random import randint
import numpy as np
import matplotlib.pyplot as plt
import csv
import xlsxwriter
import sys


colors = [(round(153/255, 2), round(153/255, 2), round(0/255, 2)),
          (round(204/255, 2), round(102/255, 2), round(0/255, 2)),
          (round(0/255, 2), round(204/255, 2), round(102/255, 2)),
          (round(0/255, 2), round(102/255, 2), round(204/255, 2)),
          (round(255/255, 2), round(102/255, 2), round(178/255, 2))]


def nucFreq_generator():
    x = []
    y = []

    for i in range(0, 50000, 500):
        x.append(i)
        y.append(randint(38, 42))

    for i in range(-10, 11):
        x.append(i*500 + 55000)
        val = ((-1) * (i**4) + 10000) * 160 / 10000 + randint(38, 42)
        y.append(val)

    for i in range(60500, 110001, 500):
        x.append(i)
        y.append(randint(38, 42))

    plt.plot(x, y)
    plt.xlabel("Assembly position (bp)")
    plt.ylabel("Sequence read depth")
    plt.savefig("../figures/nucFreq.png")


def kmer_distribution_generator():
    mu, sigma = 60, 10
    s = np.random.normal(mu, sigma, 1000)
    num_bins = 100

    n, bins, patches = plt.hist(s, num_bins,
                                density=1,
                                color='blue',
                                alpha=0.7)

    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2))

    plt.plot(bins, y, '--', color='red')
    #plt.xlabel('# of occurrences')
    #plt.ylabel('# of k-mers')
    plt.savefig('../figures/kmer_distribution.png')


def quast_parser(infile):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    filtered_data = dict()
    #print(data[0][1:])
    missing = False
    insert_at = -1
    if "canu.contigs" not in data[0][1:]:
        missing, insert_at = True, 3
    elif "verkko" not in data[0][1:]:
        missing, insert_at = True, 4
    metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp']
    for row in data:
        if row[0] in metrics:
            if not missing:
                filtered_data[row[0]] = row[1:]
            else:
                if insert_at == 4:
                    filtered_data[row[0]] = row[1:] + ["-"]
                elif insert_at == 3:
                    filtered_data[row[0]] = row[1:3] + ["-"] + row[3:]
    return filtered_data


def get_quast_reports(repeat_sizes, copies, snps, depths):
    quast_data = dict()
    path = "/Users/sakshar5068/Desktop/repeat_assembler/quast2/"
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    file = path + repeat_size + "_" + copy + "_" + snp + "/" + depth + "/report.tsv"
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    #print(id)
                    quast_data[id] = quast_parser(file)
    return quast_data


def plot_function_depth_fixed(quast_data, repeat_size, copies, snps, depth, metric):
    groups = list()
    RAmbler, Hifiasm, HiCANU, Verkko = list(), list(), list(), list()
    for copy in copies:
        for snp in snps:
            groups.append(copy + "-" + snp)
            id = repeat_size + "_" + copy + "_" + snp + "_" + depth
            RAmbler.append(quast_data[id][metric][0])
            #RAmblerWithoutOLC.append(quast_data[id][metric][1])
            Hifiasm.append(quast_data[id][metric][1])
            HiCANU.append(quast_data[id][metric][2])
            Verkko.append(quast_data[id][metric][3])
    N = len(copies) * len(snps)
    filename = "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/fixed_repeat_size_n_depth/" + metric + ".vs.Copies-SNP_rates_RepeatSize_" + repeat_size + "_Read_coverage_" + depth
    workbook = xlsxwriter.Workbook(filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    for i in range(N):
        worksheet.write("A" + str(i+1), groups[i])
        worksheet.write("B" + str(i+1), RAmbler[i])
        worksheet.write("C" + str(i+1), Hifiasm[i])
        worksheet.write("D" + str(i+1), HiCANU[i])
        worksheet.write("E" + str(i+1), Verkko[i])
    workbook.close()
    """
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, RAmbler, width, color=colors[0])
    #bar2 = plt.bar(ind + width, RAmblerWithoutOLC, width, color=colors[1])
    bar2 = plt.bar(ind + width, Hifiasm, width, color=colors[1])
    bar3 = plt.bar(ind + width * 2, HiCANU, width, color=colors[2])
    bar4 = plt.bar(ind + width * 3, Verkko, width, color=colors[3])
    xlabel = "Copies-SNP_rates for repeat size " + repeat_size + " and read coverage " + depth + "x"
    plt.xlabel(xlabel)
    plt.ylabel(metric)
    # plt.ylim([98000, 107000])
    plt.xticks(ind + width * 1.5, groups, rotation=45)
    plt.legend((bar1, bar2, bar3, bar4),
                   ('RAmbler', 'Hifiasm', 'HiCANU', 'Verkko'), loc='upper center',
                   bbox_to_anchor=(0.5, 1.2), ncol=4)
    plt.tight_layout()
    plt.savefig(
            "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/fixed_repeat_size_n_depth/" + metric + ".vs.Copies-SNP_rates_RepeatSize_" + repeat_size + "_Read_coverage_" + depth + ".png")
    """


def quast_plotter_depth_fixed(quast_data, repeat_sizes, copies, snps, depths, metric):
    for repeat_in_k in repeat_sizes:
        repeat_size = repeat_in_k + "000"
        for depth in depths:
            plot_function_depth_fixed(quast_data, repeat_size, copies, snps, depth, metric)


def plot_function_snp_fixed(quast_data, repeat_size, copies, snp, depths, metric):
    groups = list()
    RAmbler, Hifiasm, HiCANU, Verkko = list(), list(), list(), list()
    for copy in copies:
        for depth in depths:
            groups.append(copy + "-" + depth)
            id = repeat_size + "_" + copy + "_" + snp + "_" + depth
            RAmbler.append(quast_data[id][metric][0])
            #RAmblerWithoutOLC.append(quast_data[id][metric][1])
            Hifiasm.append(quast_data[id][metric][1])
            HiCANU.append(quast_data[id][metric][2])
            Verkko.append(quast_data[id][metric][3])
    N = len(copies) * len(depths)
    filename = "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/fixed_repeat_size_n_snp_rate/" + metric + ".vs.Copies-Read_Coverages_RepeatSize_" + repeat_size + "_SNP_" + snp
    workbook = xlsxwriter.Workbook(filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    for i in range(N):
        worksheet.write("A" + str(i+1), groups[i])
        worksheet.write("B" + str(i+1), RAmbler[i])
        worksheet.write("C" + str(i+1), Hifiasm[i])
        worksheet.write("D" + str(i+1), HiCANU[i])
        worksheet.write("E" + str(i+1), Verkko[i])
    workbook.close()
    """
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, RAmbler, width, color=colors[0])
    #bar2 = plt.bar(ind + width, RAmblerWithoutOLC, width, color=colors[1])
    bar2 = plt.bar(ind + width, Hifiasm, width, color=colors[1])
    bar3 = plt.bar(ind + width * 2, HiCANU, width, color=colors[2])
    bar4 = plt.bar(ind + width * 3, Verkko, width, color=colors[3])
    xlabel = "Copies-Read_Coverages for repeat size " + repeat_size + " with mutations per " + snp + " bp"
    plt.xlabel(xlabel)
    plt.ylabel(metric)
    # plt.ylim([98000, 107000])
    plt.xticks(ind + width * 1.5, groups, rotation=45)
    plt.legend((bar1, bar2, bar3, bar4),
                   ('RAmbler', 'Hifiasm', 'HiCANU', 'Verkko'), loc='upper center',
                   bbox_to_anchor=(0.5, 1.2), ncol=4)
    plt.tight_layout()
    plt.savefig(
            "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/fixed_repeat_size_n_snp_rate/" + metric + ".vs.Copies-Read_Coverages_RepeatSize_" + repeat_size + "_SNP_" + snp + ".png")
    """


def quast_plotter_snp_fixed(quast_data, repeat_sizes, copies, snps, depths, metric):
    for repeat_in_k in repeat_sizes:
        repeat_size = repeat_in_k + "000"
        for snp in snps:
            plot_function_snp_fixed(quast_data, repeat_size, copies, snp, depths, metric)


def plot_RAmbler_depth_fixed(quast_data, repeat_sizes, copies, snps, depth, metric):
    groups = list()
    metric_vals = [[], [], [], [], []]
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            groups.append(repeat_size + "-" + copy)
            i = 0
            for snp in snps:
                id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                metric_vals[i].append(quast_data[id][metric][0])
                i += 1
    N = len(repeat_sizes) * len(copies)
    filename = "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/RAmbler_fixed_depth/" + metric + ".vs.Repeat_sizes-Copies for read coverage " + depth
    workbook = xlsxwriter.Workbook(filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    for i in range(N):
        worksheet.write("A" + str(i+1), groups[i])
        worksheet.write("B" + str(i+1), metric_vals[0][i])
        worksheet.write("C" + str(i+1), metric_vals[1][i])
        worksheet.write("D" + str(i+1), metric_vals[2][i])
        worksheet.write("E" + str(i+1), metric_vals[3][i])
        worksheet.write("F" + str(i+1), metric_vals[4][i])
    workbook.close()
    """
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, metric_vals[0], width, color=colors[0])
    bar2 = plt.bar(ind + width, metric_vals[1], width, color=colors[1])
    bar3 = plt.bar(ind + width * 2, metric_vals[2], width, color=colors[2])
    bar4 = plt.bar(ind + width * 3, metric_vals[3], width, color=colors[3])
    bar5 = plt.bar(ind + width * 4, metric_vals[4], width, color=colors[4])
    xlabel = "Repeat_sizes-Copies for read coverage " + depth + "x"
    plt.xlabel(xlabel)
    plt.ylabel(metric)
    # plt.ylim([98000, 107000])
    plt.xticks(ind + width * 2, groups, rotation=45)
    plt.legend((bar1, bar2, bar3, bar4, bar5),
                   ('per 100 bp', 'per 250 bp', 'per 500 bp', 'per 1000 bp', 'per 2000 bp'), loc='upper center',
                   bbox_to_anchor=(0.5, 1.2), ncol=3)
    plt.tight_layout()
    plt.savefig(
            "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/RAmbler_fixed_depth/" + metric + ".vs.Repeat_sizes-Copies for read coverage " + depth + ".png")
    """


def plot_RAmbler_snp_fixed(quast_data, repeat_sizes, copies, snp, depths, metric):
    groups = list()
    metric_vals = [[], [], []]
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            groups.append(repeat_size + "-" + copy)
            i = 0
            for depth in depths:
                id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                metric_vals[i].append(quast_data[id][metric][0])
                i += 1
    N = len(repeat_sizes) * len(copies)
    filename = "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/RAmbler_fixed_snp/" + metric + ".vs.Repeat_sizes-Copies for SNP "+ snp
    workbook = xlsxwriter.Workbook(filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    for i in range(N):
        worksheet.write("A" + str(i+1), groups[i])
        worksheet.write("B" + str(i+1), metric_vals[0][i])
        worksheet.write("C" + str(i+1), metric_vals[1][i])
        worksheet.write("D" + str(i+1), metric_vals[2][i])
    workbook.close()
    """
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, metric_vals[0], width, color=colors[1])
    bar2 = plt.bar(ind + width, metric_vals[1], width, color=colors[2])
    bar3 = plt.bar(ind + width * 2, metric_vals[2], width, color=colors[3])
    xlabel = "Repeat_sizes-Copies for mutations per " + snp + " bp"
    plt.xlabel(xlabel)
    plt.ylabel(metric)
    # plt.ylim([98000, 107000])
    plt.xticks(ind + width, groups, rotation=45)
    plt.legend((bar1, bar2, bar3),
                   ('20x', '30x', '40x'), loc='upper center',
                   bbox_to_anchor=(0.5, 1.2), ncol=3)
    plt.tight_layout()
    plt.savefig(
            "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmbler_v2/RAmbler_fixed_snp/" + metric + ".vs.Repeat_sizes-Copies for SNP "+ snp + ".png")
    """


def plot_RAmbler_depth_snp_fixed(quast_data, repeat_sizes, copies, snp, depth, metric):
    groups = list()
    metric_vals = [[], [], [], []]
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            groups.append(repeat_size + "-" + copy)
            i = 0
            for shared_unikmer in [5, 10, 15, 20]:
                id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                metric_vals[i].append(quast_data[id][metric][i])
                i += 1
    N = len(repeat_sizes) * len(copies)
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, metric_vals[0], width, color=colors[0])
    bar2 = plt.bar(ind + width, metric_vals[1], width, color=colors[1])
    bar3 = plt.bar(ind + width * 2, metric_vals[2], width, color=colors[2])
    bar4 = plt.bar(ind + width * 3, metric_vals[3], width, color=colors[3])
    xlabel = "Repeat_sizes-Copies for read coverage " + depth + "x and mutations per " + snp + " bp"
    plt.xlabel(xlabel)
    plt.ylabel(metric)
    # plt.ylim([98000, 107000])
    plt.xticks(ind + width * 1.5, groups, rotation=45)
    plt.legend((bar1, bar2, bar3, bar4),
                   ('5', '10', '15', '20'), loc='upper center',
                   bbox_to_anchor=(0.5, 1.2), ncol=4)
    plt.tight_layout()
    plt.savefig(
            "/Users/sakshar5068/Desktop/repeat_assembler/figures/RAmblerWithoutOLC_shared_unikmers/" + metric + ".vs.Repeat_sizes-Copies for read coverage " + depth + " and mutations per " + snp + " bp.png")


def unikmers_vs_snp_plotter(file):
    kmer_file = open(file, "r")
    snp2unikmers_map = dict()
    for i in range(125):
        id = kmer_file.readline().strip().split('_')
        count = int(kmer_file.readline().strip().split(' ')[0])
        if (id[0], id[1]) not in snp2unikmers_map:
            snp2unikmers_map[(id[0], id[1])] = dict()
        snp2unikmers_map[(id[0], id[1])][int(id[2])] = count
    kmer_file.close()
    snps = [[], [], [], [], []]
    groups = []
    for repeat_size in ["5000", "10000", "15000", "20000"]:
        for copies in ["2", "5", "10", "20"]:
            temp_repeat_size = str(int(repeat_size)//1000)+"k"
            groups.append(temp_repeat_size+"-"+copies)
            i = 0
            for snp in [100, 250, 500, 1000, 2000]:
                snps[i].append(snp2unikmers_map[(repeat_size, copies)][snp])
                i += 1
    N = 16
    ind = np.arange(N)
    width = 0.15
    bar1 = plt.bar(ind, snps[0], width, color=colors[0])
    bar2 = plt.bar(ind+width, snps[1], width, color=colors[1])
    bar3 = plt.bar(ind+width*2, snps[2], width, color=colors[2])
    bar4 = plt.bar(ind+width*3, snps[3], width, color=colors[3])
    bar5 = plt.bar(ind+width*4, snps[4], width, color=colors[4])
    plt.xlabel("Repeat_size-Copies")
    plt.ylabel("# of true unikmers")
    plt.ylim([98000, 107000])
    plt.xticks(ind + width*2, groups, rotation=45)
    plt.legend((bar1, bar2, bar3, bar4, bar5), ('per 100 bp', 'per 250 bp', 'per 500 bp', 'per 1000 bp', 'per 2000 bp'))
    plt.tight_layout()
    plt.savefig("../figures/unikmers_vs_snp.png")
    return snp2unikmers_map


def preprocess_quast_data(quast_data):
    modified_quast_data = dict()
    for key in quast_data.keys():
        # print(key)
        modified_quast_data[key] = dict()
        for m in quast_data[key].keys():
            modified_quast_data[key][m] = []
            vals = quast_data[key][m]
            for i in range(len(vals)):
                if vals[i] == '-':
                    vals[i] = '0'
            #if len(vals) == 3:
            #    vals.append('0')
            if m in ['# contigs', 'NG50']:
                new_vals = [int(j) for j in vals]
            else:
                new_vals = [float(j) for j in vals]
            modified_quast_data[key][m] += new_vals
    return modified_quast_data


repeat_sizes = ["5", "10", "15", "20"]
copies = ["2", "5", "10"]
snps = ["100", "250", "500", "1000", "2000"]
depths = ["20", "30", "40"]
#metrics = ['# contigs'] #'NG50'] #'# contigs'] #, 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp']
quast_data = get_quast_reports(repeat_sizes, copies, snps, depths)
#print(quast_data["20000_2_1000_20"])
modified_quast_data = preprocess_quast_data(quast_data)
#print(modified_quast_data["20000_2_1000_20"])
#for id in modified_quast_data.keys():
#    print(id, modified_quast_data[id])
#for key in modified_quast_data.keys():
#    print(key)
#    for m in modified_quast_data[key].keys():
#        print(m, modified_quast_data[key][m])
metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp']
#for metric in metrics:
#quast_plotter_depth_fixed(modified_quast_data, repeat_sizes, copies, snps, depths, metrics[3])
#for metric in metrics:
quast_plotter_snp_fixed(modified_quast_data, repeat_sizes, copies, snps, depths, metrics[3])
#print(quast_parser("../sim5_depth50_max20000/quast_full/report.tsv"))
#snp2unikmers_map = unikmers_vs_snp_plotter("../true_unikmers.txt")
#plot_RAmbler_depth_snp_fixed(modified_quast_data, repeat_sizes, copies, sys.argv[1], sys.argv[2], metric)
#plot_RAmbler_depth_fixed(modified_quast_data, repeat_sizes, copies, snps, sys.argv[1], metrics[3])
#plot_RAmbler_snp_fixed(modified_quast_data, repeat_sizes, copies, sys.argv[1], depths, metrics[3])