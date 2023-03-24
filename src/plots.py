from random import randint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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
    metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp', '# misassemblies']
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
    path = "/Users/sakshar5068/Desktop/repeat_assembler/quast_reproduced/default/"
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
            if m in ['# contigs', 'NG50', '# misassemblies']:
                new_vals = [int(j) for j in vals]
            else:
                new_vals = [float(j) for j in vals]
            modified_quast_data[key][m] += new_vals
    return modified_quast_data


def get_bar_chart_for_misassemblies(quast_data, repeat_sizes, copies, snps, depths):
    misassemblies_map = [{'N/A': 0}, {'N/A': 0}, {'N/A': 0}, {'N/A': 0}]
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for i in range(4):
                        misassembly = quast_data[id][metrics[4]][i]
                        contig = quast_data[id][metrics[0]][i]
                        if contig == 0:
                            misassemblies_map[i]['N/A'] += 1
                        elif misassembly not in misassemblies_map[i].keys():
                            misassemblies_map[i][misassembly] = 1
                        else:
                            misassemblies_map[i][misassembly] += 1
    for i in range(4):
        print(misassemblies_map[i])


def get_box_plot_for_ng50_wrt_ref_size(quast_data, repeat_sizes, copies, snps, depths):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    ng50_wrt_ref_size = [np.zeros(experiment_no), np.zeros(experiment_no), np.zeros(experiment_no),
                                  np.zeros(experiment_no)]
    i = 0
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            ref_size = 100000 + int(repeat_size) * int(copy)
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(4):
                        #gf, contigs_no = quast_data[id][metrics[2]][j], quast_data[id][metrics[0]][j]
                        #if contigs_no != 0:
                        ng50_wrt_ref_size[j][i] = quast_data[id][metrics[1]][j] / ref_size
                    i += 1
    print(i, experiment_no)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(ng50_wrt_ref_size)

    # x-axis labels
    ax.set_xticklabels(['RAmbler', 'Hifiasm', 'HiCANU', 'Verkko'])

    # Adding title
    plt.title("box plot of NG50 w.r.t. reference genome size for different assemblers")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # show plot
    plt.savefig("../figures/box_plots/ng50_wrt_ref_size/5K-20K_2-10_100-2000_20-40.png")


def get_box_plot_for_genome_fraction_per_contig(quast_data, repeat_sizes, copies, snps, depths):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    genome_fraction_per_contig = [np.zeros(experiment_no), np.zeros(experiment_no), np.zeros(experiment_no), np.zeros(experiment_no)]
    i = 0
    #for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(4):
                        gf, contigs_no = quast_data[id][metrics[2]][j], quast_data[id][metrics[0]][j]
                        if contigs_no != 0:
                            genome_fraction_per_contig[j][i] = gf / contigs_no
                    i += 1
    print(i, experiment_no)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(genome_fraction_per_contig)

    # x-axis labels
    ax.set_xticklabels(['RAmbler', 'Hifiasm', 'HiCANU', 'Verkko'])

    # Adding title
    plt.title("box plot of genome fraction per contig for different assemblers")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # show plot
    plt.savefig("../figures/box_plots/gf_per_contig/10K-20K_5-10_100-2000_30.png")


def get_sub_plots_for_ng50_vs_misassemblies(quast_data, repeat_sizes, copies, snps, depths, xlim, ylim):
    copy_no, snp_no = len(copies), len(snps)
    ng50s, misassemblies = np.zeros((copy_no, snp_no, 4)), np.zeros((copy_no, snp_no, 4))
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        i = 0
        for copy in copies:
            j = 0
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(4):
                        ng50, misassembly = quast_data[id][metrics[1]][k] / 1000, quast_data[id][metrics[4]][k]
                        ng50s[i][j][k], misassemblies[i][j][k] = ng50, misassembly
                j += 1
            i += 1
    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green', 'violet']
    labels = ['RAmbler', 'Hifiasm', 'HiCANU', 'Verkko']
    handles = []
    for i in range(4):
        handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[i], markerfacecolor='None',
                                     linestyle='None', label=labels[i]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("misassemblies vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    #fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('number of misassemblies')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=4)
    #plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(4):
                axs[i, j].plot(ng50s[i][j][k], misassemblies[i][j][k], marker=shapes[k], markeredgecolor=face_colors[k], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig("../figures/sub_plots/misassemblies.vs.ng50_" + repeat_sizes[0] + "K_5-10_100-2000_" + depths[0] + ".png")


def get_sub_plots_for_ng50_vs_gf_per_contig(quast_data, repeat_sizes, copies, snps, depths, xlim, ylim):
    copy_no, snp_no = len(copies), len(snps)
    ng50s, gf_per_contig = np.zeros((copy_no, snp_no, 4)), np.zeros((copy_no, snp_no, 4))
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        i = 0
        for copy in copies:
            j = 0
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(4):
                        ng50 = quast_data[id][metrics[1]][k] / 1000
                        gf, contigs_no = quast_data[id][metrics[2]][k], quast_data[id][metrics[0]][k]
                        if contigs_no != 0:
                            gf_per_contig[i][j][k] = gf / contigs_no
                        ng50s[i][j][k] = ng50
                j += 1
            i += 1

    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green', 'violet']
    labels = ['RAmbler', 'Hifiasm', 'HiCANU', 'Verkko']
    handles = []
    for i in range(4):
        handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[i], markerfacecolor='None',
                                     linestyle='None', label=labels[i]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("gf/contig vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('Genome fraction per contig')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=4)
    plt.tight_layout()
    #plt.rcParams["figure.autolayout"] = True
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(4):
                axs[i, j].plot(ng50s[i][j][k], gf_per_contig[i][j][k], marker=shapes[k], markeredgecolor=face_colors[k], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig("../figures/sub_plots/gf_per_contig.vs.ng50_" + repeat_sizes[0] + "K_5-10_100-2000_" + depths[0] + ".png")


def get_sub_plots_for_ng50_vs_contig_no(quast_data, repeat_sizes, copies, snps, depths, xlim, ylim):
    copy_no, snp_no = len(copies), len(snps)
    ng50s, contig_nos = np.zeros((copy_no, snp_no, 4)), np.zeros((copy_no, snp_no, 4))
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        i = 0
        for copy in copies:
            j = 0
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(4):
                        ng50 = quast_data[id][metrics[1]][k] / 1000
                        contigs_no = quast_data[id][metrics[0]][k]
                        ng50s[i][j][k], contig_nos[i][j][k] = ng50, contigs_no
                j += 1
            i += 1

    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green', 'violet']
    labels = ['RAmbler', 'Hifiasm', 'HiCANU', 'Verkko']
    handles = []
    for i in range(4):
        handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[i], markerfacecolor='None', linestyle='None', label=labels[i]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("Contig no. vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    #fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('number of contigs')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=4)
    #plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(4):
                axs[i, j].plot(ng50s[i][j][k], contig_nos[i][j][k], marker=shapes[k], markeredgecolor=face_colors[k], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig("../figures/sub_plots/contig.vs.ng50_" + repeat_sizes[0] + "K_5-10_100-2000_" + depths[0] + ".png")


def subplotter(quast_data, repeat_sizes, copies, snps, depths):
    x_lim = [200, 250, 300]
    y_lim = (0, 20)
    for r in range(1, len(repeat_sizes)+1):
        repeat = repeat_sizes[r-1:r]
        for d in range(1, len(depths)+1):
            dep = depths[d-1:d]
            get_sub_plots_for_ng50_vs_contig_no(quast_data, repeat, copies, snps, dep, x_lim[r-1], y_lim)


repeat_sizes = ["10", "15", "20"] #, "15", "20"] #["5", "10", "15", "20"]
copies = ["5", "10"] #["2", "5", "10"]
snps = ["100", "250", "500", "1000", "2000"] #["100", "250", "500", "1000", "2000"]
depths = ["30"] #["20", "30", "40"]
#experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
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
metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp', '# misassemblies']
#for metric in metrics:
#quast_plotter_depth_fixed(modified_quast_data, repeat_sizes, copies, snps, depths, metrics[3])
#for metric in metrics:
#quast_plotter_snp_fixed(modified_quast_data, repeat_sizes, copies, snps, depths, metrics[3])
#print(quast_parser("../sim5_depth50_max20000/quast_full/report.tsv"))
#snp2unikmers_map = unikmers_vs_snp_plotter("../true_unikmers.txt")
#plot_RAmbler_depth_snp_fixed(modified_quast_data, repeat_sizes, copies, sys.argv[1], sys.argv[2], metric)
#plot_RAmbler_depth_fixed(modified_quast_data, repeat_sizes, copies, snps, sys.argv[1], metrics[3])
#plot_RAmbler_snp_fixed(modified_quast_data, repeat_sizes, copies, sys.argv[1], depths, metrics[3])
get_box_plot_for_genome_fraction_per_contig(modified_quast_data, repeat_sizes, copies, snps, depths)
#get_box_plot_for_ng50_wrt_ref_size(modified_quast_data, repeat_sizes, copies, snps, depths)
#get_sub_plots_for_ng50_vs_misassemblies(modified_quast_data, repeat_sizes, copies, snps, depths)
#get_sub_plots_for_ng50_vs_gf_per_contig(modified_quast_data, repeat_sizes, copies, snps, depths)
#subplotter(modified_quast_data, repeat_sizes, copies, snps, depths)
#get_bar_chart_for_misassemblies(modified_quast_data, repeat_sizes, copies, snps, depths)