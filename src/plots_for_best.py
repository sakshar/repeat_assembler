from random import randint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import csv
import seaborn as sns
import pandas as pd
import xlsxwriter
import sys


colors = [(round(153/255, 2), round(153/255, 2), round(0/255, 2)),
          (round(204/255, 2), round(102/255, 2), round(0/255, 2)),
          (round(0/255, 2), round(204/255, 2), round(102/255, 2)),
          (round(0/255, 2), round(102/255, 2), round(204/255, 2)),
          (round(255/255, 2), round(102/255, 2), round(178/255, 2))]


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
    metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# misassemblies']
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


def get_quast_reports(path, repeat_sizes, copies, snps, depths):
    quast_data = dict()
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    file = path + "/" + repeat_size + "_" + copy + "_" + snp + "/" + depth + "/report.tsv"
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    quast_data[id] = quast_parser(file)
    return quast_data


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


def get_box_plots(out_path, quast_data, repeat_sizes, copies, snps, depths, assemblers, method):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    gf_per_contig = np.zeros((experiment_no, len(assemblers)))
    contigs = np.zeros((experiment_no, len(assemblers)))
    misassemblies = np.zeros((experiment_no, len(assemblers)))
    ng50_wrt_ref = np.zeros((experiment_no, len(assemblers)))
    effective_gf_per_contig = np.zeros((experiment_no, len(assemblers)))
    ComAcCon = np.zeros((experiment_no, len(assemblers)))
    i = 0
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            ref_size = 100000 + int(repeat_size) * int(copy)
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(len(assemblers)):
                        gf, contigs_no = quast_data[id][metrics[2]][k], quast_data[id][metrics[0]][k]
                        ng50, misassembly = quast_data[id][metrics[1]][k], quast_data[id][metrics[3]][k]
                        ng50_wrt_ref[i][k] = ng50 / ref_size
                        misassemblies[i][k] = misassembly
                        if contigs_no != 0:
                            contigs[i][k] = contigs_no
                            gf_per_contig[i][k] = (gf / 100.0) / contigs_no
                            effective_gf_per_contig[i][k] = (gf / 100.0) / (contigs_no + misassembly)
                            ComAcCon[i][k] = 2 * effective_gf_per_contig[i][k] * ng50_wrt_ref[i][k] / (effective_gf_per_contig[i][k] + ng50_wrt_ref[i][k])
                    i += 1
    print(i, experiment_no)
    plot_tag = box_plot[method]

    sns.set_theme()

    plt.figure(figsize=(10, 7))
    #ax = fig.add_subplot(111)
    plt.title("box plot of " + plot_tag + " for different assemblers")
    # need to convert the numpy array into a dataframe
    #           thresholds misassemblies
    # 0             5           1
    # 1             10          0
    df = pd.DataFrame()
    df_assemblers, df_vals = [], []
    for j in range(experiment_no):
        for k in range(len(assemblers)):
            df_assemblers.append(assemblers[k])
            if method == 0:
                df_vals.append(ComAcCon[j][k])
            elif method == 1:
                df_vals.append(gf_per_contig[j][k])
            elif method == 2:
                df_vals.append(ng50_wrt_ref[j][k])
    df['assemblers'], df[plot_tag] = df_assemblers, df_vals

    sns.boxplot(df, x='assemblers', y=plot_tag)
    sns.color_palette("deep")

    # show plot
    plt.savefig(out_path+plot_tag+"_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_violin_plots(out_path, quast_data, repeat_sizes, copies, snps, depths, assemblers, method):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    contigs = np.zeros((experiment_no, len(assemblers)))
    misassemblies = np.zeros((experiment_no, len(assemblers)))
    i = 0
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            ref_size = 100000 + int(repeat_size) * int(copy)
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(len(assemblers)):
                        gf, contigs_no = quast_data[id][metrics[2]][k], quast_data[id][metrics[0]][k]
                        ng50, misassembly = quast_data[id][metrics[1]][k], quast_data[id][metrics[3]][k]
                        misassemblies[i][k] = misassembly
                        contigs[i][k] = contigs_no
                    i += 1
    print(i, experiment_no)
    plot_tag = violin_plot[method]

    sns.set_theme()

    plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111)
    plt.title("violin plot of " + plot_tag + " for different assemblers")
    # need to convert the numpy array into a dataframe
    #           thresholds misassemblies
    # 0             5           1
    # 1             10          0
    df = pd.DataFrame()
    df_assemblers, df_vals = [], []
    for j in range(experiment_no):
        for k in range(len(assemblers)):
            df_assemblers.append(assemblers[k])
            if method == 0:
                df_vals.append(contigs[j][k])
            elif method == 1:
                df_vals.append(misassemblies[j][k])

    df['assemblers'], df[plot_tag] = df_assemblers, df_vals

    sns.violinplot(df, x='assemblers', y=plot_tag)
    sns.color_palette("deep")
    # show plot
    plt.savefig(out_path+plot_tag+"_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


# need to modify the sub-plot methods
def get_sub_plots(out_path, quast_data, repeat_sizes, copies, snps, depths, xlim, ylim, assemblers, method):
    copy_no, snp_no = len(copies), len(snps)
    comparison_no = len(assemblers)
    ng50s, misassemblies = np.zeros((copy_no, snp_no, comparison_no)), np.zeros((copy_no, snp_no, comparison_no))
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        i = 0
        for copy in copies:
            j = 0
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for k in range(comparison_no):
                        ng50, misassembly = quast_data[id][metrics[1]][k] / 1000, quast_data[id][metrics[3]][k]
                        ng50s[i][j][k], misassemblies[i][j][k] = ng50, misassembly
                j += 1
            i += 1
    """
    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green']
    handles = []
    for i in range(len(assemblers)):
        handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[j], markerfacecolor='None',
                                     linestyle='None', label=xticklabels[i*len(thresholds)+j]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("misassemblies vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    #fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('number of misassemblies')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    #plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(comparison_no):
                axs[i, j].plot(ng50s[i][j][k], misassemblies[i][j][k], marker=shapes[k % len(thresholds)], markeredgecolor=face_colors[k // len(thresholds)], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(out_path + "sub_plots_all/misassemblies.vs.ng50_" + repeat_sizes[0] + "K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1] + "_" + depths[0] + ".png")
    """


def subplotter(out_path, quast_data, repeat_sizes, copies, snps, depths, assemblers, sub_plotter_method):
    x_lim = [200, 250, 300]
    for r in range(1, len(repeat_sizes)+1):
        repeat = repeat_sizes[r-1:r]
        for d in range(1, len(depths)+1):
            dep = depths[d-1:d]
            y_lim = sub_plot_ylimits[sub_plot[sub_plotter_method]]
            get_sub_plots(out_path, quast_data, repeat, copies, snps, dep, x_lim[r-1], y_lim, assemblers, sub_plotter_method)


assemblers = ["RAmbler", "Hifiasm", "HiCANU", "Verkko"]
repeat_sizes = ["15", "20"] #, "15", "20"] #["5", "10", "15", "20"]
copies = ["5", "10"] #["2", "5", "10"]
snps = ["250", "500", "1000"] #["100", "250", "500", "1000", "2000"]
depths = ["20", "30", "40"] #["20", "30", "40"]
sub_plot = {0: "ng50_vs_contig_no", 1: "ng50_vs_misassemblies", 2: "ng50_vs_gf_per_contig"}
box_plot = {0: "com-ac-con", 1: "gf_per_contig", 2: "ng50_wrt_ref"}
violin_plot = {0: "contigs", 1: "misassemblies"}
sub_plot_ylimits = {sub_plot[0]: (0, 7),
                   sub_plot[1]: (-1, 5),
                   sub_plot[2]: (0, 110)}
#sub_plotter_method = 0
#box_plot_method = 2
#violin_plot_method = 1
best_hyperparameters = "15_15"
input_path = "/Users/sakshar5068/Desktop/repeat_assembler/best/" + best_hyperparameters + "/quast/"
figure_path = "/Users/sakshar5068/Desktop/repeat_assembler/best/" + best_hyperparameters + "/figures/"
quast_data = get_quast_reports(input_path, repeat_sizes, copies, snps, depths)

modified_quast_data = preprocess_quast_data(quast_data)
metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# misassemblies']
for box in range(3):
    get_box_plots(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, assemblers, box)
for violin in range(2):
    get_violin_plots(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, assemblers, violin)
#subplotter(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds, sub_plotter_method)