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


def quast_parser_for_hyperparameters(infile):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    filtered_data = dict()
    #print(data[0][1:])
    metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# misassemblies']
    for row in data:
        if row[0] in metrics:
            filtered_data[row[0]] = row[1:]
    return filtered_data


def get_quast_reports(path, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    quast_data = dict()
    for tolerance in tolerances:
        quast_data[tolerance] = dict()
        for repeat_size in repeat_sizes:
            repeat_size += "000"
            for copy in copies:
                for snp in snps:
                    for depth in depths:
                        file = path + tolerance + "/" + repeat_size + "_" + copy + "_" + snp + "/" + depth + "/report.tsv"
                        id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                        quast_data[tolerance][id] = quast_parser_for_hyperparameters(file)
    return quast_data


def preprocess_quast_data(quast_data):
    modified_quast_data = dict()
    for tolerance in quast_data.keys():
        modified_quast_data[tolerance] = dict()
        for key in quast_data[tolerance].keys():
            # print(key)
            modified_quast_data[tolerance][key] = dict()
            for m in quast_data[tolerance][key].keys():
                modified_quast_data[tolerance][key][m] = []
                vals = quast_data[tolerance][key][m]
                for i in range(len(vals)):
                    if vals[i] == '-':
                        vals[i] = '0'
                #if len(vals) == 3:
                #    vals.append('0')
                if m in ['# contigs', 'NG50', '# misassemblies']:
                    new_vals = [int(j) for j in vals]
                else:
                    new_vals = [float(j) for j in vals]
                modified_quast_data[tolerance][key][m] += new_vals
    return modified_quast_data


def get_box_plot_for_ComAcCon_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    #comparison_no = len(tolerances) * len(thresholds)
    #xticklabels = []
    #for to in tolerances:
        #for th in thresholds:
            #xticklabels.append(to+"_"+th)
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    effective_gf_per_contig = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    ng50_wrt_ref = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    i = 0
    #for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            ref_size = 100000 + int(repeat_size) * int(copy)
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            gf, contigs_no = quast_data[tolerances[j]][id][metrics[2]][k], quast_data[tolerances[j]][id][metrics[0]][k]
                            ng50, misassembly = quast_data[tolerances[j]][id][metrics[1]][k], quast_data[tolerances[j]][id][metrics[3]][k]
                            ng50_wrt_ref[j][i][k] = ng50 / ref_size
                            if contigs_no != 0:
                                effective_gf_per_contig[j][i][k] = (gf / 100.0) / (contigs_no + misassembly)
                    i += 1
    ComAcCon = 2 * effective_gf_per_contig * ng50_wrt_ref / (effective_gf_per_contig + ng50_wrt_ref)
    print(i, experiment_no)

    sns.set_theme()

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    #fig.suptitle("box plot of Com-Ac-Con for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.rcParams.update({'font.size': 22})
    plt.xlabel('thresholds')
    plt.ylabel('assembly-score')
    #plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        # need to convert the numpy array into a dataframe
        #           thresholds misassemblies
        # 0             5           1
        # 1             10          0
        df = pd.DataFrame()
        df_thresholds, df_comaccons = [], []
        for j in range(experiment_no):
            for k in range(len(thresholds)):
                df_thresholds.append(int(thresholds[k]))
                df_comaccons.append(ComAcCon[i][j][k])
        df['thresholds'], df['assembly-score'] = df_thresholds, df_comaccons

        sns.boxplot(df, ax=axs[i], x='thresholds', y='assembly-score')
        sns.color_palette("deep")
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    """
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(genome_fraction_per_contig)

    # x-axis labels
    ax.set_xticklabels(xticklabels, rotation=90)

    # Adding title
    plt.title("box plot of genome fraction per contig for different hyperparameters (to_th)")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    """
    # show plot
    plt.savefig(out_path + "ComAcCon_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_box_plot_for_genome_fraction_per_contig_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    #comparison_no = len(tolerances) * len(thresholds)
    #xticklabels = []
    #for to in tolerances:
        #for th in thresholds:
            #xticklabels.append(to+"_"+th)
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    genome_fraction_per_contig = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    i = 0
    #for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            gf, contigs_no = quast_data[tolerances[j]][id][metrics[2]][k], quast_data[tolerances[j]][id][metrics[0]][k]
                            misassembly = quast_data[tolerances[j]][id][metrics[3]][k]
                            if contigs_no != 0:
                                genome_fraction_per_contig[j][i][k] = (gf / 100.0) / (contigs_no + misassembly)
                    i += 1
    print(i, experiment_no)

    sns.set_theme()

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    #fig.suptitle("box plot of genome fraction per contig for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.rcParams.update({'font.size': 22})
    plt.xlabel('thresholds')
    plt.ylabel("effective genome fraction per contig (" + r'$\zeta$)')
    #plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        # need to convert the numpy array into a dataframe
        #           thresholds misassemblies
        # 0             5           1
        # 1             10          0
        df = pd.DataFrame()
        df_thresholds, df_gfs = [], []
        for j in range(experiment_no):
            for k in range(len(thresholds)):
                df_thresholds.append(int(thresholds[k]))
                df_gfs.append(genome_fraction_per_contig[i][j][k])
        df['thresholds'], df['gf_per_contig'] = df_thresholds, df_gfs

        sns.boxplot(df, ax=axs[i], x='thresholds', y='gf_per_contig')
        sns.color_palette("deep")
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    """
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(genome_fraction_per_contig)

    # x-axis labels
    ax.set_xticklabels(xticklabels, rotation=90)

    # Adding title
    plt.title("box plot of genome fraction per contig for different hyperparameters (to_th)")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    """
    # show plot
    plt.savefig(out_path + "eff_gf_per_contig_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_box_plot_for_ng50_wrt_ref_size_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    #comparison_no = len(tolerances) * len(thresholds)
    #xticklabels = []
    #for to in tolerances:
        #for th in thresholds:
            #xticklabels.append(to+"_"+th)
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    ng50_wrt_ref_size = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    i = 0
    # for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            ref_size = 100000 + int(repeat_size) * int(copy)
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            ng50_wrt_ref_size[j][i][k] = quast_data[tolerances[j]][id][metrics[1]][k] / ref_size
                    i += 1
    print(i, experiment_no)


    sns.set_theme()

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    #fig.suptitle("box plot of NG50 w.r.t. reference genome size for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.rcParams.update({'font.size': 22})
    plt.xlabel('thresholds')
    plt.ylabel("NG50 w.r.t reference genome (" + r'$\eta$)')
    # plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        # need to convert the numpy array into a dataframe
        #           thresholds misassemblies
        # 0             5           1
        # 1             10          0
        df = pd.DataFrame()
        df_thresholds, df_ng50s = [], []
        for j in range(experiment_no):
            for k in range(len(thresholds)):
                df_thresholds.append(int(thresholds[k]))
                df_ng50s.append(ng50_wrt_ref_size[i][j][k])
        df['thresholds'], df['NG50_wrt_ref'] = df_thresholds, df_ng50s

        sns.boxplot(df, ax=axs[i], x='thresholds', y='NG50_wrt_ref')
        sns.color_palette("deep")
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    """
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(ng50_wrt_ref_size)

    # x-axis labels
    ax.set_xticklabels(xticklabels, rotation=90)

    # Adding title
    plt.title("box plot of NG50 w.r.t. reference genome size for different hyperparameters (to_th)")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    """
    # show plot
    plt.savefig(out_path + "ng50_wrt_ref_size_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_box_plot_for_misassemblies_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    #comparison_no = len(tolerances) * len(thresholds)
    #xticklabels = []
    #for to in tolerances:
        #for th in thresholds:
            #xticklabels.append(to+"_"+th)
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    misassemblies = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    #misassemblies_freq = np.zeros((5, comparison_no))
    i = 0
    #for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            misassembly = quast_data[tolerances[j]][id][metrics[3]][k]
                            misassemblies[j][i][k] = misassembly
                            #misassemblies_freq[misassembly][j * len(thresholds) + k] += 1
                    i += 1
    print(i, experiment_no)

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    fig.suptitle(
        "box plot of misassemblies for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('thresholds')
    plt.ylabel('# of misassemblies')
    # plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        axs[i].boxplot(misassemblies[i])
        axs[i].set_xticklabels(thresholds)
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    """
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(misassemblies)

    # x-axis labels
    ax.set_xticklabels(xticklabels, rotation=90)

    # Adding title
    plt.title("box plot of misassemblies for different hyperparameters (to_th)")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    """
    # show plot
    #np.savetxt(out_path + "box_plots_all/misassemblies_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".csv", misassemblies_freq.transpose(), delimiter=",")
    plt.savefig(out_path + "box_plots_all/misassemblies_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_box_plot_for_contigs_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    #comparison_no = len(tolerances) * len(thresholds)
    #xticklabels = []
    #for to in tolerances:
        #for th in thresholds:
            #xticklabels.append(to+"_"+th)
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    contigs = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    #misassemblies_freq = np.zeros((5, comparison_no))
    i = 0
    #for i in range(experiment_no):
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            contig = quast_data[tolerances[j]][id][metrics[0]][k]
                            contigs[j][i][k] = contig
                            #misassemblies_freq[misassembly][j * len(thresholds) + k] += 1
                    i += 1
    print(i, experiment_no)

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    fig.suptitle(
        "box plot of contigs for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('thresholds')
    plt.ylabel('# of contigs')
    # plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        axs[i].boxplot(contigs[i])
        axs[i].set_xticklabels(thresholds)
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    """
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # Creating axes instance
    bp = ax.boxplot(misassemblies)

    # x-axis labels
    ax.set_xticklabels(xticklabels, rotation=90)

    # Adding title
    plt.title("box plot of misassemblies for different hyperparameters (to_th)")

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    """
    # show plot
    #np.savetxt(out_path + "box_plots_all/misassemblies_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".csv", misassemblies_freq.transpose(), delimiter=",")
    plt.savefig(out_path + "box_plots_all/contigs_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_violin_plots_for_misassemblies_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    misassemblies = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    i = 0
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            misassembly = quast_data[tolerances[j]][id][metrics[3]][k]
                            misassemblies[j][i][k] = misassembly
                    i += 1
    print(i, experiment_no)

    sns.set_theme()

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    #fig.suptitle("boxen plot of misassemblies for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.rcParams.update({'font.size': 22})
    plt.xlabel('thresholds')
    plt.ylabel('# of mis-assemblies')
    # plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        # need to convert the numpy array into a dataframe
        #           thresholds misassemblies
        # 0             5           1
        # 1             10          0
        df = pd.DataFrame()
        df_thresholds, df_misassemblies = [], []
        for j in range(experiment_no):
            for k in range(len(thresholds)):
                df_thresholds.append(int(thresholds[k]))
                df_misassemblies.append(misassemblies[i][j][k])
        df['thresholds'], df['misassemblies'] = df_thresholds, df_misassemblies

        sns.boxenplot(df, ax=axs[i], x='thresholds', y='misassemblies')
        sns.color_palette("deep")
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.savefig(out_path + "misassemblies_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")


def get_violin_plots_for_contigs_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds):
    experiment_no = len(repeat_sizes) * len(copies) * len(snps) * len(depths)
    contigs = np.zeros((len(tolerances), experiment_no, len(thresholds)))
    i = 0
    for repeat_size in repeat_sizes:
        repeat_size += "000"
        for copy in copies:
            for snp in snps:
                for depth in depths:
                    id = repeat_size + "_" + copy + "_" + snp + "_" + depth
                    for j in range(len(tolerances)):
                        for k in range(len(thresholds)):
                            contig = quast_data[tolerances[j]][id][metrics[0]][k]
                            contigs[j][i][k] = contig
                    i += 1
    print(i, experiment_no)

    sns.set_theme()

    fig, axs = plt.subplots(len(tolerances), figsize=(20, 15))
    #fig.suptitle("boxen plot of contigs for different hyperparameters (to_th)")
    # fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.rcParams.update({'font.size': 22})
    plt.xlabel('thresholds')
    plt.ylabel('# of contigs')
    # plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    # plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(len(tolerances)):
        # need to convert the numpy array into a dataframe
        #           thresholds misassemblies
        # 0             5           1
        # 1             10          0
        df = pd.DataFrame()
        df_thresholds, df_contigs = [], []
        for j in range(experiment_no):
            for k in range(len(thresholds)):
                df_thresholds.append(int(thresholds[k]))
                df_contigs.append(contigs[i][j][k])
        df['thresholds'], df['contigs'] = df_thresholds, df_contigs

        sns.boxenplot(df, ax=axs[i], x='thresholds', y='contigs')
        sns.color_palette("deep")
        axs[i].set_ylabel(tolerances[i])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.savefig(out_path + "contigs_"+repeat_sizes[0]+"K-"+repeat_sizes[-1]+"K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1]+"_"+depths[0]+"-"+depths[-1]+".png")



# need to modify the sub-plot methods
def get_sub_plots_for_ng50_vs_misassemblies_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, xlim, ylim, tolerances, thresholds):
    copy_no, snp_no = len(copies), len(snps)
    comparison_no = len(tolerances) * len(thresholds)
    xticklabels = []
    for to in tolerances:
        for th in thresholds:
            xticklabels.append(to + "_" + th)
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
    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green']
    handles = []
    for i in range(len(tolerances)):
        for j in range(len(thresholds)):
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


def get_sub_plots_for_ng50_vs_gf_per_contig_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, xlim, ylim, tolerances, thresholds):
    copy_no, snp_no = len(copies), len(snps)
    comparison_no = len(tolerances) * len(thresholds)
    xticklabels = []
    for to in tolerances:
        for th in thresholds:
            xticklabels.append(to + "_" + th)
    ng50s, gf_per_contig = np.zeros((copy_no, snp_no, comparison_no)), np.zeros((copy_no, snp_no, comparison_no))
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
                        ng50 = quast_data[id][metrics[1]][k] / 1000
                        gf, contigs_no = quast_data[id][metrics[2]][k], quast_data[id][metrics[0]][k]
                        if contigs_no != 0:
                            gf_per_contig[i][j][k] = gf / contigs_no
                        ng50s[i][j][k] = ng50
                j += 1
            i += 1

    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green']
    handles = []
    for i in range(len(tolerances)):
        for j in range(len(thresholds)):
            handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[j], markerfacecolor='None',
                                     linestyle='None', label=xticklabels[i*len(thresholds)+j]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("gf/contig vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('Genome fraction per contig')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    plt.tight_layout()
    #plt.rcParams["figure.autolayout"] = True
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(comparison_no):
                axs[i, j].plot(ng50s[i][j][k], gf_per_contig[i][j][k], marker=shapes[k % len(thresholds)], markeredgecolor=face_colors[k // len(thresholds)], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(out_path + "sub_plots_all/gf_per_contig.vs.ng50_" + repeat_sizes[0] + "K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1] + "_" + depths[0] + ".png")


def get_sub_plots_for_ng50_vs_contig_no_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, xlim, ylim, tolerances, thresholds):
    copy_no, snp_no = len(copies), len(snps)
    comparison_no = len(tolerances) * len(thresholds)
    xticklabels = []
    for to in tolerances:
        for th in thresholds:
            xticklabels.append(to + "_" + th)
    ng50s, contig_nos = np.zeros((copy_no, snp_no, comparison_no)), np.zeros((copy_no, snp_no, comparison_no))
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
                        ng50 = quast_data[id][metrics[1]][k] / 1000
                        contigs_no = quast_data[id][metrics[0]][k]
                        ng50s[i][j][k], contig_nos[i][j][k] = ng50, contigs_no
                j += 1
            i += 1

    shapes = ['^', 'o', 's', 'd']
    face_colors = ['red', 'blue', 'green']
    handles = []
    for i in range(len(tolerances)):
        for j in range(len(thresholds)):
            handles.append(mlines.Line2D([], [], markeredgecolor=face_colors[i], marker=shapes[j], markerfacecolor='None',
                                     linestyle='None', label=xticklabels[i*len(thresholds)+j]))
    fig, axs = plt.subplots(copy_no, snp_no)
    fig.suptitle("Contig no. vs. NG50 for repeats of " + repeat_sizes[0] + " Kbp with read coverage " + depths[0] + "x")
    #fig.tight_layout()
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('NG50 (Kbp)')
    plt.ylabel('number of contigs')
    plt.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=6)
    #plt.rcParams["figure.autolayout"] = True
    plt.tight_layout()
    for i in range(copy_no):
        for j in range(snp_no):
            for k in range(comparison_no):
                axs[i, j].plot(ng50s[i][j][k], contig_nos[i][j][k], marker=shapes[k % len(thresholds)], markeredgecolor=face_colors[k // len(thresholds)], markerfacecolor='none')
            axs[i, j].set_title(copies[i]+', '+snps[j]+' bp', fontsize=8)
            axs[i, j].set_xlim(left=0, right=xlim)
            axs[i, j].set_ylim(bottom=ylim[0], top=ylim[1])

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(out_path + "sub_plots_all/contig.vs.ng50_" + repeat_sizes[0] + "K_"+copies[0]+"-"+copies[-1]+"_"+snps[0]+"-"+snps[-1] + "_" + depths[0] + ".png")


def subplotter_hyperparameters(out_path, quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds, sub_plotter_method):
    x_lim = [200, 250, 300]
    for r in range(1, len(repeat_sizes)+1):
        repeat = repeat_sizes[r-1:r]
        for d in range(1, len(depths)+1):
            dep = depths[d-1:d]
            y_lim = sub_plot_ylimits[sub_plot[sub_plotter_method]]
            if sub_plotter_method == 0:
                get_sub_plots_for_ng50_vs_contig_no_hyperparameters(out_path, quast_data, repeat, copies, snps, dep, x_lim[r-1], y_lim, tolerances, thresholds)
            elif sub_plotter_method == 1:
                get_sub_plots_for_ng50_vs_misassemblies_hyperparameters(out_path, quast_data, repeat, copies, snps, dep, x_lim[r-1], y_lim, tolerances, thresholds)
            elif sub_plotter_method == 2:
                get_sub_plots_for_ng50_vs_gf_per_contig_hyperparameters(out_path, quast_data, repeat, copies, snps, dep, x_lim[r-1], y_lim, tolerances, thresholds)


tolerances = ["1", "5", "10", "15", "20"]
thresholds = ["5", "10", "15", "20", "25", "30", "35", "40", "45", "50"]
repeat_sizes = ["10", "15", "20"] #, "15", "20"] #["5", "10", "15", "20"]
copies = ["2", "5", "10"] #["2", "5", "10"]
snps = ["100", "250", "500", "1000", "2000"] #["100", "250", "500", "1000", "2000"]
depths = ["20", "30", "40"] #["20", "30", "40"]
sub_plot = {0: "ng50_vs_contig_no", 1: "ng50_vs_misassemblies", 2: "ng50_vs_gf_per_contig"}
sub_plot_ylimits = {sub_plot[0]: (0, 7),
                   sub_plot[1]: (-1, 5),
                   sub_plot[2]: (0, 110)}
sub_plotter_method = 0
input_path = "/Users/sakshar5068/Desktop/repeat_assembler/hyperparameters/quast/"
figure_path = "/Users/sakshar5068/Desktop/repeat_assembler/figures_paper/hyperparameters/"
quast_data = get_quast_reports(input_path, repeat_sizes, copies, snps, depths, tolerances, thresholds)

modified_quast_data = preprocess_quast_data(quast_data)
#metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# mismatches per 100 kbp', '# misassemblies']
metrics = ['# contigs', 'NG50', 'Genome fraction (%)', '# misassemblies']
#print(modified_quast_data)
#get_box_plot_for_genome_fraction_per_contig_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_box_plot_for_ng50_wrt_ref_size_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_box_plot_for_misassemblies_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_box_plot_for_contigs_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_box_plot_for_ComAcCon_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_violin_plots_for_misassemblies_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#get_violin_plots_for_contigs_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds)
#subplotter_hyperparameters(figure_path, modified_quast_data, repeat_sizes, copies, snps, depths, tolerances, thresholds, sub_plotter_method)