# Dependencies
import pandas as pd
import scipy.stats as stats
import numpy as np
import matplotlib.pylab as plt
import matplotlib.pyplot as pypl
import matplotlib.colors
from adjustText import adjust_text
import gseapy as gp
import seaborn as sns
from matplotlib.lines import Line2D

# Python built-in dependencies
import itertools
import functools
import math
import multiprocessing
import warnings
import os


class RNAseq:
    def __init__(self):
        print("Welcome to BroMynn's Bioinformatics Suite for RNA-sequencing data.")
        self.__isPV_calculated = False
        self.__isAnnotDB_given = False
        self.__isGMT_given = False
        self.__isCorrel_calculated = False
        self.__correl_genes_included = set()
        self.__correl_indices_included = set()
        self.__isGSEA_performed = False


    # ------------------------------------------------------------------------------------------------------------#
    def __removeduplicate(self, data):
        countdict = {}
        for element in data:
            if element in countdict.keys():

                # increasing the count if the key(or element)
                # is already in the dictionary
                countdict[element] += 1
            else:
                # inserting the element as key  with count = 1
                countdict[element] = 1
        cal = []
        for key in countdict.keys():
            cal.append(key)
        return cal

    def __label_indices(self, i):
        return [index for index in range(len(self.labels_received)) if self.labels_received[index] == i]

    def __comb_naming_rule(self, i, j):
        return str(i) + '-' + str(j)

    def __pv_fdr_correction(self, p_vals):
        ranked_p_values = stats.rankdata(p_vals)
        fdr = p_vals * len(p_vals) / ranked_p_values
        fdr[fdr > 1] = 1
        return fdr

    # def __convert_ensembl_to_hgnc_df(self, data, keep_which_duplicate='maximum'):
    #     for i in range(data.shape[-1]):
    #         # Convert
    #         cal = data.iloc[:, [i]].merge(self.__annotation_db,
    #                          how="outer",
    #                          on="Ensembl_id"
    #                          )
    #
    #         if keep_which_duplicate != 'maximum':
    #             raise NotImplementedError(
    #                 "Multiple Ensembl IDs sometimes map to a single HGNC gene symbol due to alternate splice and etc.\
    #                 The convention is to keep the Ensembl ID with the maximum value, although other methods are possible.\
    #                 In this version of LT2M, duplicate handling is limited to the maximum policy."
    #             )
    #         dropped = cal[['HGNC_id', data.iloc[:, [i]].keys()[0]]].dropna().sort_values(ascending=False,
    #                                                                         by=data.iloc[:, [i]].keys()[0]).drop_duplicates(
    #             keep='first')
    #         # final = dropped[['HGNC_id', data.keys()[0]]]
    #         if i == 0:
    #             final = dropped.reset_index(drop=True).set_index('HGNC_id')
    #         else:
    #             final = final.merge(dropped,
    #                         how="outer",
    #                         on="HGNC_id")
    #     else:
    #         if 'HGNC_id' in final.index:
    #             final = final.set_index('HGNC_id')
    #         # final = final.loc[(final != 0).any(axis=1)]
    #     return final

    def __convert_ensembl_to_hgnc_df(self, data, keep_which_duplicate='maximum'):
        cal = data.rename({self.__annotation_db.index[i] : self.__annotation_db.iloc[i, 0] for i in
                           range(self.__annotation_db.shape[0])})
        dropped = cal.dropna()
        dropped['cal_sum'] = dropped.sum(axis=1)
        dropped['original_index'] = [i for i in range(dropped.shape[0])]
        dropped['indices'] = [i for i in dropped.index]
        final = dropped.sort_values(by='cal_sum', ascending=False).drop_duplicates(subset='indices', keep='first').sort_values(
            'original_index').iloc[:, :-3].rename_axis('HGNC_id')
        return final

    def __convert_ensembl_to_hgnc_single(self, ensembl_id):
        if ensembl_id in self.__annotation_db.index:
            return self.__annotation_db.loc[ensembl_id][0]
        else:
            return ensembl_id

    # ------------------------------------------------------------------------------------------------------------#

    def __normalize_with_TPM(path, read_all_files_in_path, reads_column, gene_length_column, ensembl_id_column,):
        if read_all_files_in_path == True:
            to_read = os.listdir(path)
        else:
            if type(read_all_files_in_path) != list:
                raise ValueError("If you're not processing all files in the path, you should designate the file in each file's name contained in list.")
            else:
                to_read = read_all_files_in_path

        for k in to_read:
            cal = pd.read_csv(path + '\\' + k, sep='\t')
            RPK = []
            for i in range(cal.shape[0]):
                RPK.append(cal.iloc[i, reads_column] / cal.iloc[i, gene_length_column])
            else:
                cal["RPK"] = RPK

            RPK_sum = cal["RPK"].sum()
            #     print(RPK_sum)

            TPM = []
            for j in range(cal.shape[0]):
                TPM.append(cal["RPK"][j] / RPK_sum * 1_000_000)
            else:
                cal["TPM"] = TPM
            if k == to_read[0]:
                total = cal[["Geneid", "TPM"]]
                total = total.rename(
                    columns={"Geneid": "Ensembl_id", "TPM": k})
                total.set_index("Ensembl_id", inplace=True)
            else:
                cal = cal[["Geneid", "TPM"]].rename(
                    columns={"Geneid": "Ensembl_id", "TPM": k})
                cal.set_index("Ensembl_id", inplace=True)
                total = total.merge(cal,
                                    how="outer",
                                    on="Ensembl_id"
                                    )
                # print(total)

        return total

    def __normalize_with_RPKM(self, path, read_all_files_in_path=True, reads_column=2, gene_length_column=1,
                             ensembl_column=0, ):
        raise NotImplementedError("RPKM is not yet implemented. Please use TPM.")

    def __normalize_with_FPKM(self, path, read_all_files_in_path=True, reads_column=2, gene_length_column=1,
                             ensembl_column=0, ):
        raise NotImplementedError("FPKM is not yet implemented. Please use TPM.")

    def __normalize_with_CPM(self, path, read_all_files_in_path=True, reads_column=2, gene_length_column=1,
                             ensembl_column=0, ):
        raise NotImplementedError("CPM is not yet implemented. Please use TPM.")


    ## Function pointers for RNA normalization
    __rna_normalization = {
        "TPM" : __normalize_with_TPM,
        "RPKM" : __normalize_with_RPKM,
        "FPKM" : __normalize_with_FPKM,
        "CPM": __normalize_with_CPM,
    }

    # ------------------------------------------------------------------------------------------------------------#

    ## Function pointers for statistical significance
    __significance_tests = {
        # Welch's t-test
        "ttest": functools.partial(stats.ttest_ind, equal_var=False),
        "mannwhitney": stats.mannwhitneyu,
    }

    __fc_scale = {
        "log": np.log2,
        "none": lambda x: x,
    }

    __averaging_methods = {
        "arithmetic": np.mean,
        "harmonic": stats.hmean,
    }

    # ------------------------------------------------------------------------------------------------------------#

    def read_feature_count(self, path, read_all_files_in_path=True, normalizer="TPM", reads_column=2, gene_length_column=1, ensembl_id_column=0, output=True):

        cal = self.__rna_normalization[normalizer](path, read_all_files_in_path, reads_column, gene_length_column, ensembl_id_column,)

        self.data = cal

        if output==True:
            cal.to_csv(path.split('\\')[-1]+'.tsv', sep='\t')


    def read_tsv(self, file_path, sep="\t", header='infer', index_col=0):
        self.data = pd.read_csv(file_path, sep=sep, header=header, index_col=index_col)

    def designate_labels(self, a_list):
        self.labels_received = a_list
        self.__unique_labels = self.__removeduplicate(a_list)

    def set_annotation_db(self, path, sep='\t'):
        self.__isAnnotDB_given = True

        self.__annotation_db_path = path
        self.__annotation_db_sep = sep

        # Load the Ensembl id to gene symbol database
        annotation_db = pd.read_csv(self.__annotation_db_path, sep=self.__annotation_db_sep)
        self.__annotation_db = annotation_db.rename(columns={"Probe Set ID": "Ensembl_id", "Gene Symbol": "HGNC_id"})[
            ["Ensembl_id", "HGNC_id"]].set_index("Ensembl_id")

    def set_gene_sets_of_interest(self, a_list):
        self.gene_sets_of_interest = a_list

    def set_GMT(self, path):
        self.GMT_path = path
        self.__isGMT_given = True

    def convert_data_to_hgnc(self):
        if self.__isAnnotDB_given != True:
            raise ValueError("Annotation DB is not specified.")

        self.data = self.__convert_ensembl_to_hgnc_df(self.data)


    # ------------------------------------------------------------------------------------------------------------#

    def extract_DEGs(self, p_threshold=0.05, fc_threshold=0.5, sig_test='ttest', fc_scale='log', fc_mean='arithmetic',
                     convert_to_hgnc=True,
                     volcano=True, volcano_figsize=(10, 6), dpi=300, dotsize=0.2, dotsize_deg=0.8, volcano_annotate=True, fontsize=5, fixed_range=(False, -10, 10, 10), max_annotations=100,
                     heatmap=True, heatmap_normalization=True, heatmap_figsize=(5, 10), heatmap_max_size=100, heatmap_fontsize=7, heatmap_fname='heatmap'):
        # Broadcast indices to calculate
        cal_pv_table = {}
        cal_fc_table = {}

        # Perform significance test and fold change calculation between every combination of label indices
        for i, j in itertools.combinations(self.__unique_labels, 2):
            if self.__isPV_calculated == True:
                break
            group_a = self.__label_indices(i)
            group_b = self.__label_indices(j)

            cal_pv_table[self.__comb_naming_rule(i, j)] = [
                self.__significance_tests[sig_test](self.data.iloc[k, group_a], self.data.iloc[k, group_b]).pvalue
                for k in range(self.data.shape[0])
            ]

            if fc_scale == 'log':
                # Note that when calculating in log scale, 1 is automatically added to prevent data from going 0.
                cal_fc_table[self.__comb_naming_rule(i, j)] = [
                    self.__fc_scale[fc_scale](self.__averaging_methods[fc_mean](self.data.iloc[k, group_a]) + 1)
                    -
                    self.__fc_scale[fc_scale](self.__averaging_methods[fc_mean](self.data.iloc[k, group_b]) + 1)
                    for k in range(self.data.shape[0])
                ]
            else:
                cal_fc_table[self.__comb_naming_rule(i, j)] = [
                    self.__fc_scale[fc_scale](self.__averaging_methods[fc_mean](self.data.iloc[k, group_a])
                                              /
                                              self.__averaging_methods[fc_mean](self.data.iloc[k, group_b]))
                    for k in range(self.data.shape[0])
                ]
        else:
            self.pv_table = pd.DataFrame(cal_pv_table).set_index(self.data.index)
            self.fc_table = pd.DataFrame(cal_fc_table).set_index(self.data.index)
            self.__isPV_calculated = True

        cal_DEGs = {}
        for i, j in itertools.combinations(self.__unique_labels, 2):
            cal_DEGs[self.__comb_naming_rule(i, j)] = list(np.intersect1d(
                self.pv_table[self.__comb_naming_rule(i, j)][
                    self.pv_table[self.__comb_naming_rule(i, j)] < p_threshold].index,
                self.fc_table[self.__comb_naming_rule(i, j)][
                    (self.fc_table[self.__comb_naming_rule(i, j)] > abs(fc_threshold))
                    | (self.fc_table[self.__comb_naming_rule(i, j)] < -abs(fc_threshold))].index,
            ))
        if convert_to_hgnc == True:
            self.__DEGs_original_index = cal_DEGs
            self.DEGs = {k: [self.__convert_ensembl_to_hgnc_single(z) for z in cal_DEGs[k]] for k in cal_DEGs}
        else:
            self.__DEGs_original_index = cal_DEGs
            self.DEGs = cal_DEGs

        if volcano == True:
            for k, l in itertools.combinations(self.__unique_labels, 2):
                i = self.__comb_naming_rule(k, l)
                cal_indices = self.pv_table[i].notna()
                x_total = self.fc_table.loc[cal_indices][i]
                y_total = -np.log10(self.pv_table.loc[cal_indices][i])
                plt.figure(figsize=volcano_figsize, dpi=dpi)
                plt.scatter(
                    x=x_total,
                    y=y_total,
                    s=dotsize,
                    label="Not significant",
                    color='gray'
                )

                upregulated_indices = np.intersect1d(
                    x_total[x_total > abs(fc_threshold)].index, y_total[y_total > -np.log10(p_threshold)].index
                )

                downregulated_indices = np.intersect1d(
                    x_total[x_total < -abs(fc_threshold)].index, y_total[y_total > -np.log10(p_threshold)].index
                )

                plt.scatter(
                    x=x_total[upregulated_indices],
                    y=y_total[upregulated_indices],
                    s=dotsize_deg,
                    label="Up in "+k,
                    color="red",
                )

                plt.scatter(
                    x=x_total[downregulated_indices],
                    y=y_total[downregulated_indices],
                    s=dotsize_deg,
                    label="Up in "+l,
                )

                if volcano_annotate == True and len(upregulated_indices) + len(downregulated_indices) != 0:
                    texts = []
                    if max_annotations < len(upregulated_indices) + len(downregulated_indices):
                        ranking_cal = list(upregulated_indices)
                        ranking_cal.extend(list(downregulated_indices))
                        # Puts priority on PV * FC value
                        selected_indices = y_total[ranking_cal][
                            (y_total[ranking_cal] * x_total[ranking_cal].apply(lambda x: abs(x))).rank(
                                ascending=False) <= max_annotations].dropna().index
                    else:
                        selected_indices = list(upregulated_indices)
                        selected_indices.extend((list(downregulated_indices)))

                    if convert_to_hgnc == True and self.__isAnnotDB_given == True:
                        for m in x_total[selected_indices].index:
                            texts.append(
                                plt.text(s=self.__convert_ensembl_to_hgnc_single(m), x=x_total[m], y=y_total[m],
                                         fontsize=fontsize))
                    else:
                        for m in x_total[selected_indices].index:
                            texts.append(plt.text(s=m, x=x_total[m], y=y_total[m], fontsize=fontsize))

                    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.1))

                plt.xlabel("log(FC)")
                plt.ylabel("-log(Pval)")
                plt.axvline(fc_threshold, color="grey", linestyle="--", linewidth=1)
                plt.axvline(-fc_threshold, color="grey", linestyle="--", linewidth=1)
                plt.axhline(-np.log10(p_threshold), color="grey", linestyle="--", linewidth=1)

                if fixed_range[0] == True:
                    plt.xlim(fixed_range[1], fixed_range[2])
                    plt.ylim(0, fixed_range[3])

                plt.legend()
                plt.title("◀ "+l+" : "+k+" ▶")
                plt.savefig(i + '_volcano_pv_' + str(p_threshold) + '_fc_' + str(fc_threshold) + '.png')

        if heatmap == True:
            norm = matplotlib.colors.Normalize(-3, 3)
            colors = [
                [norm(-3.0), (0, 0, 255 / 255)],
                [norm(-0.5), "white"],
                [norm(0.5), "white"],
                [norm(3.0), (255 / 255, 0, 0)],
            ]
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

            cal_DEGs_indices = []
            for i in self.__DEGs_original_index.keys():
                cal_DEGs_indices.extend(self.__DEGs_original_index[i])
            else:
                cal_DEGs_indices = list(set(cal_DEGs_indices))

            if len(cal_DEGs_indices) >= heatmap_max_size:
                # Find common DEGs first
                final = []
                common_DEGs = []
                for i in self.__DEGs_original_index.keys():
                    if common_DEGs != []:
                        common_DEGs = self.__DEGs_original_index[i]
                    else:
                        common_DEGs = [ab for ab in common_DEGs if ab in self.__DEGs_original_index[i]]
                else:
                    final.extend(common_DEGs)

                # Puts priority on PV * FC value again
                indices_minus_common = [ac for ac in cal_DEGs_indices if ac not in common_DEGs]
                cal = (self.pv_table.loc[indices_minus_common].apply(lambda x:np.log10(x)) * self.fc_table.loc[indices_minus_common].apply(lambda x:abs(x))).sum(axis=1)
                cal_DEGs_indices = cal[cal.rank(ascending=False) <= heatmap_max_size - len(common_DEGs)].index



            if convert_to_hgnc == True:
                conversion_table = {z:self.__convert_ensembl_to_hgnc_single(z) for z in cal_DEGs_indices}
                if heatmap_normalization == True:
                    heatmap_matrix = stats.zscore(self.data, axis=1).loc[cal_DEGs_indices].rename(index=conversion_table)
                else:
                    heatmap_matrix = self.data.loc[cal_DEGs_indices].rename(index=conversion_table)
            else:
                if heatmap_normalization == True:
                    heatmap_matrix = stats.zscore(self.data, axis=1).loc[cal_DEGs_indices]
                else:
                    heatmap_matrix = self.data.loc[cal_DEGs_indices]

            the_heatmap = sns.clustermap(
                heatmap_matrix,
                cmap=cmap,
                # zscore was calculated in the previous heatmap_matrix
                # z_score=0,
                col_cluster=True,
                row_cluster=True,
                vmin=-3,
                vmax=3,
                xticklabels=1,
                yticklabels=1,
                figsize=heatmap_figsize,
                fmt='.2f',
                cbar_pos=(-0.1, 0.8, 0.05, 0.18),
            )
            the_heatmap.ax_heatmap.set_yticklabels(the_heatmap.ax_heatmap.get_ymajorticklabels(), fontsize=heatmap_fontsize)
            the_heatmap.fig.suptitle(heatmap_fname, y=1.01).get_figure().savefig(heatmap_fname, bbox_inches='tight', dpi=dpi)


    def analyze_enrichr(self, gene_list=None, gene_sets_to_select=None, organism='Human', drop_non_hgnc=True, cutoff=0.05, top_term=10,
                        dotsize=15, figsize=(6, 20), dpi=300, save_table=True):
        if gene_sets_to_select == None:
            cal_gene_sets = self.gene_sets_of_interest
        elif type(gene_sets_to_select) == list:
            cal_gene_sets = gene_sets_to_select

        if gene_list == None:
            gene_lists_dict = self.DEGs
        elif type(gene_list) == dict:
            gene_lists_dict = gene_list
        elif type(gene_list) in [list, tuple]:
            gene_lists_dict = {'Custom_gene_list': gene_list}
        else:
            raise ValueError("Neither DEGs are previously extracted nor gene list or gene lists dictionary is passed.")

        for i in gene_lists_dict:
            if drop_non_hgnc == True:
                cal_DEGs = [k for k in self.DEGs[i] if k in self.__annotation_db.values]
            else:
                cal_DEGs = self.DEGs

            try:
                enr = gp.enrichr(
                    gene_list=cal_DEGs,
                    gene_sets=cal_gene_sets,
                    organism=organism,
                )

                if save_table == True:
                    enr.results.to_csv(i + '_enrichr.tsv', sep='\t')

                gp.dotplot(
                    enr.results,
                    column="Adjusted P-value",
                    cutoff=cutoff,
                    x='Gene Sets',
                    top_term=top_term,
                    figsize=figsize,
                    size=dotsize,
                    title="Enrichr",
                    xticklabels_rot=45,  # rotate xtick labels
                    show_ring=True,  # set to False to remove outer ring
                    marker='o',
                    dpi=dpi,
                    ofname=i + '_enrichr.png'
                )

                gp.dotplot(
                    enr.results,
                    column="Adjusted P-value",
                    cutoff=cutoff,
                    x='Gene Sets',
                    top_term=top_term,
                    figsize=figsize,
                    size=dotsize,
                    title="Enrichr",
                    xticklabels_rot=45,  # rotate xtick labels
                    show_ring=True,  # set to False to remove outer ring
                    marker='o',
                    dpi=dpi,
                )
            except:
                print("An error occurred for " + i)

    def analyze_custom_ORA(self, input_gmt=None, drop_non_hgnc=False, cutoff=0.05, top_term=10,
                           dotsize=15, figsize=(6, 20), dpi=300, save_table=True):

        if self.__isGMT_given == True:
            cal_gene_sets = self.GMT_path
        elif self.__isGMT_given == False and type(input_gmt) in [str, pd.DataFrame, dict]:
            cal_gene_sets = input_gmt
        else:
            ValueError("No gene sets are given.")

        for i in self.DEGs:
            if drop_non_hgnc == True:
                cal_DEGs = [k for k in self.DEGs[i] if k in self.__annotation_db.values]
            else:
                cal_DEGs = self.DEGs

            try:
                enr = gp.enrich(
                    gene_list=cal_DEGs,
                    gene_sets=cal_gene_sets,
                )

                if save_table == True:
                    enr.results.to_csv(i + '_enrichr.tsv', sep='\t')

                gp.dotplot(
                    enr.results,
                    column="Adjusted P-value",
                    cutoff=cutoff,
                    x='Gene Sets',
                    top_term=top_term,
                    figsize=figsize,
                    size=dotsize,
                    title="Enrichr",
                    xticklabels_rot=45,  # rotate xtick labels
                    show_ring=True,  # set to False to remove outer ring
                    marker='o',
                    dpi=dpi,
                    ofname=i + '_enrichr.png'
                )

                gp.dotplot(
                    enr.results,
                    column="Adjusted P-value",
                    cutoff=cutoff,
                    x='Gene Sets',
                    top_term=top_term,
                    figsize=figsize,
                    size=dotsize,
                    title="Enrichr",
                    xticklabels_rot=45,  # rotate xtick labels
                    show_ring=True,  # set to False to revmove outer ring
                    marker='o',
                    dpi=dpi,
                )
            except:
                print("An error occurred for " + i)

    def analyze_correlation(self, genes_of_interest=None, pv_cutoff=0.001, correl_cutoff=0.8, labels_to_include=None, corregulation_subtyping=True):

        # Checking the integrity of genes of interest
        if genes_of_interest == None:
            genes_to_investigate = list(self.data.index)
        elif type(genes_of_interest) == list:
            genes_to_investigate = genes_of_interest
        else:
            raise ValueError("An inappropriate gene list is given. Must be a list with index numbers of the group.")

        # Checking if these genes have previously calculated
        if self.__correl_genes_included != set(genes_to_investigate):
            self.__isCorrel_calculated = False

        # Checking the integrity to groups to include indices
        included_indices = []
        if labels_to_include == None:
            for i in self.__removeduplicate(self.labels_received):
                included_indices.extend(self.__label_indices(i))

        elif type(labels_to_include) == list:
            for i in labels_to_include:
                included_indices.extend(self.__label_indices(i))
        else:
            raise ValueError("An inappropriate gene list is given. Must be a list with index numbers of the group.")

        # Checking if the given parameter has previously calculated
        if self.__correl_indices_included != set(included_indices):
            self.__isCorrel_calculated = False

        # Saving the current genes & indices combination
        self.__correl_indices_included = set(included_indices)
        self.__correl_genes_included = set(genes_to_investigate)

        warnings.filterwarnings(action='ignore')
        correlation_dict = {}
        unabbreviated_correl_dict = {}
        if self.__isCorrel_calculated == False:
            for i in genes_to_investigate:
                # first empty list for gene name, second empty list for correlation coeff, third empty list for pvalue
                correlation_dict[i] = [[], [], []]
                unabbreviated_correl_dict[i] = [[], [], []]

                for k in range(self.data.shape[0]):
                    cal = stats.pearsonr(
                        self.data.loc[i].iloc[included_indices],
                        self.data.iloc[k, included_indices],
                    )
                    unabbreviated_correl_dict[i][0].append(self.data.iloc[k, :].name)  # gene name
                    unabbreviated_correl_dict[i][1].append(cal.correlation)  # correlation
                    unabbreviated_correl_dict[i][2].append(cal.pvalue)  # pvalue
                    if cal.pvalue < pv_cutoff and abs(cal.correlation) > correl_cutoff:
                        correlation_dict[i][0].append(self.data.iloc[k, :].name)  # gene name
                        correlation_dict[i][1].append(cal.correlation)  # correlation
                        correlation_dict[i][2].append(cal.pvalue)  # pvalue
            else:
                self.__all_correl = unabbreviated_correl_dict
                self.correl = correlation_dict
                self.__isCorrel_calculated = True
        else:
            self.correl = {}
            for i in self.__all_correl.keys():
                correlation_dict[i] = [[], [], []]
                for j in range(len(self.__all_correl[i][0])):
                    if self.__all_correl[i][2][j] < pv_cutoff and abs(self.__all_correl[i][1][j]) > correl_cutoff:
                        correlation_dict[i][0].append(self.__all_correl[i][0][j])  # gene name
                        correlation_dict[i][1].append(self.__all_correl[i][1][j])  # correlation
                        correlation_dict[i][2].append(self.__all_correl[i][2][j])  # pvalue
            else:
                self.correl = correlation_dict
                self.__isCorrel_calculated = True
        warnings.filterwarnings(action='default')

        if corregulation_subtyping == True:
            # Find commonly significant genes
            common_genes_with_significant_correlation = None
            for i in genes_to_investigate:
                if common_genes_with_significant_correlation == None:
                    common_genes_with_significant_correlation = set(self.correl[i][0])
                else:
                    common_genes_with_significant_correlation.intersection_update(self.correl[i][0])

            # Container for DataFrame-ized Correl Values
            data_frame_collection = [pd.DataFrame(self.correl[i]).transpose().set_index(0) for i in
                                     genes_to_investigate]

            # All possible correlation sign signatures except 0
            x = [1 for _ in range(len(genes_to_investigate))]
            x.extend(-1 for _ in range(len(genes_to_investigate)))
            signature = sorted(list(set(itertools.permutations(x, len(genes_to_investigate)))),reverse=True)

            corregulated_dict = {}
            for sig_index in range(len(signature)):
                corregulated_dict[signature[sig_index]] = []
                for i in common_genes_with_significant_correlation:
                    if (np.sign([j.loc[i][1] for j in data_frame_collection]) == signature[sig_index]).sum() == len(
                            signature[sig_index]):
                        corregulated_dict[signature[sig_index]].append(i)
            self.corregulated = corregulated_dict



    def analyze_GSEA(self, label1, label2, analysis_name, convert_hgnc=True,
                     ranking_method='signal_to_noise', permutation=1000, permutation_type='phenotype', override='ask',
                     fdr_threshold=0.10, gene_set_min=10, gene_set_max=500,
                     top_term=False, fontsize=14, fig_size=(16,32),
                     fig_save=True,):

        if self.__isGSEA_performed == False:
            self.GSEA_results = {}


        indices_to_include = []
        class_vector = []
        indices_to_include.extend(self.__label_indices(label1))
        class_vector.extend([label1 for i in range(len(self.__label_indices(label1)))])
        indices_to_include.extend(self.__label_indices(label2))
        class_vector.extend([label2 for i in range(len(self.__label_indices(label2)))])

        if convert_hgnc == True:
            gene_exp = self.__convert_ensembl_to_hgnc_df(self.data.iloc[:, indices_to_include])
        else:
            gene_exp = self.data.iloc[:, indices_to_include]

        gene_sets = self.GMT_path

        new_gsea_flag = True

        if self.__isGSEA_performed == True and self.__comb_naming_rule(label1, label2) in self.GSEA_results.keys():
            if override == 'ask':
                the_answer = input("You have performed GSEA with the same class vectors. Do you wish to override the previous results? [y/n]")
                if the_answer in ['Y','y']:
                    new_gsea_flag = True
                elif the_answer in ['N','n']:
                    new_gsea_flag = False
                else:
                    raise ValueError("Invalid answer")
            if override == True:
                new_gsea_flag = True


        if new_gsea_flag == True:
            print("GSEA initiated for ", analysis_name)
            gs_res = gp.gsea(data=gene_exp,
                             gene_sets=gene_sets,
                             cls=class_vector,
                             # set permutation_type to phenotype if samples >=15
                             permutation_type=permutation_type,
                             permutation_num=permutation,  # reduce number to speed up test
                             outdir=None,  # do not write output to disk
                             method=ranking_method,
                             min_size=gene_set_min,
                             max_size=gene_set_max,
                             threads=multiprocessing.cpu_count(),
                             seed=7)
            print("Saving results.")
            gs_res.res2d.to_csv(analysis_name + ".tsv", sep="\t")
            terms = gs_res.res2d.Term

            self.__isGSEA_performed = True
            self.GSEA_results[self.__comb_naming_rule(label1, label2)] = gs_res



        else:
            gs_res = self.GSEA_results[self.__comb_naming_rule(label1, label2)]
            terms = gs_res.res2d.Term

        try:
            ax = gp.dotplot(gs_res.res2d,
                            column="FDR q-val",
                            title="",
                            cmap=pypl.cm.viridis,
                            size=10,
                            top_term=len(gs_res.res2d) if top_term == False else top_term,
                            figsize=fig_size,
                            cutoff=fdr_threshold,
                            ofname=analysis_name + "_summary.png",
                            y_fontsize=fontsize,
                            dpi=300
                            )
        except:
            pass

        if fig_save == True:
            for i in range(len(gs_res.res2d)):
                if gs_res.results[terms[i]]['fdr'] >= fdr_threshold:
                    continue
                gp.gseaplot(rank_metric=gs_res.ranking, term=terms[i],
                            ofname=analysis_name + '_' + str(i) + '.png', **gs_res.results[terms[i]])



    def iterate_GSEA_for_all_labels(self, convert_hgnc=True,
                     ranking_method='signal_to_noise', permutation=1000, permutation_type='phenotype', override='ask',
                     fdr_threshold=0.10, gene_set_min=10, gene_set_max=500,
                     top_term=False, fontsize=14, fig_size=(16,32),
                     fig_save=True,):
        for i, j in itertools.combinations(self.__unique_labels, 2):
            self.analyze_GSEA(label1=i, label2=j, analysis_name=self.__comb_naming_rule(i,j), convert_hgnc=convert_hgnc,
                     ranking_method=ranking_method, permutation=permutation, permutation_type=permutation_type, override=override,
                     fdr_threshold=fdr_threshold, gene_set_min=gene_set_min, gene_set_max=gene_set_max,
                     top_term=top_term, fontsize=fontsize, fig_size=fig_size,
                     fig_save=fig_save)

    def analyze_GSEA_prerank(self, label1, label2, analysis_name, prerank_mode='fcsign_log10pv', permutation=1000,  gene_set_min=10, gene_set_max=500,
                             fig_save_mode=True, fdr_threshold=0.1, top_term=False, fontsize=14, fig_size=(16,32)):
        if prerank_mode == 'fcsign_log10pv':
            if self.__isPV_calculated == False:
                self.extract_DEGs(volcano=False, heatmap=False)
            naming = self.__comb_naming_rule(label1, label2)
            prerank_data = self.__convert_ensembl_to_hgnc_df(
                pd.DataFrame(
                    (
                            self.pv_table[naming].apply(lambda x:-np.log10(x)) * self.fc_table[naming].apply(lambda x:np.sign(x))
                    ).fillna(0).sort_values(ascending=False)
                )
            )

        gene_sets = self.GMT_path

        pre_res = gp.prerank(rnk=prerank_data,
                             gene_sets=gene_sets,
                             threads=multiprocessing.cpu_count(),
                             min_size=gene_set_min,
                             max_size=gene_set_max,
                             permutation_num=permutation,  # reduce number to speed up testing
                             outdir=None,  # don't write to disk
                             seed=6,
                             verbose=True,  # see what's going on behind the scenes
                             )

        pre_res.res2d.to_csv(analysis_name + ".tsv", sep="\t")

        try:
            ax = gp.dotplot(pre_res.res2d,
                         column="FDR q-val",
                         title="",
                         cmap=plt.cm.viridis,
                         size=10,
                         top_term=len(pre_res.res2d) if top_term == False else top_term,
                         figsize=fig_size,
                         cutoff=fdr_threshold,
                         ofname=naming + "_summary.png",
                         y_fontsize=fontsize,
                         dpi=300
                         )
        except:
            pass


        terms = pre_res.res2d.Term

        for i in range(len(pre_res.res2d)):
            if fig_save_mode != True:
                continue

            if pre_res.results[terms[i]]['fdr'] >= fdr_threshold:
                continue

            # save figure
            gp.gseaplot(rank_metric=pre_res.ranking, term=terms[i],
                     ofname=analysis_name+ '_' + str(i) + '.png', **pre_res.results[terms[i]])

    def analyze_vertical_perceptron(self, label1, label2, accuracy=1, ivns_threshold=1, sort_by='IVNS', convert_hgnc=True, save_fig=True, label1_color='k', label2_color='g', dot_size=10, dot_transparency=0.6, verbose=True):
        # First, gather all indices to use, orderly.
        y = self.__label_indices(label1)
        y.extend(self.__label_indices(label2))

        # Second, select data to use with the gathered indices. Note that the data is transposed.
        x = self.data.iloc[:, y].transpose()
        # Then, use length of the first label to assign y-labels (first label = 0, second label = 1)
        y = [0 if i < len(self.__label_indices(label1)) else 1 for i in range(len(y))]


        target_accuracy = accuracy

        # Select the sequence of indices for y-label 0 and y-label 1
        zeros = [i for i in range(y.__len__()) if y[i] == 0]
        ones = [i for i in range(y.__len__()) if y[i] == 1]

        # self.x = x
        # self.y = y
        # self.zeros = zeros
        # self.ones = ones

        # Now comes the training to find the separating line.
        success = {}
        # For accuracy of 100%, I have a very fast and exact version of the algorithm.
        if accuracy == 1:
            # Iterate for each gene.
            for i in range(x.shape[1]):
                # There is only two versions of separation in accuracy 100%, Label 1 being on left side or Label 1 being on the right side.
                decision_boundary_1 = (x.iloc[ones, i].max() + x.iloc[zeros, i].min()) / 2
                decision_boundary_2 = (x.iloc[zeros, i].max() + x.iloc[ones, i].min()) / 2

                if ((x.iloc[ones, i] < decision_boundary_1).sum() + (x.iloc[zeros, i] > decision_boundary_1).sum()) / x.shape[0] - ((x.iloc[zeros, i] < decision_boundary_2).sum() + (x.iloc[ones, i] > decision_boundary_2).sum()) / x.shape[0] >= 0:
                    # if the_left_is_zero = False:
                    if ((x.iloc[ones, i] < decision_boundary_1).sum() + (x.iloc[zeros, i] > decision_boundary_1).sum()) / x.shape[0] >= target_accuracy:
                        # The dictionary has gene names for keys, and for each row represents the following.
                        success[x.iloc[zeros, i].name] = [
                            # 1st row: decision boundary
                            decision_boundary_1,
                            # 2nd row: accuracy
                            ((x.iloc[ones, i] < decision_boundary_1).sum() + (
                                        x.iloc[zeros, i] > decision_boundary_1).sum()) / x.shape[0],
                            # 3rd row: distance from the decision boundary
                            decision_boundary_1 - x.iloc[ones, i].max(),
                            # 4th row: the mean of each group divided by their standard deviation
                            x.iloc[zeros, i].mean() / x.iloc[zeros, i].std() - x.iloc[ones, i].mean() / x.iloc[
                                ones, i].std(),
                            # 5th row: Inter-group Variance Normalized Score (IVNS)
                            (decision_boundary_1 - x.iloc[ones, i].max()) / (
                                        (x.iloc[zeros, i].std() + x.iloc[ones, i].std()) / 2),
                        ]
                else:
                    # if the_left_is_zero = True:
                    if ((x.iloc[zeros, i] < decision_boundary_2).sum() + (x.iloc[ones, i] > decision_boundary_2).sum()) / x.shape[0] >= target_accuracy:
                        success[x.iloc[zeros, i].name] = [
                            decision_boundary_2,
                            ((x.iloc[zeros, i] < decision_boundary_2).sum() + (
                                        x.iloc[ones, i] > decision_boundary_2).sum()) / x.shape[0],
                            decision_boundary_2 - x.iloc[zeros, i].max(),
                            x.iloc[ones, i].mean() / x.iloc[ones, i].std() - x.iloc[zeros, i].mean() / x.iloc[
                                zeros, i].std(),
                            (decision_boundary_2 - x.iloc[zeros, i].max()) / (
                                        (x.iloc[ones, i].std() + x.iloc[zeros, i].std()) / 2),]
        else:
            for i in range(x.shape[1]):
                accuracy_dictionary = {}
                for zero,one in itertools.product(zeros, ones):
                    cal_decision_boundary = (x.iloc[zero, i] + x.iloc[one,i])/2
                    # This version cannot tell which label is left, so this version tries a rule-of-thumb: the one with the higher accuracy must be correct.
                    candidate1 = ((x.iloc[ones, i] < cal_decision_boundary).sum() + (x.iloc[zeros, i] > cal_decision_boundary).sum()) / x.shape[0]
                    candidate2 = ((x.iloc[ones, i] > cal_decision_boundary).sum() + (x.iloc[zeros, i] < cal_decision_boundary).sum()) / x.shape[0]
                    if candidate1 >= candidate2:
                        cal_accuracy = candidate1
                    else:
                        cal_accuracy = candidate2

                    if cal_accuracy >= accuracy:
                        if cal_accuracy not in accuracy_dictionary.keys():
                            accuracy_dictionary[cal_accuracy] = [(zero, one)]
                        else:
                            accuracy_dictionary[cal_accuracy].append((zero, one))
                else:
                    if accuracy_dictionary == {}:
                        continue
                    else:
                        if len(accuracy_dictionary[max(accuracy_dictionary)]) > 1:
                            if verbose == True:
                                print("Warning: "+str(x.iloc[:, i].name)+" has more than one decision boundary. Selecting where maximum seperation is the greatest.")
                            pre_maximum_cal = {}
                            for t1, t2 in accuracy_dictionary[max(accuracy_dictionary)]:
                                cal_decision_boundary = (x.iloc[t1, i] + x.iloc[t2, i])/2
                                cal_maximum = abs(cal_decision_boundary - x.iloc[t1, i])

                                if pre_maximum_cal == {}:
                                    pre_maximum_cal[cal_maximum] = [t1, t2]
                                else:
                                    if cal_maximum >= list(pre_maximum_cal.keys())[0]:
                                        pre_maximum_cal[cal_maximum] = [t1, t2]

                            the_zero = pre_maximum_cal[list(pre_maximum_cal.keys())[0]][0]
                            the_one = pre_maximum_cal[list(pre_maximum_cal.keys())[0]][1]
                        else:
                            the_zero = accuracy_dictionary[max(accuracy_dictionary)][0][0]
                            the_one = accuracy_dictionary[max(accuracy_dictionary)][0][1]

                        final_decision_boundary = (x.iloc[the_zero, i] + x.iloc[the_one, i])/2,
                        success[x.iloc[:, i].name] = [
                            # 1st row: decision boundary
                            final_decision_boundary,
                            # 2nd row: accuracy
                            max(accuracy_dictionary),
                            # 3rd row: distance from the decision boundary
                            abs(final_decision_boundary - x.iloc[the_zero, i]),
                            # 4th row: the mean of each group divided by their standard deviation
                            x.iloc[zeros, i].mean() / x.iloc[zeros, i].std() - x.iloc[ones, i].mean() / x.iloc[
                                ones, i].std(),
                            # 5th row: Inter-group Variance Normalized Score (IVNS)
                            abs(final_decision_boundary - x.iloc[the_zero, i]) / (
                                    (x.iloc[zeros, i].std() + x.iloc[ones, i].std()) / 2),
                        ]
        self.perceptron_success = success

        if success == {}:
            print("No genes separate accordingly to the designated accuracy.")
            return None

        sorting_rules = {"IVNS":4, "accuracy":1}

        sorting_column = sorting_rules[sort_by]

        perceptron_results = pd.DataFrame(success).transpose().sort_values(by=sorting_column, ascending=False)
        perceptron_results.index.name = 'Ensembl_id'

        interesting_genes = perceptron_results[perceptron_results[4] >= ivns_threshold].sort_values(by=sorting_column, ascending=False)

        if convert_hgnc == True:
            self.perceptron_results = self.__convert_ensembl_to_hgnc_df(perceptron_results).sort_values(by=sorting_column, ascending=False)
            self.perceptron_interesting = self.__convert_ensembl_to_hgnc_df(interesting_genes).sort_values(by=sorting_column, ascending=False)
        else:
            self.perceptron_results = perceptron_results
            self.perceptron_interesting = interesting_genes

        if save_fig == True:
            counter = 0
            rank_str_formatter = '0' + str(int(math.log(len(interesting_genes.index), 10) + 1))
            for i in interesting_genes.index:
                counter += 1
                gene_id = i

                fig, ax = plt.subplots()
                ax.scatter(np.array(self.data.iloc[:, self.__label_indices(label2)].loc[gene_id]), [1 for _ in self.__label_indices(label2)], color=label2_color, s=dot_size*8, alpha=dot_transparency)
                ax.scatter(np.array(self.data.iloc[:, self.__label_indices(label1)].loc[gene_id]), [0 for _ in self.__label_indices(label1)], color=label1_color, s=dot_size*8, alpha=dot_transparency)


                plt.axvline(x=interesting_genes.loc[i][0], color='r', label='axvline - full height', linestyle='--')

                title_name = self.__convert_ensembl_to_hgnc_single(i) if convert_hgnc == True else i

                if accuracy == 1:
                    plt.title(title_name + ' (#' + str(counter) + ')' + '\n' + 'accuracy = ' + str(
                        round(interesting_genes.loc[i][1],
                              2)) + '\n' + 'max separation = ' + str(
                        round(interesting_genes.loc[i][2][0],
                              2)) + '\n' + 'groups variance normalized separation = ' + str(
                        round(interesting_genes.loc[i][4][0], 2)))
                else:
                    plt.title(title_name + ' (#' + str(counter) + ')' + '\n'  + 'accuracy = ' + str(
                        round(interesting_genes.loc[i][1],
                              2)) + '\n'  + 'max separation = ' + str(
                        round(interesting_genes.loc[i][2][0],
                              2)) + '\n' + 'groups variance normalized separation = ' + str(
                        round(interesting_genes.loc[i][4][0], 2)))

                plt.ylim((-0.1, 1.1))
                yax = ax.axes.get_yaxis()
                yax = yax.set_visible(False)

                custom_legends = [
                    Line2D([0], [0], marker='o', color='w', label=label2, markerfacecolor=label2_color, markersize=dot_size),
                    Line2D([0], [0], marker='o', color='w', label=label1, markerfacecolor=label1_color, markersize=dot_size)
                ]

                plt.legend(handles=custom_legends)

                plt.savefig('#' + format(counter, rank_str_formatter) + ' ' + title_name + '.png',
                            dpi=300, bbox_inches='tight')
                plt.show()
                plt.close()
