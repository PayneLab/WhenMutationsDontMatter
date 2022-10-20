#!/usr/bin/env python
# coding: utf-8

# # Part 1: Create database set
# By inputting a target gene, a nested dictionary will be created, with each key being the name of a cancer and each value being a dictionary containing the genomic, proteomic, transcriptomic, and phosphoproteomic data. In other words, the dictionary will be:
# 
# {
#     cancer : {   
#     
#             "gen" : <cancer's genomic data>,
#             "prot" : <cancer's proteomic data>,
#             "trans": <cancer's transcriptomic data>,
#             "phos": <cancer's phosphoproteomic data>
#         }
#     ... for each cancer
# }

import cptac
import numpy as np
import cptac.utils as ut
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import warnings
from math import log
from sys import argv
warnings.filterwarnings("ignore")


try:
    target_gene = argv[1]
    if argv[2] in ['All', 'None']:
        path_indices = argv[2]
    else:
        path_indices = argv[2:]
except IndexError:
    raise ValueError(f"Script requires at least two arguments."
                    " arg1 = Name of gene in question."
                    " All other args specify the indices of the pathways in question."
                    " You may also pass \"All\" or \"None\" after the gene name."
                    "\n\n"
                    "{ut.get_pathways_with_proteins(proteins = target_gene, database = 'wikipathways')}"
                    )
## Create datasets
dataset_names = ['Brca', 'Ccrcc', 'Coad', 'Gbm', 'Hnscc', 'Lscc', 'Luad', 'Ov', 'Ucec', 'Pdac']
omics_key = {"prot":"proteomics", "trans":"transcriptomics","phos":"phosphoproteomics"}
datasets = {}
for name in dataset_names:
    print(f"Loading {name}...")
    exec(f"cancer = cptac.{name}()")
    genotype_data = cancer.get_genotype_all_vars(mutations_gene=target_gene, omics_source = "washu", mutations_source="washu")
    genotype_data.astype({'Mutation':'category', 'Mutation_Status':'category'})
    print('   Loaded genomics...')
    
    proteomic_data = cancer.get_proteomics(source = 'umich')
    proteomic_data = proteomic_data.droplevel("Database_ID", axis = 1)
    proteomic_data = proteomic_data.astype({col: 'float32' for col in proteomic_data.columns})
    print('   Loaded proteomics...')
    
    phospho_data = cancer.get_phosphoproteomics(source = 'umich')
    phospho_data = phospho_data.droplevel("Database_ID", axis = 1).droplevel("Peptide", axis = 1)
    phospho_data = phospho_data.astype({col: 'float32' for col in phospho_data.columns})
    print('   Loaded phosphoproteomics...')
          
    
    transcript_data = cancer.get_transcriptomics(source = "washu")
    transcript_data = transcript_data.droplevel("Database_ID", axis = 1)
    transcript_data = transcript_data.astype({col: 'float32' for col in transcript_data.columns})
    print('   Loaded transcriptomics...')

    
    datasets[name] = {
                        "gen" : genotype_data,
                        "prot": proteomic_data,
                        "trans": transcript_data,
                        "phos" : phospho_data
                         }    

# # Part 2: Select mutations
# Determine which cancers have mutation types with rates above 20%.
# To facilitate that, the mutation types from the get_genotype_all_vars must be grouped into functional types.
mutation_types = ["Deletion", "Amplification", "Truncation", "Indel", "Missense"]
# Grouping the mutation types for each cancer
for name, dataset in datasets.items():
    genotypes = dataset["gen"]
    print(genotypes)
    genotypes.replace(["Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del"], "Truncation", inplace = True)
    genotypes.replace(["In_Frame_Del", "In_Frame_Ins"], "Indel", inplace = True)
    genotypes.replace(["Missense_Mutation", "Wildtype_Tumor"], ["Missense", "Wildtype"], inplace = True)
    genotypes = genotypes.loc[~genotypes["Mutation"].isin(["Intron","Silent"])]
    datasets[name]["gen"] = genotypes
# Create a bargraph displaying the mutation types
plt.rcParams["figure.figsize"] = [15,21]
mutation_data = pd.DataFrame()
for cancer, dataset in datasets.items():
    mutation_data[cancer]= datasets[cancer]["gen"].groupby("Mutation").count()["Location"]
mutation_data = mutation_data.fillna(0)
mutation_data = mutation_data.transpose()
mutation_data.drop(columns="Wildtype").plot.bar(rot=0, width=0.8)

plt.xlabel("Cancer Type", fontsize=16)
plt.ylabel("Count of Mutation", fontsize=16)
plt.title("Count of Mutation Type by Cancer", fontsize=20)
# Determine which cancers have mutation rates over 20%
variants_to_analyze = {}
for name, dataset in datasets.items():
    mut_type_freq = dataset["gen"]["Mutation"].value_counts() / len(dataset["gen"])
    
    mut_type_freq = mut_type_freq.loc[mut_type_freq >= 0.2]
    variants_to_analyze[name] = set(mut_type_freq.index.tolist())
    variants_to_analyze[name].discard("Wildtype")
print(str(variants_to_analyze))
mutation_data = mutation_data.apply(lambda row: row / row.sum(), axis=1)
ax = mutation_data.drop(columns="Wildtype").plot.bar(rot=0, width=0.8)
plt.axhline(y=0.2, xmin=0, xmax=1, linestyle="dotted", color="grey")

plt.xlabel("Cancer Type", fontsize=16)
plt.ylabel("Percent of Samples Containing Mutation Type", fontsize=16)
plt.title("Percent of Samples Contaning Mutation Type by Cancer", fontsize=20)
ax.set_yticklabels(["{:,.0%}".format(val) for val in ax.get_yticks()])
ax.get_yticklabels()[2].set_weight("bold")
plt.show()
# plt.figure1_plot_mutations(datasets.values(), list(datasets.keys()), save_to_path = "Test_Figure.png")


# # Part 3: Cis effects
# Determines which genes that, given a mutation, will change the expression in their own proteomics, transcriptomics, or phosphoproteomics.This analysis is repeated for every cancer, and every mutation type.
# 
# Writes out a file containing the infomation for each mutation and cancer, and creates boxplot summarizing the data
'''
Parameters
pval_symbol: String. The symbol associated with the significance of the pvalue (NS, *, **, or ***)
x1, x2: Float. x coordinates for the start and end of the annotation line.
line_start: Float. y coordinate for the horizontal line (Default = .05).
line_height: Float. The height of the vertical lines at the start and end of the horizontal line (Default = .05).

Returns a formated annotation line showing the significance symbol of each pvalue.
'''
def format_pval_annotation(pval_symbol, x1, x2, line_start = .05, line_height=.05):
    # for manual adjustment to pval annotations
    
    y, h = line_start, line_height
    plt.plot([x1, x1, x2, x2], #draw horizontal line
             [y, y+h, y+h, y], #vertical line
             lw=1.5, color= '.3')
    plt.text((x1+x2)*.5, # half between x coord
             y+h, pval_symbol, horizontalalignment='center', verticalalignment='bottom', color = "black")
def generate_cis_plot(boxplot_data, omics_type, pval_df):
    plt.rcParams['figure.figsize'] = (15,10) # size of plot
#     if len(gene_list) > 15:
#         plt.rcParams['figure.figsize']=(20,10)
    sns.set(font_scale = 1.2, style = "white")

    for site in boxplot_data["Site"].unique():
        mutations = sorted(boxplot_data["Mutation"].unique(), reverse = True)
        boxplot = sns.boxplot(x = "cancer_type", y = "P_Value", hue = "Mutation", showfliers = False, 
                              data = boxplot_data[boxplot_data["Site"] == site], hue_order = mutations)
        boxplot = sns.stripplot(x = "cancer_type", y = "P_Value", hue = "Mutation", jitter = True, hue_order = mutations,
                               color = ".3", dodge = True, data = boxplot_data[boxplot_data["Site"] == site])
#         gene_site = site if site == target_gene else target_gene+'.'+site
        gene_site = site
        boxplot.set(xlabel = "", ylabel = f"{gene_site}_{omics_key[omics_type]}")
        boxplot.set(title = f"Pancan Analysis of {gene_site} {omics_key[omics_type].capitalize()} Mutations")
#         if len(gene_list) > 15:
#             boxplot.set_xticklabels(boxplot.get_xticklabels(),rotation=50)
        
        num = len(boxplot_data["Mutation"].unique())
        handles, labels = boxplot.get_legend_handles_labels()
        legend = plt.legend(handles[0:num], labels[0:num], loc = "lower right")
        
        ## CREATE PVAL ANNOTATIONS
        symbols = {}
        start = -0.3
        cancer_list = boxplot.get_xticklabels() # get a list of cancer types actually plotted 
        first_mutation = mutations[1]
        bottom,top = plt.ylim()
        annot_size = top + 0.5
        annot_scale = (top - bottom) * 0.05 # sets annotation height between mutations based on y-axis scale
        for cancer in cancer_list:
            cancer = cancer.get_text()
            mutation_list = sorted(boxplot_data.loc[boxplot_data["cancer_type"] == cancer]["Mutation"].unique().tolist(),
                                   reverse = True)
            mutation_list.remove("Wildtype")
            offset = 0.4
            annot_offset = annot_scale
            for mutation in mutation_list:
                pval = pval_df.loc[(pval_df["Cancer"] == cancer)&(pval_df["Mutation"] == mutation)&
                                   (pval_df["Site"] == site)]['P_Value']
                if pval.tolist():
                    pval = pval.tolist()[0]
                    if pval  <= 0.001:
                         symbols[cancer] = '***'
                    elif pval  <= 0.01:
                        symbols[cancer] = '**'  
                    elif pval <= 0.05:
                        symbols[cancer] = '*'
                    else:
                        symbols[cancer] = 'ns'
                    if mutation == first_mutation:
                        format_pval_annotation(symbols[cancer], start, start + 0.2, annot_size)
                    else:
                        format_pval_annotation(symbols[cancer], start, start + offset, annot_size + annot_offset)
                        offset += 0.2
                        annot_offset += annot_scale
            start += 1
        plt.show()
        plt.clf()
        plt.close()
            
#         fig = boxplot.get_figure()
#         fig.savefig(f"Cis_Effect_output/Test_{gene_site}_{omics_key[omics_type]}.png")
#         fig.clf()
def analyze_cis_mutation(omics_type, genes_to_analyze = [target_gene], write_out = True, file_name = None, output_boxplot = True, write_insignificant=False):
    """Analyzes the cis effects for a given mutation type across all cancers.
    
    Parameters:
    omics_type (str): The omics to analyze. Possible values are in list ["prot", "trans", "phos"]
    genes_to_analyze (list): A list of strings representing genes to analyze. For cis effects, this is always [target_gene]
    write_out (bool, default=True): Whether to write out the data to a separate file.
    output_boxplot (bool, default=True): Outputs the data as a boxplot given by the generate_cis_plot function
    write_insignificant (bool, default=True): Returns data only when the p-value <= 0.05. If write_out, it only writes out the data under that condition.
    
    Returns a dataframe containing output of scipy.ttest_ind() between "Wildtype" and mutation_type tumors.
    """
    if write_out:# Sets file name and clears target files
        import os
        if file_name is None:
            effect_type = "Cis" if genes_to_analyze == [target_gene] else "Trans"
            file_name = f"{effect_type}_Effect_output/{omics_key[omics_type].capitalize()}_Analysis.txt" 
            if '/' in file_name:
                directory = file_name[:file_name.rfind('/')]
                if not os.path.exists(directory):
                    os.makedirs(directory)
        with open(file_name, 'w') as out_file: pass
                
           
    result = pd.DataFrame()
    for cancer_type, dataset in datasets.items():
        ## PREPARE DATA ##
        genotype_and_omic = dataset[omics_type].copy()
        if isinstance(genotype_and_omic.columns, pd.MultiIndex):
            genotype_and_omic.columns = ['_'.join(col) for col in genotype_and_omic.columns.values]
        if genes_to_analyze is not None:
            genotype_and_omic = genotype_and_omic.filter(regex='(' + '|'.join(genes_to_analyze) + ')(_.*)?$', axis=1)
        else:
            genotype_and_omic.drop(target_gene, axis = 1, inplace=True)
        genotype_and_omic["Mutation"] = dataset["gen"]["Mutation"].copy()
        genotype_and_omic = genotype_and_omic.loc[:,~genotype_and_omic.columns.duplicated()]
        
        for mutation_type in mutation_types:
            if mutation_type not in variants_to_analyze[cancer_type]: 
                continue
                    
            print(f"   {cancer_type} {mutation_type}s") #debug temp
            target_mut_data = genotype_and_omic.loc[genotype_and_omic["Mutation"].isin(["Wildtype", mutation_type])]
            ## ANALYZE P-VALUES ##
            # DIFF b/w trans effects: loops through each column to avoid correction
            return genotype_and_omic
            for site in target_mut_data.drop("Mutation", axis=1).columns:
                mutation_analysis = ut.wrap_ttest(target_mut_data, label_column="Mutation",
                                    correction_method = "fdr_bh", return_all=True)
                mutation_analysis.rename({'Comparison':'Site'}, axis = 1, inplace = True)
                mutation_analysis['Cancer'], mutation_analysis['Mutation'] = cancer_type, mutation_type
            
                result = pd.concat([result, mutation_analysis])

            #END DIFF
            
                ## GATHER DATA FOR BOXPLOT ##
                if output_boxplot:
                    new_data = genotype_and_omic[(genotype_and_omic["Mutation"] == mutation_type)|
                                                (genotype_and_omic["Mutation"] == "Wildtype")]
                    new_data["Site"], new_data["cancer_type"] = site, cancer_type
                    new_data.rename({site:"P_Value"}, axis=1, inplace=True)
                    try:
                        boxplot_data = pd.concat([boxplot_data, new_data], axis=0)
                    except NameError:
                        boxplot_data = new_data
                        
    result.drop_duplicates(inplace=True)

            ## WRITE OUT RESULT ##
    if write_out:
        with open(file_name, 'a') as out_file:
            output = result.copy()
            if not write_insignificant:
                output = output.loc[output['P_Value'] <= 0.05]
            
            for cancer in output["Cancer"].unique():
                for mutation in output["Mutation"].unique():
                    mut_data = output.loc[(output["Cancer"] == cancer)&(output["Mutation"] == mutation)]
                    out_file.write(f"{cancer.upper()} {mutation.upper()}: ")
                    out_file.write(f"{len(mut_data.loc[mut_data['P_Value'] <= 0.05])} significant entries\n")
                    for row in mut_data.iterrows():
                        data = row[1]
                        out_file.write(f"   {data['Site']}:".ljust(15) + f"{data['P_Value']}\n")
                    out_file.write('\n')

    ## CREATE BOXPLOT ##
    if output_boxplot:
        generate_cis_plot(boxplot_data, omics_type, result) 

    result.reset_index(drop=True, inplace=True)
    result["Omics"] = omics_type
    return result
cis_effects = pd.concat([analyze_cis_mutation(omics_type, genes_to_analyze=[target_gene]) for omics_type in ["prot", "trans", "phos"]])
cis_effects.reset_index(drop=True, inplace=True)
cis_effects


# # Part 4: Trans Effects
# Determines which genes that, given a mutation, will change the expression in another protein's proteomics, transcriptomics, or phosphoproteomics.This analysis is repeated for every cancer, and every mutation type.
# 
# To reduce the computational intensity and the corrections needed, the scope will be restricted to genes that share specific pathways with the target gene. These pathways are hypothesized to be significant in the gene's function, and **must be specificed beforehand**, either at the beginning of this notebook or in the following cells.
# 
# The analysis will then be repeated for all genes in the dataset. As the number of comparisons is greater, the correction applied will be more substantial.
# 
# Writes out a file containing the infomation for each mutation and cancer, and creates boxplot summarizing the data
# if path_indices == []:
#     print()
interacting_proteins = ut.get_proteins_in_pathways(paths["pathway"].loc[path_indices], database = 'wikipathways')
interacting_proteins = set(interacting_proteins["member"])

# genes_measured = [gene for dataset in datasets.values() for gene in dataset["prot"].columns]
genes_measured = set([gene for gene in dataset["prot"].columns for dataset in datasets.values()])
genes_to_analyze = sorted([protein for protein in interacting_proteins if protein in genes_measured])
no_analysis = sorted([protein for protein in interacting_proteins if protein not in genes_measured])# debug temp

if target_gene in genes_to_analyze:
    genes_to_analyze.remove(target_gene)
print(f"Interacting genes: {len(interacting_proteins)}")
print(f"Not measured in cptac: {len(no_analysis)}") #debug temp
print(f"Genes to analyze: {len(genes_to_analyze)}")
def analyze_trans_mutation(omics_type, genes_to_analyze = None, write_out = True, file_name = None, output_boxplot = True, write_insignificant=False):
    """Analyzes the cis effects for a given mutation type across all cancers.
    
    Parameters:
    omics_type (str): The omics to analyze. Possible values are in list ["prot", "trans", "phos"]
    genes_to_analyze (list, default=None): The list of genes to consider for trans effects. Defaults to all genes in dataset, excluding target_gene
    write_out (bool, default=True): Whether to write out the data to a separate file.
    output_boxplot (bool, default=True): Outputs the data as a boxplot given by the generate_trans_plot function
    write_insignificant (bool, default=True): Returns data only when the p-value <= 0.05. If write_out, it only writes out the data under that condition.
    
    Returns the output of scipy.ttest_ind() between "Wildtype" and mutation_type tumors. If analyzing phosphoproteomics, returns a list of such outputs
    """
    if write_out:# Sets file name and clears target files
        import os
        if file_name is None:
            effect_type = "Cis" if genes_to_analyze == [target_gene] else "Trans"
            file_name = f"{effect_type}_Effect_output/{omics_key[omics_type].capitalize()}_Analysis_trans.txt"
        if '/' in file_name:
            directory = file_name[:file_name.rfind('/')]
            if not os.path.exists(directory):
                os.makedirs(directory)
        with open(file_name, 'w') as out_file: pass
                
                
    for cancer_type, dataset in datasets.items():
        ## PREPARE DATA ##
        genotype_and_omic = dataset[omics_type].copy()
        if isinstance(genotype_and_omic.columns, pd.MultiIndex):
            genotype_and_omic.columns = ['_'.join(col) for col in genotype_and_omic.columns.values]
        if genes_to_analyze is not None:
            genotype_and_omic = genotype_and_omic.filter(regex='(' + '|'.join(genes_to_analyze) + ')(_.*)?$', axis=1)
        else:
            genotype_and_omic.drop(genotype_and_omic.filter(regex = 'Test').columns, axis = 1, inplace = True)
        genotype_and_omic["Mutation"] = dataset["gen"]["Mutation"].copy()
        genotype_and_omic = genotype_and_omic.loc[:,~genotype_and_omic.columns.duplicated()]
        
        for mutation_type in mutation_types:
            if mutation_type not in variants_to_analyze[cancer_type]: 
                continue
            print(f"   {cancer_type} {mutation_type}s")

            target_mut_data = genotype_and_omic.loc[genotype_and_omic["Mutation"].isin(["Wildtype", mutation_type])]
   
            ## ANALYZE P-VALUES ##
            mutation_analysis = ut.wrap_ttest(target_mut_data, label_column="Mutation",
                                correction_method = "fdr_bh", return_all=True)
            mutation_analysis.rename({'Comparison':'Site'}, axis = 1, inplace = True)
            mutation_analysis['Cancer'], mutation_analysis['Mutation'] = cancer_type, mutation_type
            
            try:
                result = pd.concat([result, mutation_analysis])
            except NameError:
                result = mutation_analysis

            ## GATHER DATA FOR BOXPLOT ##
            if output_boxplot:
                new_data = genotype_and_omic[(genotype_and_omic["Mutation"] == mutation_type)|
                                            (genotype_and_omic["Mutation"] == "Wildtype")]
                new_data["cancer_type"] = cancer_type
                try:
                    boxplot_data = pd.concat([boxplot_data, new_data], axis=0)
                except NameError:
                    boxplot_data = new_data
    
    result.drop_duplicates(inplace=True)

    ## WRITE OUT RESULT ##
    if write_out:
        with open(file_name, 'a') as out_file:
            output = result.copy()
            for cancer in output["Cancer"].unique():
                for mutation in output["Mutation"].unique():
                    mut_data = output.loc[(output["Cancer"] == cancer)&(output["Mutation"] == mutation)]
                    if len(mut_data) == 0: 
                        continue
                    if not write_insignificant:
                        mut_data = mut_data.loc[mut_data['P_Value'] <= 0.05]
                    out_file.write(f"{cancer.upper()} {mutation.upper()}: ")
                    out_file.write(f"{len(mut_data.loc[mut_data['P_Value'] <= 0.05])} significant entries\n")
                    for row in mut_data.iterrows():
                        data = row[1]
                        out_file.write(f"   {data['Site']}:".ljust(15) + f"{data['P_Value']}\n")
                    out_file.write('\n')

    ## CREATE BOXPLOT ##
#     if output_boxplot:
#         pval_df = pd.DataFrame(pval_dict)
#         generate_cis_plot(boxplot_data, omics_type, pval_df) 
        
    return result
trans_effects = None
for omics_type in ["prot", "trans", "phos"]:
    print(omics_type.upper() + ':')
    data = analyze_trans_mutation(omics_type, genes_to_analyze=genes_to_analyze, write_insignificant=False)
    data["Omics"] = omics_type
    if trans_effects is not None:
        trans_effects = pd.concat([trans_effects, data])
    else:
        trans_effects = data
    print()
trans_effects.reset_index(drop=True, inplace=True)
trans_effects
trans_effects_all = None
for omics_type in ["prot", "trans", "phos"]:
    print(omics_type.upper() + ':')
    data = analyze_trans_mutation(omics_type, genes_to_analyze=None, write_insignificant=False, file_name=f"Trans_Effect_output/{omics_key[omics_type].capitalize()}_Analysis_trans_ALL")
    data["Omics"] = omics_type
    if trans_effects_all is not None:
        trans_effects_all = pd.concat([trans_effects_all, data])
    else:
        trans_effects_all = data
    print()
trans_effects_all.reset_index(drop=True, inplace=True)
trans_effects_all


# # Part 5: Discussion
# Looks for ptaterns in the analysis that will indicate meaningful conclusions

# ## 5.1: Overlapping Genes
# A pan-cancer effect is, by definition, one that is present across many types of cancers. We will look at the most significant effects in each cancer for each mutation, and find the overlap with all other cancers and mutations.
def find_common_effects(effects, omics_type, slice_by="percentile", slice_val=25, file_name = None):
    most_significant = []
    num_entries = 0
    for cancer in effects["Cancer"].unique():
        for mutation in effects["Mutation"].unique():
            data = effects.loc[(effects["Cancer"] == cancer) & (effects["Mutation"] == mutation) & (effects["Omics"] == omics_type)]
            # This slices the data based on what input parameter was put in
            if slice_by == "p_value":
                most_significant.append(data[data["P_Value"] < slice_val])
            elif slice_by == "rank":
                most_significant.append(data.head(slice_val))
            elif slice_by=="percentile":
                most_significant.append(data.head(int(len(data) * slice_val / 100)))

            #Not all combinations of cancer+mutation exist in dataframe; only counting those that do
            if len(data) > 0:
                num_entries += 1
    
    #Transforming into dataframe
    most_significant = pd.concat(most_significant)
    most_significant = most_significant["Site"].value_counts()
    most_significant = most_significant[most_significant > 1]
    
    #Everything below is just for outputting data
    import os
    if file_name is None:
        effect_type = "Cis" if genes_to_analyze == [target_gene] else "Trans"
        file_name = f"{effect_type}_Effect_output/{omics_key[omics_type].capitalize()}_Common_Effects_by_{slice_by}.txt"
    if '/' in file_name:
        directory = file_name[:file_name.rfind('/')]
        if not os.path.exists(directory):
            os.makedirs(directory)
    #Writing out file
    with open(file_name, 'w') as out_file:
        out_file.write(f"{target_gene} {omics_key[omics_type].capitalize()}\n")
        out_file.write(f"Total number of cancer+mutation types measured: {num_entries}\n")
        if slice_by=="p_value":
            out_file.write(f"Taking the genes where p < {slice_val}\n")
        elif slice_by=="rank":
            out_file.write(f"Taking the {slice_val} most significant genes from each\n")
        elif slice_by=="percentile":
            out_file.write(f"Taking the top {slice_val} percent most significant genes from each\n")
        out_file.write(most_significant.to_string())
    
    #Creating graph
    data = most_significant.copy()
    if len(data) > 40:
        data=data[:40]
    if len(data) <= 5:
        return most_significant # Graphs break when only one or fewer data points
    plt.rcParams["figure.figsize"]=[12,len(data)/1.5]
    data.plot.barh(legend=None)
    plt.gca().invert_yaxis()
    plt.xlim(0,num_entries)
    for idx, val in enumerate(data.values):
        plt.text(val-plt.xlim()[1]*.035, idx+0.15, str(val).rjust(2), fontsize=14, fontweight="bold", color="white")
    plt.suptitle(f"    Common {omics_key[omics_type].title()[:-1]} Trans Effects", fontsize=20, fontweight="bold", y=1.2-0.2*(log(len(data)/5))**.5)
    if slice_by=="p_value":
        chart_title=f"trans effects where p < {slice_val}"
    elif slice_by=="rank":
        chart_title=f"the top {slice_val} most significant trans effects"
    elif slice_by=="percentile":
        chart_title=f"the top {slice_val}% most significant trans effects"
    plt.title(f"(Considering {chart_title})", fontsize=18)
    plt.xlabel("Number of Cancer+Mutation Combinations in which Trans Effect was Observed", fontweight="bold", fontsize=14, y=1)
    plt.savefig(file_name.replace(".txt", ".png"), dpi=75)
    plt.show()
    
    return most_significant

for omics in ["prot", "trans", "phos"]:
    find_common_effects(trans_effects, omics, slice_by="p_value", slice_val=0.05)
    find_common_effects(trans_effects, omics, slice_by="percentile", slice_val=25)
    find_common_effects(trans_effects_all, omics, file_name=f"Trans_Effect_output/{omics_key[omics].capitalize()}_Common_Effects_All_by_p_value.txt", slice_by="p_value", slice_val=0.05)
    find_common_effects(trans_effects_all, omics, file_name=f"Trans_Effect_output/{omics_key[omics].capitalize()}_Common_Effects_All_by_percentile.txt", slice_by="percentile", slice_val=25)

