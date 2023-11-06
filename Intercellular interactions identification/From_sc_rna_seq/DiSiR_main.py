#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import pandas as pd
import sys
import subprocess
import csv
import numpy as np
###############################################################################
################################ Set user INPUTS ##############################
###############################################################################
# Input directory with spring data:"matrix.mtx, categorical_coloring_data.json, genes.txt":
# 1. matrix.mtx: matrix of single-cell gene expression data where its rows are associated with cells and its columns are associated with genes
# 2. categorical_coloring_data.json: meta data file which includes cell type labels of cells and also other information such as: disease, sub-cell-states and etc.
# 3. genes.txt: gene names of gene expression matrix
indir='Large-scale multicellular modeling of the arthritic joint/From_sc_rna_seq/'


# Input directory with function scripts:
indir_fuctions="Large-scale multicellular modeling of the arthritic joint/From_sc_rna_seq/" 

# Output directory to save DiSiR results:
outdir="Large-scale multicellular modeling of the arthritic joint/From_sc_rna_seq/"

# Comma-separated interaction file: one interaction per line, for example: "IL1A | IL1R1,IL1A | IL1RAP":
subunit_interactions_path=""

# Select metadata track for analysis, default="CellStates":
select_track = "CellStatesID"

# Select groups in the metadata track for exclusion, default="Unclassified, Other":
exclude_tracks = "Unclassified, Other"

# Number of iterations for permutating data in statistical test: 
iteration = 100

# Threshold on number of cells expressed each ligand or receptor per cell type:
threshold_numbers = 0

# Threshold on scaled (max-normalized) average expression of each ligand or receptor within a cell type:
threshold_expressions = 0

# Threshold on p-value for filtering non-significant LR interactions:
threshold = 0.05

###############################################################################
##################################### Read inputs #############################
###############################################################################
subset = True

# Path to input gene expression matrix
scRNA_path=indir + "total_RA_data_mono_fibro_TH1_Bioturing.txt"
#scRNA_path = indir + '/matrix.mtx'

# Path to names of all genes in gene expression matrix
gene_names_all_path = indir + "genes_names_fibro_mono_TH1.txt"
#gene_names_all_path = indir + '/genes.txt'

# Path to metadata
metadata_path = indir + "output_mono_fibro_TH1.json"
#metadata_path = indir + '/categorical_coloring_data.json'

###################################################################################################
###################################### RUN MOD CELLPHONEDB ########################################
###################################################################################################
import sys
sys.path.append(indir_fuctions)
import DiSiR_functions as ds

np.random.seed(0)

scRNA_array, cell_type_labels = ds.filter_out_unannotated_cells(scRNA_path,
                             metadata_path, select_track, 
                             exclude_tracks, subset)

scRNA_array=scRNA_array.apply(pd.to_numeric, errors='coerce', axis=1)
expressions_scRNA, gene_names = ds.gene_expressions_with_selected_genes(scRNA_array,
                                            gene_names_all_path,
                                            subunit_interactions_path)


average_expression_df, number_expression_df, fraction_expression_df, unique_cell_type_labels, cell_type_numbers, cell_type_specific_numbers, cell_type_specific_expressions = ds.calculate_celltype_average_expressions(expressions_scRNA,
                                           gene_names,
                                           cell_type_labels)

##############################################################################
######################### Transfer to rank space #############################
##############################################################################
# average_expression_df = ds.transform_to_rank(average_expression_df)

##############################################################################
interactions_prefiltered, interactions_gene_names, interactions_class_names = ds.calculate_LR_interactions_matrix_prefiltered(average_expression_df,
                            gene_names,
                            unique_cell_type_labels,
                            cell_type_specific_numbers,
                            cell_type_specific_expressions,
                            threshold_numbers,
                            threshold_expressions)


p_values, p_values_gene_names, p_values_celltype_names, expression_celltype_matrix = ds.calculate_LR_interactions_pvalues(interactions_prefiltered,
                                   interactions_gene_names,
                                   interactions_class_names,
                                   expressions_scRNA,
                                   cell_type_numbers,
                                   gene_names,
                                   unique_cell_type_labels,
                                   iteration)


results_table=pd.concat([p_values_gene_names,p_values_celltype_names,expression_celltype_matrix,p_values],axis=1)
p_values.to_csv(r'DISIR_pvalues_fibro_macro_TH1.txt', header=None, index=None, sep=' ', mode='a')
p_values_celltype_names.to_csv(r'DISIR_celltype_names_fibro_macro_TH1.txt', header=None, index=None, sep=' ', mode='a')
expression_celltype_matrix.to_csv(r'DISIR_expression_scores_fibro_macro_TH1.txt', header=None, index=None, sep=' ', mode='a')
p_values_gene_names.to_csv(r'DISIR_gene_names_fibro_macro_TH1.txt', header=None, index=None, sep=' ', mode='a',na_rep='NULL')




subunit_interactions_data = pd.read_csv(subunit_interactions_path, header = None, index_col = None)

if len(subunit_interactions_data.columns) > 1:
    # SUBUNIT INTERACTIONS
    graph_data_nodes, graph_data_links = ds.build_interactions_graph_subunit(p_values_matrix,
                                 p_values_gene_names,
                                 unique_cell_type_labels,
                                 expression_celltype_matrix,
                                 subunit_interactions_data,
                                 threshold)
 
else:
    # NO SUBUNIT INTERACTIONS
    selected_inreactions = subunit_interactions_data
    graph_data_nodes, graph_data_links = ds.build_interactions_graph(p_values_matrix,
                                 p_values_gene_names,
                                 unique_cell_type_labels,
                                 selected_inreactions,
                                 expression_celltype_matrix,
                                 threshold)

heatmap_dict = ds.calculate_celltype_interactions_heatmap_new(graph_data_links,
                                            unique_cell_type_labels,
                                            average_expression_df)

###############################################################################
#################### Define and save outputs and plots #######################
###############################################################################

if not os.path.isdir(f'{outdir}'):
    os.mkdir(outdir)
    #subprocess.call(f'mkdir {outdir}', shell=True)
    
for relation in heatmap_dict:
    
    if not os.path.isdir(f'{outdir}/{relation}'):
        relation_path=os.path.join(outdir,relation)
        os.mkdir(relation_path)
        #subprocess.call(f'mkdir {outdir}/{relation}', shell=True)
    heatmap_dict[f'{relation}']['heatmap'].to_csv(relation_path + "/Heatmap.csv")
    heatmap_dict[f'{relation}']['heatmap_all_interactions'].to_csv(relation_path +'/Heatmap_all_interactions.csv')


if not os.path.isdir(f'{outdir}/graph_data'):
    graph_path=os.path.join(outdir,"graph_data")
    os.mkdir(graph_path)
    #subprocess.call(f'mkdir {outdir}/graph_data', shell=True)
    
    
ds.plot_cell_dist(average_expression_df, 
               number_expression_df,
               fraction_expression_df,
               outdir)


with open(f'{outdir}/graph_data/Links.csv','w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('from', 'to', 'weight', 'name') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(graph_data_links.shape[0]):
         thewriter.writerow([graph_data_links.iloc[row,0], graph_data_links.iloc[row,1], graph_data_links.iloc[row,2], graph_data_links.iloc[row,3]])
print()
             
with open(f'{outdir}/graph_data/Nodes.csv','w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('id', 'name', 'type', 'y', 'x') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(graph_data_nodes.shape[0]):
         thewriter.writerow([graph_data_nodes.iloc[row,0], graph_data_nodes.iloc[row,1], graph_data_nodes.iloc[row,2], graph_data_nodes.iloc[row,3], graph_data_nodes.iloc[row,4]])
print()


###############################################################################
###############################################################################
###############################################################################
