#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:34:41 2019
Edited December 11 2019

calculate_LOE_scores.py

A script to calculate multiple Lines of Evidence (LOE) scores for every gene in the genome.
Originally developed in Perl by Deborah Weighill. 
    Weighill, D., Jones, P., Shah, M., Ranjan, P., Muchero, W., Schmutz, J., et al. (2018). 
    Pleiotropic and Epistatic Network-Based Discovery: Integrated Networks for Target Gene Discovery. 
    Front. Energy Res. 6, 267997. doi:10.3389/fenrg.2018.00030.

Adapted for Python and public release by Anna Furches.
    Furches, A., Kainer, D., Weighill, D., Large, A., Jones, P., Walker, A.M., Romero, J., Gazolla, JGFM, 
    Joubert, W., Shah, M., Streich, J., Ranjan, P., Schmutz, J., Sreedasyam, A., Macaya-Sanz, D., Zhao, N., 
    Martin, M.Z., Rao, X., Dixon, R.A., DiFazio, S., Tschaplinski, T.J., Chen, J-C., Tuskan, G.A., Jacobson, D.
    Finding New Cell Wall Regulatory Genes in Populus trichocarpa Using Multiple Lines of Evidence.
    Front. Plant Sci. (in review).



WARNINGS: 
This script was written on a Mac/Unix system; behavior with other systems has not been investigated.
This script uses the Python libraries sys, collections, and pandas. Otherwise, dependencies have not been investigated.



USAGE:
> python3 calculate_LOE_scores.py <input_file_list>
    - "input_file_list" is a text file of paths to the following input files: (1) genes to be scored, (2) anchor genes and phenotypes, (3...n) networks from which LOE scores are to be calculated.
    - each file path should be on a separate line

    

INPUT
The following file naming conventions are required, but files may be listed in any order:
    - GENES TO BE SCORED: path/to/file/genes_<other_text_of_arbitrary_length>.txt
    - ANCHORS: path/to/file/anchors_<other_text_of_arbitrary_length>.txt
    - NETWORK (each on a separate line): path/to/file/<datatype>_<other_text_of_arbitrary_length>.sif
        - "datatype" should be replaced with a single word to identify the data layer in the scores table; 
        - EXAMPLES: "coex", "cometh", "snpcor","metabGWAS"
    - an underscore "_" must be present after the first word in the file name for proper delimitation
       
File formatting requirements:
    - files must be tab-delimited
    - any file ending may be used (i.e., "txt","sif","tsv"), python ignores these
    - networks must be in >= 3-column format with 1 edge per line, with nodes and edges in the following order: 
        node1     edge     node2
    - if >3 columns are included in network files, only the first 3 columns will be used

Gene and phenotype id formatting recommendations:
    - make sure gene and phenotype id's are consistent across all files



OUTPUT
A tab-delimited table with total breadth score, total depth score, and layer-specific depth scores for 
each gene id in the "GENES TO BE SCORED" file.
    - columns will be labeled as such:
        gene_id     breadth_tot     depth_tot    layername1     layername2 ... layername(n)
    - "layername" columns contain layer-specific depth scores
    - "layername" is replaced by the user specified data type provided in the input network file name
        - EXAMPLE: path/to/file/coex_extrawords_somedate.sif ==> translates to column label "depth_coex"
        - EXAMPLE: path/to/file/metab-GWAS_extrawords_somedate.sif ==> translates to column label "depth_metab-GWAS"
    - output file (scores table) will be named as follows: genes_<other_text_of_arbitrary_length>_LOEscores.txt
    - self loops and reciprocal edges are NOT included in LOE scores 
"""



#LOAD LIBRARIES
import sys
from collections import defaultdict
import pandas as pd



# LOAD INPUT FILES
filepaths = sys.argv[1]

with open(filepaths,'r') as fp:
    filelist = [i.strip('\n') for i in fp]
print('\nINPUT FILES:\n  ' + '\n  '.join(filelist) + '\n')



# PARSE FILE TYPES
filedict = dict.fromkeys(filelist)
for f in filelist:
    filetype = f.split('/')[-1].split('_')[0]
    filedict[f] = filetype
filetype_list = set([filedict[k] for k in filedict.keys()])
print('FILE TYPES:\n  ' + '\n  '.join(filetype_list) + '\n')



# ASSIGN FILE TYPES
layer_list = []
for k in filedict.keys():
    if filedict[k] == 'genes':
        genes = k
    elif filedict[k] == 'anchors':
        anc_file = k
    elif filedict[k] != 'genes' and filedict[k] != 'anchors':
        layer_list.append(k)
print('Genes to score: ' + genes.split('/')[-1])
print('Anchor genes and phenotypes: ' + anc_file.split('/')[-1])
print('Networks to search:\n  ' + '\n  '.join([i.split('/')[-1] for i in layer_list]) + '\n')



# CREATE BREADTH AND DEPTH DICTIONARIES
with open(genes,'r') as gf: 
    gene_list = [line.strip('\n') for line in gf]
print('Number of genes to be scored: ' + str(len(gene_list)) + '\n')
print('Gene list preview:\n  ' + '\n  '.join(gene_list[0:3]) + '\n')

breadth_dict = defaultdict(lambda: defaultdict(int))
for g in gene_list:
    for l in layer_list:
        breadth_dict[g][filedict[l]] = 0
#print(str(breadth_dict))        
print('Breadth dictionary preview:\n' + str(dict(list(breadth_dict.items())[0:1])))

depth_dict = defaultdict(lambda: defaultdict(int))
for g in gene_list:
    for l in layer_list:
        depth_dict[g][filedict[l]] = 0
print('Depth dictionary preview:\n' + str(dict(list(depth_dict.items())[0:1])) + '\n')



# LOAD ANCHOR GENES AND METABS
with open(anc_file,'r') as af:
    anchors = set([line.strip('\n') for line in af])
print('Number of anchor genes and phenos: ' + str(len(anchors)) + '\n')



# CALCULATE LOE SCORES
print('** SCORING LAYERS **')
for layer in layer_list:
    current_layer = filedict[layer]
    print('Scoring ' + current_layer + ' layer.')
    #keep track of edges to avoid scoring reciprocal edges
    recip_dict = defaultdict(lambda: defaultdict(str))
    with open(layer,'r') as lf:
        for line in lf:
            node1 = line.strip('\n').split('\t')[0]
            edge = line.strip('\n').split('\t')[1]
            node2 = line.strip('\n').split('\t')[2]
            nodes = sorted([node1,node2])
            #don't score self loops
            if node1 != node2:
                if node1 or node2 in anchors:
                    if edge not in recip_dict[nodes[0]][nodes[1]]:
                        recip_dict[nodes[0]][nodes[1]] = edge
                        if node1 in anchors and node2 not in anchors:
                            if node2 in gene_list:
                            #score node2
                                if breadth_dict[node2][current_layer] == 0:
                                    breadth_dict[node2][current_layer] += 1
                                depth_dict[node2][current_layer] += 1
                        elif node1 not in anchors and node2 in anchors:
                            if node1 in gene_list:
                            #score node1
                                if breadth_dict[node1][current_layer] == 0:
                                    breadth_dict[node1][current_layer] += 1
                                depth_dict[node1][current_layer] += 1
                        elif node1 in anchors and node2 in anchors:
                            #score both nodes
                            if node1 in gene_list:
                                if breadth_dict[node1][current_layer] == 0:
                                    breadth_dict[node1][current_layer] += 1
                                depth_dict[node1][current_layer] += 1
                            if node2 in gene_list:
                                if breadth_dict[node2][current_layer] == 0:
                                    breadth_dict[node2][current_layer] += 1
                                depth_dict[node2][current_layer] += 1



# CREATE LOE SCORES TABLE
print('Tallying final scores.\n')
            
# convert breadth scores dictionary to pandas dataframe
br_df = pd.DataFrame.from_dict(breadth_dict,orient='index',dtype=int)
br_df['breadth_tot'] = br_df.sum(axis=1)
br_df.reset_index(inplace=True)
br_df.rename(columns={'index':'geneID'},inplace=True)

# convert depth scores dictionary to pandas dataframe
dep_df = pd.DataFrame.from_dict(depth_dict,orient='index',dtype=int)
dep_df['depth_tot'] = dep_df.sum(axis=1)
dep_df.reset_index(inplace=True)
dep_df.rename(columns={'index':'geneID'},inplace=True)

# merge breadth and depth scores dataframes
LOE_scores = br_df[['geneID','breadth_tot']].merge(dep_df,how='outer',on='geneID')
LOE_scores.set_index('geneID',inplace=True)
LOE_scores = LOE_scores.fillna(0).astype(int)
print('Number of genes in scores table: ' + str(len(LOE_scores)))
print('Preview of scores table:\n' + str(LOE_scores.head(3)) + '\n')



# WRITE OUT SCORES FILE
print("Writing to: " + (genes[0:-4] + '_LOEscores.txt'))
LOE_scores.to_csv(genes[0:-4] + '_LOEscores.txt',sep='\t',header=True,index=True)
