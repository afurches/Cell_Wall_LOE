#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#script to calculate intersect scores with a weighted GO-term functional network
#created by Anna Furches, as used in the manuscript:
#   Furches, A., Kainer, D., Weighill, D., Large, A., Jones, P., Walker, A.M., Romero, J., Gazolla, JGFM,
#   Joubert, W., Shah, M., Streich, J., Ranjan, P., Schmutz, J., Sreedasyam, A., Macaya-Sanz, D., Zhao, N.,
#   Martin, M.Z., Rao, X., Dixon, R.A., DiFazio, S., Tschaplinski, T.J., Chen, J-C., Tuskan, G.A., Jacobson, D.
#   Finding New Cell Wall Regulatory Genes in Populus trichocarpa Using Multiple Lines of Evidence.
#   Front. Plant Sci. (in review).




# USAGE: python3 calc_intersect_scores.py <json_of_GOterm_network> <filepaths_of_networks_to_be_intersected.txt> <desired_output_filepath.json>

import sys
import json



#json file of GO-term functional network converted to dictionary format
GOjson = sys.argv[1]

#list of network filepaths to be intersected with the GO-term network
LOE_randlist = sys.argv[2]

#desired filepath of output file
int_scores_file = sys.argv[3]



#create list of networks to be intersect with the GO-term network
LOEfiles = []
with open(LOE_randlist,'r') as randlist:
    for line in randlist:
        LOEfiles.append(line.strip('\n'))

print('Number of files to be intersected: ' + str(len(LOEfiles)))



#load the GO-term json into memory
with open(GOjson,encoding='utf-8') as gd:
    go_dict = json.loads(gd.read())

print('Length of GO dictionary: ' + str(len(go_dict)))



# intersect and record all scores for each file in network list
all_int_scores_dict = {}
for z in LOEfiles:
    with open(z,'r') as temploe:
        print('Intersecting LOE file: ' + z.split('/')[-1])

        #make network dict
        loe_dict = {}
        for line in temploe:
            nodeA = line.strip('\n').split('\t')[0][:-5]  #remove .v3.0
            edge = line.strip('\n').split('\t')[1]
            nodeB = line.strip('\n').split('\t')[2][:-5]
            #sort edges to avoid writing reciprocal edges (avoid dups)
            nodes = sorted([nodeA,nodeB])  
            node1 = nodes[0]
            node2 = nodes[1]
            if node1 != node2:  #don't write self-loops
                if node1 not in loe_dict:
                    loe_dict[node1] = {}
                if node2 not in loe_dict[node1]:
                    loe_dict[node1][node2] = []
                if edge not in loe_dict[node1][node2]:  #only write node1-node2 edge relationships 1x per edge
                    loe_dict[node1][node2].append(edge)
        #write out network dict
        dict_path = ('/'.join(z.split('/')[:-1]) + '/' + z.split('/')[-1][:-4] + '_dict.json')
        print('Writing out LOE dict: ' + dict_path.split('/')[-1])
        with open(dict_path,'w') as outdict_loe1:
            outdict_loe1.write(json.dumps(loe_dict))
        
        
        #create dict of edge multipliers
        loe_totals = {}
        for i in loe_dict.keys():
            for j in loe_dict[i]:
                if i not in loe_totals.keys():
                    loe_totals[i] = {}
                if j not in loe_totals[i]:
                    loe_totals[i][j] = len(loe_dict[i][j])
        #write out totals dict
        totals_path = ('/'.join(z.split('/')[:-1]) + '/' + z.split('/')[-1][:-4] + '_dict-edgecounts.json')
        print('Writing LOE edge counts to dict: ' + totals_path.split('/')[-1])
        with open(totals_path,'w') as outdict_loe2:
            outdict_loe2.write(json.dumps(loe_totals))
            
        
        #calculate GO-v-network intersect score for each node-node pair
        int_dict = {}
        for i in loe_totals.keys():
            for j in loe_totals[i]:
                try:
                    if go_dict[i][j]:
                        if i not in int_dict.keys():
                            int_dict[i] = {}
                        if j not in int_dict[i]:
                            int_dict[i][j] = go_dict[i][j] * loe_totals[i][j]
                except KeyError:
                    pass
        #write out
        int_dict_path = ('/'.join(z.split('/')[:-1]) + '/' + z.split('/')[-1][:-4] + '_dict-edgeintersects.json')
        print('Writing LOE edge intersect scores to dict: ' + int_dict_path.split('/')[-1])
        with open(int_dict_path,'w') as outdict_loe3:
             outdict_loe3.write(json.dumps(int_dict))
 
 
        #calculate total score for this network file
        total_int_score = 0
        for k in int_dict.keys():
            for m in int_dict[k]:
                total_int_score += int_dict[k][m]
        #write to dict
        print('Total intersect score: ' + str(total_int_score))
        tempkey = z.split('/')[-1]  #LOE file name
        all_int_scores_dict[tempkey] = total_int_score



#write intersect scores dict to file
print('Intersections complete. Writing all intersect scores to file.')
with open(int_scores_file,'w') as allscores:
     allscores.write(json.dumps(all_int_scores_dict))
