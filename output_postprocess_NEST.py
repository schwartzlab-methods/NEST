import numpy as np
import csv
import pickle
from scipy import sparse
import scipy.io as sio
import scanpy as sc 
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import stlearn as st
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
from typing import List
import qnorm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict
import pandas as pd
import gzip
from kneed import KneeLocator
import copy 
import altairThemes
import altair as alt
import argparse
import gc


##########################################################
data_name = 'PDAC_64630' # 'PDAC_130355_D1' #'PDAC_140694' #'PDAC_130355_B1' #'V1_Human_Lymph_Node_spatial' #'LUAD_GSM5702473_TD1' #'PDAC_64630' #LUAD_GSM5702473_TD1
current_directory = '/cluster/projects/schwartzgroup/fatema/find_ccc/'
##########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset')
    args = parser.parse_args()
    filter_min_cell = 5
    threshold_expression = 98

elif data_name == 'PDAC_130355_B1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/exp2_B1/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_130355_B1', help='The name of dataset')
    args = parser.parse_args()
    filter_min_cell = 5
    threshold_expression = 98.1
    
elif data_name == 'PDAC_130355_A1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/exp2_A1/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_130355_A1', help='The name of dataset')
    args = parser.parse_args()
    filter_min_cell = 5
    threshold_expression = 98.7
	
elif data_name == 'PDAC_130355_D1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/exp1/exp1_D1/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_130355_D1', help='The name of dataset')
    args = parser.parse_args()
    filter_min_cell = 5
    threshold_expression = 98
	
elif data_name == 'PDAC_140694':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/V10M25-60_C1_PDA_140694_Pa_P_Spatial10x/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_140694', help='The name of dataset')
    #parser.add_argument( '--model_name', type=str, default='gat_r1_2attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()
    filter_min_cell = 1
    threshold_expression = 98
####### get the gene id, cell barcode, cell coordinates ######

if data_name == 'LUAD_GSM5702473_TD1':
    # read the mtx file
    temp = sc.read_10x_mtx(args.data_path)
    print(temp)
    #sc.pp.log1p(temp)
    sc.pp.filter_genes(temp, min_cells=filter_min_cell)
    print(temp)
    #sc.pp.highly_variable_genes(temp) #3952
    #temp = temp[:, temp.var['highly_variable']]
    #print(temp)
    
    gene_ids = list(temp.var_names) 
    cell_barcode = np.array(temp.obs.index)
    
    # now read the tissue position file. It has the format: 
    #df = pd.read_csv('/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/spatial/tissue_positions_list.csv', sep=",",header=None)   # read dummy .tsv file into memory
    df = pd.read_csv('/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv', sep=",",header=None)   # read dummy .tsv file into memory
    tissue_position = df.values
    barcode_vs_xy = dict() # record the x and y coord for each spot
    for i in range (0, tissue_position.shape[0]):
        barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][5], tissue_position[i][4]] #for some weird reason, in the .h5 format, the x and y are swapped
        #barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][4], tissue_position[i][5]] 
    
    coordinates = np.zeros((cell_barcode.shape[0], 2)) # insert the coordinates in the order of cell_barcodes
    for i in range (0, cell_barcode.shape[0]):
        coordinates[i,0] = barcode_vs_xy[cell_barcode[i]][0]
        coordinates[i,1] = barcode_vs_xy[cell_barcode[i]][1]
    

else:
    adata_h5 = st.Read10X(path=args.data_path, count_file='filtered_feature_bc_matrix.h5') #count_file=args.data_name+'_filtered_feature_bc_matrix.h5' )
    print(adata_h5)
    sc.pp.filter_genes(adata_h5, min_cells=filter_min_cell)
    print(adata_h5)
    gene_ids = list(adata_h5.var_names)
    coordinates = adata_h5.obsm['spatial']
    cell_barcode = np.array(adata_h5.obs.index)
    

##################### get metadata: barcode_info ###################################
 

with gzip.open(args.metadata_to + args.data_name +'_barcode_info', 'rb') as fp:  #b, a:[0:5]   _filtered
    barcode_info = pickle.load(fp)




############# load output graph #################################################

filename = ["r1_", "r2_", "r3_", "r4_", "r5_", "r6_", "r7_", "r8_", "r9_", "r10_"]
total_runs = 5
start_index = 0 # 0 #5 if pdac 64630

distribution_rank = []
all_edge_sorted_by_rank = []
for layer in range (0, 2):
    distribution_rank.append([])
    all_edge_sorted_by_rank.append([])

layer = -1
percentage_value = 0

for l in [2,3]: #, 3]: # 2 = layer 2, 3 = layer 1 
    layer = layer + 1
    csv_record_dict = defaultdict(list)
    for run_time in range (start_index, start_index+total_runs):
        gc.collect()
        #run_time = 2
        run = run_time
        print('run %d'%run)

        attention_scores = []
        for i in range (0, datapoint_size):
            attention_scores.append([])   
            for j in range (0, datapoint_size):	
                attention_scores[i].append([])   
                attention_scores[i][j] = []

        distribution = []
        ##############################################
        #X_attention_filename = args.embedding_data_path + args.data_name + '/' + args.data_name + '_cellchat_nichenet_threshold_distance_bothAbove_cell98th_tanh_3dim_'+filename[run_time]+'attention_l1.npy'
        X_attention_bundle = np.load(X_attention_filename, allow_pickle=True) #_withFeature
        for index in range (0, X_attention_bundle[0].shape[1]):
            i = X_attention_bundle[0][0][index]
            j = X_attention_bundle[0][1][index]
            distribution.append(X_attention_bundle[l][index][0])

        ################# scaling the attention scores so that layer 1 and 2 will be comparable ##############################        
        min_attention_score = 1000
        max_value = np.max(distribution)
        min_value = np.min(distribution)
        distribution = []
        for index in range (0, X_attention_bundle[0].shape[1]):
            i = X_attention_bundle[0][0][index]
            j = X_attention_bundle[0][1][index]
            scaled_score = (X_attention_bundle[l][index][0]-min_value)/(max_value-min_value)
            attention_scores[i][j].append(scaled_score) #X_attention_bundle[2][index][0]
            if min_attention_score > scaled_score:
                min_attention_score = scaled_score
            distribution.append(scaled_score)
            
            
        if min_attention_score<0:
            min_attention_score = -min_attention_score
        else: 
            min_attention_score = 0
        
        print('min attention score %g, total edges %d'%(min_attention_score, len(distribution)))
        
        ccc_index_dict = dict()
        threshold_down =  np.percentile(sorted(distribution), 0)
        threshold_up =  np.percentile(sorted(distribution), 100)
        connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))
        for j in range (0, datapoint_size):
            #threshold =  np.percentile(sorted(attention_scores[:,j]), 97) #
            for i in range (0, datapoint_size):
                atn_score_list = attention_scores[i][j]
                #print(len(atn_score_list))
                #s = min(0,len(atn_score_list)-1)
                for k in range (0, len(atn_score_list)):
                    if attention_scores[i][j][k] >= threshold_down and attention_scores[i][j][k] <= threshold_up: #np.percentile(sorted(distribution), 50):
                        connecting_edges[i][j] = 1
                        ccc_index_dict[i] = ''
                        ccc_index_dict[j] = ''
    
    
    
        graph = csr_matrix(connecting_edges)
        n_components, labels = connected_components(csgraph=graph,directed=True, connection = 'weak',  return_labels=True) #
        print('number of component %d'%n_components)
    
        count_points_component = np.zeros((n_components))
        for i in range (0, len(labels)):
             count_points_component[labels[i]] = count_points_component[labels[i]] + 1
    
        print(count_points_component)
    
        id_label = 2 # initially all are zero. =1 those who have self edge but above threshold. >= 2 who belong to some component
        index_dict = dict()
        for i in range (0, count_points_component.shape[0]):
            if count_points_component[i]>1:
                index_dict[i] = id_label
                id_label = id_label+1
    
        print(id_label)
    
    
        for i in range (0, len(barcode_info)):
        #    if barcode_info[i][0] in barcode_label:
            if count_points_component[labels[i]] > 1:
                barcode_info[i][3] = index_dict[labels[i]] #2
            elif connecting_edges[i][i] == 1 and len(lig_rec_dict[i][i])>0: 
                barcode_info[i][3] = 1
            else:
                barcode_info[i][3] = 0
    
        #######################
    
    
 
        ###############
        csv_record = []
        csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'attention_score', 'component', 'from_id', 'to_id'])
        for j in range (0, len(barcode_info)):
            for i in range (0, len(barcode_info)):
                
                if i==j:
                    if len(lig_rec_dict[i][j])==0:
                        continue
                 
                atn_score_list = attention_scores[i][j]
                for k in range (0, len(atn_score_list)):
                    if attention_scores[i][j][k] >= threshold_down and attention_scores[i][j][k] <= threshold_up: 
                        if barcode_info[i][3]==0:
                            print('error')
                        elif barcode_info[i][3]==1:
                            csv_record.append([barcode_info[i][0], barcode_info[j][0], lig_rec_dict[i][j][k][0], lig_rec_dict[i][j][k][1], min_attention_score + attention_scores[i][j][k], '0-single', i, j])
                        else:
                            csv_record.append([barcode_info[i][0], barcode_info[j][0], lig_rec_dict[i][j][k][0], lig_rec_dict[i][j][k][1], min_attention_score + attention_scores[i][j][k], barcode_info[i][3], i, j])
 
        ###########	
        #run = 1
        #csv_record_dict = defaultdict(list)
        print('records found %d'%len(csv_record))
        for i in range (1, len(csv_record)): 
            key_value = str(csv_record[i][6]) +'-'+ str(csv_record[i][7]) + '-' + csv_record[i][2] + '-' + csv_record[i][3]# + '-'  + str( csv_record[i][5])
            csv_record_dict[key_value].append([csv_record[i][4], run])
        
    for key_value in csv_record_dict.keys():
        run_dict = defaultdict(list)
        for scores in csv_record_dict[key_value]:
            run_dict[scores[1]].append(scores[0])
        
        for runs in run_dict.keys():
            run_dict[runs] = np.mean(run_dict[runs]) 
        
 
        csv_record_dict[key_value] = []
        for runs in run_dict.keys():
            csv_record_dict[key_value].append([run_dict[runs],runs])
  
    #######################################
    
    all_edge_list = []
    for key_value in csv_record_dict.keys():
        edge_score_runs = []
        edge_score_runs.append(key_value)
        for runs in csv_record_dict[key_value]:
            edge_score_runs.append(runs[0]) #
        
        all_edge_list.append(edge_score_runs) # [key_value, score_by_run1, score_by_run2, etc.]

    ## Find the rank product #####################################################################
    ## all_edge_list has all the edges along with their scores for different runs in following format: 
    ## [edge_1_info, score_by_run1, score_by_run2, etc.], [edge_2_info, score_by_run1, score_by_run2, etc.], ..., [edge_N_info, score_by_run1, score_by_run2, etc.]
    edge_rank_dictionary = defaultdict(list)
    # sort the all_edge_list by each run's rank 
    print('total runs %d'%total_runs)
    for runs in range (0, total_runs):
        sorted_list_temp = sorted(all_edge_list, key = lambda x: x[runs+1], reverse=True) # sort based on attention score by current run: large to small
        for rank in range (0, len(sorted_list_temp)):
            edge_rank_dictionary[sorted_list_temp[rank][0]].append(rank+1) # small rank being high attention
            
    max_weight = len(sorted_list_temp)
    all_edge_vs_rank = []
    for key_val in edge_rank_dictionary.keys():
        rank_product = 1
        attention_score_list = csv_record_dict[key_value]
        avg_score = 0 #[]
        total_weight = 0
        for i in range (0, len(edge_rank_dictionary[key_val])):
            rank_product = rank_product * edge_rank_dictionary[key_val][i]
            weight_by_run = max_weight - edge_rank_dictionary[key_val][i]
            avg_score = avg_score + attention_score_list[i][0] * weight_by_run
            #avg_score.append(attention_score_list[i][0])  
            total_weight = total_weight + weight_by_run
            
        avg_score = avg_score/total_weight # lower weight being higher attention np.max(avg_score) #
        all_edge_vs_rank.append([key_val, rank_product**(1/total_runs), avg_score])  # small rank being high attention
        distribution_rank[layer].append(rank_product**(1/total_runs))
        
    all_edge_sorted_by_rank[layer] = sorted(all_edge_vs_rank, key = lambda x: x[1]) # small rank being high attention 

#############################################################################################################################################

################################ or ###############################################################################################################
percentage_value = 20 ##100 #20 # top 20th percentile rank, low rank means higher attention score
csv_record_intersect_dict = defaultdict(list)
edge_score_intersect_dict = defaultdict(list)
for layer in range (0, 2):
    threshold_up = np.percentile(distribution_rank[layer], percentage_value) #np.round(np.percentile(distribution_rank[layer], percentage_value),2)
    for i in range (0, len(all_edge_sorted_by_rank[layer])):
        if all_edge_sorted_by_rank[layer][i][1] <= threshold_up:
            csv_record_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(i)
            edge_score_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][2])
###########################################################################################################################################
## get the aggregated rank for all the edges ##
distribution_temp = []
for key_value in csv_record_intersect_dict.keys():  
    arg_index = np.argmin(csv_record_intersect_dict[key_value])
    csv_record_intersect_dict[key_value] = np.min(csv_record_intersect_dict[key_value]) # smaller rank being the higher attention
    edge_score_intersect_dict[key_value] = edge_score_intersect_dict[key_value][arg_index]
    distribution_temp.append(csv_record_intersect_dict[key_value]) 

#################

################################################################################
csv_record_dict = copy.deepcopy(csv_record_intersect_dict)
combined_score_distribution = []
csv_record = []
csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'edge_rank', 'component', 'from_id', 'to_id', 'attention_score'])
for key_value in csv_record_dict.keys():
    item = key_value.split('-')
    i = int(item[0])
    j = int(item[1])
    ligand = item[2]
    receptor = item[3]        
    edge_rank = csv_record_dict[key_value]        
    score = edge_score_intersect_dict[key_value] # weighted average attention score, where weight is the rank, lower rank being higher attention score
    label = -1 
    csv_record.append([barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score])
    combined_score_distribution.append(score)

        
print('common LR count %d'%len(csv_record))

##### scale the attention scores from 0 to 1 : high score representing higher attention ########
score_distribution = []
for k in range (1, len(csv_record)):
    score_distribution.append(csv_record[k][8])

min_score = np.min(score_distribution)
max_score = np.max(score_distribution)
for k in range (1, len(csv_record)):
    scaled_score = (csv_record[k][8]-min_score)/(max_score-min_score) 
    csv_record[k][8] = scaled_score


##### save the file for downstream analysis ########
csv_record_final = []
csv_record_final.append(csv_record[0])
for k in range (1, len(csv_record)):
    ligand = csv_record[k][2]
    receptor = csv_record[k][3]
    #if ligand =='CCL19' and receptor == 'CCR7':
    csv_record_final.append(csv_record[k])



    
df = pd.DataFrame(csv_record_final) # output 4
df.to_csv('/cluster/home/t116508uhn/64630/NEST_combined_rank_product_output_'+args.data_name+'_split_top20percent.csv', index=False, header=False)
#df.to_csv('/cluster/home/t116508uhn/64630/NEST_combined_rank_product_output_'+args.data_name+'_h2048_top20percent.csv', index=False, header=False)
#df.to_csv('/cluster/home/t116508uhn/64630/NEST_combined_rank_product_output_'+args.data_name+'_top20percent.csv', index=False, header=False)
#df.to_csv('/cluster/home/t116508uhn/64630/NEST_combined_rank_product_output_'+args.data_name+'_all.csv', index=False, header=False)
#df.to_csv('/cluster/home/t116508uhn/64630/NEST_combined_rank_product_output_'+args.data_name+'_all_TcellZone.csv', index=False, header=False)

