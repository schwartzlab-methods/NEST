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

spot_diameter = 89.43 #pixels
data_name = 'PDAC_130355_D1'


if data_name == 'LUAD_GSM5702473_TD1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='LUAD_GSM5702473_TD1', help='The name of dataset')
    parser.add_argument( '--model_name', type=str, default='gat_r1_3attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()  	
#############################################################   
elif data_name == 'V1_Human_Lymph_Node_spatial':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='V1_Human_Lymph_Node_spatial', help='The name of dataset')
    args = parser.parse_args()


elif data_name == 'PDAC_64630':
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


elif data_name == 'PDAC_130355_D1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/exp1/exp1_D1/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_130355_D1', help='The name of dataset')
    #parser.add_argument( '--model_name', type=str, default='gat_r1_2attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()
    filter_min_cell = 5
    threshold_expression = 98
####### get the gene id, cell barcode, cell coordinates ######

if data_name == 'LUAD_GSM5702473_TD1':
    # read the mtx file
    temp = sc.read_10x_mtx(args.data_path)
    print(temp)
    sc.pp.filter_genes(temp, min_cells=filter_min_cell)
    print(temp)
    
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
    temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(adata_h5.X)))  
    adata_X = np.transpose(temp)  
    cell_vs_gene = copy.deepcopy(adata_X)
    print('min gene count after quantile transformation %g'%np.min(cell_vs_gene))    

##################### make cell metadata: barcode_info ###################################
i=0
barcode_serial = dict()
for cell_code in cell_barcode:
    barcode_serial[cell_code]=i
    i=i+1
    
i=0
barcode_info=[]
for cell_code in cell_barcode:
    barcode_info.append([cell_code, coordinates[i,0],coordinates[i,1], 0]) # last entry will hold the component number later
    i=i+1



#### needed if split data is used ##############
i=0
node_id_sorted_xy=[]
for cell_code in cell_barcode:
    node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
    i=i+1
	
node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + args.data_name+'_'+'node_id_sorted_xy', 'wb') as fp:  #b, a:[0:5]   
	pickle.dump(node_id_sorted_xy, fp)
################################################

gene_info=dict()
for gene in gene_ids:
    gene_info[gene]=''

gene_index=dict()    
i = 0
for gene in gene_ids: 
    gene_index[gene] = i
    i = i+1
################# for running Niches ###############################
'''
gene_vs_cell = np.transpose(cell_vs_gene)  
np.save("/cluster/projects/schwartzgroup/fatema/find_ccc/gene_vs_cell_quantile_transformed_"+args.data_name, gene_vs_cell)
df = pd.DataFrame(gene_ids)
df.to_csv('/cluster/projects/schwartzgroup/fatema/find_ccc/gene_ids_'+args.data_name+'.csv', index=False, header=False)
df = pd.DataFrame(cell_barcode)
df.to_csv('/cluster/projects/schwartzgroup/fatema/find_ccc/cell_barcode_'+args.data_name+'.csv', index=False, header=False)
'''   
####################
''' 
for i in range (0, cell_vs_gene.shape[0]):
    max_value = np.max(cell_vs_gene[i][:])
    min_value = np.min(cell_vs_gene[i][:])
    for j in range (0, cell_vs_gene.shape[1]):
        cell_vs_gene[i][j] = (cell_vs_gene[i][j]-min_value)/(max_value-min_value)

with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'cell_vs_gene_quantile_transformed_scaled', 'wb') as fp:  #b, a:[0:5]   
	pickle.dump(cell_vs_gene, fp)
'''

ligand_dict_dataset = defaultdict(list)
cell_cell_contact = dict() 
   
cell_chat_file = '/cluster/home/t116508uhn/Human-2020-Jin-LR-pairs_cellchat.csv'
df = pd.read_csv(cell_chat_file)
for i in range (0, df["ligand_symbol"].shape[0]):
    ligand = df["ligand_symbol"][i]
    #if ligand not in gene_marker_ids:
    if ligand not in gene_info:
        continue
        
    if df["annotation"][i] == 'ECM-Receptor':    
        continue
        
    receptor_symbol_list = df["receptor_symbol"][i]
    receptor_symbol_list = receptor_symbol_list.split("&")
    for receptor in receptor_symbol_list:
        if receptor in gene_info:
        #if receptor in gene_marker_ids:
            ligand_dict_dataset[ligand].append(receptor)
            #######
            if df["annotation"][i] == 'Cell-Cell Contact':
                cell_cell_contact[receptor] = ''
            #######                
            
print(len(ligand_dict_dataset.keys()))

nichetalk_file = '/cluster/home/t116508uhn/NicheNet-LR-pairs.csv'   
df = pd.read_csv(nichetalk_file)
for i in range (0, df["from"].shape[0]):
    ligand = df["from"][i]
    #if ligand not in gene_marker_ids:
    if ligand not in gene_info:
        continue
    receptor = df["to"][i]
    #if receptor not in gene_marker_ids:
    if receptor not in gene_info:
        continue
    ligand_dict_dataset[ligand].append(receptor)
    
##############################################################
print('number of ligands %d '%len(ligand_dict_dataset.keys()))
count_pair = 0
for gene in list(ligand_dict_dataset.keys()): 
    ligand_dict_dataset[gene]=list(set(ligand_dict_dataset[gene]))
    gene_info[gene] = 'included'
    for receptor_gene in ligand_dict_dataset[gene]:
        gene_info[receptor_gene] = 'included'
        count_pair = count_pair + 1
        
print('number of pairs %d '%count_pair)       

count = 0
included_gene=[]
for gene in gene_info.keys(): 
    if gene_info[gene] == 'included':
        count = count + 1
        included_gene.append(gene)
        
print('number of affected genes %d '%count)
affected_gene_count = count
######################################

lr_gene_index = []
for gene in gene_info.keys(): 
    if gene_info[gene] == 'included':
        lr_gene_index.append(gene_index[gene])

lr_gene_index = sorted(lr_gene_index)
cell_vs_lrgene = cell_vs_gene[:, lr_gene_index]
with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'cell_vs_lrgene_quantile_transformed_'+args.data_name, 'wb') as fp:  #b, a:[0:5]   
	pickle.dump(cell_vs_lrgene, fp)
	
''''''
######################################

ligand_list = list(ligand_dict_dataset.keys())  
print('len ligand_list %d'%len(ligand_list))
total_relation = 0
l_r_pair = dict()
count = 0
lr_id = 0
for gene in list(ligand_dict_dataset.keys()): 
    ligand_dict_dataset[gene]=list(set(ligand_dict_dataset[gene]))
    l_r_pair[gene] = dict()
    for receptor_gene in ligand_dict_dataset[gene]:
        l_r_pair[gene][receptor_gene] = lr_id 
        lr_id  = lr_id  + 1
        
print('total type of l-r pairs found: %d'%lr_id )


from sklearn.metrics.pairwise import euclidean_distances
distance_matrix = euclidean_distances(coordinates, coordinates)

dist_X = np.zeros((distance_matrix.shape[0], distance_matrix.shape[1]))
for j in range(0, distance_matrix.shape[1]):
    max_value=np.max(distance_matrix[:,j])
    min_value=np.min(distance_matrix[:,j])
    for i in range(distance_matrix.shape[0]):
        dist_X[i,j] = 1-(distance_matrix[i,j]-min_value)/(max_value-min_value)
        	
    #list_indx = list(np.argsort(dist_X[:,j]))
    #k_higher = list_indx[len(list_indx)-k_nn:len(list_indx)]
    for i in range(0, distance_matrix.shape[0]):
        if distance_matrix[i,j] > spot_diameter*4: #i not in k_higher:
            dist_X[i,j] = 0 #-1
            
cell_rec_count = np.zeros((cell_vs_gene.shape[0]))

cell_percentile = []
for i in range (0, cell_vs_gene.shape[0]):
    y = sorted(cell_vs_gene[i]) # sort each row/cell
    x = range(1, len(y)+1)
    kn = KneeLocator(x, y, curve='convex', direction='increasing')
    kn_value = y[kn.knee-1]
    cell_percentile.append([np.percentile(y, 10), np.percentile(y, 20),np.percentile(y, 90), np.percentile(y, threshold_expression), kn_value])

##############################################################################
count_total_edges = 0
activated_cell_index = dict()

cells_ligand_vs_receptor = []
for i in range (0, cell_vs_gene.shape[0]):
    cells_ligand_vs_receptor.append([])
    
for i in range (0, cell_vs_gene.shape[0]):
    for j in range (0, cell_vs_gene.shape[0]):
        cells_ligand_vs_receptor[i].append([])
        cells_ligand_vs_receptor[i][j] = []
start_index = 0 #args.slice
end_index = len(ligand_list) #min(len(ligand_list), start_index+100)
included_LR = defaultdict(dict)
for g in range(start_index, end_index): 
    gene = ligand_list[g]
    for i in range (0, cell_vs_gene.shape[0]): # ligand
        count_rec = 0    
        if cell_vs_gene[i][gene_index[gene]] < cell_percentile[i][3]:
            continue
        
        for j in range (0, cell_vs_gene.shape[0]): # receptor
            if distance_matrix[i,j] > spot_diameter*4:
                continue

            for gene_rec in ligand_dict_dataset[gene]:
                if cell_vs_gene[j][gene_index[gene_rec]] >= cell_percentile[j][3]: # or cell_vs_gene[i][gene_index[gene]] >= cell_percentile[i][4] :#gene_list_percentile[gene_rec][1]: #global_percentile: #
                    if gene_rec in cell_cell_contact and distance_matrix[i,j] > spot_diameter:
                        continue

                    communication_score = cell_vs_gene[i][gene_index[gene]] * cell_vs_gene[j][gene_index[gene_rec]]
                    relation_id = l_r_pair[gene][gene_rec]
                    #print("%s - %s "%(gene, gene_rec))
                    if communication_score<=0:
                        print('zero valued ccc score found')
                        continue	
                    cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])
                    included_LR[gene][gene_rec] = ''
                    count_rec = count_rec + 1
                    count_total_edges = count_total_edges + 1
                    activated_cell_index[i] = ''
                    activated_cell_index[j] = ''

                            
        cell_rec_count[i] =  count_rec   
        #print("%d - %d "%(i, count_rec))
        #print("%d - %d , max %g and min %g "%(i, count_rec, max_score, min_score))
    
    print(g)
    
print('total number of edges in the input graph %d '%count_total_edges)

 
################################################################################
ccc_index_dict = dict()
row_col = []
edge_weight = []
lig_rec = []
count_edge = 0
max_local = 0
#local_list = np.zeros((102))
for i in range (0, len(cells_ligand_vs_receptor)):
    #ccc_j = []
    for j in range (0, len(cells_ligand_vs_receptor)):
        if distance_matrix[i][j] <= spot_diameter*4: 
            count_local = 0
            if len(cells_ligand_vs_receptor[i][j])>0:
                for k in range (0, len(cells_ligand_vs_receptor[i][j])):
                    gene = cells_ligand_vs_receptor[i][j][k][0]
                    gene_rec = cells_ligand_vs_receptor[i][j][k][1]
                    count_edge = count_edge + 1
                    count_local = count_local + 1
                    #print(count_edge)                      
                    mean_ccc = cells_ligand_vs_receptor[i][j][k][2]
                    row_col.append([i,j])
                    edge_weight.append([dist_X[i,j], mean_ccc,cells_ligand_vs_receptor[i][j][k][3]])
                    lig_rec.append([gene, gene_rec])                      
                
                if max_local < count_local:
                    max_local = count_local


print('len row col %d'%len(row_col))
print('count local %d'%max_local) 



##########
#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" +args.data_name+'_adjacency_records_GAT_selective_lr_STnCCC_separate_'+'bothAbove_cell99th', 'wb') as fp:  #b, a:[0:5]   
with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" +args.data_name+'_adjacency_records_GAT_selective_lr_STnCCC_separate_'+'bothAbove_cell98th_3d', 'wb') as fp:  #b, a:[0:5]  _filtered 
    pickle.dump([row_col, edge_weight, lig_rec], fp)
             
with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" +args.data_name+'_cell_vs_gene_quantile_transformed', 'wb') as fp:  #b, a:[0:5]   _filtered
	pickle.dump(cell_vs_gene, fp)




	
