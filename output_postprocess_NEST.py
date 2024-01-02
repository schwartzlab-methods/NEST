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
	
i=0
node_id_sorted_xy=[]
for cell_code in cell_barcode:
    node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
    i=i+1
	
node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))



####### load annotations ##############################################
if data_name == 'LUAD_GSM5702473_TD1':
    '''
    ccc_too_many_cells_LUAD = pd.read_csv('/cluster/projects/schwartzgroup/fatema/CCST/exp2_D1_ccc_toomanycells_cluster.csv')
    ccc_too_many_cells_LUAD_dict = dict()
    for i in range(0, len(ccc_too_many_cells_LUAD)):
        ccc_too_many_cells_LUAD_dict[ccc_too_many_cells_LUAD['cell'][i]] = int(ccc_too_many_cells_LUAD['cluster'][i])
    
    for i in range(0, len(barcode_info)):
        barcode_info[i][3] = ccc_too_many_cells_LUAD_dict[barcode_info[i][0]]
    	
    '''
    barcode_type=dict()
    for i in range (0, len(barcode_info)):
        barcode_type[barcode_info[i][0]] = 0 
    
    
elif data_name == 'V1_Human_Lymph_Node_spatial':
    '''
    pathologist_label_file='/cluster/home/t116508uhn/human_lymphnode_Spatial10X_manual_annotations.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        if pathologist_label[i][1] == 'GC': 
            barcode_type[pathologist_label[i][0]] = 1 
        #elif pathologist_label[i][1] =='':
        #    barcode_type[pathologist_label[i][0]] = 0 #'stroma_deserted'
        else:
            barcode_type[pathologist_label[i][0]] = 0
    '''
    pathologist_label_file=current_directory + '/spot_vs_type_dataframe_V1_HumanLympNode.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    	
    spot_label = []
    for i in range (1, len(pathologist_label)):
        spot_label.append([pathologist_label[i][0], float(pathologist_label[i][1]), float(pathologist_label[i][2]), float(pathologist_label[i][3])])
        
    spot_label = sorted(spot_label, key = lambda x: x[3], reverse=True) # descending order of 
    barcode_Tcell = []
    barcode_B = []
    barcode_GC = []
    for i in range (0, len(spot_label)):
        if spot_label[i][1] >= (spot_label[i][2] + spot_label[i][3])*2:
            barcode_Tcell.append(spot_label[i][0])
        elif spot_label[i][2] >= (spot_label[i][1] + spot_label[i][3])*2:
            barcode_B.append(spot_label[i][0])
        elif spot_label[i][3] >= (spot_label[i][1] + spot_label[i][2])*2:
            barcode_GC.append(spot_label[i][0])
           
    barcode_type=dict()
    for i in range (0, len(barcode_Tcell)):
        barcode_type[barcode_Tcell[i]] = 1 # tcell
        
    for i in range (0, len(barcode_B)):
        barcode_type[barcode_B[i]] = 2
      
    for i in range (0, len(barcode_GC)):
        barcode_type[barcode_GC[i]] = 3
        
    for i in range (0, len(spot_label)):
        if spot_label[i][0] not in barcode_type:
            barcode_type[spot_label[i][0]] = 0
    #############################################################################
    '''
    data_list=dict()
    data_list['pathology_label']=[]
    data_list['X']=[]
    data_list['Y']=[]
    
    for i in range (0, len(barcode_info)):
        data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(barcode_info[i][2])
    
    
    data_list_pd = pd.DataFrame(data_list)
    set1 = altairThemes.get_colour_scheme("Set1", 4)
    set1[0] = '#000000'
    
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
        color=alt.Color('pathology_label:N', scale=alt.Scale(range=set1)),
        tooltip=['pathology_label']
    )
    
    save_path = '/cluster/home/t116508uhn/64630/'
    chart.save(save_path+'V1_humanLymphNode.html') #   
    '''            
elif data_name == 'PDAC_64630':
    pathologist_label_file='/cluster/home/t116508uhn/IX_annotation_artifacts.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)	
    	
    barcode_type=dict() # record the type (annotation) of each spot (barcode)
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]
        
    '''
    ########## sabrina ###########################################	
    pathologist_label_file='/cluster/projects/schwartzgroup/fatema/find_ccc/singleR_spot_annotation_Sabrina.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
    
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        if pathologist_label[i][1] == 'tumor': #'Tumour':
            barcode_type[pathologist_label[i][0]] = 1 #'tumor'
        elif pathologist_label[i][1] =='stroma_deserted':
            barcode_type[pathologist_label[i][0]] = 0 #'stroma_deserted'
        elif pathologist_label[i][1] =='acinar_reactive':
            barcode_type[pathologist_label[i][0]] = 2 #'acinar_reactive'
        else:
            barcode_type[pathologist_label[i][0]] = 0 #'zero' 
    '''	
    #################################################################		
elif data_name == 'PDAC_140694':
    spot_type = []
    pathologist_label_file='/cluster/home/t116508uhn/V10M25-060_C1_T_140694_Histology_annotation_IX.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
            spot_type.append(line[1])
            
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1] 
	    
elif data_name == 'PDAC_130355_B1':
    spot_type = []
    pathologist_label_file='/cluster/home/t116508uhn/V10M26-61_B1_IX_annotation_pathology.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
            spot_type.append(line[1])
            
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1] 
	    
elif data_name == 'PDAC_130355_A1':
    spot_type = []
    pathologist_label_file='/cluster/home/t116508uhn/V10M25-61_A1_N_130355_Histology_annotation_IX.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
            spot_type.append(line[1])
            
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1] 
	    
elif data_name == 'PDAC_130355_D1':
    spot_type = []
    pathologist_label_file='/cluster/home/t116508uhn/v10M25-060_D1_N_130355_Histology_annotation_IX.csv' #IX_annotation_artifacts.csv' #
    pathologist_label=[]
    with open(pathologist_label_file) as file:
        csv_file = csv.reader(file, delimiter=",")
        for line in csv_file:
            pathologist_label.append(line)
            spot_type.append(line[1])
            
    barcode_type=dict()
    for i in range (1, len(pathologist_label)):
        barcode_type[pathologist_label[i][0]] = pathologist_label[i][1] 

##### load input graph ##########################

#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'adjacency_records_GAT_selective_lr_STnCCC_c_'+'all_avg', 'rb') as fp:  #b, a:[0:5]           
#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'adjacency_records_GAT_synthetic_region1_onlyccc_70', 'wb') as fp:
#    row_col, edge_weight = pickle.load(fp)
with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'adjacency_records_GAT_selective_lr_STnCCC_separate_'+'bothAbove_cell98th_3d', 'rb') as fp:  #b, a:[0:5]   
#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'adjacency_records_GAT_selective_lr_STnCCC_separate_'+'all_kneepoint_woBlankedge', 'rb') as fp:  #b, a:[0:5]   
#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" + 'adjacency_records_GAT_omniPath_separate_'+'threshold_distance_density_kneepoint', 'rb') as fp:  #b, a:[0:5]   
#with gzip.open("/cluster/projects/schwartzgroup/fatema/find_ccc/" +args.data_name+ '_adjacency_records_GAT_selective_lr_STnCCC_separate_'+'bothAbove_cell98th_3d', 'rb') as fp: 
    row_col, edge_weight, lig_rec = pickle.load(fp) # density_

datapoint_size = len(barcode_info)
lig_rec_dict = []
self_loop_found = defaultdict(dict)
for i in range (0, datapoint_size):
    lig_rec_dict.append([])  
    for j in range (0, datapoint_size):	
        lig_rec_dict[i].append([])   
        lig_rec_dict[i][j] = []

for index in range (0, len(row_col)):
        i = row_col[index][0]
        j = row_col[index][1]
        lig_rec_dict[i][j].append(lig_rec[index])  
        self_loop_found[i][j] = ''

with gzip.open('/cluster/home/t116508uhn/64630/'+'self_loop_record_'+args.data_name, 'wb') as fp:  #b, a:[0:5]   _filtered
	pickle.dump(self_loop_found, fp)

############################################################################
'''
attention_scores = []
datapoint_size = len(barcode_info)
for i in range (0, datapoint_size):
    attention_scores.append([])   
    for j in range (0, datapoint_size):	
        attention_scores[i].append([])   
        attention_scores[i][j] = []

distribution = []
ccc_index_dict = dict()
combined_score_distribution_ccl19_ccr7 = []
combined_score_distribution = []
for index in range (0, len(row_col)):
    i = row_col[index][0]
    j = row_col[index][1]
    if i==j:
        if len(lig_rec_dict[i][j])==0:
            continue 
    if (barcode_type[cell_barcode[i]]==1 and barcode_type[cell_barcode[j]]==1) != True: #i in spot_interest_list and j in spot_interest_list:
    #if lig_rec[index][0]!='CCL19' or lig_rec[index][1] != "CCR7": #lig_rec[index][0]=='IL21' and lig_rec[index][1] == "IL21R": #
        continue
    if edge_weight[index][1]>0:
        ligand = lig_rec[index][0]
        receptor = lig_rec[index][1]
        if ligand =='CCL19' and receptor == 'CCR7':
            combined_score_distribution_ccl19_ccr7.append(edge_weight[index][1])
        
        combined_score_distribution.append(edge_weight[index][1])
        attention_scores[i][j].append(edge_weight[index][1]) # * edge_weight[index][0]) # * edge_weight[index][2])
        distribution.append(edge_weight[index][1]) # * edge_weight[index][0]) # * edge_weight[index][2])
        ccc_index_dict[i] = ''
        ccc_index_dict[j] = ''   
        #lig_rec_dict[i][j].append(lig_rec[index])          

        
some_dict = dict(A=combined_score_distribution, B=combined_score_distribution_ccl19_ccr7)
df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in some_dict.items()]))
df = df.rename(columns={'A': 'all_pairs', 'B': 'CCL19_CCR7'})
source = df
#########################################################################
chart = alt.Chart(source).transform_fold(
    ['all_pairs',
     'CCL19_CCR7'],
    as_ = ['distribution_type', 'value']
).transform_density(
    density = 'value',
    groupby=['distribution_type'],        
    steps=0.5
).mark_area(opacity=0.7).encode(
    alt.X('value:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('distribution_type:N')
)   
chart.save('/cluster/home/t116508uhn/64630/'+'region_of_interest_cccScore_distribution.html')

'''

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
        ##X_attention_filename = args.embedding_data_path + args.data_name + '/' + 'PDAC_cellchat_nichenet_threshold_distance_bothAbove_cell98th_tanh_3dim_'+filename[run_time]+'attention_l1.npy' #a
        #X_attention_filename = args.embedding_data_path + args.data_name + '/' + 'PDAC_cellchat_nichenet_threshold_distance_bothAbove_cell98th_tanh_3dim_h2048_'+filename[run_time]+'attention_l1.npy' #a 
        X_attention_bundle = np.load(X_attention_filename, allow_pickle=True) #_withFeature
        for index in range (0, X_attention_bundle[0].shape[1]):
            i = X_attention_bundle[0][0][index]
            j = X_attention_bundle[0][1][index]
            #if barcode_type[barcode_info[i][0]] != 1 or barcode_type[barcode_info[j][0]] != 1:
            #    continue
            distribution.append(X_attention_bundle[l][index][0])

        ################# scaling the attention scores so that layer 1 and 2 will be comparable ##############################        
        min_attention_score = 1000
        max_value = np.max(distribution)
        min_value = np.min(distribution)
        distribution = []
        for index in range (0, X_attention_bundle[0].shape[1]):
            i = X_attention_bundle[0][0][index]
            j = X_attention_bundle[0][1][index]
            #if barcode_type[barcode_info[i][0]] != 1 or barcode_type[barcode_info[j][0]] != 1:
            #    continue
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

            
        '''
        ####################following if split###########################
        for set_id in range(0, len(edge_list)):
            print('subgraph %d'%set_id)
            ##############
            set1_exist_dict = defaultdict(dict)
            for i in range (0, datapoint_size):  
                for j in range (0, datapoint_size):	
                    set1_exist_dict[i][j]=-1
    
            for edge in edge_list[set_id]:
                row_col = edge[0]
                new_i = row_col[0]
                new_j = row_col[1]
                i = id_map_new_old[set_id][new_i] 
                j = id_map_new_old[set_id][new_j] 
                set1_exist_dict[i][j] = 1
            
            ############
            #X_attention_filename = args.embedding_data_path + args.data_name + '/' + args.data_name + '_cellchat_nichenet_threshold_distance_bothAbove_cell98th_tanh_3dim_split_'+filename[run_time]+'attention_l1_'+str(set_id+1)+'.npy' #_h1024
            X_attention_filename = args.embedding_data_path + args.data_name + '/' + 'PDAC_cellchat_nichenet_threshold_distance_bothAbove_cell98th_tanh_3dim_split_'+filename[run_time]+'attention_l1_'+str(set_id+1)+'.npy'
            X_attention_bundle = np.load(X_attention_filename, allow_pickle=True) #_withFeature args.data_name + 
    
            for index in range (0, X_attention_bundle[0].shape[1]):
                new_i = X_attention_bundle[0][0][index]
                new_j = X_attention_bundle[0][1][index] 
                # these returned i and j are new
                i = id_map_new_old[set_id][new_i] 
                j = id_map_new_old[set_id][new_j]           
                            
                if i in set1_exist_dict and j in set1_exist_dict[i] and set1_exist_dict[i][j]==1:
                ###################################
                    attention_scores[i][j].append(X_attention_bundle[l][index][0]) 
                    distribution.append(X_attention_bundle[l][index][0])

        min_attention_score = 1000
        max_value = np.max(distribution)
        min_value = np.min(distribution)
        distribution = []    
        for i in range (0, datapoint_size):  
            for j in range (0, datapoint_size):    
                if len(attention_scores[i][j])>0:
                    for k in range (0, len(attention_scores[i][j])):
                        scaled_score = (attention_scores[i][j][k]-min_value)/(max_value-min_value)
                        attention_scores[i][j][k] = scaled_score 
                        distribution.append(attention_scores[i][j][k])
                        if min_attention_score > scaled_score:
                            min_attention_score = scaled_score                    
                        
        
        ##################################
        if min_attention_score<0:
            min_attention_score = -min_attention_score
        else: 
            min_attention_score = 0
            
        '''
        ####################################################################################



        
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
    
    
        '''
        data_list=dict()
        data_list['pathology_label']=[]
        data_list['component_label']=[]
        data_list['X']=[]
        data_list['Y']=[]
    
        for i in range (0, len(barcode_info)):
            #if barcode_type[barcode_info[i][0]] == 'zero':
            #    continue
            data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
            data_list['component_label'].append(barcode_info[i][3])
            data_list['X'].append(barcode_info[i][1])
            data_list['Y'].append(barcode_info[i][2])
    
    
        data_list_pd = pd.DataFrame(data_list)
        #set1 = altairThemes.get_colour_scheme("Set1", len(data_list_pd["component_label"].unique()))
        set1 = altairThemes.get_colour_scheme("Set1", id_label)
        set1[0] = '#000000'
    
        chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
            alt.X('X', scale=alt.Scale(zero=False)),
            alt.Y('Y', scale=alt.Scale(zero=False)),
            shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
            color=alt.Color('component_label:N', scale=alt.Scale(range=set1)),
            tooltip=['component_label']
        )#.configure_legend(labelFontSize=6, symbolLimit=50)
        # output 2
        save_path = '/cluster/home/t116508uhn/64630/'
        chart.save(save_path+args.data_name+'_altair_plot_bothAbove98_3dim_tanh_3heads_l2attention_th95_'+filename[run_time]+'.html')
        '''
        ##############
        '''
        region_list =[2, 3, 9, 11, 4, 5, 7]
        
        spot_interest_list = []
        for i in range (0, len(barcode_info)):
            if data_list['component_label'][i] in region_list:
                
                spot_interest_list.append(i)
        '''
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
        '''
        df = pd.DataFrame(csv_record)
        df.to_csv('/cluster/home/t116508uhn/64630/input_test.csv', index=False, header=False)
        ############
        alt.themes.register("publishTheme", altairThemes.publishTheme)
        # enable the newly registered theme
        alt.themes.enable("publishTheme")
        inFile = '/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]
        df = readCsv(inFile)
        df = preprocessDf(df)
        outPathRoot = inFile.split('.')[0]
        p = plot(df)
        #outPath = '/cluster/home/t116508uhn/64630/test_hist_'+args.data_name+'_'+filename[run_time]+'_th99p7_h512_l2attention_'+str(len(csv_record))+'edges.html' #filteredl2attention__ l2attention_
        outPath = '/cluster/home/t116508uhn/64630/test_hist_'+args.data_name+'_'+filename[run_time]+'_selective_only_Tcellzone_th90_h512_'+str(len(csv_record))+'edges.html' #filteredl2attention__ l2attention_
        p.save(outPath)	# output 3
        '''
        ###########	
        #run = 1
        #csv_record_dict = defaultdict(list)
        print('records found %d'%len(csv_record))
        for i in range (1, len(csv_record)):
            # insert if both of them are in tcell zone
            #if (barcode_type[barcode_info[csv_record[i][6]][0]] == 1 and barcode_type[barcode_info[csv_record[i][7]][0]] == 1)!=True:
            #    continue
                
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
'''
csv_record_intersect_dict = defaultdict(list)
for layer in range (0, 2):
    for i in range (0, len(all_edge_sorted_by_rank[layer])):
        csv_record_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(i)
'''
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
'''
keyval = '2219-502-CCL21-CXCR4' #'1965-1868-CCL21-CXCR4'
for i in range (0, len( all_edge_sorted_by_rank)):
    if  all_edge_sorted_by_rank[0][i][0]==keyval:
        print(all_edge_sorted_by_rank[0][i])
        break
'''
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

############################################### IGNORE the rest ###########################################################################
############################################### ONLY for Human Lymph Node #################################################################
'''
combined_score_distribution_ccl19_ccr7 = []
for k in range (1, len(csv_record)):
    i = csv_record[k][6]
    j = csv_record[k][7]
    ligand = csv_record[k][2]
    receptor = csv_record[k][3]
    if ligand =='CCL19' and receptor == 'CCR7':
        combined_score_distribution_ccl19_ccr7.append(csv_record[k][4])
        
some_dict = dict(A=combined_score_distribution, B=combined_score_distribution_ccl19_ccr7)
'''
'''
     
    df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in some_dict.items()]))

    df = df.rename(columns={'A': 'all_pairs', 'B': 'CCL19_CCR7'})

    source = df

    chart = alt.Chart(source).transform_fold(
        ['all_pairs',
         'CCL19_CCR7'],
        as_ = ['distribution_type', 'value']
    ).transform_density(
        density = 'value',
        bandwidth=0.3,
        groupby=['distribution_type'],        
        counts = True,
        steps=100
    ).mark_area(opacity=0.5).encode(
        alt.X('value:Q'),
        alt.Y('density:Q', stack='zero' ),
        alt.Color('distribution_type:N')
    )#.properties(width=400, height=100)


chart = alt.Chart(source).transform_fold(
    ['all_pairs', 'CCL19_CCR7'],
    as_=['Distribution Type', 'Attention Score']
).mark_bar(
    opacity=0.5,
    binSpacing=0
).encode(
    alt.X('Attention Score:Q', bin=alt.Bin(maxbins=100)),
    alt.Y('count()', stack=None),
    alt.Color('Distribution Type:N')
)

chart.save(save_path+'region_of_interest_filtered_combined_attention_distribution.html')
'''
############################################################################################

connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))  
for k in range (1, len(csv_record_final)):
    i = csv_record_final[k][6]
    j = csv_record_final[k][7]
    connecting_edges[i][j]=1
        
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
    if count_points_component[labels[i]] > 1:
        barcode_info[i][3] = index_dict[labels[i]] #2
    elif connecting_edges[i][i] == 1 and len(lig_rec_dict[i][i])>0: 
        barcode_info[i][3] = 1
    else: 
        barcode_info[i][3] = 0

# update the label based on new component numbers
#max opacity
for record in range (1, len(csv_record_final)):
    i = csv_record_final[record][6]
    label = barcode_info[i][3]
    csv_record_final[record][5] = label


######################################################################################################################## 

###########	list those spots who are participating in CCC ##################
filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
inFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]
df = pd.read_csv(inFile, sep=",")
csv_record_final = df.values.tolist()
i=0
j=0
csv_record_final.append([barcode_info[i][0], barcode_info[j][0], 'no-ligand', 'no-receptor', 0, 0, i, j]) # dummy for histogram

df_column_names = list(df.columns)
csv_record_final = [df_column_names] + csv_record_final

active_spot = defaultdict(list)
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    record = csv_record_final[record_idx]
    i = record[6]
    pathology_label = barcode_type[barcode_info[i][0]]
    component_label = record[5]
    X = barcode_info[i][1]
    Y = -barcode_info[i][2]
    opacity = record[4]
    active_spot[i].append([pathology_label, component_label, X, Y, opacity])
    
    j = record[7]
    pathology_label = barcode_type[barcode_info[j][0]]
    component_label = record[5]
    X = barcode_info[j][1]
    Y = -barcode_info[j][2]
    opacity = record[4]   
    active_spot[j].append([pathology_label, component_label, X, Y, opacity])
    ''''''
    
######### color the spots in the plot with opacity = attention score #################
opacity_list = []
for i in active_spot:
    sum_opacity = []
    for edges in active_spot[i]:
        sum_opacity.append(edges[4])
        
    avg_opacity = np.max(sum_opacity) #np.mean(sum_opacity)
    opacity_list.append(avg_opacity) 
    active_spot[i]=[active_spot[i][0][0], active_spot[i][0][1], active_spot[i][0][2], active_spot[i][0][3], avg_opacity]

min_opacity = np.min(opacity_list)
max_opacity = np.max(opacity_list)
min_opacity = min_opacity - 5
#######################################################################################
data_list=dict()
data_list['pathology_label']=[]
data_list['component_label']=[]
data_list['X']=[]
data_list['Y']=[]   
data_list['opacity']=[] 

for i in range (0, len(barcode_info)):        
    if i in active_spot:
        data_list['pathology_label'].append(active_spot[i][0])
        data_list['component_label'].append(active_spot[i][1])
        data_list['X'].append(active_spot[i][2])
        data_list['Y'].append(active_spot[i][3])
        data_list['opacity'].append((active_spot[i][4]-min_opacity)/(max_opacity-min_opacity))
        
    else:
        data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
        data_list['component_label'].append(0)
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(-barcode_info[i][2])
        data_list['opacity'].append(0.1)


#id_label= len(list(set(data_list['component_label'])))
data_list_pd = pd.DataFrame(data_list)
set1 = altairThemes.get_colour_scheme("Set1", id_label)
set1[0] = '#000000'
chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
    alt.X('X', scale=alt.Scale(zero=False)),
    alt.Y('Y', scale=alt.Scale(zero=False)),
    shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
    color=alt.Color('component_label:N', scale=alt.Scale(range=set1)),
    #opacity=alt.Opacity('opacity:N'), #"opacity",
    tooltip=['component_label'] #,'opacity'
)#.configure_legend(labelFontSize=6, symbolLimit=50)

# output 6
save_path = '/cluster/home/t116508uhn/64630/'
chart.save(save_path+'altair_plot_test.html')
####################################################################################################################
filename_str = 'NEST_combined_output_'+args.data_name+'.csv'
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")
inFile = '/cluster/home/t116508uhn/64630/'+filename_str #'/cluster/home/t116508uhn/64630/input_test.csv' #sys.argv[1]
df = readCsv(inFile)
df = preprocessDf(df)
outPathRoot = inFile.split('.')[0]
p = plot(df)
outPath = '/cluster/home/t116508uhn/64630/test_hist_temp.html'
p.save(outPath)	# output 5
#####################################################################################################################




##################################################
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
set1 = altairThemes.get_colour_scheme("Set1", id_label)
colors = set1
colors[0] = '#000000'
ids = []
x_index=[]
y_index=[]
colors_point = []
for i in range (0, len(barcode_info)):    
    ids.append(i)
    x_index.append(barcode_info[i][1])
    y_index.append(barcode_info[i][2])    
    colors_point.append(colors[barcode_info[i][3]]) 
  
max_x = np.max(x_index)
max_y = np.max(y_index)


from pyvis.network import Network
import networkx as nx

barcode_type=dict()
for i in range (1, len(pathologist_label)):
    if 'tumor'in pathologist_label[i][1]: #'Tumour':
        barcode_type[pathologist_label[i][0]] = 1
    else:
        barcode_type[pathologist_label[i][0]] = 0
    '''
    elif pathologist_label[i][1] == 'stroma_deserted':
        barcode_type[pathologist_label[i][0]] = 0
    elif pathologist_label[i][1] =='acinar_reactive':
        barcode_type[pathologist_label[i][0]] = 2
    else:
        barcode_type[pathologist_label[i][0]] = 'zero' #0
    '''
g = nx.MultiDiGraph(directed=True) #nx.Graph()
for i in range (0, len(barcode_info)):
    label_str =  str(i)+'_c:'+str(barcode_info[i][3])+'_'
    #if barcode_type[barcode_info[i][0]] == 'zero':
    #    continue
    if barcode_type[barcode_info[i][0]] == 0: #stroma
        marker_size = 'circle'
        label_str = label_str + 'stroma'
    elif barcode_type[barcode_info[i][0]] == 1: #tumor
        marker_size = 'box'
        label_str = label_str + 'tumor'
    else:
        marker_size = 'ellipse'
        label_str = label_str + 'acinar_reactive'
	
    g.add_node(int(ids[i]), x=int(x_index[i]), y=int(y_index[i]), label = label_str, pos = str(x_index[i])+","+str(-y_index[i])+" !", physics=False, shape = marker_size, color=matplotlib.colors.rgb2hex(colors_point[i]))    

nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True

count_edges = 0
for k in range (1, len(csv_record_final)):
    i = csv_record_final[k][6]
    j = csv_record_final[k][7]    
    ligand = csv_record_final[k][2]
    receptor = csv_record_final[k][3]
    title_str =  "L:"+ligand+", R:"+receptor
    edge_score = csv_record_final[k][4]
    g.add_edge(int(i), int(j), label = title_str, value=np.float64(edge_score), color=colors_point[i] ) 
    count_edges = count_edges + 1
     
nt.from_nx(g)
nt.show('mygraph.html')
cp mygraph.html /cluster/home/t116508uhn/64630/mygraph.html


from networkx.drawing.nx_agraph import write_dot
write_dot(g, "/cluster/home/t116508uhn/64630/test_interactive.dot")

