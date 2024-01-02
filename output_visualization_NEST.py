import numpy as np
import csv
import pickle
from scipy import sparse
import scanpy as sc
import matplotlib
matplotlib.use('Agg') 
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import stlearn as st
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
#from typing import List
import qnorm
from scipy.sparse import csr_matrix
#from scipy.sparse.csgraph import connected_components
from collections import defaultdict
import pandas as pd
import gzip
#import copy 
import argparse
import os
from scipy.sparse.csgraph import connected_components

#sys.path.append("/home/gw/code/utility/altairThemes/")
#if True:  # In order to bypass isort when saving
#    import altairThemes
import altairThemes
import altair as alt

alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")

spot_diameter = 89.43 #pixels
current_directory = '/cluster/home/t116508uhn/64630/'
top_edge_count = 1300
##########################################################
# readCsv, preprocessDf, plot: these three functions are taken from GW's repository                                                                                                                                                                     /mnt/data0/gw/research/notta_pancreatic_cancer_visium/plots/fatema_signaling/hist.py                                                                                                                                                                                         
import scipy.stats


def readCsv(x):
  """Parse file."""
  #colNames = ["method", "benchmark", "start", "end", "time", "memory"]
  df = pd.read_csv(x, sep=",")
 
  return df

def preprocessDf(df):
  """Transform ligand and receptor columns."""
  df["ligand-receptor"] = df["ligand"] + '-' + df["receptor"]
  df["component"] = df["component"] #.astype(str).str.zfill(2)

  return df


def plot(df):
  set1 = altairThemes.get_colour_scheme("Set1", len(df["component"].unique()))
  set1[0] = '#000000'
  base = alt.Chart(df).mark_bar().encode(
            x=alt.X("ligand-receptor:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y=alt.Y("count()"),
            color=alt.Color("component:N", scale = alt.Scale(range=set1)),
            order=alt.Order("component:N", sort="ascending"),
            tooltip=["component"]
        )
  p = base

  return p

####################### Set the name of the sample you want to visualize ###################################

data_name = 'PDAC_64630' #'PDAC_130355_D1' #'PDAC_140694' #'PDAC_130355_B1' #'PDAC_64630' #'LUAD_GSM5702473_TD1' #LUAD_GSM5702473_TD1



##########################################################
if data_name == 'LUAD_GSM5702473_TD1':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='LUAD_GSM5702473_TD1', help='The name of dataset')
    parser.add_argument( '--model_name', type=str, default='gat_r1_3attr', help='model name')
    
    args = parser.parse_args()  	
#############################################################   
elif data_name == 'V1_Human_Lymph_Node_spatial':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='V1_Human_Lymph_Node_spatial', help='The name of dataset')
    parser.add_argument( '--model_name', type=str, default='gat_r1_2attr', help='model name')
    parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()


elif data_name == 'PDAC_64630':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_64630', help='The name of dataset')
    #parser.add_argument( '--model_name', type=str, default='gat_2attr', help='model name')
    #parser.add_argument( '--slice', type=int, default=0, help='starting index of ligand')
    args = parser.parse_args()


elif data_name == 'PDAC_140694':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_path', type=str, default='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/V10M25-60_C1_PDA_140694_Pa_P_Spatial10x/outs/' , help='The path to dataset') 
    parser.add_argument( '--embedding_data_path', type=str, default='new_alignment/Embedding_data_ccc_rgcn/' , help='The path to attention') #'/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/'
    parser.add_argument( '--data_name', type=str, default='PDAC_140694', help='The name of dataset')
    args = parser.parse_args()
    filter_min_cell = 1
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


####### get the gene id, cell barcode, cell coordinates ######

if data_name == 'LUAD_GSM5702473_TD1':
    # read the mtx file
    temp = sc.read_10x_mtx(args.data_path)
    print(temp)
    sc.pp.filter_genes(temp, min_cells=1)
    print(temp)
    
    #gene_ids = list(temp.var_names) 
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
    #sc.pp.filter_genes(adata_h5, min_cells=filter_min_cell)
    #print(adata_h5)
    #gene_ids = list(adata_h5.var_names)
    coordinates = adata_h5.obsm['spatial']
    cell_barcode = np.array(adata_h5.obs.index)
    #temp = adata_h5.X


######################### get the cell vs gene matrix ##################
'''
temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(temp)))  
cell_vs_gene = np.transpose(temp)  
'''
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

 

####### load annotations ##############################################
if data_name == 'LUAD_GSM5702473_TD1':
    barcode_type=dict()
    for i in range (0, len(barcode_info)):
        barcode_type[barcode_info[i][0]] = 0 
    
    
elif data_name == 'V1_Human_Lymph_Node_spatial':
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


###############################  read which spots have self loops ################################################################
with gzip.open(current_directory +'self_loop_record_'+data_name, 'rb') as fp:  #'/cluster/projects/schwartzgroup/fatema/find_ccc/'+
    self_loop_found = pickle.load(fp)
    

######################### read the NEST output in csv format ####################################################

filename_str = 'NEST_combined_rank_product_output_'+args.data_name+'_top20percent.csv' #h2048_
inFile = current_directory +filename_str 
df = pd.read_csv(inFile, sep=",")
csv_record = df.values.tolist()

## sort the edges based on their rank (column 4 = score) column, low to high, low being higher attention score
csv_record = sorted(csv_record, key = lambda x: x[4])

## add the column names and take first top_edge_count edges
# columns are: from_cell, to_cell, ligand_gene, receptor_gene, attention_score, component, from_id, to_id
df_column_names = list(df.columns)
csv_record_final = [df_column_names] + csv_record[0:top_edge_count]

## add a dummy row at the end for the convenience of histogram preparation (to keep the color same as altair plot)
i=0
j=0
csv_record_final.append([barcode_info[i][0], barcode_info[j][0], 'no-ligand', 'no-receptor', 0, 0, i, j]) # dummy for histogram

######################## unique ligand-receptor count ##############################################################
unique_ligand_receptor = dict()
for k in range (1, len(csv_record_final)-1): 
    ligand = csv_record_final[k][2]
    receptor = csv_record_final[k][3]
    unique_ligand_receptor[ligand+'-'+receptor] = 1

print('unique ligand-receptor count is %d'%(len(unique_ligand_receptor.keys())))
####################### pattern finding ##########################################################################
# make a dictionary to keep record of all the outgoing edges [to_node, ligand, receptor] for each node
each_node_outgoing = defaultdict(list)
for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
    i = csv_record_final[k][6]
    j = csv_record_final[k][7]
    if i == j:
        continue        
    ligand = csv_record_final[k][2]
    receptor = csv_record_final[k][3]
    each_node_outgoing[i].append([j, ligand, receptor])

# all possible 2 hop pattern count
pattern_distribution = defaultdict(list)
# pattern_distribution['ligand-receptor to ligand-receptor']=[1,1,1,1, ...]
for i in each_node_outgoing:
    for tupple in each_node_outgoing[i]: # first hop
        j = tupple[0]
        lig_rec_1 = tupple[1]+'-'+tupple[2]
        if j in each_node_outgoing:
            for tupple_next in each_node_outgoing[j]: # second hop
                k = tupple_next[0]
                if k == i or k == j:
                    continue
                lig_rec_2 = tupple_next[1]+'-'+tupple_next[2]
                pattern_distribution[lig_rec_1 + ' to ' + lig_rec_2].append(1)

'''                
                if lig_rec_1 in pattern_distribution:
                    pattern_distribution[lig_rec_1][lig_rec_2].append(1)
                else:
                    pattern_distribution[lig_rec_1]=defaultdict(list)
                    pattern_distribution[lig_rec_1][lig_rec_2].append(1)
'''
two_hop_pattern_distribution = []
same_count = 0
for key in pattern_distribution:
    count = len(pattern_distribution[key])
    two_hop_pattern_distribution.append([key, count]) 
    #if lig_rec_1 == lig_rec_2:
    #    same_count = same_count + 1

two_hop_pattern_distribution = sorted(two_hop_pattern_distribution, key = lambda x: x[1], reverse=True) 

data_list=dict()
data_list['X']=[]
data_list['Y']=[]   
for key in pattern_distribution:
    for i in range (0, len(pattern_distribution[key])):
        data_list['X'].append(key)
        data_list['Y'].append(1)
        
data_list_pd = pd.DataFrame(data_list)
chart= alt.Chart(data_list_pd).mark_bar().encode(
            x=alt.X("X:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y=alt.Y("count()"),
        )
chart.save(current_directory +'pattern_distribution.html')




######################## connected component finding #################################

connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))  
for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
    i = csv_record_final[k][6]
    j = csv_record_final[k][7]
    connecting_edges[i][j]=1
        
graph = csr_matrix(connecting_edges)
n_components, labels = connected_components(csgraph=graph,directed=True, connection = 'weak',  return_labels=True) # It assigns each SPOT to a component based on what pair it belongs to
print('number of connected components %d'%n_components) 

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
    elif connecting_edges[i][i] == 1 and (i in self_loop_found and i in self_loop_found[i]): # that is: self_loop_found[i][i] do exist 
        barcode_info[i][3] = 1
    else: 
        barcode_info[i][3] = 0

# update the label based on found component numbers
#max opacity
for record in range (1, len(csv_record_final)-1):
    i = csv_record_final[record][6]
    label = barcode_info[i][3]
    csv_record_final[record][5] = label
 
############################################### Optional filtering ########################################################
'''
## change the csv_record_final here if you want histogram for specific components/regions only. e.g., if you want to plot only stroma region, or tumor-stroma regions etc.    ##
#region_of_interest = [...] 
csv_record_final_temp = []
csv_record_final_temp.append(csv_record_final[0])
component_dictionary_dummy = dict()
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    # if at least one spot of the pair is tumor, then plot it
    if (barcode_type[csv_record_final[record_idx][0]] == 'tumor' or barcode_type[csv_record_final[record_idx][1]] == 'tumor'): #((barcode_type[csv_record_final[record_idx][0]] == 'tumor' and barcode_type[csv_record_final[record_idx][1]] == 'tumor') or (barcode_type[csv_record_final[record_idx][0]] != 'tumor' and barcode_type[csv_record_final[record_idx][1]] != 'tumor')):
        csv_record_final_temp.append(csv_record_final[record_idx])
    if csv_record_final[record_idx][5] not in component_dictionary_dummy:
        component_dictionary_dummy[csv_record_final[record_idx][5]] = csv_record_final[record_idx]
        
# insert just one record from each other components so that the color scheme does not change in the altair scatter plot and histogram :-(
for component_id in component_dictionary_dummy:
    csv_record_final_temp.append(component_dictionary_dummy[component_id])
 
csv_record_final_temp.append(csv_record_final[len(csv_record_final)-1])
csv_record_final = copy.deepcopy(csv_record_final_temp)
'''

#####################################
component_list = dict()
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    record = csv_record_final[record_idx]
    i = record[6]
    j = record[7]
    component_label = record[5]
    barcode_info[i][3] = component_label #?
    barcode_info[j][3] = component_label #?
    component_list[component_label] = ''

component_list[0] = ''
unique_component_count = max(len(component_list.keys()), id_label)


##################################### Altair Plot ##################################################################
## dictionary of those spots who are participating in CCC ##
active_spot = defaultdict(list)
for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
    record = csv_record_final[record_idx]
    i = record[6]
    pathology_label = barcode_type[barcode_info[i][0]]
    component_label = record[5]
    X = barcode_info[i][1]
    Y = -barcode_info[i][2]
    opacity = np.float(record[8])
    active_spot[i].append([pathology_label, component_label, X, Y, opacity])
    
    j = record[7]
    pathology_label = barcode_type[barcode_info[j][0]]
    component_label = record[5]
    X = barcode_info[j][1]
    Y = -barcode_info[j][2]
    opacity = np.float(record[8])   
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
#min_opacity = min_opacity - 5

#### making dictionary for converting to pandas dataframe to draw altair plot ###########
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
        data_list['component_label'].append(0) # make it zero so it is black
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(-barcode_info[i][2])
        data_list['opacity'].append(0.1)
        # barcode_info[i][3] = 0



# converting to pandas dataframe

data_list_pd = pd.DataFrame(data_list)
id_label = len(list(set(data_list['component_label'])))#unique_component_count
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

chart.save(current_directory +'altair_plot_test.html')
###################################  Histogram plotting #################################################################################

df = pd.DataFrame(csv_record_final)
df.to_csv(current_directory+'temp_csv.csv', index=False, header=False)
df = readCsv(current_directory+'temp_csv.csv')
os.remove(current_directory+'temp_csv.csv') # delete the intermediate file
df = preprocessDf(df)
p = plot(df)
outPath = current_directory+'histogram_test.html'
p.save(outPath)	


############################  Network Plot ######################
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
set1 = altairThemes.get_colour_scheme("Set1", unique_component_count)
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

g = nx.MultiDiGraph(directed=True) #nx.Graph()
for i in range (0, len(barcode_info)):
    label_str =  str(i)+'_c:'+str(barcode_info[i][3])+'_' # label of the node or spot is consists of: spot id, component number, type of the spot 
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



# scale the edge scores [0 to 1] to make plot work
score_list = []
for k in range (1, len(csv_record_final)-1):
    score_list.append(csv_record_final[k][8])

min_score = np.min(score_list)
max_score = np.max(score_list)

count_edges = 0
for k in range (1, len(csv_record_final)-1):
    i = csv_record_final[k][6]
    j = csv_record_final[k][7] 
    #if colors_point[i]==:
    #    continue
    ligand = csv_record_final[k][2]
    receptor = csv_record_final[k][3]
    edge_score = csv_record_final[k][8]
    edge_score = (max_score-edge_score)/(max_score-min_score)   
    title_str =  "L:" + ligand + ", R:" + receptor+ ", "+ str(edge_score) #+
    g.add_edge(int(i), int(j), label = title_str, color=colors_point[i], value=np.float64(edge_score)) #
    count_edges = count_edges + 1

print("total edges plotted: %d"%count_edges)

nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True     
nt.from_nx(g)
nt.save_graph('mygraph.html')
os.system('cp mygraph.html /cluster/home/t116508uhn/64630/mygraph.html')

################## The rest can be ignored ##############################################################################
 

# convert it to dot file to be able to convert it to pdf or svg format for inserting into the paper
from networkx.drawing.nx_agraph import write_dot
write_dot(g, "/cluster/home/t116508uhn/64630/test_interactive.dot")

'''
#These commands are to be executed in the linux terminal to convert the .dot file to pdf/svg:

cat test_interactive.dot.dot  | sed 's/ellipse/triangle/g'   | sed 's/tumor",/tumor",style="filled",/g'   | sed 's/L:\([^ ]\+\), R:/\1-/g'   | sed 's/label="[0-9][^"]*"/label=""/g' | awk -F'=' '{ if ($1 == "penwidth") {print $1 "=" ($2 ^ 6) ","} else {print $0 }}'   | tr '\n' ' '   | sed "s/;/\n/g"  > tmp

cat tmp   | dot -Kneato -n -y -Tpdf -Efontname="Arial" -Nlabel="" -Nwidth=1.5 -Nheight=1.5 -Npenwidth=8	> test.pdf

cat tmp   | dot -Kneato -n -y -Tsvg -Efontname="Arial" -Nlabel="" -Nwidth=1.5 -Nheight=1.5 -Npenwidth=8	> test.svg
'''
