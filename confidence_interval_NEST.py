print('package loading')
import numpy as np
import csv
import pickle
import statistics
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import median_abs_deviation
from scipy.stats import skew
from collections import defaultdict
import pandas as pd
from random import choices
import gzip
#from kneed import KneeLocator
import copy 
import argparse
import gc
import os


##########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--ccc_list_path', type=str, default='output/', help='CCC list path') # default='PDAC_64630',
    parser.add_argument( '--model_name', type=str, help='Name of the trained model', required=True)
    parser.add_argument( '--data_name', type=str, help='Name of the data.', required=True)
    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--top20', type=int, default=-1, help='set to 1 to output the CCC having the attention score within the 95th percent confidence interval of top20 cutoff')
    parser.add_argument( '--std', type=int, default=-1, help='set to 1 to output the CCC having the attention score within the 95th percent confidence interval of standard deviation')

    
    args = parser.parse_args()
    if args.output_path=='output/':
        args.output_path = args.output_path + args.data_name + '/'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    if args.ccc_list_path=='output/':
        args.ccc_list_path = args.ccc_list_path + args.data_name + '/'
        
    ccc_list_filename = args.ccc_list_path + args.model_name+'_allCCC.csv'
##################### get metadata: barcode_info ###################################
    df = pd.read_csv(ccc_list_filename, sep=",")
    csv_record = df.values.tolist()
    df_column_names = list(df.columns)
    # from_cell, to_cell, ligand_gene, receptor_gene, rank, component, from_id, to_id,  attention_score
    total_edge = len(csv_record)
    edge_index = []
    attention_score_original = []
    for i in range (0, len(csv_record)):
        edge_index.append(i)
        attention_score_original.append(csv_record[i][8])
        
    new_sets = []
    for i in range (0, num_new_set):
        temp_set = choices(edge_index, k=(2*total_edge)/3)
        new_sets.append(new_sets)

    # replace new_sets[i] with attention scores
    for i in range (0, num_new_set):
        for j in range(0, len(new_sets[i])):
            new_sets[i][j] =  csv_record[new_sets[i][j]][8]

    # now new_sets[i] is a set of attention scores.
    if args.top20 == 1:
        top20_list = []
        for i in range (0, len(new_sets)):
            top20_score = np.percentile(new_sets[i],80) # top 20% means above 80th percentile since the list is in 
                                                        # ascending order and higher value means more strong
            top20_list.append(top20_score)

        lower_limit = np.percentile(top20_list,2.7)
        upper_limit = np.percentile(top20_list, 97)
        print('95th percent confidence interval is: %g to %g'%(lower_limit, upper_limit))
        # output only ccc which are between this confidence interval
        csv_record_final = []
        for i in range (0, len(csv_record)):
            if lower_limit <= csv_record[i][8] and csv_record[i][8] <= upper_limit:
                csv_record_final.append(csv_record[i])
                
        csv_record_final = [df_column_names] + csv_record_final
        df = pd.DataFrame(csv_record_final) # output 4
        df.to_csv(args.output_path + args.model_name+'_CCC_list_top20_confidence.csv', index=False, header=False)
        
    if args.std == 1:
        std_list = []
        for i in range (0, len(new_sets)):
            std_score = np.std(new_sets[i]) # top 20% means above 80th percentile since the list is in 
                                            # ascending order and higher value means more strong
            std_list.append(std_score)

        lower_limit = np.percentile(std_list,2.7)
        upper_limit = np.percentile(std_list, 97)
        print('95th percent confidence interval is: %g to %g'%(lower_limit, upper_limit))
        # output only ccc which are between this confidence interval
        csv_record_final = []
        mean_score = np.mean(attention_score_original)
        for i in range (0, len(csv_record)):
            if lower_limit <= np.abs(mean_score-csv_record[i][8]) and np.abs(mean_score-csv_record[i][8]) <= upper_limit:
                csv_record_final.append(csv_record[i])
                
        csv_record_final = [df_column_names] + csv_record_final
        df = pd.DataFrame(csv_record_final) # output 4
        df.to_csv(args.output_path + args.model_name+'_CCC_list_std_confidence.csv', index=False, header=False)
    
    ###########################################################################################################################################

