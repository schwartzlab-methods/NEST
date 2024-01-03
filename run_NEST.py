import os
import sys
import numpy as np
import torch
from datetime import datetime 
import time
import random
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # =========================== must be provided ===============================
    parser.add_argument( '--data_name', type=str, help='Name of the dataset') #default='PDAC_64630', 
    parser.add_argument( '--model_name', type=str, help='Provide a model name')
    #parser.add_argument( '--num_cells', type=int, help='Number of cells or spots in the dataset')
    #=========================== default is set ======================================
    parser.add_argument( '--num_epoch', type=int, default=50000, help='Number of epochs or iterations for model training')
    parser.add_argument( '--model_path', type=str, default='model', help='Path to save the model state') # We do not need this for output generation  
    parser.add_argument( '--embedding_path', type=str, default='embedding_data', help='Path to save the node embedding and attention scores') 
    parser.add_argument( '--hidden', type=int, default=512, help='Hidden layer dimension (dimension of node embedding)')
    parser.add_argument( '--training_data', type=str, default='input_graph/', help='Path to input graph. ')
    parser.add_argument( '--heads', type=int, default=1, help='Number of heads in the attention model')
    parser.add_argument( '--dropout', type=float, default=0)
    parser.add_argument( '--lr_rate', type=float, default=0.00001)
    parser.add_argument( '--manual_seed', type=str, default='no')
    parser.add_argument( '--seed', type=int )
    #=========================== optional ======================================
    parser.add_argument( '--load_init', type=int, default=0, help='Load initial model state for the given model_name')  
    parser.add_argument( '--retrain', type=int, default=0 , help='Load last model state to retrain for the given model_name')

    args = parser.parse_args() 

    #parser.add_argument( '--options', type=str)
    #parser.add_argument( '--withFeature', type=str, default='r1') 
    #parser.add_argument( '--workflow_v', type=int, default=1)
    #parser.add_argument( '--datatype', type=str)


    args.training_data = args.training_data + args.data_name + '/'
    args.embedding_path = args.embedding_path + args.data_name +'/'
    args.model_path = args.model_path + args.data_name +'/'


    print(args.data_name+', '+str(args.heads)+', '+args.training_data+', '+str(args.hidden) )

    if args.manual_seed == 'yes':
        torch.manual_seed(args.seed)
        random.seed(args.seed)
        np.random.seed(args.seed)


    if not os.path.exists(args.embedding_data_path):
        os.makedirs(args.embedding_data_path) 
    if not os.path.exists(args.model_path):
        os.makedirs(args.model_path) 

    print ('------------------------Model and Training Details--------------------------')
    print(args) 

    
    from train_CCC_gat import CCC_on_ST
    start_time = time.time()
    CCC_on_ST(args)
    end_time = time.time() - start_time
    print('time elapsed %g min'%(end_time/60))
    


    
