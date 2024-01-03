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
    # =========================== args ===============================
    parser.add_argument( '--data_name', type=str, default='V1_Breast_Cancer_Block_A_Section_1', help="'MERFISH' or 'V1_Breast_Cancer_Block_A_Section_1") 
    parser.add_argument( '--data_path', type=str, default='generated_data/', help='data path')
    parser.add_argument( '--model_path', type=str, default='model') 
    parser.add_argument( '--embedding_data_path', type=str, default='Embedding_data') 
    #parser.add_argument( '--result_path', type=str, default='results') 
    parser.add_argument( '--load', type=int, default=0, help='Load pretrained DGI model')
    parser.add_argument( '--num_epoch', type=int, default=5000, help='numebr of epoch in training DGI')
    parser.add_argument( '--hidden', type=int, default=256, help='hidden channels in DGI') 
    parser.add_argument( '--retrain', type=int, default=0 , help='Run which clustering at the end')
    parser.add_argument( '--model_load_path', type=str, default='model')
    parser.add_argument( '--model_name', type=str, default='r1')
    parser.add_argument( '--training_data', type=str, default='provide please')
    parser.add_argument( '--heads', type=int, default=1)
    parser.add_argument( '--num_cells', type=int, default=1)
    parser.add_argument( '--options', type=str)
    parser.add_argument( '--withFeature', type=str, default='r1') 
    parser.add_argument( '--workflow_v', type=int, default=1)
    parser.add_argument( '--datatype', type=str)
    parser.add_argument( '--dropout', type=float, default=0)
    parser.add_argument( '--lr_rate', type=float, default=0.00001)
    parser.add_argument( '--manual_seed', type=str, default='no')
    parser.add_argument( '--seed', type=int)
    args = parser.parse_args() 

    args.embedding_data_path = args.embedding_data_path + args.data_name +'/'
    args.model_path = args.model_path + args.data_name +'/'
    #args.result_path = args.result_path +'/'+ args.data_name +'/'
    args.model_load_path = args.model_load_path + args.data_name +'/'

    print(args.model_name+', '+str(args.heads)+', '+args.training_data+', '+str(args.hidden) )

    if args.manual_seed == 'yes':
        torch.manual_seed(args.seed)
        random.seed(args.seed)
        np.random.seed(args.seed)


    start_time = time.time()
    if not os.path.exists(args.embedding_data_path):
        os.makedirs(args.embedding_data_path) 
    if not os.path.exists(args.model_path):
        os.makedirs(args.model_path) 
    #args.result_path = args.result_path+'/'
    #if not os.path.exists(args.result_path):
    #    os.makedirs(args.result_path) 
    print ('------------------------Model and Training Details--------------------------')
    print(args) 

    
    from train_CCC_gat import CCC_on_ST
    CCC_on_ST(args)
    end_time = time.time() - start_time
    print('time elapsed %g min'%(end_time/60))
    


    
