# NO NEED NOW ## This file will be removed
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from sklearn import metrics
from scipy import interp
from sklearn.metrics import roc_curve, auc, roc_auc_score
import gzip
import numpy as np
from scipy import sparse
import pickle
import pandas as pd
import scanpy as sc
import anndata as ad

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import DataLoader


rootPath = os.path.dirname(sys.path[0])
os.chdir(rootPath+'/CCC_project')

def CCC_on_ST(args):
    from CCC_gat import get_graph, train_DGI
    # Parameters
    batch_size = 1  # Batch size
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)

    print("-----------Deep Graph Infomax-------------")
    data_list, num_feature = get_graph(args.data_path+args.training_data)
    data_loader = DataLoader(data_list, batch_size=batch_size)

    # train the model
    DGI_model = train_DGI(args, data_loader=data_loader, in_channels=num_feature)

    # training done. Now just save the final result here 
    for data in data_loader:
        data.to(device)
        X_embedding, _, _ = DGI_model(data)
        X_embedding = X_embedding.cpu().detach().numpy()
        X_embedding_filename =  args.embedding_data_path + args.model_name + '_Embed_X.npy'
        np.save(X_embedding_filename, X_embedding)
          
        X_attention_index = DGI_model.encoder.attention_scores_mine[0]
        X_attention_index = X_attention_index.cpu().detach().numpy()

        X_attention_score_normalized_l1 = DGI_model.encoder.attention_scores_mine_l1[1]
        X_attention_score_normalized_l1 = X_attention_score_normalized_l1.cpu().detach().numpy()

        X_attention_score_normalized = DGI_model.encoder.attention_scores_mine[1]
        X_attention_score_normalized = X_attention_score_normalized.cpu().detach().numpy()

        X_attention_score_unnormalized = DGI_model.encoder.attention_scores_mine_unnormalized
        X_attention_score_unnormalized = X_attention_score_unnormalized.cpu().detach().numpy()

        X_attention_score_unnormalized_l1 = DGI_model.encoder.attention_scores_mine_unnormalized_l1
        X_attention_score_unnormalized_l1 = X_attention_score_unnormalized_l1.cpu().detach().numpy()

#        X_attention_bundle = [X_attention_index, X_attention_score, X_attention_score_unnormalized]
#        X_attention_filename =  args.embedding_data_path + args.model_name + '_attention.npy'
#        np.save(X_attention_filename, X_attention_bundle)

        print('making the bundle to save')
        X_attention_bundle = [X_attention_index, X_attention_score_normalized_l1, X_attention_score_unnormalized, X_attention_score_unnormalized_l1, X_attention_score_normalized]
        X_attention_filename =  args.embedding_data_path + args.model_name + '_attention_l1.npy'
        np.save(X_attention_filename, X_attention_bundle)

        print('both attention score saved')
        print(X_attention_bundle[3][0:10]) 

            
            
    print("DGI is finished")
