from scipy import sparse
import pickle
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import Linear, DeepGraphInfomax
from torch_geometric.data import Data, DataLoader
import gzip

from GATv2Conv_NEST import GATv2Conv

def get_graph(X, training_data_name):

    # create sparse matrix
    f = gzip.open(training_data_name , 'rb')
    row_col, edge_weight, lig_rec = pickle.load(f)


    print("row_col %d"%len(row_col))
    #print(edge_weight)
    ###########

    edge_index = torch.tensor(np.array(row_col), dtype=torch.long).T
    edge_attr = torch.tensor(np.array(edge_weight), dtype=torch.float)
    print('X shape ')
    print(X.shape)
    graph_bags = []
    graph = Data(x=torch.tensor(X, dtype=torch.float), edge_index=edge_index, edge_attr=edge_attr)

    graph_bags.append(graph)
    print('get graph done')
    return graph_bags


class Encoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, heads, dropout):
        super(Encoder, self).__init__()
        print('incoming channel %d'%in_channels)

        heads = heads
        self.conv =  GATv2Conv(in_channels, hidden_channels, edge_dim=3, heads=heads, concat = False)#, dropout=dropout)
        self.conv_2 =  GATv2Conv(hidden_channels, hidden_channels, edge_dim=3, heads=heads, concat = False)#, dropout=0)

        self.attention_scores_mine_l1 = 'attention_l1'
        self.attention_scores_mine_unnormalized_l1 = 'attention_unnormalized_l1'

        self.attention_scores_mine = 'attention'
        self.attention_scores_mine_unnormalized = 'attention_unnormalized'

        #self.prelu = nn.Tanh(hidden_channels)
        self.prelu = nn.PReLU(hidden_channels)


    def forward(self, data):

        # layer 1
        x, attention_scores, attention_scores_unnormalized = self.conv(data.x, data.edge_index, edge_attr=data.edge_attr, return_attention_weights = True)
        self.attention_scores_mine_l1 = attention_scores
        self.attention_scores_mine_unnormalized_l1 = attention_scores_unnormalized
#        x = F.dropout(x, p=0.5, training=self.training)
#        x = F.elu(x)

        # layer 2
        x, attention_scores, attention_scores_unnormalized  = self.conv_2(x, data.edge_index, edge_attr=data.edge_attr, return_attention_weights = True)  # <---- ***
        self.attention_scores_mine = attention_scores #self.attention_scores_mine_l1 #attention_scores
        self.attention_scores_mine_unnormalized = attention_scores_unnormalized #self.attention_scores_mine_unnormalized_l1 #attention_scores_unnormalized

        ###############################
        x = self.prelu(x)

        return x #, attention_scores

class my_data():
    def __init__(self, x, edge_index, edge_attr):
        self.x = x
        self.edge_index = edge_index
        self.edge_attr = edge_attr


def corruption(data):
    #print('inside corruption function')
    x = data.x[torch.randperm(data.x.size(0))]
    return my_data(x, data.edge_index, data.edge_attr)


def train_DGI(args, data_loader, in_channels):

    loss_curve = np.zeros((args.num_epoch//500+1))
    loss_curve_counter = 0

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    DGI_model = DeepGraphInfomax(
        hidden_channels=args.hidden,
        encoder=Encoder(in_channels=in_channels, hidden_channels=args.hidden, heads=args.heads, dropout = args.dropout),
        summary=lambda z, *args, **kwargs: torch.sigmoid(z.mean(dim=0)),
        corruption=corruption).to(device)
    #print('initialized DGI model')
    #print(DGI_model.encoder.attention_scores_mine)
    #DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=0.005, weight_decay=5e-4)
    DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=args.lr_rate) #1e-5)#5 #6
    #DGI_optimizer = torch.optim.RMSprop(DGI_model.parameters(), lr=1e-5)
    DGI_filename = args.model_path+'DGI'+ args.model_name  +'.pth.tar'
    if args.load:
        DGI_model.load_state_dict(torch.load(DGI_filename))
    else:
	import datetime
        start_time = datetime.datetime.now()
        min_loss=10000
        if args.retrain==1:
            DGI_load_path = args.model_load_path+'DGI'+ args.model_name+'.pth.tar'
            DGI_model.load_state_dict(torch.load(DGI_load_path))
            DGI_optimizer.load_state_dict(torch.load(args.model_path+'DGI_optimizer_'+ args.model_name  +'.pth.tar'))

        #print('Saving init model state ...')
        torch.save(DGI_model.state_dict(), args.model_path+'DGI_init'+ args.model_name  + '.pth.tar')
        torch.save(DGI_optimizer.state_dict(), args.model_path+'DGI_optimizer_init'+ args.model_name  + '.pth.tar')
        #print('training starts ...')
        for epoch in range(args.num_epoch):
            DGI_model.train()
            DGI_optimizer.zero_grad()
            DGI_all_loss = []

            for data in data_loader:
                data = data.to(device)
                '''x_old, attention_score = DGI_model(data=data)
                print('model output')
                print(x_old)
                pos_z = x_old[0]
                neg_z = x_old[1]
                summary = x_old[2]'''

                pos_z, neg_z, summary = DGI_model(data=data)
                #print('epoch %d '%epoch)
                #print(DGI_model.encoder.attention_scores_mine)
                DGI_loss = DGI_model.loss(pos_z, neg_z, summary)
                DGI_loss.backward()
                DGI_all_loss.append(DGI_loss.item())
                DGI_optimizer.step()

#                loss_curve[loss_curve_counter] = np.mean(DGI_all_loss)
#                loss_curve_counter = loss_curve_counter + 1

            if ((epoch)%500) == 0:
                print('Epoch: {:03d}, Loss: {:.4f}'.format(epoch+1, np.mean(DGI_all_loss)))
                loss_curve[loss_curve_counter] = np.mean(DGI_all_loss)
                loss_curve_counter = loss_curve_counter + 1

                if np.mean(DGI_all_loss)<min_loss:
                    min_loss=np.mean(DGI_all_loss)
                    torch.save(DGI_model.state_dict(), DGI_filename)
                    torch.save(DGI_optimizer.state_dict(), args.model_path+'DGI_optimizer_'+ args.model_name  +'.pth.tar')
                    save_tupple=[pos_z, neg_z, summary]
                    saved_attention = DGI_model.encoder.attention_scores_mine
                    saved_attention_unnormalized = DGI_model.encoder.attention_scores_mine_unnormalized
                    saved_attention_l1 = DGI_model.encoder.attention_scores_mine_l1
                    saved_attention_unnormalized_l1 = DGI_model.encoder.attention_scores_mine_unnormalized_l1
                    print(DGI_model.encoder.attention_scores_mine_unnormalized_l1[0:10])











#            if ((epoch)%60000) == 0:
#                DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=1e-6)  #5 #6

        end_time = datetime.datetime.now()

#        torch.save(DGI_model.state_dict(), DGI_filename)
        print('Training time in seconds: ', (end_time-start_time).seconds)
        DGI_model.load_state_dict(torch.load(DGI_filename))
        print("debug loss")
        DGI_loss = DGI_model.loss(pos_z, neg_z, summary)
        print("debug loss latest tupple %g"%DGI_loss.item())
        DGI_loss = DGI_model.loss(save_tupple[0], save_tupple[1], save_tupple[2])
        print("debug loss min loss tupple %g"%DGI_loss.item())

        DGI_model.encoder.attention_scores_mine = saved_attention
        DGI_model.encoder.attention_scores_mine_unnormalized = saved_attention_unnormalized

        DGI_model.encoder.attention_scores_mine_l1 = saved_attention_l1
        DGI_model.encoder.attention_scores_mine_unnormalized_l1 = saved_attention_unnormalized_l1


        logfile=open(args.model_path+'DGI'+ args.model_name+'_loss_curve.csv', 'wb')
        np.savetxt(logfile,loss_curve, delimiter=',')
        logfile.close()

    return DGI_model

