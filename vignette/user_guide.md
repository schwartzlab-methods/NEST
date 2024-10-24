## Data Preprocess 
### Arguments
--data_name = Name of the dataset. Type is String. Required

--data_from = Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix. Type is String. Required.

--data_to = Path to save the input graph (to be passed to GAT). Type is String. Default is 'input_graph/'.

--metadata_to = Path to save the metadata. Type is String. default='metadata/'

--filter_min_cell = Minimum number of cells for gene filtering. Type is Int. default=1. 

--threshold_gene_exp = Threshold percentile for gene expression. Genes above this percentile are considered active. Type is float. default=98

--tissue_position_file = If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file. Type is String. default='None'

--spot_diameter = Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel). Type is float. default=89.43

--split = How many split sections? Type is Int. default=0. 
--neighborhood_threshold = Set neighborhood threshold distance in terms of same unit as spot diameter. Type is float. If not set, then it is set to four times of the spot diameter.

--database_path = Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database. Type is String. default='database/NEST_database.csv'

## Model Run ###
### Arguments
    --data_name = Name of the dataset. Type is String. Required.  
    --model_name = Provide a model name. Type is String. Required. 
    --run_id = Please provide a running ID, for example: 0, 1, 2, etc. Five runs are recommended. Type is Int.
    #=========================== default is set ======================================
    --num_epoch = Number of epochs or iterations for model training. Number of epochs or iterations for model training. Type is Int. default=60000
    
    --model_path = Path to save the model state. default='model/'  
    --embedding_path = Path to save the node embedding and attention scores. default='embedding_data/'
    --hidden = Hidden layer dimension or dimension of node embedding. default=512
    --training_data = Path to input graph. default='input_graph/'
    --heads = Number of heads in the attention model.  default=1
    --dropout = Dropout value for model training. default=0.
    --lr_rate = Model training parameter. default=0.00001
    --manual_seed = Set 'Yes' if you want to input your own seed for random number generation. type=str, default='no')
    --seed = Input the seed if --manual_seed=Yes.  , type=int )
    --split = Set 1 to split the input graph to accommodate the memory. type=int, default=0)
    --total_subgraphs = If --split = 1, then input the number of desired subgraph.  type=int, default=1)
    --metadata_to = Path to load the metadata. , type=str, default='metadata/'
    #=========================== optional ======================================
    --load', type=int, default=0, help='Load a previously saved model state')  
    --load_model_name', type=str, default='None' , help='Provide the model name that you want to reload')


## Result Postprocess 
### Arguments

## Output Visualize 
### Arguments
