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

## Model Run 
### Arguments
    --data_name = Name of the dataset. Type is String. Required.  
    --model_name = Provide a model name. Type is String. Required. 
    --run_id = Please provide a running ID, for example: 0, 1, 2, etc. Five runs are recommended. Type is Int.
    #=========================== default is set ======================================
    --num_epoch = Number of epochs or iterations for model training. Number of epochs or iterations for model training. Type is Int. default=60000
    
    --model_path = Path to save the model state. Type is String. default='model/'  
    --embedding_path = Path to save the node embedding and attention scores. Type is String. default='embedding_data/'
    --hidden = Hidden layer dimension or dimension of node embedding. Type is Int. default=512
    --training_data = Path to input graph. Type is String. default='input_graph/'
    --heads = Number of heads in the attention model. Type is Int. default=1
    --dropout = Dropout value for model training. Type is Float. default=0.
    --lr_rate = Model training parameter. Type is Float. default=0.00001
    --manual_seed = Set 'Yes' if you want to input your own seed for random number generation. Type is String. default='no'. 
    --seed = Input the seed if --manual_seed=Yes. Type is Int. 
    --split = Set 1 to split the input graph to accommodate the memory. Type is Int. default=0. 
    --total_subgraphs = If --split = 1, then input the number of desired subgraph. Type is Int. default=1
    --metadata_to = Path to load the metadata. Type is String. default='metadata/'
    #=========================== optional ======================================
    --load = Set 1 to load a previously saved model state. Type is Int. default=0.  
    --load_model_name = Provide the model name that you want to reload. Type is String. default='None'

## Model Prediction Ensemble 
### Arguments
    --data_namec = The name of dataset. Type is String. Required.
    --model_name = Name of the trained model.  Type is String. Required.
    --total_runs = How many runs for ensemble (at least 2 are preferred). Type is Int. Required.
    #######################################################################################################
    --embedding_path = Path to grab the attention scores from.  Type is String. default='embedding_data/'
    --metadata_from =  Path to grab the metadata.  Type is String. default='metadata/' 
    --data_from = Path to grab the input graph from (to be passed to GAT).  Type is String. default='input_graph/'
    --output_path = Path to save the visualization results, e.g., histograms, graph etc.  Type is String. default='output/'
    --top_percent = Top N percentage communications to pick. Type is Float. default=20.
    

## Output Visualize 
### Arguments
    --data_name = The name of dataset. Type is String.  Required.
    --model_name = Name of the trained model. Type is String. Required.
    --top_edge_count = Number of the top communications to plot. To plot all insert -1. Type is Int. default=1500.
    --top_percent = Top N percentage communications to pick. Type is Int. default=20    
    --metadata_from = Path to grab the metadata. Type is String. default='metadata/' 
    --output_path = Path to save the visualization results, e.g., histograms, graph etc. Type is String. default='output/'
    --barcode_info_file = Path to load the barcode information file produced during data preprocessing step. Type is String. default=''
    --annotation_file_path = Path to load the annotation file in csv format (if available). Type is String. default=''
    --selfloop_info_file = Path to load the selfloop information file produced during data preprocessing step. Type is String. default=''
    --top_ccc_file = Path to load the selected top CCC file produced during data postprocessing step. Type is String. default=''
    --output_name = Output file name prefix according to user's choice. Type is String.  default=''
    --filter = Set --filter=1 if you want to filter the CCC. Type is Int. default=0
    --filter_by_ligand_receptor = Set ligand-receptor pair, e.g., --filter_by_ligand_receptor="CCL19-CCR7" if you want to filter the CCC by LR pair. Type is String. default=''
    --filter_by_annotation = Set cell or spot type, e.g., --filter_by_annotation="T-cell" if you want to filter the CCC. Type is String. default=''
    --filter_by_component = Set component id, e.g., --filter_by_component=9 if you want to filter by component id. Type is String. Type is Int.  default=-1
    --histogram_attention_score = Set --histogram_attention_score=1 if you want to sort the histograms of CCC by attention score. Type is String. default=-1
