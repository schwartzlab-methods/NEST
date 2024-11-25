We provide the arguments for running various steps of CellNEST model along with the preferred values as default. 

## Data Preprocess 
### Arguments
    --data_name = Name of the dataset. Type is String. Required
    --data_from = Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix. Type is String. Required.
    --data_to = Path to save the input graph (to be passed to GAT). Type is String. Default = 'input_graph/'.
    --metadata_to = Path to save the metadata. Type is String. default='metadata/'
    --filter_min_cell = Minimum number of cells for gene filtering. Type is Int. Default=1. 
    --threshold_gene_exp = Threshold percentile for gene expression. Genes above this percentile are considered active. Type is float. Default=98
    --tissue_position_file = If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file. Type is String. Default='None'
    --spot_diameter = Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel). Type is float. Default=89.43
    --split = How many split sections? Type is Int. default=0. 
    --neighborhood_threshold = Set neighborhood threshold distance in terms of same unit as spot diameter. Type is float. If not set, then it is set to four times of the spot diameter.
    --database_path = Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database. Type is String. Default='database/NEST_database.csv'


We have the following optional arguments for integrating intracellular pathways with the ligand-receptor coexpression:
### Arguments for integrating intracellular signaling with inter-cellular signals
    --intra_database_path = Provide your desired ligand-receptor database path here. Default database is a postprocessed NicheNet database as explained in the paper. Type is String. Default='database/nichenet_pathways_NEST.csv'. 
    --add_intra = Set it to 1 for intracellular signaling pathway. Type is Int. Default=1
    --num_hops = Maximum number of hops for intracellular signaling pathway search. Type is Int. Default=10
    --threshold_gene_exp_intra = Threshold percentile for gene expression. Genes above this percentile are considered active. Type is float. Default=20

## Model Run 
### Arguments
    --data_name = Name of the dataset. Type is String. Required.  
    --model_name = Provide a model name. Type is String. Required. 
    --run_id = Please provide a running ID, for example: 0, 1, 2, etc. Five runs are recommended. Type is Int.
    #=========================== default is set ======================================
    --num_epoch = Number of epochs or iterations for model training. Number of epochs or iterations for model training. Type is Int. Default=60000
    --model_path = Path to save the model state. Type is String. Default='model/'  
    --embedding_path = Path to save the node embedding and attention scores. Type is String. Default='embedding_data/'
    --hidden = Hidden layer dimension or dimension of node embedding. Type is Int. Default=512
    --training_data = Path to input graph. Type is String. Default='input_graph/'
    --heads = Number of heads in the attention model. Type is Int. Default=1
    --dropout = Dropout value for model training. Type is Float. Default=0.
    --lr_rate = Model training parameter. Type is Float. Default=0.00001
    --manual_seed = Set 'Yes' if you want to input your own seed for random number generation. Type is String. Default='no'. 
    --seed = Input the seed if --manual_seed=Yes. Type is Int. 
    --split = Set 1 to split the input graph to accommodate the memory. Type is Int. Default=0. 
    --total_subgraphs = If --split = 1, then input the number of desired subgraph. Type is Int. Default=1
    --metadata_to = Path to load the metadata. Type is String. Default='metadata/'
    #=========================== optional ======================================
    --load = Set 1 to load a previously saved model state. Type is Int. Default=0.  
    --load_model_name = Provide the model name that you want to reload. Type is String. Default='None'

## Model Prediction Ensemble 
### Arguments
    --data_name = The name of dataset. Type is String. Required.
    --model_name = Name of the trained model.  Type is String. Required.
    --total_runs = How many runs for ensemble (at least 2 are preferred). Type is Int. Required.
    #######################################################################################################
    --embedding_path = Path to grab the attention scores from.  Type is String. Default='embedding_data/'
    --metadata_from =  Path to grab the metadata.  Type is String. Default='metadata/' 
    --data_from = Path to grab the input graph from (to be passed to GAT).  Type is String. Default='input_graph/'
    --output_path = Path to save the visualization results, e.g., histograms, graph etc.  Type is String. Default='output/'
    --top_percent = Top N percentage communications to pick. Type is Float. Default=20.
    

## Output Visualize 
### Arguments
    --data_name = The name of dataset. Type is String.  Required.
    --model_name = Name of the trained model. Type is String. Required.
    --top_edge_count = Number of the top communications to plot. To plot all insert -1. Type is Int. Default=1500.
    --top_percent = Top N percentage communications to pick. Type is Int. Default=20    
    --metadata_from = Path to grab the metadata. Type is String. Default='metadata/' 
    --output_path = Path to save the visualization results, e.g., histograms, graph etc. Type is String. Default='output/'
    --barcode_info_file = Path to load the barcode information file produced during data preprocessing step. Type is String. Default path points to the metadata folder.
    --annotation_file_path = Path to load the annotation file in csv format (if available). Type is String. Default path points to the metadata folder.
    --selfloop_info_file = Path to load the selfloop information file produced during data preprocessing step. Type is String. Default path points to the metadata folder.
    --top_ccc_file = Path to load the selected top CCC file produced during data postprocessing step. Type is String. Default path points to the output folder.
    --output_name = Output file name prefix according to user's choice. Type is String.  Default path points to the metadata folder.
    --filter = Set --filter=1 if you want to filter the CCC. Type is Int. Default=0
    --filter_by_ligand_receptor = Set ligand-receptor pair, e.g., --filter_by_ligand_receptor="CCL19-CCR7" if you want to filter the CCC by LR pair. Type is String. 
    --filter_by_annotation = Set cell or spot type, e.g., --filter_by_annotation="T-cell" if you want to filter the CCC. Type is String.
    --filter_by_component = Set component id, e.g., --filter_by_component=9 if you want to filter by component id. Type is String. Type is Int.  Default=-1
    --histogram_attention_score = Set --histogram_attention_score=1 if you want to sort the histograms of CCC by attention score. Type is String. Default=-1