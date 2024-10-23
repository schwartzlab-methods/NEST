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


## Result Postprocess 
### Arguments

## Output Visualize 
### Arguments
