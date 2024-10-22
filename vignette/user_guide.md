## Data Preprocess 
### Arguments
--data_name = 
--data_from =
--data_to = 
--data_name = 
--data_from = 
--data_to =
--metadata_to =
--filter_min_cell = 
--threshold_gene_exp = 
--tissue_position_file = 
--spot_diameter = 
--split =
--neighborhood_threshold =
--database_path = 
    parser.add_argument( '--data_name', type=str, help='Name of the dataset', required=True)  
    parser.add_argument( '--data_from', type=str, required=True, help='Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix.')
    ################# default is set ################################################################
    parser.add_argument( '--data_to', type=str, default='input_graph/', help='Path to save the input graph (to be passed to GAT)')
    parser.add_argument( '--metadata_to', type=str, default='metadata/', help='Path to save the metadata')
    parser.add_argument( '--filter_min_cell', type=int, default=1 , help='Minimum number of cells for gene filtering') 
    parser.add_argument( '--threshold_gene_exp', type=float, default=98, help='Threshold percentile for gene expression. Genes above this percentile are considered active.')
    parser.add_argument( '--tissue_position_file', type=str, default='None', help='If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file.')
    parser.add_argument( '--spot_diameter', type=float, default=89.43, help='Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel).')
    parser.add_argument( '--split', type=int, default=0 , help='How many split sections?') 
    parser.add_argument( '--neighborhood_threshold', type=float, default=0 , help='Set neighborhood threshold distance in terms of same unit as spot diameter') 
    parser.add_argument( '--database_path', type=str, default='database/NEST_database.csv' , help='Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database.') 

## Model Run ###
### Arguments


## Result Postprocess 
### Arguments

## Output Visualize 
### Arguments
