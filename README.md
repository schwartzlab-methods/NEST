
# NEST
## Requirements:

###   Used python packages:

python == 3.7.2

numpy == 1.21.6

pytorch (with GPU support) == 3.7.2

torch == 1.12.0

torch_geometric== '2.1.0'

torch_sparse== '0.6.15'

torch_scatter== '2.0.9'

pickle == 4.0

scipy == 1.7.3

stlearn == '0.4.8'

qnorm == '0.8.1'

pandas == '1.3.5'

scanpy == '1.9.1'

altair == '4.2.0'

(Please download this additional package as well : https://github.com/schwartzlab-methods/altair-themes)

csv == '1.0'

matplotlib ==  '3.5.2'

gc == 

gzip ==

collections ==

###   System requirements: 
This model is tested on CentOS 7 and GPU servers with versions: Nvidia P100 and V100. This model is expected to run on any Linux server, e.g., Compute Canada as well.  
  
### Make the system recognize 'nest' command to run the model:

Download the NEST repository at your desired location and change your current working directory to NEST. Make the bash script 'nest' executable by executing following command:  
````
chmod +x nest
````

Then copy this executable bash script to your '/usr/local/bin/' or '/usr/bin/' directory so that your system can recognize this command:
````
cp nest /usr/local/bin/
````

(Instructions above are to be executed once only)

## Instruction to run NEST:

We use publicly available Visium sample on human lymph node (https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0) for the demonstration purpose. Please download the following two files:

a. The filtered feature matrix from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5

b. The spatial imaging data from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz (please unzip the spatial imaging data)

Both should be kept under the same directory, e.g., data/V1_Human_Lymph_Node_spatial/ directory. We have provided a default ligand-receptor database by merging the records from CellChat and NicheNet database. This is kept under 'database/' directory and will be used by NEST unless some other database is referred by the user.   

Change your current working directory to the downloaded NEST repository. Then execute following commands to run NEST on the human lymph node sample. 
   
1. NEST takes two main inputs: spatial transcriptomics dataset and a ligand-receptor database. Assuming that the spatial dataset is in "data/V1_Human_Lymph_Node_spatial/" directory and the ligand-receptor database is in 'database/NEST_database.csv', data preprocessing for input graph generation can be done as follows:
````
nest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/'
````
It will create two folders in the current working directories: "input_graph/V1_Human_Lymph_Node_spatial/" and "metadata/V1_Human_Lymph_Node_spatial/" to save the preprocessed input data. Please use the argument --help to see all available input parameters.  

2. To train a NEST model on the preprocessed 'V1_Human_Lymph_Node_spatial' data use following command with preferred model name. If the same experiment if repeated multiple times for model ensemble, each time a different run_id should be used and the run_id is expected to be consecutive. For example, if it is run five times then the run_id for the five runs should be 1, 2, 3, 4, and 5 respectively. By default the model will be trained for 80,000 epochs. Please use the argument --help to see all available input parameters. Please note that, the script will use GPU if avaiable, otherwise it will use CPU. Since this is a time consuming step, we suggest to run this step in the background and print the outputs in a separate log file as follows:

````
nohup nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &
````

  It will save trained model state with minimum loss in 'model/V1_Human_Lymph_Node_spatial/' and the corresponding attention scores and node embedding in 'embedding_data/V1_Human_Lymph_Node_spatial/'.   

3. To postprocess the model output, i.e., ensemble of multiple runs (through rank of product) and producing list of top 20% highly ranked communications we have to run following commands:

````
nest postprocess --dataname='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial' --total_runs=5 
````

  In the command, we use --total_runs=5 assuming that the model is run five times. The top 20% highly ranked communications are saved in a file named as 'V1_Human_Lymph_Node_spatial_top20percent.csv' in "output/V1_Human_Lymph_Node_spatial/".  

4. To visualize the output graph, i.e., finding connected components and ploting them, we run following command:

````
nest visualize --dataname='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial'
````

  This will generate output in four formats: altair plot, histogram plot, networkx plot, and dot file for pdf generation. 

## Vignette
For detail explanation on the available parameters and their usage, please see the vignettes:

1. [Generate active CCC lists given a spatial transcriptomics data](vignette/workflow.md)
2. [Filter CCC list for specific region / cell type / specific ligand-receptor pair](vignette/filter_ccc_list_for_type_region.html)
   
    
