
# NEST
## Requirements:
###   System requirements: 
[to be]
###   Required packages in python: 
[to be] 
  
## Instruction to run NEST:
We use publicly available Visium sample on human lymph node (https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0) for the demonstration purpose. 

Please download the following two files:

a. The filtered feature matrix from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5

b. The spatial imaging data from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz
Both should be kept under the same directory, e.g., data/V1_Human_Lymph_Node_spatial/ directory.
   
1. Assuming that the spatial dataset is in "data/V1_Human_Lymph_Node_spatial/" directory, data preprocessing for input graph generation can be done as follows:
````
    python data_preprocess_NEST.py --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/'
````
It will create two folders in the current working directories: "input_graph/V1_Human_Lymph_Node_spatial/" and "metadata/V1_Human_Lymph_Node_spatial/" to save the preprocessed input data. Please use the argument --help to see all available input parameters.  

2. To train a NEST model on the preprocessed 'V1_Human_Lymph_Node_spatial' data use following command with preferred model name. If the same experiment if repeated multiple times for model ensemble, each time a different run_id should be used and the run_id is expected to be consecutive. For example, if it is run five times then the run_id for the five runs should be 1, 2, 3, 4, and 5 respectively. By default the model will be trained for 80,000 epochs. Please use the argument --help to see all available input parameters. Please note that, the script will use GPU if avaiable, otherwise it will use CPU. 

````
    nohup python -u run_NEST.py  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
    nohup python -u run_NEST.py  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
    nohup python -u run_NEST.py  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
    nohup python -u run_NEST.py  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
    nohup python -u run_NEST.py  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &

````

  It will save trained model state with minimum loss in 'model/V1_Human_Lymph_Node_spatial/' and the corresponding attention scores and node embedding in 'embedding_data/V1_Human_Lymph_Node_spatial/'.   

3. To postprocess the model output, i.e., ensemble of multiple runs (through rank of product) and producing list of top 20% highly ranked communications we have to run following commands:

````
    python output_postprocess_NEST.py --dataname='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial' --total_runs=5 
````

  In the command, we use --total_runs=5 assuming that the model is run five times. The top 20% highly ranked communications are saved in a file named as 'V1_Human_Lymph_Node_spatial_top20percent.csv' in "output/V1_Human_Lymph_Node_spatial/".  

4. To visualize the output graph, i.e., finding connected components and ploting them, we run following command:

````
    nohup python -u output_visualization_NEST.py --dataname='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial'  > output.log &
````

  This will generate output in four formats: altair plot, histogram plot, networkx plot, and dot file for pdf generation. 
