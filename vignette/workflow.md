
This workflow will use the Visium sample on human lymph node as a use case (https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0) for the demonstration purpose. Please download the following two files:

a. The filtered feature matrix from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5

b. The spatial imaging data from here: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz (please unzip the spatial imaging data)

Both should be kept under the same directory, e.g., data/V1_Human_Lymph_Node_spatial/ directory. We have provided a default ligand-receptor database by merging the records from CellChat and NicheNet database. This is kept under 'database/' directory and will be used by NEST unless some other database is referred by the user.   

Change your current working directory to the downloaded NEST repository to run our model.

## Data preprocessing 

We first preprocess the data before passing it to NEST. It takes two main inputs: spatial transcriptomics dataset and a ligand-receptor database. Assuming that the spatial dataset is in "data/V1_Human_Lymph_Node_spatial/" directory and the ligand-receptor database is in 'database/NEST_database.csv', data preprocessing for input graph generation can simply be done as follows:
````
nest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/'
````
This method applies Quantile normalization on the gene expression matrix and generates an input graph where each spot in the ST data becomes a vertex in the graph and each vertex is connected with its neighbouring vertices. The neighborhood is decided based on the --neighborhood_threshold parameter with the default value: spot_diameter*4 (--spot_diameter=89.43 is the default value), i.e., a vertex will be directly connected with all other vertices who are positioned within that distance. Each connection represents a neighbourhood relation (corresponding to a ligand-receptor pair from the database) and number of total connections in an input graph depends on two more parameters:  --threshold_gene_exp and --filter_min_cell. The default values are --threshold_gene_exp=98 (for each cell, genes having expression above 98th percentile are considered active) and --filter_min_cell=5 (gene will be kept if it is expressed in at least 5 spots). 

Lower values for --threshold_gene_exp and --filter_min_cell will generate more connections and higher values will generate less number of connections in the input graph which largly decides how much GPU memory will the model use. We try to generate as many connections as we can to predict more CCC at the end. For example, the results presented in our paper was generated using this preprocessing command:
````
nest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/' --filter_min_cell=1 
````

The --data_name parameter is used to decide the target directories to save the processed data. For example, above command creates two folders in the current working directories: 
1. "input_graph/V1_Human_Lymph_Node_spatial/": Contains
   - V1_Human_Lymph_Node_spatial_adjacency_records: The input graph
   - V1_Human_Lymph_Node_spatial_cell_vs_gene_quantile_transformed: The quantile normalized gene expression matrix
2. "metadata/V1_Human_Lymph_Node_spatial/": Contains
   - V1_Human_Lymph_Node_spatial_barcode_info: A list having the information on barcodes and their coordinates.
   - V1_Human_Lymph_Node_spatial_self_loop_record: A dictionary object saving the information on barcodes having autocrine and juxtacrine (in case of spot based data) information. Used later for efficient visualization.      
  
Please use the argument --help to see all available input parameters.  

## Run NEST to generate CCC list

We recommend running the model at least 5 times with different seeds and then ensemble the outputs to get more consistent result. We can run the following commands in the terminal: 
```
nohup nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=2 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=3 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=4 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=5 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &
```
If you have enough GPU memory you can start running all of them in parallel. Model running takes a couple of hours to finish, so running the model in the background is recommended. If you are using Compute Canada servers, then a sample script to submit the GPU job can look like this: 
```
to be added
```


## Postprocessing output to generate a list of strong CCC

To post-process the model output, i.e., ensemble of multiple runs (through the product of ranks of predicted CCCs) and produce a list of top 20% highly ranked communications, we have to run the following commands:

````
nest postprocess --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --total_runs=5 
````

  In the command, we use --total_runs=5 assuming that the model is run five times (if you run the model just once, use --total_runs=1). The top 20% highly ranked communications are saved in a file: "output/V1_Human_Lymph_Node_spatial/NEST_V1_Human_Lymph_Node_spatial_top20percent.csv". This step uses --top_percent=20 by default. If you would prefer different percentage, e.g., 10%, please pass the parameter --top_percent=10 while running the command.


## Visualize the list of stronger CCC in different formats

We can use the CCC list "output/V1_Human_Lymph_Node_spatial/NEST_V1_Human_Lymph_Node_spatial_top20percent.csv" for various downstream analysis. Here we show how to visualize those CCC on the tissue surface. Instead of plotting all the CCC (which is 436,103 and highly memory consuming) on the tissue surface, we plot top 40,000 CCC as follows:

````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=40000
````

This step looks for the top 20% CCC list by default, but if you used different percentage, for example, top 10% CCC list in the previous postprocessing step, then please pass the parameter --top_percent=10 while running the command. 
This step generates the following four files under the directory 'output/V1_Human_Lymph_Node_spatial/': 
1. NEST_V1_Human_Lymph_Node_spatial_component_plot.html (in Altair) 
2. NEST_V1_Human_Lymph_Node_spatial_histogram_test.html (in Altair)
3. NEST_V1_Human_Lymph_Node_spatial_mygraph.html (in NetworkX)
4. NEST_V1_Human_Lymph_Node_spatial_test_interactive.dot

Although the NetworkX plot shows the appealing view of CCC, it can be very big and memory-consuming to open in the web-browser and inconvenient to share. Therefore we prefer to convert the corresponding *.dot file to a *.pdf and *.svg file by executing the following command (takes input the path of *.dot file as an argument): 

```
nest output_graph_picture output/V1_Human_Lymph_Node_spatial/NEST_V1_Human_Lymph_Node_spatial_test_interactive.dot
```
It will generate two files: edge_graph.svg and edge_graph.pdf in the current working directory, which are easy to view and share. 
 
The screenshots of the component plot and histograms are provided below (you can find the original files inside the vignette directory). 
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top40000.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top40000.png)

Please note that component 1 is dedicated to only those spots that are singleton, i.e., they only have self-loops but do not communicate with neighboring spots. In the plot above, we do not see component 1 because all the spots are talking to neighbors or are inactive (black), but there is no singleton.  

### Seeing comparatively stronger CCC by changing --top_edge_count parameter from high to low

If we want to pinpoint the location of particular communications on the tissue or which regions of the tissue involve which set of communications, we can gradually reduce the value of '--top_edge_count' to see more disjoint subgraphs as follows:   

````
nest visualize --dataname='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=10000
````
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top10000.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top10000.png)

````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000
````
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top3000.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top3000.png)

We see that the TGFB1 signaling is more prominent in the middle region with green components. 

### Supplying annotations to use different shapes for the cell types in the altair plot

Supplying the annotation file path in a *.csv format can assign different shape to different spot types as shown below: 

````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000 --annotation_file_path='data/V1_Human_Lymph_Node_spatial_annotation.csv'
````

![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top3000_annotated.png)

Next, we run the following command to convert the *.dot file to *.pdf and *.png file: 

```
nest output_graph_picture output/V1_Human_Lymph_Node_spatial/NEST_V1_Human_Lymph_Node_spatial_test_interactive.dot
```

It will generate two files: edge_graph.svg and edge_graph.pdf. The screenshot of the edge_graph.pdf is provided bellow:

![png file of the generated altair plot for top 40000 CCC](../images/edge_graph.pdf)

### Data preprocess when spatial transcriptomic data comes in non-Visium format

If data comes in a non-Visium format, we suggest to supply two files separately to the data preprocessing step: 
1. A *.mtx.gz file (having the Gene expression, Gene ID, and Barcodes)
2. Tissue position file in a *.csv format.

As an example, we request to download the Lung Adenocarcinoma (LUAD) sample from the Gene Expression Omnibus under accession number GSE189487 and keep it under the 'data/' directory: 

Then we can run the preprocess as follows:
```
nest preprocess --data_name='LUAD_TD1' --data_from='data/LUAD_GSM5702473_TD1/' --tissue_position_file='data/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv'
```

Similar approach can be followed for single-cell resolution, e.g., MERFISH data. 


   




