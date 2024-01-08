
# NEST
## Requirements:
###   System requirements: 
[to be]
###   Required packages in python: 
[to be]
  
## Instruction to run NEST:
   
1. Assuming that the spatial dataset is in "data/PDAC_64630/" directory, data preprocessing for input graph generation can be done as follows:  

``    nohup python -u data_preprocess_NEST.py --dataname='PDAC_64630' > output.log & ``    

It will create two folders in the current working directories: "input_graph/PDAC_64630/" and "metadata/PDAC_64630/" to save the preprocessed input data. Please use the argument --help to see all available input parameters.  

2. To train a NEST model on the preprocessed 'PDAC_64630' data use following command with preferred model name. If the same experiment if repeated multiple times for model ensemble, each time a different run_id should be used and the run_id is expected to be consecutive. For example, if it is run five times then the run_id for the five runs should be 1, 2, 3, 4, and 5 respectively. By default the model will be trained for 60,000 epochs. Please use the argument --help to see all available input parameters. Please note that, the script will use GPU if avaiable, otherwise it will use CPU. 

``nohup python -u run_NEST.py --dataname='PDAC_64630' --model_name="PDAC_64630_NEST" --run_id=1 > output.log & ``

  It will save trained model state with minimum loss in 'model/PDAC_64630/' and the corresponding attention scores and node embedding in 'embedding_data/PDAC_64630/'.   

3. To postprocess the model output, i.e., ensemble of multiple runs (through rank of product) and producing list of top 20% highly ranked communications we have to run following commands:

``nohup python -u output_postprocess_NEST.py --dataname='PDAC_64630' --total_runs=5 > output.log & ``

  In the command, we use --total_runs=5 assuming that the model is run five times. The top 20% highly ranked communications are saved in a file named as 'PDAC_64630_top20percent.csv' in "output/PDAC_64630/".  

4. To visualize the output graph, i.e., finding connected components and ploting them, we run following command:

``nohup python -u output_visualization_NEST.py --dataname='PDAC_64630' > output.log & ``

  This will generate output in four formats: altair plot, histogram plot, networkx plot, and dot file for pdf generation. 
