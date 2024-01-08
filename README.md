
# NEST
Requirements:
   System requirements: [to be]
   Required packages in python: [to be]
  
Instruction to run NEST:
   
1. Assuming that the spatial dataset is in "data/PDAC_64630/" directory, data preprocessing for input graph generation can be done as follows:  

nohup python -u data_preprocess_NEST.py --dataname='PDAC_64630' > output.log &

It will create two folders in the current working directories: "input_graph/PDAC_64630/" and "metadata/PDAC_64630/" to save the preprocessed input data. Please use the argument --help to see all available input parameters.  

2. To train a NEST model on the preprocessed 'PDAC_64630' data use following command with preferred model name. If the same experiment if repeated multiple times for model ensemble, each time a different run_id should be used and the run_id is expected to be consecutive. For example, if it is run five times then the run_id for the five runs should be 1, 2, 3, 4, and 5 respectively.    

nohup python -u run_NEST.py --dataname='PDAC_64630' --model_name="PDAC_64630_NEST" --run_id=1 > output.log &

3. 

