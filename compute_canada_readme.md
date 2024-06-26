# Python package installation

Run following commands to create virtual environment and activate it before installing python libraries:
```
module load python/3.10
virtualenv --no-download /project/[group_name]/[user_name]/nest_venv
source /project/[group_name]/[user_name]/nest_venv/bin/activate
```
Please replace [group_name] and [user_name] with your assigned group and name.

Upgrade pip3:
```
pip3 install --no-index --upgrade pip
```

Install required libraries as follows:  
```
pip3 install ipython
pip3 install pandas
pip3 install altair vega_datasets
pip3 install flake8 pytest jinja2 sphinx m2r docutils
pip3 install collection
pip3 install stlearn
pip3 install scanpy
pip3 install qnorm
pip3 install csv pickle gzip matplotlib scipy sklearn 
```

Next, see which pytorch and CUDA versions are available:
```
avail_wheels "torch*"
module load cuda [press tab to see the available versions]
```

Next, based on the Pytorch and CUDA version, run with wheels links for installing. Here we are using pytorch 1.13.1 with CUDA 11.7:

```
pip3 install torch-scatter -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-sparse -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-cluster -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-geometric -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
```

# Submit GPU jobs

There is a sample shell script written for submitting GPU job as follows:
```
sbatch gpu_job_submit_compute_canada.sh
```
Please see the 'gpu_job_submit_compute_canada.sh' in the NEST repository to understand the requested GPU resources. This will execute just one run of NEST. If enough GPU is available, multiple NEST runs can be executed in parallel by replacing the contents of line 35 with the following:
```
nohup nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=1 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=2 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=3 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=4 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
nohup nest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=5 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &
```




