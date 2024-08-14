### Pulling the nest_image.sif

We can pull the singularity image as follows:

```
mkdir nest_container 
cd nest_container
singularity pull nest_image.sif library://fatema/collection/nest_image.sif:latest
```
[** NOTE: Instructions for pulling the image on Digital Alliance are provided at the end.](https://github.com/schwartzlab-methods/NEST/blob/main/vignette/running_NEST_singularity_container.md#running-nest-container-on-digital-alliance-using-apptainer)

### Running NEST model through the downloaded singularity image

For the demonstration purpose, let us assume that the paths of NEST repository, NEST image, and data are as follows:
1. NEST repository: /cluster/projects/prof-group/fatema/NEST
2. NEST image: /cluster/projects/prof-group/fatema/nest_container/nest_image.sif
3. Data: /cluster/projects/prof-group/fatema/NEST/data/

I will be referring to these while running the commands. First, we navigate to the NEST repository:
```
cd /cluster/projects/prof-group/fatema/NEST
```
Then we run NEST container. We can use "singularity run" to directly execute the NEST commands, or 
open the container shell/terminal using "singularity shell" command and interactively run the NEST commands there. 
Here, I am opening the shell as follows:

```
singularity shell --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif
```
Then, on the shell, I execute the usual NEST preprocessing command:
```
Singularity> bash nest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='/cluster/projects/prof-group/fatema/NEST/data/V1_Human_Lymph_Node_spatial/'
```


Once it is done, you can find the preprocessed data in the /metadata and input_graph/ directory. 
```
Singularity> ls metadata/V1_Human_Lymph_Node_spatial/
Singularity> ls input_graph/V1_Human_Lymph_Node_spatial/
```



Now, you can start the model training on the shell. But you have to ensure that it has the access to the GPU. 
That is why, we have to exit from the current shell and restart the shell using --nv command as follows:
```
$ singularity shell --nv --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif
```
Here, --nv option is used for GPU accessibility. You can see the allocated GPU in the shell as below. 
```
Singularity> nvidia-smi
```

Then you can start training as below: 
```
Singularity> bash nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 500 --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1
```


However, running the model training like this is never feasible due to the long time required (couple of hours to couple of days based on the size of the input data). 
If you get disconnected from the servers for some reason, your training will also stop. 
That is why we have two options:

1. We can run the job in the background using “singularity run” command as follows:
```
nohup singularity run --nv --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif bash nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 500 --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_NEST_V1_Human_Lymph_Node_spatial_run1.log &
```
2. If we are using high performance computing (HPC) service with Slurm for cluster management and job scheduling, then we have to submit the job to the cluster as follows (the parameters might vary depending on the system):
``` 
sbatch -A prof-group_gpu -p gpu --gres=gpu:1 --constraint gpu32g gpu_job_container_NEST_human_lymph_node_run1.sh
```

The contents of gpu_job_container_NEST_human_lymph_node_run1.sh:
```
#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a job on a cluster.
# ---------------------------------------------------------------------
#SBATCH -c 16
#SBATCH --mem=30GB
#SBATCH --time=72:00:00
#SBATCH --job-name=NEST_lymph
#SBATCH --output=some_name-%j.out
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------


# module load
module load singularity/3.11.0

# navigate to the directory that has NEST repository
cd /cluster/projects/prof-group/fatema/NEST/
echo "Current working directory: `pwd`"

# nest_image.sif is here: /cluster/projects/prof-group/fatema/nest_container/nest_image.sif

# checking the available gpu memory before starting the training
nvidia-smi

# run your python script with parameters using singularity
echo "lymph. Going to start process: run 1"

singularity run --nv --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif bash nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 60000 --manual_seed='yes' --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1


echo "Job finished with exit code $? at: `date`"
```

After the training is finished, we run postprocessing. Usually, we run the training 5 times with five DIFFERENT seeds and --run_id starting from 1 to 5, and then postprocess those 5 runs. For running postprocessing, we can open the singularity shell as before, or just use the “singularity run” command as below: 
```
$ singularity run --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif bash nest postprocess --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --total_runs=5
```

After the postprocessing is done, we can find the list of top 20% communications (strength or attention score is over 80th percentile) inside the output/ directory:
```
$ ls output/V1_Human_Lymph_Node_spatial/
NEST_V1_Human_Lymph_Node_spatial_top20percent.csv
```

Then we run visualization script for this list as below:
```
$ singularity run --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif bash nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000
```
We will find the visualization files in the same output/ directory.
```
$ ls output/V1_Human_Lymph_Node_spatial/
```


Then we can convert the *.dot file to pdf and svg as follows:
```
$ singularity run --home=/cluster/projects/prof-group/fatema/nest_container/ /cluster/projects/prof-group/fatema/nest_container/nest_image.sif bash nest output_graph_picture output/V1_Human_Lymph_Node_spatial/NEST_V1_Human_Lymph_Node_spatial_test_interactive.dot 
```
It will generate edge_graph.pdf and edge_graph.svg in the current working directory.


### Running NEST container on Digital Alliance using Apptainer
#### Pulling the Singularity image on Digital Alliance using Apptainer
```
module load apptainer
apptainer pull nest_image.sif library://fatema/collection/nest_image.sif:latest
```
The command might show the following error:

FATAL:   Unable to get library client configuration: remote has no library client (see https://apptainer.org/docs/user/latest/endpoint.html#no-default-remote)

In that case, execute the following commands to resolve the error (details are provided on the official site: https://apptainer.org/docs/user/1.0/endpoint.html#restoring-pre-apptainer-library-behavior):
```
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud
apptainer remote list
```

After that we are able to pull the image using the following command:
```
apptainer pull nest_image.sif library://fatema/collection/nest_image.sif:latest
```
All the Singularity commands mentioned above also work with Apptainer if the 'singularity' term is replaced with 'apptainer'. Few things to note:
1. We may have to move the NEST repository and data to the NEST container so that those are visible to the container while running. 
```
cd /project/def-prof-group/fatema/
mv /project/def-prof-group/fatema/NEST /cluster/projects/prof-group/fatema/nest_container/
```
Therefore, the paths of NEST repository, NEST image, and data are as follows:

NEST repository: /project/def-prof-group/fatema/nest_container/NEST

NEST image: /project/def-prof-group/fatema/nest_container/nest_image.sif

Data: /project/def-prof-group/fatema/nest_container/NEST/data/
   
2. Additionally, we may need to set the home directory to the container path as follows:
```
apptainer shell --home=/project/def-prof-group/fatema/nest_container/ /project/def-prof-group/fatema/nest_container/nest_image.sif
```
It will open the Singularity shell as the image is a Singularity image. 
If we don't use --home, it may incorrectly look at the system paths for the Python packages. 

#### NEST Preprocessing:
```
cd /project/def-prof-group/fatema/nest_container/NEST

apptainer run --home=/project/def-prof-group/fatema/nest_container/ /project/def-prof-group/fatema/nest_container/nest_image.sif bash nest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/'
```

#### Running NEST model: 
We write a shell script `gpu_job_digital_alliance_container.sh’ as follows:
```
#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a GPU job on Graham at Digital Alliance cluster.
# ---------------------------------------------------------------------
#SBATCH --account=def-prof-group
#SBATCH --gres=gpu:v100:1    	# Request 2 v100 GPUs
#SBATCH --constraint=cascade,v100 # They will be cascaded to offer 32 GB
#SBATCH --cpus-per-task=16  # CPU Cores
#SBATCH --mem=32000M    	# CPU Memory
#SBATCH --time=72:00:00
#SBATCH --job-name=human_lymph_r1
#SBATCH --output=human_lymph_r1.out
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------

echo "available gpu:"
nvidia-smi

echo "Current working directory: `pwd`"
# change current working directory to NEST
cd /project/def-prof-group/fatema/nest_container/NEST/
echo "Current working directory: `pwd`"

# load necessary modules
module load apptainer

# run your python script with parameters using singularity
echo "lymph. Going to start process: run 1"

apptainer run --nv --home=/project/def-prof-group/fatema/nest_container/ /project/def-prof-group/fatema/nest_container/nest_image.sif bash nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 60000  --manual_seed='yes' --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1 

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"

```
Please change the parameters and paths to the NEST repository as required. We may save this script here: /project/prof-group/fatema/gpu_job_digital_alliance_container.sh

Then navigate to that directory and submit the job as follows:
```
cd /project/def-prof-group/fatema/
sbatch gpu_job_digital_alliance_container.sh
```
Job status can be viewed as follows:
```
squeue -u fatema
```
This is for one run. Usually, we run the training 5 times with five DIFFERENT seeds and --run_id starting from 1 to 5, and then postprocess those 5 runs. And then proceed for postprocessing.

#### NEST postprocessing:
```
apptainer run --home=/project/def-prof-group/fatema/nest_container/ /project/def-prof-group/fatema/nest_container/nest_image.sif bash nest postprocess --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --total_runs=5
```

#### NEST visualization:
```
apptainer run --home=/project/def-prof-group/fatema/nest_container/ /project/def-prof-group/fatema/nest_container/nest_image.sif bash nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000
```

