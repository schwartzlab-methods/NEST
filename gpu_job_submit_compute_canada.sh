#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a gpu job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --account=dummy_group
#SBATCH --gres=gpu:1        # Request GPU "generic resources"
#SBATCH --cpus-per-task=16  # Cores proportional to GPUs: 6 on Cedar, 16 on Graham.
#SBATCH --mem=63500M        # Memory proportional to GPUs: 31500 Cedar, 63500 Graham.
#SBATCH --time=24:00:00
#SBATCH --job-name=test_job
#SBATCH --output=V1_human_lymph_log-%j.out
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------
# activate your virtual environment
source /home/[user_name]/ENV/bin/activate

# load necessary modules
module load python/3.10

echo 'Available gpu: '
nvidia-smi

# Change current working directory to the downloaded NEST repository
cd /project/[group_name]/[user_name]/NEST

# run your python script with parameters
echo 'Running NEST'
nohup nest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=1 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &

echo 'To track the training progress see the output_human_lymph_node_run1.log under the working directory'

echo 'Going to sleep for 24 hours. If the process is finished before that, please kill this job to free up GPU resources. '
sleep 86400
# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
