#!/bin/bash                                                                                                                     
#SBATCH -p core                                                                                                                 
#SBATCH -t 3-00:00:00                                                                                                           
#SBATCH -A snic0000-000-00                                                                                                      
#SBATCH --array=0-500                                                                                                         
#SBATCH -e DIR_error/%A_%a.error                                                                                                

python msp_sims_NULL.py $SLURM_ARRAY_TASK_ID &
python msp_sims_ALT1.py $SLURM_ARRAY_TASK_ID &
python msp_sims_ALT2.py $SLURM_ARRAY_TASK_ID &
python msp_sims_ALT3.py $SLURM_ARRAY_TASK_ID &
python msp_sims_ALT4.py $SLURM_ARRAY_TASK_ID &

wait
exit 0

########################################         
 
# This SLURM array submission script iterates through 500 jobs (replicate simulations) per model. It implements the necessary simulation script for each model, e.g. msp_sims_ALT1.py.
# The resulting 500 tree sequence files per model are compressed and outputted to "tree_seq_files/"
# WARNING: careful as even though output is compressed, this will require A LOT of memory.
# Submit jobs as: "sbatch submit_msp_sims_arrayjob.sh".
 

