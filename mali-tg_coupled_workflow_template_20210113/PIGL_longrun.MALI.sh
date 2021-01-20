#!/bin/bash -l

##SBATCH --qos=debug
#SBATCH --qos=regular
##SBATCH --qos=premium
#SBATCH --constraint=knl
#SBATCH -t 4:00:00
#SBATCH -J PIGL
#SBATCH -N 5
#SBATCH --tasks-per-node=68
#SBATCH -A m1795

#source /global/project/projectdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh

echo "Starting MALI at time:"
date
#time srun -n 340 --cpu-bind=cores --hint=nomultithread /global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/MPAS-Model/landice_model
time srun -n 340 --cpu-bind=cores --hint=nomultithread /global/cscratch1/sd/hoffman2/THWAITES_1km_2021/landice_model_intel_20210115
echo "Finished MALI at time:"
date


