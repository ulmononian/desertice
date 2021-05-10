#!/bin/bash -l

#SBATCH --qos=debug
##SBATCH --qos=regular
##SBATCH --qos=premium
#SBATCH --constraint=knl
#SBATCH -t 00:30:00
#SBATCH -J mali-tg-wft
#SBATCH -n 68
#SBATCH --tasks-per-node=68
##SBATCH -A m1795
#SBATCH -A m1041

# Script to alternately run MALI and GIA in a data-coupled fashion.
# There are assumptions of a one year coupling interval.
#
# Workflow is like:
#
# 0       1       2
# | MALI  |
# |------>|
# |       |
# |  GIA  |
# |------>|
#         |
#         | MALI  |
#         |------>|
#         |       |
#         |  GIA  |
#         |------>|


# ===================
# Set these locations and vars
MALI=./landice_model
MALI_INPUT=thwaites.4km.cleaned.nc

source /global/project/projectdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
meshvars="latCell,lonCell,xCell,yCell,zCell,indexToCellID,latEdge,lonEdge,xEdge,yEdge,zEdge,indexToEdgeID,latVertex,lonVertex,xVertex,yVertex,zVertex,indexToVertexID,cellsOnEdge,nEdgesOnCell,nEdgesOnEdge,edgesOnCell,edgesOnEdge,weightsOnEdge,dvEdge,dcEdge,angleEdge,areaCell,areaTriangle,cellsOnCell,verticesOnCell,verticesOnEdge,edgesOnVertex,cellsOnVertex,kiteAreasOnVertex"

#  First run MALI
echo "Starting MALI at time:"
date
#srun -n 36 $MALI
#srun -n 68 $MALI
time srun -n 68 --cpu-bind=cores --hint=nomultithread $MALI
echo "Finished MALI at time:"
date
