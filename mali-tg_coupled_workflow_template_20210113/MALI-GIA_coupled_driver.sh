#!/bin/bash -l

##SBATCH --qos=debug
#SBATCH --qos=regular
##SBATCH --qos=premium
#SBATCH --constraint=knl
#SBATCH -t 1:00:00
#SBATCH -J TG1km
#SBATCH -N 5
#SBATCH --tasks-per-node=68
#SBATCH -A m1795

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
GIAPATH=.
MALI=./landice_model_intel_20210115
MALI_INPUT=Thwaites_1to8km_r02_20210119.nc
MALI_OUTPUT=output-cpl.nc
MALI_NL=namelist.landice
MALI_STREAMS=streams.landice
START_YEAR=2015
RUN_DURATION=300
CPL_DT=1.0

RESTART_SCRIPT=1 # should be 0 or 1
# ==================

# Other things you could change
GIAGRID=gia_grid.nc
MALILOAD=mali_load.nc
GIAOUTPUT=uplift_GIA.nc


export HDF5_USE_FILE_LOCKING=FALSE
#source /users/mhoffman/setup_badger_mods.20181206.sh
source /global/project/projectdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
meshvars="latCell,lonCell,xCell,yCell,zCell,indexToCellID,latEdge,lonEdge,xEdge,yEdge,zEdge,indexToEdgeID,latVertex,lonVertex,xVertex,yVertex,zVertex,indexToVertexID,cellsOnEdge,nEdgesOnCell,nEdgesOnEdge,edgesOnCell,edgesOnEdge,weightsOnEdge,dvEdge,dcEdge,angleEdge,areaCell,areaTriangle,cellsOnCell,verticesOnCell,verticesOnEdge,edgesOnVertex,cellsOnVertex,kiteAreasOnVertex"

END_ITER=`python -c "import math; end_iter=int(math.ceil($RUN_DURATION / $CPL_DT)); print(end_iter)"`
echo END_ITER=$END_ITER

# Match output interval in MALI streams to END_ITER & ensure format matches that needed by MALI
output_int=`python -c "output_int=int($CPL_DT); print('{0:04d}'.format(output_int))"`
sed -i.SEDBACKUP -e "/output-cpl.nc/,/<\/streams>/ s/output_interval.*/output_interval=\"$output_int-00-00_00:00:00\">/" $MALI_STREAMS
echo output_int=$output_int

# Match run duration in MALI namelist to coupling interval & ensure format matches that needed by MALI
cpl_dt_formatted=`python -c "cpl_dt_int=int($CPL_DT); print('{0:04d}'.format(cpl_dt_int))"`
sed -i.SEDBACKUP "s/config_run_duration.*/config_run_duration = '$cpl_dt_formatted-00-00_00:00:00'/" $MALI_NL
echo cpl_dt=$cpl_dt_formatted


if [ $RESTART_SCRIPT -eq 1 ]; then
  start_ind=`cat coupler_restart.txt`
  #python -c "s=int($CPL_DT*$start_ind); print('{0:04d}-01-01_00:00:00'.format(s))" > restart_timestamp
  python -c "s=int($CPL_DT*$start_ind)+${START_YEAR}; print('{0:04d}-01-01_00:00:00'.format(s))" > restart_timestamp
  echo "new restart_timestamp value: " `cat restart_timestamp`
else
  start_ind=0
fi
echo start_ind=$start_ind

# if [ $RESTART_SCRIPT -eq 1 ]; then
#   # REDO THE FINAL ITERATION
#   RSTTIME=`head -c 20 restart_timestamp | tail -c 19`
#   RSTYR=`echo $RSTTIME|head -c 4`
#   echo RSTYR=$RSTYR
#   startyear=`python -c "newyr=int('$RSTYR'); print '{0:04d}'.format(newyr-$CPL_DT)"`
#   echo startyear=$startyear
#   echo " $startyear-01-01_00:00:00" > restart_timestamp
# else
#   startyear=0
# fi

mkdir iteration_archive

# for i in $(seq $startyear $END_ITER); do
for i in $(seq $start_ind $END_ITER); do

   echo ""; echo ""
   echo "=================================================="
   echo "Starting iteration $i"; echo ""; echo ""

   # Check if initial run or restart
   if [ $i -eq 0 ] && [ $RESTART_SCRIPT -eq 0 ]; then
      echo "This is the first iteration of a new simulation: Preparing new run."

      # Set up GIA mesh
      echo "Setting up GIA mesh"
      $GIAPATH/create_GIA_domain.py -m $MALI_INPUT -g $GIAGRID
      rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

      # Set restart flag to false, to be safe
      echo "Setting restart flags in namelist"
      sed -i.SEDBACKUP "s/config_do_restart.*/config_do_restart = .false./" $MALI_NL
      #sed -i.SEDBACKUP "s/config_start_time.*/config_start_time = '0000-01-01_00:00:00'/" $MALI_NL
      #sed -i.SEDBACKUP "s/config_start_time.*/config_start_time = '2015-01-01_00:00:00'/" $MALI_NL
      sed -i.SEDBACKUP "s/config_start_time.*/config_start_time = '${START_YEAR}-01-01_00:00:00'/" $MALI_NL
      
      # Set GIA model args for a cold run
      GIAARGS="-d $CPL_DT"
   else # this is a restart
      echo "This iteration is a restart for MALI.  Preparing MALI restart run."
      # Set restart flag to restart (will be done every time, but that's ok)
      sed -i.SEDBACKUP "s/config_do_restart.*/config_do_restart = .true./" $MALI_NL
      sed -i.SEDBACKUP "s/config_start_time.*/config_start_time = 'file'/" $MALI_NL

      # Set GIA model args for a restart
      GIAARGS="-r -d $CPL_DT"
   fi

   mkdir iteration_archive/iter_$i

   # First run MALI
   echo "Starting MALI at time:"
   date
   #srun -n 36 $MALI
   #srun -n 68 $MALI
   #time srun -n 68 --cpu-bind=cores --hint=nomultithread $MALI
   time srun -n 340 --cpu-bind=cores --hint=nomultithread $MALI    
   echo "Finished MALI at time:"
   date

   mv log.landice* iteration_archive/iter_${i}
   mv log.albany* iteration_archive/iter_${i}

   # interpolate ice load to GIA grid
   # copy second to last time of the needed fields.  This represents the prior year.
   echo "Copying thickness and bed topography fields from $MALI_OUTPUT into $MALILOAD."
   time ncks -A -d Time,-2 -v thickness,bedTopography,$meshvars $MALI_OUTPUT $MALILOAD
   echo "Finished copying thickness and bed topography fields into $MALILOAD."
   rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
   cp $MALILOAD iteration_archive/iter_${i}
   
   echo "Interpolating $MALILOAD onto the GIA grid."
   time $GIAPATH/interp_MALI-GIA.py -d g -m $MALILOAD -g $GIAGRID
   echo "Finished interpolating $MALILOAD onto the GIA grid."
   rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
   cp iceload.nc iteration_archive/iter_${i}

   # Run GIA model
   echo "Starting GIA model"
   time $GIAPATH/mali-gia-driver.py $GIAARGS
   echo "Finished GIA model"
   rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
   cp $GIAOUTPUT iteration_archive/iter_${i}

   # interpolate bed topo to MALI grid
   echo "Interpolating GIA uplift field onto the MALI grid."
   time $GIAPATH/interp_MALI-GIA.py -d m -m $MALI_INPUT -g $GIAOUTPUT
   echo "Finished interpolating GIA uplift field onto MALI grid." 
   rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
   cp bedtopo_update_mpas.nc iteration_archive/iter_${i}

   # Stick new bed topo into restart file
   # (could also input it as a forcing file... not sure which is better)
   #RSTTIME=`head -c 20 restart_timestamp | tail -c 19 | tr : .`
   RSTTIME=`head -c 20 restart_timestamp | tail -c 20 | tr : .`
   RSTFILE=restart.$RSTTIME.nc
   echo restart time=$RSTTIME
   echo restart filename=$RSTFILE

   if [ $RESTART_SCRIPT -eq 1 ]; then
#     let jj=${i}+1
#     cp $RSTFILE $RSTFILE.iter${jj}.bak  # back up first (maybe remove later)
     cp $RSTFILE iteration_archive/iter_${i}
   else
     cp $RSTFILE iteration_archive/iter_${i}
   fi

   # ncks -A -v bedTopography bedtopo_update_mpas.nc $RSTFILE
   # rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
   echo "Copying upliftRate field to the restart file."
   time ncks -A -v upliftRate bedtopo_update_mpas.nc $RSTFILE
   echo "Finished copying upliftRate field to the restart file."
   rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi



   echo "Finished iteration $i"

   mid_iter=$i
   echo $mid_iter>coupler_restart.txt

done;
