#!/bin/bash

##########################
# A20 particle tracking SLURM MPI job #
##########################

#SBATCH --job-name=A20interpolation

#SBATCH --ntasks-per-node=1
#SBATCH --nodes=4
#SBATCH -A nn9297k
#SBATCH --mem-per-cpu=64G
# run for  minutes
#              d-hh:mm:ss
#SBATCH --time=1-23:00:00

# short partition should do it
#SBATCH --partition=bigmem

# turn on all mail notification
#SBATCH --mail-type=ALL

SCRATCH_DIRECTORY=/cluster/projects/nn9412k/A20/DELTA/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# Activate Anaconda work environment for OpenDrift
source /cluster/home/${USER}/.bashrc
conda activate OpenDrift

# we copy everything we need to the scratch directory
# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
scp ${SLURM_SUBMIT_DIR}/interpolateNORESM_using_ESMF.py ${SCRATCH_DIRECTORY}


# we execute the job and time it
python interpolateNORESM_using_ESMF.py &> a20.output

# after the job is done we copy our output back to $SLURM_SUBMIT_DIR
#cp ${SCRATCH_DIRECTORY}/a20.output ${SLURM_SUBMIT_DIR}

# happy end
exit 0