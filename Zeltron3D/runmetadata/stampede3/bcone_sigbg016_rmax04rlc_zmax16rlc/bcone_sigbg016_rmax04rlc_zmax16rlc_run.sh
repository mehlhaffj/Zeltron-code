#!/bin/bash

# Choose number of nodes to run on
##SBATCH -A TG-AST030031N

#SBATCH -A TG-PHY140041      # Allocation name (recon)
##SBATCH -A phy160032        # Allocation name (turb)
#SBATCH --nodes=16            # Total # of nodes
#SBATCH --ntasks-per-node=32  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=1     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=12:00:00       # Total run time limit (hh:mm:ss)
#SBATCH -J zeltron            # Job name
#SBATCH -o zeltron.out        # Name of stdout output file
#SBATCH -e zeltron.err        # Name of stderr error file
##SBATCH -p debug              # Queue (partition) name
#SBATCH -p skx               # Queue (partition) name
#SBATCH --mail-user=john.mehlhaff@univ-grenoble-alpes.fr
#SBATCH --mail-type=all      # Send email at begin and end of job

##########
#
# Variables
#
##########
RUNNAME=bcone_sigbg016_rmax04rlc_zmax16rlc
contact=john.mehlhaff@univ-grenoble-alpes.fr

OUTFILE="zeltron.out"
ERRFILE="zeltron.err"

# Launch MPI code
mpirun -np $SLURM_NTASKS ./zeltron.exe

# #MSUB -r bcone_sigbg016_rmax04rlc_zmax16rlc_run
# #MSUB -n 512            # Number of tasks to use
# #MSUB -T 43200          # Walltime in seconds
# #MSUB -A gen7669
# #MSUB -q rome
# #MSUB -m scratch
# #MSUB -Q normal
# ##MSUB -Q long
# ##MSUB -Q test
# 
# ccc_mprun ./zeltron.exe > zeltron.out 2> zeltron.err
