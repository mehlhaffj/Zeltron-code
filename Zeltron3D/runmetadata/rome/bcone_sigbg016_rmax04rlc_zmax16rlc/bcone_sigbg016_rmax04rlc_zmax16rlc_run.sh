#!/bin/bash

#MSUB -r bcone_sigbg016_rmax04rlc_zmax16rlc_run
#MSUB -n 512            # Number of tasks to use
#MSUB -T 43200          # Walltime in seconds
#MSUB -A gen7669
#MSUB -q rome
#MSUB -m scratch
#MSUB -Q normal
##MSUB -Q long
##MSUB -Q test

ccc_mprun ./zeltron.exe > zeltron.out 2> zeltron.err
