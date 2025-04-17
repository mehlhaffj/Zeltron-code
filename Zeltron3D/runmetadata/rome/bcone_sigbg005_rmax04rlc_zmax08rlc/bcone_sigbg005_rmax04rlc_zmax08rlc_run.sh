#!/bin/bash

#MSUB -r bcone_sigbg005_rmax04rlc_zmax08rlc_run
#MSUB -n 128            # Number of tasks to use
#MSUB -T 1800           # Walltime in seconds
#MSUB -A gen7669
#MSUB -q rome
#MSUB -m scratch
##MSUB -Q normal
##MSUB -Q long
#MSUB -Q test

ccc_mprun ./zeltron.exe > zeltron.out 2> zeltron.err
