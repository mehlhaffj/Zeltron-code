#!/bin/bash

#MSUB -r bcone_sigbg016_rmax07rlc_zmax28rlc_run
#MSUB -n 2048           # Number of tasks to use
#MSUB -T 259200         # Walltime in seconds
#MSUB -A gen7669
#MSUB -q rome
#MSUB -m scratch
##MSUB -Q normal
#MSUB -Q long
##MSUB -Q test

ccc_mprun ./zeltron.exe > zeltron.out 2> zeltron.err
