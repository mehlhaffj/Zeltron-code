#!/bin/bash

#MSUB -r test_nozzle_inj_bgvac_icmonopole_run
#MSUB -n 16             # Number of tasks to use
#MSUB -T 1800           # Walltime in seconds
#MSUB -A gen7669
#MSUB -q rome
#MSUB -m scratch
##MSUB -Q normal
##MSUB -Q long
#MSUB -Q test

ccc_mprun ./zeltron.exe > zeltron.out 2> zeltron.err
