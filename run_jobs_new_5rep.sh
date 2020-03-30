#!/bin/bash
#
# ALICE JOB SCRIPT
# ARRAY JOB FOR MULTIPLE N-body RUNS
# THOMAS COX - FEB 2020
#
#PBS -N 3-body
#PBS -j oe
#PBS -m bea
#PBS -M tbc7@student.le.ac.uk
#PBS -l pvmem=100mb
#PBS -l walltime=02:00:00
#PBS -t 170-240
#PBS -l nodes=16:ppn=1
cd /lustre/alice3/scratch/spectre/t/tbc7
ap=`echo $PBS_ARRAYID / 100.0 | bc -l`
./restricted_3body_rep_1.e -q 1.0 -e 0.1 -r $ap

./restricted_3body_rep_2.e -q 1.0 -e 0.1 -r $ap

./restricted_3body_rep_3.e -q 1.0 -e 0.1 -r $ap

./restricted_3body_rep_4.e -q 1.0 -e 0.1 -r $ap

./restricted_3body_rep_5.e -q 1.0 -e 0.1 -r $ap


