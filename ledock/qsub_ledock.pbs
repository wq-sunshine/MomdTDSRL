#!/bin/bash
#PBS -N ledock-mpi-pbs-$1
#PBS -l nodes=1:ppn=28
#PBS -l walltime=1:00:00
#PBS -q q_csywz
#PBS -V
#PBS -S /bin/bash
#PBS -o /home/csywz/deep/ledock/mpi_test_log_$2.log

date
mpiexec -n 28 python /home/csywz/deep/Modules/mpi_ledock.py $3
date