#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --nodelist=node01
#SBATCH -o log/sbatch.calTKE.log.o
#SBATCH -e log/sbatch.calTKE.log.e

python calTKE.py > calTKE.out
echo "DONE"
