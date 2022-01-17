#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --job-name=Lotus			# Job name
#SBATCH --time=60:00:00		 		# Walltime
#SBATCH --partition=batch

module load openmpi/3.0.0/gcc-6.4.0
module load gcc/6.4.0
module load python/2.7.14
module load boost/1.76.0

source ~/.profile

./Allrun -parallel 1
