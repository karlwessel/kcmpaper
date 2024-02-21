#!/bin/bash

#SBATCH --ntasks-per-node=80
#SBATCH -N 1
#SBATCH -n 80
#SBATCH --mem=200G
#SBATCH -p medium
#SBATCH -t 2:00:0
#SBATCH -C scratch
echo "produceall.sh"
module load julia
export JULIA_PROJECT=.
julia src/produceall.jl
