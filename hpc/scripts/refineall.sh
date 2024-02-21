#!/bin/bash

#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --mem-per-cpu=80G
#SBATCH -p medium
#SBATCH --qos short
#SBATCH -t 0:60:0
echo "refineall.sh"
module load julia
export JULIA_PROJECT=.
julia src/refineall.jl
