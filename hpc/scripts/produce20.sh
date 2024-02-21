#!/bin/bash

#SBATCH --ntasks-per-node=50
#SBATCH -n 300
#SBATCH -N 6
#SBATCH --mem=320G
#SBATCH -p medium
#SBATCH -t 48:00:0
echo "produce20.sh"
module load julia
export JULIA_PROJECT=.
julia src/produce20.jl
