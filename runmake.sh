#!/bin/sh -l

ncore=$1

#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 00:10:00

#SBATCH -A perc-long
#SBATCH -p cpu
#SBATCH -q normal

project=$2

cd ~/projects/$project

echo accessing...
pwd

make -j $ncore
echo Running Make at $project

