#!/bin/bash
#SBATCH -J simrun_v1        
#SBATCH -c 32                   
#SBATCH --mem=8G                
#SBATCH -p small_cpu
#SBATCH --tmp=5G                             
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de

srun R --save < 2026_03_07_sim_v2.R