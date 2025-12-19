#!/bin/bash
#SBATCH -J simucalculable        
#SBATCH -c 32                   
#SBATCH --mem=8G                
#SBATCH -p small_cpu
#SBATCH --tmp=5G                             
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de

srun R --save < Simulationcalculable.R