#!/bin/bash

#SBATCH --account=apam # Account name is needed
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 24
#SBATCH --mem-per-cpu 4G
#SBATCH --time 0-8:00:00
#SBATCH --job-name pdfmix

echo "Python Interpreter:" $(which python)
python gen_mix.py
