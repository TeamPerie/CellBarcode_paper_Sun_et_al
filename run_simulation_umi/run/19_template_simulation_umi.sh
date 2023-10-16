#!/bin/bash
#PBS -l walltime=2:0:00
#PBS -l mem=50gb
#PBS -l nodes=1:ppn=3
#PBS -q batch
#PBS -N 19 
#PBS -o 19.log
#PBS -e 19.err


cd /data/kdi_prod/project_result/1056/18.00/analysis_4/run_simulation_umi/
source ~/.bashrc
#eval `modulecmd bash load igor`
export PS1=
# conda activate R
mamba activate /data/users/wsun/anaconda3/envs/R
#echo "Subject: [Curie Bioinfo Cluster] barcode_simulation is done" | sendmail sunwjie@gmail.com

simu_number=19

Rscript bin_simulation_umi.R \
    ./tmp/barcode_simulation_${simu_number} \
    ./umi_simu_design_matrix.tsv 
