#!/bin/bash
#PBS -l walltime=2:0:00
#PBS -l mem=30gb
#PBS -l nodes=1:ppn=3
#PBS -q batch
#PBS -N 2 
#PBS -o 2.log
#PBS -e 2.err


cd /data/kdi_prod/project_result/1056/18.00/analysis_4/run_simulation_no_umi/
source ~/.bashrc
#eval `modulecmd bash load igor`
export PS1=
# conda activate R
mamba activate /data/users/wsun/anaconda3/envs/R
#echo "Subject: [Curie Bioinfo Cluster] barcode_simulation is done" | sendmail sunwjie@gmail.com

simu_number=2

Rscript bin_simulation_no_umi.R \
    ./tmp/barcode_simulation_${simu_number} \
    ./non_umi_simu_design_matrix.tsv 
