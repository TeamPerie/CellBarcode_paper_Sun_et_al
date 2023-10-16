#!/bin/bash
#PBS -l walltime=20:0:00
#PBS -l mem=70gb
#PBS -l nodes=1:ppn=5
#PBS -q batch
#PBS -N no_umi_simulation 
#PBS -o no_umi_simulation.log
#PBS -e no_umi_simulation.err


cd /data/kdi_prod/project_result/1056/18.00/analysis_4/run_preprocess_simulation/
source ~/.bashrc
#eval `modulecmd bash load igor`
export PS1=
mamba activate /data/users/wsun/anaconda3/envs/R
Rscript run_preprocess_no_umi_simulation.R
#echo "Subject: [Curie Bioinfo Cluster] {{job_name}} is done" | sendmail sunwjie@gmail.com
