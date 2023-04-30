#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		         ##number of requested CPUs
#SBATCH --time=10:00:00             ##time requested
#SBATCH --job-name qiime_import ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd GFvsConv_Data_NoLabControls

pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"
module load Qiime/2-2022.8

#step 3, Run the code you want

#assuming metadata and manifest files are already in scratch

echo "Data Import"
#Make the sequences into a qiime recognizable format

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ./metas/manifestGFvsConv_nolabcontrols.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path demux.qza

ls 

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

#Check demux.qzv for decision on Trim and Trunc

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch pipeline_1_import.slurm
