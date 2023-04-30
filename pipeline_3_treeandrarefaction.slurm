#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		         ##number of requested CPUs
#SBATCH --time=10:00:00             ##time requested
#SBATCH --job-name TreeandRarefaction  ##name of job
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

#step 3, Philo tree and rarefaction pre-diversity

echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

echo "Rarefaction"
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 53000 \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization alpha-rarefaction.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch pipeline_3_treeandrarefaction.slurm
