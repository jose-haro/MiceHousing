#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		         ##number of requested CPUs
#SBATCH --time=10:00:00             ##time requested
#SBATCH --job-name qiime_tax_ancom ##name of job
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

echo "Make a stacked barplot"
#qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization taxa-bar-plots.qzv

echo "Filter away some sample types and do ancom"

echo "Comparing GF vs Conv only in week 7"

#qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --p-where "[week]='7'" \
  --o-filtered-table week7-table.qza
  
#qiime composition add-pseudocount \
  --i-table week7-table.qza \
  --o-composition-table comp-week7-table.qza
  
#qiime composition ancom \
  --i-table comp-week7-table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization ancom-facility-w7.qzv

echo "Comparing GF vs Conv only in week 5"

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --p-where "[week]='5'" \
  --o-filtered-table week5-table.qza

qiime composition add-pseudocount \
  --i-table week5-table.qza \
  --o-composition-table comp-week5-table.qza

qiime composition ancom \
  --i-table comp-week5-table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization ancom-facility-w5.qzv

echo "Comparing GF vs Conv only in week 3"

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --p-where "[week]='3'" \
  --o-filtered-table week3-table.qza

qiime composition add-pseudocount \
  --i-table week3-table.qza \
  --o-composition-table comp-week3-table.qza

qiime composition ancom \
  --i-table comp-week3-table.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization ancom-facility-w3.qzv


echo "Collapse levels and ancom again."

echo "For Facilities in week 7"

#qiime taxa collapse \
  --i-table week7-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table week7-table-l6.qza

#qiime composition add-pseudocount \
  --i-table week7-table-l6.qza \
  --o-composition-table comp-week7-table-l6.qza

#qiime composition ancom \
  --i-table comp-week7-table-l6.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization l6-ancom-facility-w7.qzv

echo "For Facilities in week 5"

qiime taxa collapse \
  --i-table week5-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table week5-table-l6.qza

qiime composition add-pseudocount \
  --i-table week5-table-l6.qza \
  --o-composition-table comp-week5-table-l6.qza

qiime composition ancom \
  --i-table comp-week5-table-l6.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization l6-ancom-facility-w5.qzv

echo "For Facilities in week 3"

qiime taxa collapse \
  --i-table week3-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table week3-table-l6.qza

qiime composition add-pseudocount \
  --i-table week3-table-l6.qza \
  --o-composition-table comp-week3-table-l6.qza

qiime composition ancom \
  --i-table comp-week3-table-l6.qza \
  --m-metadata-file metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization l6-ancom-facility-w3.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch pipeline_5_taxonomic.slurm
