#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		         ##number of requested CPUs
#SBATCH --time=10:00:00             ##time requested
#SBATCH --job-name qiime_diversity ##name of job
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

#step 3, Alpha and Beta diversity

#Calculate diversity metrics. You will need to designate a sampling depth. You might want to consult the rarifaction output prior to this step.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 30855 \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --output-dir core-metrics-results

#Adding Alpha non phylo
qiime diversity alpha \
  --i-table table.qza \
  --p-metric pielou_e \
  --o-alpha-diversity core-metrics-results/pielou_e_vector.qza

qiime diversity alpha \
  --i-table table.qza \
  --p-metric chao1 \
  --o-alpha-diversity core-metrics-results/chao1_vector.qza


#Statistical testing
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/pielou_e_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/pielou_e_vector-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/chao1_vector.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --o-visualization core-metrics-results/chao1_vector-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization core-metrics-results/unweighted-unifrac-facility-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column facility \
  --o-visualization core-metrics-results/bray_curtis-facility-significance.qzv \
  --p-pairwise

# For Diversity Mantel
qiime metadata distance-matrix \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --m-metadata-column week \
  --o-distance-matrix week_dm \

qiime diversity mantel \
  --i-dm1 week_dm.qza \
  --i-dm2 core-metrics-results/unweighted_unifrac_distance_matrix.qza\
  --p-label1 'Week' \
  --p-label2 'Unweighted_unifrac' \
  --o-visualization core-metrics-results/mantel_unweighted-unifrac-week-significance.qzv

#Make ordination
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --p-custom-axes week \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-week.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --p-custom-axes week \
  --o-visualization core-metrics-results/bray-curtis-emperor-week.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file ./metas/metadataGFvsConv_nolabcontrols.txt \
  --p-custom-axes week \
  --o-visualization core-metrics-results/weighted_unifrac-emperor-week.qzv


echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch pipeline_4_diversity.slurm
