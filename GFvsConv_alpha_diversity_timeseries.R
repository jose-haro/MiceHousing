##############################################################
# title: "Alpha diversity in R - qiime2 output"
# author: "ANSC595"
# date: "March 16, 2021"
##############################################################

#The first step is very important. You need to set your working 
#directory. Just like in unix we have to `cd` to where our data is, the 
#same thing is true for R.
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir diversity_results/
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/*.qz* .
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/diversity_results/*.qz* diversity_results/
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/*.tsv .
##############################################

#So where `~/Desktop/ANSC595/moving-pictures` is in the code below, you 
#need to enter the path to where you saved the tutorial or your data.

#setwd('~/Desktop/ANSC595/2022/moving_pictures/from_bell/')
setwd("C:/jose/Molecular_microbiome_exercises/R_exercises/ANSC516-repo/ANSC516/data/moving-pictures")
list.files()

# Modified from the original online version available at 
# http://rpubs.com/dillmcfarlan/R_microbiotaSOP

# and Tutorial: Integrating QIIME2 and R for data visualization 
# and analysis using qiime2R
# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

##Goal
# The goal of this tutorial is to demonstrate basic analyses of microbiota data to determine if and how communities differ by variables of interest. In general, this pipeline can be used for any microbiota data set that has been clustered into operational taxonomic units (OTUs).
#
# This tutorial assumes some basic statistical knowledge. Please consider if your data fit the assumptions of each test (normality? equal sampling? Etc.). If you are not familiar with statistics at this level, we strongly recommend collaborating with someone who is. The incorrect use of statistics is a pervasive and serious problem in the sciences so don't become part of the problem! That said, this is an introductory tutorial and there are many, many further analyses that can be done with microbiota data. Hopefully, this is just the start for your data!

##Data
# The data used here are from the qiime2 moving pictures tutorial. 
# Please see their online tutorial for an explanation of the dataset.

##Files
# We will use the following files created using the qiime2 moving pictures tutorial.

# diversity_results/evenness_vector.qza (alpha diversity)
# diversity_results/faith_pd_vector.qza (alpha diversity)
# diversity_results/observed_features_vector.qza (alpha diversity)
# diversity_results/shannon_vector.qza (alpha diversity)
# metadataGFvsConv_nolabcontrols_wcat.txt (metadata)


# Data manipulation
## Load Packages

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","

# Load data

meta<-read_q2metadata("metadataGFvsConv_nolabcontrols.txt")
str(meta)
colnames(meta)[3] <- "facility"
colnames(meta)[4] <- "week"
colnames(meta)[5] <- "pool"
str(meta)

evenness = read_qza("diversity_results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("diversity_results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("diversity_results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("diversity_results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd2 <- faith_pd [,-1]
colnames(faith_pd2)[1] <- "SampleID"
colnames(faith_pd2)[2] <- "faith_pd"
str(faith_pd2)
faith_pd <- faith_pd2

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(meta)
observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)



###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
meta = meta[,-1]
str(meta)

meta1 <- subset(meta, meta$facility!="GF")
meta1 <- subset(meta1, meta1$facility!="None")

meta2 <- subset(meta, meta$facility!="Conv")
meta2 <- subset(meta2, meta2$facility!="None")


#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or usenon-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(meta2$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta2$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta2$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta2$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta2$shannon_entropy, title = "Shannon")
ggqqplot(meta2$faith_pd, title = "Faith PD")
ggqqplot(meta2$pielou_evenness, title = "Evenness")
ggqqplot(meta2$observed_features, title = "Observed Features")



# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta2$shannon)
shapiro.test(meta2$faith_pd)
shapiro.test(meta2$pielou_evenness)
shapiro.test(meta2$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

# We see that, as expected from the graphs, shannon and evenness 
# are normally distributed.


#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest are

# body site: gut, tongue, right palm, left palm
# subject: 1 and 2
# month-year: 10-2008, 1-2009, 2-2009, 3-2009, 4-2009

## Categorical variables
# Now that we know which tests can be used, let's run them. 

levels(meta2$facility)
#Re-order the groups because the default is alphabetical order
meta2$facility.ord = factor(meta2$facility, c("None", "Conv"))
levels(meta2$facility.ord)


##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since week is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

# Does week impact evenness of the microbiota?

glm.faith.week = glm(faith_pd ~ week, data=meta2)
summary(glm.faith.week)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(faith_pd ~ week, data=meta2)
#Add the glm best fit line
plot(faith_pd ~ week, data=meta2) + abline(glm.faith.week)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.

gaussian.faith.time = glm(faith_pd ~ week, data=meta2, family="gaussian")
plot(gaussian.faith.time, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.time = glm(faith_pd ~ week, data=meta2, family="quasipoisson")
plot(qp.faith.time, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.time)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and week.

#Plot
plot(log(faith_pd) ~ week, data=meta2, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ week, data=meta2, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.time)
