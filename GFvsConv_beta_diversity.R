
#Load the packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")


library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#diversity_results/bray_curtis_pcoa_results.qza
#diversity_results/weighted_unifrac_pcoa_results.qza
#diversity_results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#diversity_results/evenness_vector.qza
#diversity_results/faith_pd_vector.qza
#diversity_results/observed_otus_vector.qza
#diversity_results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir diversity_results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/diversity_results/* diversity_results/.
##############################################


###Set your working directory
#setwd("path/to/ANSC516/ANSC-repo/ANSC516/data/moving-pictures")
setwd("C:/jose/Molecular_microbiome_exercises/R_exercises/ANSC516-repo/ANSC516/data/moving-pictures")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("metadataGFvsConv_nolabcontrols_wcat.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("metadataGFvsConv_nolabcontrols_wcat.txt")
str(metadata)
levels(metadata$`facility`)
colnames(metadata)[3] <- "facility"
colnames(metadata)[4] <- "week"
colnames(metadata)[5] <- "pool"
str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("diversity_results/bray_curtis_pcoa_results.qza")
wUF <- read_qza("diversity_results/weighted_unifrac_pcoa_results.qza")
UnUF_PCoA <- read_qza("diversity_results/unweighted_unifrac_pcoa_results.qza")
Jac_PCoA <- read_qza("diversity_results/jaccard_pcoa_results.qza")


facility_colors <- c("Red", "Blue", "Green")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
#alpha controls transparency and helps when points are overlapping
ggplot(bc_meta, aes(x=PC1, y=PC2, color=facility)) +
  theme_q2r() +
  xlab("PC1 (40.10%)") +
  ylab("PC2 (31.43%)") +
  scale_color_manual(values=c("Red", "Blue", "Green"), name = "facility")

# Now we are going to make our code a little more re-usable
my_column <- "facility"
#my_column <- "DietTreatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~week) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=6, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "facility"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape=week)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/BC-ellipse-", my_column,"-week.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


##############################################
## SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("diversity_results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= week), size = 2) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/Wuni-facility", my_column,"-week.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

wUF_dist_mat<-read_qza("diversity_results/weighted_unifrac_distance_matrix.qza")
wUF_dm <- as.matrix(wUF_dist_mat$data) 
rownames(wUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(wUF_dm),metadata$SampleID),]
rownames(wUF_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(wUF_dm ~ facility, data = metadata_sub)

write.table(PERMANOVA_out,"output/Facility_wUF_Adonis_overall.csv",sep=",", row.names = TRUE) 


###################################
#special weighted Unifrac for weeks (didn't result cause few points !!!!)
###################################

metadata2 <- read.delim("metadataGFvsConv_nolabcontrols_onlyconv.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("metadataGFvsConv_nolabcontrols_onlyconv.txt")
str(metadata)
levels(metadata$`week`)
colnames(metadata)[3] <- "facility"
colnames(metadata)[4] <- "week"
colnames(metadata)[5] <- "pool"
str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)



my_column <- "week"
week_colors <- c("Red", "Blue", "Green", "Orange")

Wuni_PCoA<-read_qza("diversity_results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=week_colors, name = "Week")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= facility), size = 2) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Week")
ggsave(paste0("output/Wuni-", my_column,"-byfacility.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

wUF_dist_mat<-read_qza("diversity_results/weighted_unifrac_distance_matrix.qza")
wUF_dm <- as.matrix(wUF_dist_mat$data) 
rownames(wUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(wUF_dm),metadata$SampleID),]
rownames(wUF_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(wUF_dm ~ facility, data = metadata_sub)

write.table(PERMANOVA_out,"output/Facility_wUF_Adonis_overall.csv",sep=",", row.names = TRUE) 



#################################################################################
## SAME thing but with UNWEIGHTED UNIFRAC
##
##

UnUF_PCoA <- read_qza("diversity_results/unweighted_unifrac_pcoa_results.qza")

UnUF_meta <- UnUF_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UnUF_meta,mean)

ggplot(UnUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UnUF_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UnUF_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/UnWuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(UnUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= week), size = 1) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UnUF_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UnUF_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/UnUF-ellipse_", my_column,"-facilitybyweek.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

UnUF_dist_mat<-read_qza("diversity_results/unweighted_unifrac_distance_matrix.qza")
UnUF_dm <- as.matrix(UnUF_dist_mat$data)
rownames(UnUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(UnUF_dm),metadata$SampleID),]
rownames(UnUF_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(UnUF_dm ~ facility, data = metadata_sub)

write.table(PERMANOVA_out,"output/UnUF_Facility_Adonis_overall.csv",sep=",", row.names = TRUE) 

##########################
#same for Jaccard

Jac_PCoA <- read_qza("diversity_results/jaccard_pcoa_results.qza")

Jac_meta <- Jac_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Jac_meta,mean)

my_column <- "facility"

ggplot(Jac_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jac_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jac_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/Jac-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Jac_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= week), size = 1) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jac_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jac_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=facility_colors, name = "Facility")
ggsave(paste0("output/Jac-ellipse_facilitybyweek", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

Jac_dist_mat<-read_qza("diversity_results/jaccard_distance_matrix.qza")
Jac_dm <- as.matrix(Jac_dist_mat$data)
rownames(Jac_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(Jac_dm),metadata$SampleID),]
rownames(Jac_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(Jac_dm ~ facility, data = metadata_sub)

write.table(PERMANOVA_out,"output/Jac_Facility_Adonis_overall.csv",sep=",", row.names = TRUE) 


######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

Facility_Pair <- pairwise.adonis2(bc_dm ~ facility, data = metadata_sub)
write.table(Facility_Pair,"output/Facility_BC_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

Facility_Pair <- pairwise.adonis2(bc_dm ~ facility, data = metadata_sub)
write.table(Facility_Pair,"output/Facility_BC_Adonis_pairwise.csv",sep=",", row.names = TRUE)

Facility_Pair_UnUF <- pairwise.adonis2(UnUF_dm ~ facility, data = metadata_sub)
write.table(Facility_Pair_UnUF,"output/Facility_UnUF_Adonis_pairwise.csv",sep=",", row.names = TRUE)
