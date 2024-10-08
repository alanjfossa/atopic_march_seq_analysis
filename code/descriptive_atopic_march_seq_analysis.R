#Purpose: The purpose of this project is to attempt to vizualize and characterize "atopic march" in HOME study children. 

#Author: Alan J. Fossa

#Load required packages----
library(tidyverse)
library(conflicted)
library(cluster) 
library(WeightedCluster)
library(devtools)
library(TraMineR) # Sequence analysis workhorse
library(TraMineRextras)
library(NbClust)
library(ggseqplot)
library(gridExtra)
library(Tmisc)

#Set working dir----
setwd("C:\\Users\\afossa\\OneDrive\\repos\\atopic_march_seq_analysis")

#Load in asthma/allergy outcomes data----
atopy<-read_csv("C:\\Users\\afossa\\Brown Dropbox\\Alan Fossa\\Data Re-Organization 2024\\Clean Data\\ALL_ECZ_WHZ.csv") 

#Define missingness threshold----
miss_thresh<-0.25

##Wheeze (any)----
whz<-atopy %>% 
  select(participant_id,contains("whz_any")) %>% 
  mutate(
    across(contains("whz_any"),factor)
  ) %>% 
  select(sort(names(.))) #Column names need to be ordered by time.

whz$pmiss<-whz %>% 
  select(contains("whz_any")) %>% 
  apply(MARGIN=1,FUN=function(x){sum(is.na(x))/ncol(.)})

whz<-whz %>% 
  dplyr::filter(pmiss<miss_thresh) %>% 
  select(-pmiss)

##Eczema----
derm<-atopy %>% 
  select(participant_id,contains("eczema")) %>% 
  mutate(
    across(contains("eczema"),factor)
  ) %>% 
  select(sort(names(.)))

derm$pmiss<-derm %>% 
  select(contains("eczema")) %>% 
  apply(MARGIN=1,FUN=function(x){sum(is.na(x))/ncol(.)})

derm<-derm %>% 
  dplyr::filter(pmiss<miss_thresh) %>% 
  select(-pmiss)

##Rhinitis----
allergy<-atopy %>% 
  select(participant_id,contains("allergy"),-food_allergy_P4) %>% 
  mutate(
    across(contains("allergy"),factor)
  ) %>% 
  select(sort(names(.)))

allergy$pmiss<-allergy %>% 
  select(contains("allergy")) %>% 
  apply(MARGIN=1,FUN=function(x){sum(is.na(x))/ncol(.)})

allergy<-allergy %>% 
  dplyr::filter(pmiss<miss_thresh)

#Sequence index plots----

##Wheezing (any)----
###Set state alphabet,labels, and colors----
whz_alphabet <- c("0","1","99")

whz_labels <- c("No wheeze","Wheeze","Missing")

no_wheeze<-"lightblue"
wheeze<-"salmon"
missing<-"beige"

###Define sequences----
whz_seq <- seqdef(
  whz, # Select data   
  var = 2:15, # Columns containing the sequences
  alphabet = whz_alphabet,
  labels = whz_labels,
  start=6, 
  missing=NA,
  left="99",
  gaps="99",
  right="99",
  void="%", 
  cpal=c(no_wheeze,wheeze,missing)
  )

View(whz_seq)

###Plot----
seqIplot(whz_seq, 
         sortv = "from.start",   # Sequence object
         with.legend = "right", # Display legend on right side of plot
         cex.legend = 0.6,  # Change size of legend
         main = "Wheezing Sequence Index Plot", # Plot title
         xlab = "Age (months)"
         ) 

##Eczema----
###Set state alphabet,labels, and colors----
derm_alphabet <- c("0","1","99")

derm_labels <- c("No eczema","Eczema","Missing")

no_eczema<-"lightblue"
eczema<-"salmon"

###Define sequences----
derm_seq <- seqdef(
  derm, # Select data   
  var = 1:14, # Columns containing the sequences
  alphabet = derm_alphabet,
  labels = derm_labels,
  start=6, 
  missing=NA,
  left="99",
  gaps="99",
  right="99",
  void="%", 
  cpal=c(no_eczema,eczema,missing)
)

###Plot----
seqIplot(derm_seq, 
         sortv = "from.start",   # Sequence object
         with.legend = "right", # Display legend on right side of plot
         cex.legend = 0.6,  # Change size of legend
         main = "Eczema Sequence Index Plot", # Plot title
         xlab = "Age (months)"
) 

##Rhinitis----
###Set state alphabet,labels, and colors----
allergy_alphabet <- c("0","1","99")

allergy_labels <- c("No rhinitis","Rhinitis","Missing")

no_rhinitis<-"lightblue"
rhinitis<-"salmon"

###Define sequences----
allergy_seq <- seqdef(
  allergy, # Select data   
  var = 1:14, # Columns containing the sequences
  alphabet = allergy_alphabet,
  labels = allergy_labels,
  start=6, 
  missing=NA,
  left="99",
  gaps="99",
  right="99",
  void="%", 
  cpal=c(no_rhinitis,rhinitis,missing)
)

###Plot----
seqIplot(allergy_seq, 
         sortv = "from.start",   # Sequence object
         with.legend = "right", # Display legend on right side of plot
         cex.legend = 0.6,  # Change size of legend
         main = "Rhinitis Sequence Index Plot", # Plot title
         xlab = "Age (months)"
) 

#Sequence analysis

##Eczema----
###Cost matrix----
derm_cost_matrix <- seqsubm(derm_seq,  # Sequence object
                      method = "TRATE",  # Method to determine costs
                      time.varying = FALSE) # Does not allow the cost to vary over time)
derm_cost_matrix

###Get costs----
####Optimal matching method----
derm_dist <- seqdist(derm_seq,
                   method = "HAM",
                   sm = derm_cost_matrix)

###Get clusters----
####PAM method----
derm_pamRange <- wcKMedRange(derm_dist, kvals=2:6) # this takes a while to run
summary(derm_pamRange, max.rank=6)

#Plot metrics
plot(derm_pamRange, stat = c("ASW","HC"), norm="zscore", lwd = 2, cex=2, col = c('#6666ff', '#cc0000'), legendpos = "topright", main = "HAM PAM Quality")
abline(v=7, col="#666666", lty="longdash", lwd = 2)

####Choose number of clusters----
derm_clust4 <- wcKMedoids(derm_dist, k=4, cluster.only=TRUE)
derm$clust4 <- wcKMedoids(derm_dist, k=4, cluster.only=TRUE)

###Visualize clusters----

####Modal plot----
derm_modal <- seqmsplot(derm_seq, group=derm_clust4, border=NA) #index plots by cluster

##Wheeze----
###Cost matrix----
whz_cost_matrix <- seqsubm(whz_seq,  # Sequence object
                            method = "TRATE",  # Method to determine costs
                            time.varying = FALSE) # Does not allow the cost to vary over time)
whz_cost_matrix

###Get costs----
####Optimal matching method----
whz_dist <- seqdist(whz_seq,
                     method = "HAM",
                     sm = whz_cost_matrix)

###Get clusters----
####PAM method----
whz_pamRange <- wcKMedRange(whz_dist, kvals=2:6) # this takes a while to run
summary(whz_pamRange, max.rank=6)

#Plot metrics
plot(whz_pamRange, stat = c("ASW","HC"), norm="zscore", lwd = 2, cex=2, col = c('#6666ff', '#cc0000'), legendpos = "topright", main = "HAM PAM Quality")
abline(v=7, col="#666666", lty="longdash", lwd = 2)

####Choose number of clusters----
whz_clust2 <- wcKMedoids(whz_dist, k=2, cluster.only=TRUE)
whz$clust2 <- wcKMedoids(whz_dist, k=2, cluster.only=TRUE)

####Modal plot----
whz_modal <- seqmsplot(whz_seq, group=whz_clust2, border=NA) #index plots by cluster

##Rhinitis----
###Cost matrix----
allergy_cost_matrix <- seqsubm(allergy_seq,  # Sequence object
                               method = "TRATE",  # Method to determine costs
                               time.varying = FALSE) # Does not allow the cost to vary over time)
allergy_cost_matrix

###Get costs----
####Optimal matching method----
allergy_dist <- seqdist(allergy_seq,
                        method = "HAM",
                        sm = allergy_cost_matrix)

###Get clusters----
####PAM method----
allergy_pamRange <- wcKMedRange(allergy_dist, kvals=2:6) # this takes a while to run
summary(allergy_pamRange, max.rank=6)

#Plot metrics
plot(allergy_pamRange, stat = c("ASW","HC"), norm="zscore", lwd = 2, cex=2, col = c('#6666ff', '#cc0000'), legendpos = "topright", main = "HAM PAM Quality")
abline(v=7, col="#666666", lty="longdash", lwd = 2)

####Choose number of clusters----
allergy_clust3 <- wcKMedoids(allergy_dist, k=4, cluster.only=TRUE)
whz$clust3 <- wcKMedoids(allergy_dist, k=4, cluster.only=TRUE)

####Modal plot----
allergy_modal <- seqmsplot(allergy_seq, group=allergy_clust3, border=NA) #index plots by cluster

