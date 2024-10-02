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
  dplyr::filter(pmiss==0) %>% 
  select(-pmiss)

##Rhinitis----
allergy<-atopy %>% 
  select(participant_id,contains("allergy"),-food_allergy_P4) %>% 
  mutate(
    across(contains("eczema"),factor)
  ) %>% 
  select(sort(names(.)))

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
         main = "Sequence Index Plot", # Plot title
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
         main = "Sequence Index Plot", # Plot title
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
         main = "Sequence Index Plot", # Plot title
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
                   method = "OM",
                   indel= 1.0,
                   sm = derm_cost_matrix)

###Get clusters----
####PAM method----
derm_pamRange <- wcKMedRange(derm_dist, kvals=2:4) # this takes a while to run
summary(derm_pamRange, max.rank=4)

#Plot metrics
plot(derm_pamRange, stat = c("ASW","HC"), norm="zscore", lwd = 2, cex=2, col = c('#6666ff', '#cc0000'), legendpos = "topright", main = "OM PAM Quality")
abline(v=7, col="#666666", lty="longdash", lwd = 2)

####Choose number of clusters----
derm_clust4 <- wcKMedoids(derm_dist, k=4, cluster.only=TRUE)
derm$clust4 <- wcKMedoids(derm_dist, k=4, cluster.only=TRUE)

###Visualize clusters----

####Modal plot----
derm_modal <- seqmsplot(derm_seq, group=derm_clust4, border=NA) #index plots by cluster
