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
whz<-read_csv("C:\\Users\\afossa\\Brown Dropbox\\Alan Fossa\\Data Re-Organization 2024\\Clean Data\\ALL_ECZ_WHZ.csv") %>% 
  select(participant_id,contains("whz_any")) %>% 
  mutate(
    across(contains("whz_any"),factor)
  )

whz$pmiss<-whz %>% 
  select(contains("whz_any")) %>% 
  apply(MARGIN=1,FUN=function(x){sum(is.na(x))/ncol(.)})

whz_small<-whz %>% 
  dplyr::filter(pmiss<0.5)

#Sequence index plot----

##Set state alphabet,labels, and colors----
state_alphabet <- c("0","1","99")

state_labels <- c("No wheeze","Wheeze","Missing")

no_wheeze<-"lightblue"
wheeze<-"salmon"
missing<-"beige"

##Define sequences----
whz_seq <- seqdef(
  whz_small, # Select data   
  var = 2:15, # Columns containing the sequences
  alphabet = state_alphabet,
  labels = state_labels,
  start=6, 
  missing=NA,
  left="99",
  gaps="99",
  right="99",
  void="%", 
  cpal=c(no_wheeze,wheeze,missing)
  )

View(whz_seq)

##Plot----
seqIplot(whz_seq, 
         sortv = "from.start",   # Sequence object
         with.legend = "right", # Display legend on right side of plot
         cex.legend = 0.6,  # Change size of legend
         main = "Sequence Index Plot", # Plot title
         xlab = "Age (months)"
         ) 


