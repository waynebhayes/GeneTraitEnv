#setwd("")
rm(list = ls())
library(stringr)
library(dplyr)
source("Shuffle.R")
source("Train_SA.R")
source("Kfold.R")
source("Validate.R")
source("Draw.R")
source("None_zero_idx_same.R")
source("Delete.R")
source("Ins_right.R")
source("Main.R")
source("Preprocessing.R")
#source("Add_synthesis.R")
source("Add_synthesis_2.R")
source("Main_synthetic.R")
source("GetpBad.R")

seed<-sample(-.Machine$integer.max:.Machine$integer.max, 1)
traits<-2
K<-4
obj_best<-0.00005

cat("Seed is ", seed, "\n")
result<-Main("Data_Matrix.csv", F, seed, traits, K, obj_best)
