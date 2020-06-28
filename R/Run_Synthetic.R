setwd("C:/Users/RW/Documents/Study/UCI/Research/Real Data/4-20 Automatic Decision of Temperature Range")

library(stringr)
library(dplyr)
library(data.table)
source("Shuffle.R")
source("Kfold.R")
source("Validate.R")
source("Draw.R")
source("None_zero_idx_same.R")
source("Delete.R")
source("Ins_right.R")
source("Main_synthetic.R")
source("Train_SA.R")
source("Preprocessing.R")
source("Add_synthesis.R")
source("Add_synthesis_2.R")
source("GetpBad_synthetic.R")
source("Train_SA_synthetic.R")
source("Main_synthetic.R")

traits<-1000
K<-4
obj_best<-0.00005
seed<-5

cat("Seed is ", seed, "\n")
result<-Main_synthetic("Data_Matrix.csv", F, seed, traits, K, obj_best)
