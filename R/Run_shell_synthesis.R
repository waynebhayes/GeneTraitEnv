library(data.table)
library(stringr)
library(dplyr)

source("../../Shuffle.R")
source("../../Kfold.R")
source("../../Validate.R")
source("../../Draw.R")
source("../../None_zero_idx_same.R")
source("../../Delete.R")
source("../../Ins_right.R")
source("../../Main_synthetic.R")
source("../../Train_SA.R")
source("../../Preprocessing.R")
source("../../Add_synthesis.R")
source("../../Add_synthesis_2.R")


args=commandArgs(T)
input<-args[1]
seed<-as.integer(args[2])
traits<-as.integer(args[3])
K<-as.integer(args[4])
result<-Main_synthetic(input, F, seed, traits, K, 0.05)

