##########################################################################
##########################################################################
# Project: MS lineage embryogenesis in C elegans
# Script purpose: predicting TF regulation activity 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct  7 17:31:46 2020
##########################################################################
##########################################################################
install.magic = FALSE
if(install.magic){
  install.packages("Rmagic")
  system('which python')
  system('python --version')
  system('pip install --user magic-impute')
  
}

##########################################
# test MAGIC and it works for the example
# original code from https://github.com/KrishnaswamyLab/MAGIC
##########################################
library(Rmagic)
library(ggplot2)
data(magic_testdata)

ss = apply(magic_testdata, 2, sum)
magic_testdata = magic_testdata[, ss>0]
MAGIC_data <- Rmagic::magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"), verbose = 1)

ggplot(MAGIC_data) +
  geom_point(aes(x=VIM, y=CDH1, color=ZEB1))


