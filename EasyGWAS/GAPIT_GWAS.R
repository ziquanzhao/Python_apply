install.packages("BiocManager")
BiocManager::install("multtest")+install.packages("gplots")+install.packages("LDheatmap")+install.packages("genetics")+install.packages("ape")+install.packages("EMMREML")+install.packages("scatterplot3d")
install.packages("bigmemory")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")


setwd("E:/R/reaslut/MYgwas")
myY <- read.table('TL222taxa.txt',header = T)
myG <- read.table('TLsample231_gene01.hmp.txt',header = F)

myY <- data.frame(myY[,1],myY[,4])
colnames(myY) <- c('Taxa','SFSN')

mygwas <- GAPIT(Y = myY,G = myG,
              kinship.algorithm = 'VanRaden',
              kinship.cluster = 'average',
              kinship.group = 'Mean',
              PCA.total = 3,  
              Model.selection = T,
              #Inter.Plot = T,
              model = c('FarmCPU','MLM','GLM'))



