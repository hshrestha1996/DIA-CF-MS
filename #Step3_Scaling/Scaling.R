######################################################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("PrInCE")
#browseVignettes("PrInCE")
#######################################################################

getwd()
setwd("./")

###################REP1_scaling###########################
pcpData <- read.csv("Normalized_and_preprocessed_common_rep1.csv", check.names = FALSE, row.names = "ID")
#head(pcpData, n = 4)
View(pcpData)
dim(pcpData)


scaled_rep1<- pcpData
scaled_rep1<- t(apply(scaled_rep1, 1, function(x)(x-min(x))/(max(x)-min(x))))
dim(scaled_rep1)
View(scaled_rep1)
write.csv(scaled_rep1, "./scaled_rep1.csv")
################################################################################

###################Repeat the above for REP2 ###########################