getwd()
setwd("./")
library(PrInCE)
library(data.table)
require(devtools)
library('CCprofiler')

###########################################################

##############Load Data#######################

scaled_rep1 <- read.csv("scaled_rep1.csv", check.names = FALSE, row.names = "ID")
scaled_rep2 <- read.csv("scaled_rep2.csv", check.names = FALSE, row.names = "ID")

####Run the complex detection for CORUM###########

set.seed(0)
######CORUM database formatted for complexes#######
complexHypotheses <- read.csv("CORUM complexes.csv")
complexHypotheses <- setDT(complexHypotheses)
complexHypotheses$subunits <- lapply(strsplit(as.character(complexHypotheses$subunits), ";"), unlist)
complexes <- t(complexHypotheses)
colnames(complexes) <- complexes[1,]
complexes<-complexes[-1,]
####################################################
detected_complex_rep1 <- detect_complexes(scaled_rep1, complexes,method = c("pearson","euclidean"),
                                          min_pairs = 10,
                                          bootstraps = 100,
                                          progress = TRUE)
# remove complexes that could not be analyzed
detected_complex_rep1 = na.omit(detected_complex_rep1)
# how many could be tested?
length(detected_complex_rep1)
# print the top complexes
detected_complex_rep1<-sort(detected_complex_rep1, decreasing = TRUE)
View(detected_complex_rep1)
# how many were significant at uncorrected, two-tailed p < 0.05?
sum(detected_complex_rep1 > 1.96)
write.csv(detected_complex_rep1,"./detected_complex_rep1.csv")

#########Repeat for rep2################
detected_complex_rep2 <- detect_complexes(scaled_rep2, complexes,method = c("pearson","euclidean"),
                                          min_pairs = 10,
                                          bootstraps = 100,
                                          progress = TRUE)
# remove complexes that could not be analyzed
detected_complex_rep2 = na.omit(detected_complex_rep2)
# how many could be tested?
length(detected_complex_rep2)
# print the top complexes
detected_complex_rep2<-sort(detected_complex_rep2, decreasing = TRUE)
View(detected_complex_rep2)
# how many were significant at uncorrected, two-tailed p < 0.05?
sum(detected_complex_rep2 > 1.96)
write.csv(detected_complex_rep2,"./detected_complex_rep2.csv")
#############################################################



#############################################################
####Run the complex detection for huMAP2###########
set.seed(0)
######Humap database formatted for complexes#######
complexHypotheses_humap <- read.csv("hu_map2_complexes.csv")
complexHypotheses_humap <- setDT(complexHypotheses_humap)
complexHypotheses_humap$subunits <- lapply(strsplit(as.character(complexHypotheses_humap$subunits), ";"), unlist)
complexes_humap <- t(complexHypotheses_humap)
colnames(complexes_humap) <- complexes_humap[1,]
complexes_humap<-complexes_humap[-1,]
####################################################
detected_complex_rep1_humap <- detect_complexes(scaled_rep1, complexes_humap,method = c("pearson","euclidean"),
                                                min_pairs = 10,
                                                bootstraps = 100,
                                                progress = TRUE)
# remove complexes that could not be analyzed
detected_complex_rep1_humap = na.omit(detected_complex_rep1_humap)
# how many could be tested?
length(detected_complex_rep1_humap)
# print the top complexes
detected_complex_rep1_humap<-sort(detected_complex_rep1_humap, decreasing = TRUE)
View(detected_complex_rep1_humap)
# how many were significant at uncorrected, two-tailed p < 0.05?
sum(detected_complex_rep1_humap > 1.96)
write.csv(detected_complex_rep1_humap,"./detected_complex_rep1_humap.csv")

###############rep2##########
detected_complex_rep2_humap <- detect_complexes(scaled_rep2, complexes_humap,method = c("pearson","euclidean"),
                                                min_pairs = 10,
                                                bootstraps = 100,
                                                progress = TRUE)

# remove complexes that could not be analyzed
detected_complex_rep2_humap = na.omit(detected_complex_rep2_humap)
# how many could be tested?
length(detected_complex_rep2_humap)
# print the top complexes
detected_complex_rep2_humap<-sort(detected_complex_rep2_humap, decreasing = TRUE)
View(detected_complex_rep2_humap)
# how many were significant at uncorrected, two-tailed p < 0.05?
sum(detected_complex_rep2_humap > 1.96)
write.csv(detected_complex_rep2_humap,"./detected_complex_rep2_humap.csv")
############ Complexes done##########
