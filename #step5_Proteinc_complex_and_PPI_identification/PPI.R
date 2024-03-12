getwd()
#setwd("../../PrInCE normalized then scaling/")
library(PrInCE)
library(data.table)
require(devtools)
library('CCprofiler')

######update the format from CORUM database#######
complexHypotheses <- read.csv("./complexHypotheses_mod.csv") 
complexHypotheses <- setDT(complexHypotheses)
complexHypotheses_mod<- complexHypotheses[complexHypotheses$complex_id %in% complexHypotheses$complex_id[duplicated(complexHypotheses$complex_id)],]
complexHypotheses_mod<- complexHypotheses_mod[complexHypotheses_mod$idno %in% complexHypotheses_mod$idno[duplicated(complexHypotheses_mod$idno)],]
colnames(complexHypotheses_mod)<- c("complex_id","complex_name","protein_id")
binaryHypotheses <- generateBinaryNetwork(complexHypotheses_mod)
CORUM <- adjacency_matrix_from_data_frame(binaryHypotheses)
dim(CORUM)
##############################

###################REP1###########################
pcpData <- read.csv("scaled_rep1.csv", check.names = FALSE, row.names = "ID")
scaled_rep1<- pcpData
################################################################################

###################Build_features_rep1###########################################
set.seed(0)
gauss <- build_gaussians(scaled_rep1, 
                          max_gaussians = 5, min_R_squared = 0.75,
                          max_iterations = 50,criterion = "AICc",
                          filter_gaussians_center = TRUE,
                          filter_gaussians_height = 0.15,
                          filter_gaussians_variance_min = 0.2,
                          filter_gaussians_variance_max = 20,min_points = 1,
                          min_consecutive = 3,impute_NA = TRUE,
                          smooth = TRUE,smooth_width = 4) #max_iterations = 50

head(gauss)

# filter profiles that were not fit
scaled_rep1 <- scaled_rep1[names(gauss), ]
dim(scaled_rep1)
feat <- calculate_features(scaled_rep1, gauss)
dim(feat)

################################################################################

###################Identify Features#####################################
set.seed(0)
# ###classifier: the type of classifier to use; one of NB, SVM, RF, LR, or ensemble, corresponding to the options described above
pred_ppi <- predict_interactions(feat, CORUM, classifier = "RF",
                                  cv_folds = 10,verbose = TRUE, trees = 100)
 
precision <- pred_ppi$precision[seq_len(4e5)]
plot(precision)

net <- threshold_precision(pred_ppi, threshold = 0.95)
nrow(net)
write.csv(net,"./Rep1_Rsq0.75_RF_cv10_trees100_netthreshold0.95.csv")

#################################################################################


###################REP2_normalization and scaling###########################
pcpData2 <- read.csv("scaled_rep2.csv", check.names = FALSE, row.names = "ID")
scaled_rep2<- pcpData2
#############################################################################

###################Build_features_rep2###########################################
set.seed(0)
gauss2 <- build_gaussians(scaled_rep2, 
                           max_gaussians = 5, min_R_squared = 0.75,
                           max_iterations = 50,criterion = "AICc",
                           filter_gaussians_center = TRUE,
                           filter_gaussians_height = 0.15,
                           filter_gaussians_variance_min = 0.2,
                           filter_gaussians_variance_max = 20,
                           min_points = 1,
                           min_consecutive = 3,impute_NA = TRUE,
                           smooth = TRUE,smooth_width = 4)
 
# filter profiles that were not fit
scaled_rep2 <- scaled_rep2[names(gauss2), ]
dim(scaled_rep2)

#feat <- calculate_features(pcpData, gauss, euclidean_distance = TRUE)
feat2 <- calculate_features(scaled_rep2, gauss2)
dim(feat2)
################################################################################

################################################################################
###################Identify Features#####################################
set.seed(0)
# ###classifier: the type of classifier to use; one of NB, SVM, RF, LR, or ensemble, corresponding to the options described above
pred_ppi2 <- predict_interactions(feat2, CORUM, classifier = "RF",
                                   cv_folds = 10,verbose = TRUE, trees = 100)
dim(pred_ppi2)
precision2 <- pred_ppi2$precision
plot(precision2)
net2 <- threshold_precision(pred_ppi2, threshold = 0.95)
nrow(net2)
write.csv(net2,"./Rep2_PPI_Rsquare0.75_RFcv10_trees100_netthreshold0.95_CORUM.csv")
# ###########################################################
# ##########################################################



###############################################################################

############################Combined analysis###################################
feat_com <- concatenate_features(list(feat, feat2))
feat_comb<-na.omit(feat_com)


# predict interactions
# set the seed to ensure reproducible output
set.seed(0)
###classifier: the type of classifier to use; one of NB, SVM, RF, LR, or ensemble, corresponding to the options described above
pred_ppi_combine <- predict_interactions(feat_comb, CORUM, classifier = "RF",
                                         cv_folds = 10,verbose = TRUE, trees = 100)

dim(pred_ppi_combine)
precision_combine <- pred_ppi_combine$precision[seq_len(4e5)]
plot(precision_combine)

net_comb <- threshold_precision(pred_ppi_combine, threshold = 0.95)
nrow(net_comb)
write.csv(net_comb,"./unscaled_Both_modgauss5_Rsquare75_RF_cvfold10_trees100_netthres0.95_CORUM.csv")
################################################################