getwd()
setwd("../Desktop/Monthly meeting/Manuscript/rebuttal/Preprocessing/")

library(PrInCE)
library(data.table)
require(devtools)


################Rep1####################################

##Load the datamatrix for preprocessing

df1= read.csv("For_BSA_normalization_rep1.csv", row.names = "Protein.Group")
gene1 <-df1[1]
proteins1 <- df1[2:ncol(df1)]

##Calculate factor based on Internal Std
BSA1 <- proteins1[grep("P02769", rownames(proteins1)), , drop = FALSE]
factor1 <- BSA1 / rowMeans(BSA1, na.rm = TRUE)
write.csv(factor1, file = "Normalization_factor_rep1.csv")

## Normalized by IS
corrected_data1 <- proteins1 / factor1[rep(1, nrow(proteins1)), ]
corrected_data_gene1 <- merge(gene1, corrected_data1, by = "row.names", all.x = TRUE)
write.csv(corrected_data_gene1, file = "Normalization_datamatrix_rep1.csv")

##Filteration for min consectutive
mydata1<-corrected_data_gene1
row.names(mydata1)<- mydata1$Row.names
mydata1 <- mydata1[3:ncol(mydata1)]
filtered1 = filter_profiles(mydata1,min_points = 1,min_consecutive = 3)

##Cleaning the protein profile
chromatograms_rep1 = clean_profiles(filtered1, smooth_width = 4) #smooth_width = 4
write.csv(chromatograms_rep1,"./Normalized_and_preprocessed_rep1.csv")


################Repeat the process for Rep2####################################
df2= read.csv("For_BSA_normalization_rep2.csv", row.names = "Protein.Group")
gene2 <-df2[1]
proteins2 <- df2[2:ncol(df2)]

##Calculate factor based on Internal Std
BSA2 <- proteins2[grep("P02769", rownames(proteins2)), , drop = FALSE]
factor2 <- BSA2 / rowMeans(BSA2, na.rm = TRUE)
write.csv(factor2, file = "Normalization_factor_rep2.csv")

## Normalized by IS
corrected_data2 <- proteins2 / factor2[rep(1, nrow(proteins2)), ]
corrected_data_gene2 <- merge(gene2, corrected_data2, by = "row.names", all.x = TRUE)
write.csv(corrected_data_gene2, file = "Normalization_datamatrix_rep2.csv")

##Filteration for min consectutive
mydata2<-corrected_data_gene2
row.names(mydata2)<- mydata2$Row.names
mydata2 <- mydata2[3:ncol(mydata2)]
filtered2 = filter_profiles(mydata2,min_points = 1,min_consecutive = 3)

##Cleaning the protein profile
chromatograms_rep2 = clean_profiles(filtered2, smooth_width = 4) #smooth_width = 4
write.csv(chromatograms_rep2,"./Normalized_and_preprocessed_rep2.csv")


######Getting the dataframe from common proteins#######

common_rep1 <- chromatograms_rep1[rownames(chromatograms_rep1) %in% rownames(chromatograms_rep2), , drop = FALSE]
common_rep1<- merge(gene1, common_rep1, by = "row.names", all.x = FALSE)
write.csv(common_rep1,"./Normalized_and_preprocessed_common_rep1.csv")


common_rep2 <- chromatograms_rep2[rownames(chromatograms_rep2) %in% rownames(chromatograms_rep1), , drop = FALSE]
common_rep2<- merge(gene2, common_rep2, by = "row.names", all.x = FALSE)
write.csv(common_rep2,"./Normalized_and_preprocessed_common_rep2.csv")

######Preprocessing complete##################