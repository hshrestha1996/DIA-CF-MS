getwd()
setwd("./")

options(java.parameters = "-Xmx8000m")
library(diann)
library(scales)
library(data.table)
library(reshape2)
library(stringr)
library(wrProteo)
library(openxlsx)

################################################################################
count_n <- function(df, from='pep', to='prot') {
  # count the number of identified peptides for a protein
  # df is diann report
  # from can be 'prec', 'pep';
  # to can be 'pep', 'prot'
  # supports prec -> pep, prec -> prot, pep -> prot
  
  # define from
  if (from == 'pep'){
    key <- paste0(df$Run, df$Stripped.Sequence)
  } else if (from == 'prec') {
    key <- paste0(df$Run, df$Precursor.Id)
  }
  x <- df[!duplicated(key),c('Run','Protein.Names','Stripped.Sequence','Precursor.Id')]
  
  # 
  if (to == 'prot'){
    y <- aggregate(x, by = list(x$Run, x$Protein.Names), FUN = length)
  } else if (to == 'pep') {
    y <- aggregate(x, by = list(x$Run, x$Stripped.Sequence), FUN = length)
  }
  
  y <- dcast(y, Group.2 ~ Group.1, value.var = 'Run')
  rownames(y) <- y$Group.2
  y$Group.2 <- NULL
  
  return(y)
}
################################################################################

# read DIA-NN report file
data <- diann_load("report.tsv")

# filter by FDR
MQ <- 0.01 # FDR cutoff
q <- MQ # precursor level q
pgq <- MQ # protein group level q
##df <- data[which(data$Q.Value <= q & data$PG.Q.Value <= pgq),]

libq <- MQ # library precursor q

#df <- data[which(data$Q.Value <= q & data$PG.Q.Value <= pgq & data$Lib.Q.Value <= libq &data$Lib.PG.Q.Value <=libq),] 

df <- data[which(data$Q.Value <= q & data$PG.Q.Value <= pgq & data$Lib.Q.Value <= libq),] 
dim(df)

# keep precursors corresponding to unique protein groups
UNIQUE <- T 
if (UNIQUE) df <- df[!grepl(';',df$Protein.Names),] 
dim(df)
#####CCprofiler input######
cc <- df[,c(2,3,14,28)]
cc<- cc[!duplicated(cc),]
dim(cc)
write.csv(cc,"./CC_inputv2.csv")

#####CCprofiler input######
######end#######


# get protein information from fasta file
prot_fasta <- readFasta2("Human_swissprot_20422_March15_and_BSA.fasta", tableOut = T) 
# replace the above fasta file with the file you used to create the library

# extract precursor information
prec_info <- df[,c("Protein.Group","Protein.Ids","Protein.Names","Genes",
                   "Stripped.Sequence","Precursor.Id","Precursor.Charge")]
prec_info <- prec_info[!duplicated(prec_info),]
prec_info <- merge(prec_info, prot_fasta[,2:3], by.x="Protein.Group", by.y="uniqueIdentifier", all.x = T)
colnames(prec_info)[c(4,8)] <- c("Gene.Symbol","Protein.Description")
# extract protein information
prot_info <- prec_info[,c(1,3,4,8)]
prot_info <- prot_info[!duplicated(prot_info),]
# extrat peptide information
pep_info <- prec_info[,c(1:5,8)]
pep_info <- pep_info[!duplicated(pep_info),]

# precursor-level quantification
prec_quan <- diann_matrix(df)
colnames(prec_quan) <- basename(colnames(prec_quan))
n_prec_quan <- colSums(!is.na(prec_quan))
print(n_prec_quan)
# add precursor info
prec_quan_anno <- as.data.frame(prec_quan)
prec_quan_anno$Precursor.Id <- rownames(prec_quan_anno)
prec_quan_anno <- merge(prec_info, prec_quan_anno, by="Precursor.Id")

# peptide-level quantification
pep_quan <- diann_maxlfq(df, sample.header = "Run",
                         group.header="Stripped.Sequence", 
                         id.header = "Precursor.Id", 
                         quantity.header = "Precursor.Normalised")
n_pep_quan <- colSums(!is.na(pep_quan))
print(n_pep_quan)
# add peptide info
pep_quan_anno <- as.data.frame(pep_quan)
pep_quan_anno$Stripped.Sequence <- rownames(pep_quan_anno)
pep_quan_anno <- merge(pep_info, pep_quan_anno, by="Stripped.Sequence")

# protein-level quantification
prot_quan <- diann_maxlfq(df, sample.header = "Run",
                          group.header="Protein.Names", 
                          id.header = "Precursor.Id", 
                          quantity.header = "Precursor.Normalised")
n_prot_quan <- colSums(!is.na(prot_quan))
print(n_prot_quan)

# add protein info
prot_quan_anno <- as.data.frame(prot_quan)
prot_quan_anno$Protein.Names <- rownames(prot_quan_anno)
prot_quan_anno <- merge(prot_info, prot_quan_anno, by="Protein.Names", all.x = F, all.y = F)

# filter protein quantification by # of peptides/precursors
prot_n_pep <- count_n(df, from='pep', to='prot')

# output 
write.xlsx(prot_quan_anno, file="prot_quan_pub1.xlsx", sheetName="Protein-level quan", rowNames=FALSE)
write.xlsx(prot_n_pep, file="prot_quan_pub2.xlsx", sheetName="# peptide per protein", append=TRUE, rowNames=TRUE)
write.xlsx(pep_quan_anno, file="pep_quan_pub3.xlsx", sheetName="Peptide-level quan", rowNames=FALSE)
write.xlsx(prec_quan_anno, file="prec_quan_pub.xlsx", sheetName="Precursor-level quan", row.names=FALSE)

