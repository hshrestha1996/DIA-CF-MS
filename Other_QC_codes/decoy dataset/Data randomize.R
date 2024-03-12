#install.packages("picante")
library('vegan')
library('picante')


df <- read.csv("Normalized_and_preprocessed_common_rep2.csv", check.names = FALSE, row.names = "Row.names")
gene <- df[1]
df <- df[,-1]
View(df)


set.seed(4)
df_shuffled <- randomizeMatrix(df,null.model = "richness",iterations = 1000)
#View(df_shuffled)
#df_shuffled
df_shuffled_gene <- merge(gene, df_shuffled, by = "row.names", all.x = TRUE)


write.csv(df_shuffled_gene,"./Decoy_192fractions_rep2.csv")
