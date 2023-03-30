library(stringr)

## 1. single cell csv store barcodes, sample_name and celltypes
df1 <- read.csv("../data/sham/singlecell.csv")
df2 <- read.csv("../data/day3/singlecell.csv")
df3 <- read.csv("../data/day10/singlecell.csv")
df <- rbind(df1, df2, df3)

df <- subset(df, select = c("X", "sample", "celltype"))
df$celltype <- gsub("/", ".", df$celltype)
colnames(df) <- c("Barcode", "Sample", "Cluster")


## 2. split df by the sample name
df_healthy_female  <- subset(df, Sample == "Healthy_Female")
df_healthy_male    <- subset(df, Sample == "Healthy_Male")
df_day3_female     <- subset(df, Sample == "Day3_Female")
df_day3_male       <- subset(df, Sample == "Day3_Male")
df_day10_female    <- subset(df, Sample == "Day10_Female")
df_day10_male      <- subset(df, Sample == "Day10_Male")


df_healthy_female$Sample  <- NULL 
df_healthy_male$Sample    <- NULL 
df_day3_female$Sample     <- NULL 
df_day3_male$Sample       <- NULL 
df_day10_female$Sample    <- NULL 
df_day10_male$Sample      <- NULL 


##
dir.create("../data/clusters")
write.table(df_healthy_female,file = "../data/clusters/Healthy_Female.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_healthy_male  ,file = "../data/clusters/Healthy_Male.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_day3_female   ,file = "../data/clusters/Day3_Female.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_day3_male     ,file = "../data/clusters/Day3_Male.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_day10_female  ,file = "../data/clusters/Day10_Female.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(df_day10_male    ,file = "../data/clusters/Day10_Male.txt", row.names = FALSE, sep = "\t", quote = FALSE)
