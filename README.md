# UCLA-AD-vs-NA-RNA-Seq-Analysis

#Package and Cran Install
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", force=TRUE)

  BiocManager::install("DESeq2", force=TRUE)
  BiocManager::install("EnhancedVolcano", force=TRUE)
  BiocManager::install("pheatmap", force=TRUE)
  BiocManager::install("RColorBrewer", force=TRUE)
  BiocManager::install("gProfileR", force=TRUE)
  BiocManager::install("apeglm", force=TRUE)
  BiocManager::install("ashr", force=TRUE)
  BiocManager::install("DEGreport", force=TRUE)
  BiocManager::install("pathwview", force = TRUE)
  BiocManager::install("clusterProfiler", force = TRUE)
  BiocManager::install("edgeR", force = TRUE)
  BiocManager::install("AnnotationDbi", force = TRUE)
  BiocManager::install("clusterProfiler", force = TRUE)
  install.packages("ggfortify")
  install.packages("tibble")
  install.packages("dplyr")
```

#Library Install
```
library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")
library("gProfileR")
library("apeglm")
library("ashr")
library("DEGreport")
library("AnnotationDbi")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggfortify")
library(dplyr)
library(tibble)
library(tidyr)
```

#Load Count Data and Meta Data
```
cts <- read.csv("Data/ALL Cell Lines Raw Counts.new.protcodeing.csv") 
rownames <- cts[ , 1]

# remove first column (gene names) 
cts <- cts[ , -1]

#convert matrix to numbers
cts <- apply(cts, 2, as.numeric)

# move gene names to rownames
rownames(cts) <- rownames

data <- as.data.frame(cts)


#Convert all counts to integers
data[,1] <- as.integer(data[,1])
data[,2] <- as.integer(data[,2])
data[,3] <- as.integer(data[,3])
data[,4] <- as.integer(data[,4])
data[,5] <- as.integer(data[,5])
data[,6] <- as.integer(data[,6])
data[,7] <- as.integer(data[,7])
data[,8] <- as.integer(data[,8])
data[,9] <- as.integer(data[,9])
data[,10] <- as.integer(data[,10])
data[,11] <- as.integer(data[,11])
data[,12] <- as.integer(data[,12])
data[,13] <- as.integer(data[,13])
data[,14] <- as.integer(data[,14])
data[,15] <- as.integer(data[,15])
data[,16] <- as.integer(data[,16])
data[,17] <- as.integer(data[,17])
data[,18] <- as.integer(data[,18])


#Import Metadata file
meta <- read.csv("Meta/ALL Cell Lines RAW Meta.new.protcodeing.csv", row.names =1, na.strings = "")
```

#Data Model Check
```
mean_counts <- apply(data[, 1:18], 1, mean)
variance_counts <- apply(data[, 1:18], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() + scale_x_log10()

#Checks column names
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/09daaae2-a6cf-44c5-bedb-aca4b96ce32b)

#Creating DESEQ2 Object for Sample Analysis
```
#Create DESEQ Object
dds.new <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Cell.Line + Condition)
dds.new$Condition <- relevel(dds.new$Condition, ref = "AD")

#Generate Normalized Counts
dds.new <- estimateSizeFactors(dds.new)
sizeFactors(dds.new)
normalized_counts <- counts(dds.new, normalized=TRUE)
```

# Sample Clustering and PCA Plot
```
#Quality Control

#Transform normalized counts using the rlog transformation
rld <- rlog(dds.new, blind=F)

# Input is a matrix of log transformed values
rld <- rlog(dds.new, blind=F)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
pcaplot <- autoplot(pca, data = meta , colour = "Condition", shape = "Cell.Line", size = 3) + theme(text=element_text(family="Times New Roman", face="bold", size=12))
pcaplot
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/bd828e59-a83a-42ed-83b9-f120a00976be)

Sample Clustering Heatmap
```
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

pheatmap(rld_cor)
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/2d1af065-959f-4a13-933a-36dbba2203bb)

Dispersion Estimates 
```
## Create DESeq object
dds.new <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Cell.Line + Condition)


## Run analysis
dds.new <- DESeq(dds.new, test = "Wald")

## Total number of raw counts per sample
colSums(counts(dds.new))

## Total number of normalized counts per sample
colSums(counts(dds.new, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds.new)
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/43075edf-4f4a-4306-8536-77e406aa0596)

Sample Comparison Contrasts and MA Plot
```
## Define contrasts, extract results table, and shrink the log2 fold changes
contrast_adhNa <- c("Condition","NA","AD")

res_tableAdhNA_unshrunken <- results(dds.new, contrast= contrast_adhNa, alpha = 0.05)
res_table.new <- lfcShrink(dds.new, contrast= contrast_adhNa, res=res_tableAdhNA_unshrunken, type = "ashr")


plotMA(res_tableAdhNA_unshrunken, ylim=c(-10,10))
plotMA(res_table.new, ylim=c(-10,10))
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/6c8c1cda-689f-406f-a060-7648e5f57353)
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/7a3565e7-0e02-4224-a3f0-c6981f798f0c)

Significant Genes
```
## Summarize results
summary(res_table.new)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

##########Creating Significant Genes Variable

res_table.new_tb <- res_table.new %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

sigGenes <- res_table.new_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sigGenes
```

Gene Expression Heat Map
Gene expression does not cluster by Gene Expression
```
norm_OEsig <- normalized_counts[,c(1,2:19)] %>% 
  dplyr::filter(gene %in% sigGenes$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### Annotate our heatmap (optional)
annotation <- adhNA_meta %>% 
  dplyr::select(samplename, Condition) %>% 
  data.frame(row.names = "samplename")

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

geneheatmap <- pheatmap(norm_OEsig, 
                        color = heat_colors, 
                        cluster_rows = T, 
                        show_rownames = F,
                        legend_labels = "Log2FoldChange\n",
                        annotation = annotation, 
                        border_color = NA, 
                        fontsize = 12, 
                        scale = "row",
                        cluster_cols = T,
                        fontfamily = "serif", angle_col = 0)
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/d7154417-f99c-4112-96ad-aecde5a803ec)

Volcano Plot
```
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table.new_tb <- res_table.new_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.5)

## Create a column to indicate which genes to label
res_table.new_tb <- res_table.new_tb %>% arrange(padj) %>% mutate(genelabels = "")
res_table.new_tb$genelabels[1:20] <- res_table.new_tb$gene[1:20]

allvolcanoplot <- EnhancedVolcano(res_table.new_tb,
                                  lab = res_table.new_tb$genelabels,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  pCutoff = 0.000000001,
                                  FCcutoff = 0.5,
                                  pointSize = 2.0,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = TRUE,
                                  colAlpha = 4/5,
                                  legendPosition = "top",
                                  legendLabSize = 12,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5,
                                  colConnectors = 'black',
                                  title = "HK336 AD vs NA DEGs",
                                  subtitle = "Differentially Expressed Genes",
                                  caption = "") + theme(axis.text = element_text(size = 12),plot.subtitle = element_text(hjust = 0.5 ),plot.title = element_text(hjust = 0.5 )) + 
  theme(text=element_text(family="Times New Roman", face="bold", size=12)) 
allvolcanoplot
```
![image](https://github.com/anubhavchandla/UCLA-AD-vs-NA-RNA-Seq-Analysis/assets/166166948/5e712f8f-a241-4b94-be5c-faeca00764b9)









