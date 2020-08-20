Packages and code to perform linear regression, to find correlation between signatures of unknown aetiology and genes, and to perform gene set enrichment analysis on these genes.


```{r}
library(limma)
library(statmod)
library(devtools)
library(MAST)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(plyr)
library(AnnotationFilter)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(formattable)
library(tidyr)
library(gridExtra)
library(xtable)
register(SerialParam())

# Linear regression to find correlation between genes and signatures extracted from SKCM


# Read in file containing TPM normalized RNA-HTSeq data
SKCM_Norm <- read.delim("/Users/maahg/Documents/SKCM_Primary_Tumor_log2_TPM.txt", header = TRUE, sep = " ")

# Read in exposures/mutational weights
SKCM_exposures <- read.delim("/Users/maahg/Documents/SKCM.snvs.exposures.tsv", header = TRUE, sep = "\t")

# Ensure normalized data and exposures have all the same samples
SKCM_exposures <- SKCM_exposures[,colnames(SKCM_exposures) %in% c(colnames(SKCM_Norm))]
SKCM_Norm <- SKCM_Norm[,colnames(SKCM_Norm) %in% c(colnames(SKCM_exposures))]


colnames(SKCM_exposures) <- gsub('.', '-', colnames(SKCM_exposures), fixed = TRUE)
colnames(SKCM_Norm) <- gsub('.', '-', colnames(SKCM_exposures), fixed = TRUE)

rownames(SKCM_exposures) <-  c("SBS49", "SBS7a(1)", "SBS6", "SBS7b", "SBS7a(2)", "SBS17b")


SKCM_Norm <- as.matrix(SKCM_Norm)
SKCM_exposures <- as.matrix(SKCM_exposures)

# Check dataframes
#dim(SKCM_Norm)
#dim(SKCM_exposures)
#head(SKCM_exposures)
#head(SKCM_Norm)

# Linear regression

reg_SKCM <- apply(SKCM_Norm, 1, function(x) lm(x ~ scale(SKCM_exposures[1,], scale = 103) + scale(SKCM_exposures[2,], scale = 103) + scale(SKCM_exposures[3,], scale = 103) + scale(SKCM_exposures[4,], scale = 103) + scale(SKCM_exposures[5,], scale = 103) + scale(SKCM_exposures[6,]) - 1)$coefficients)

# Row six contains signature SBS17b, we perform linear regression on this row

SKCM_SBS17b <- reg_SKCM[6,]

# Checking to see if there is a normal distribution of log2 fold change
plot(density(SKCM_SBS17b))

SKCM_SBS17b <- as.data.frame(SKCM_SBS17b)
SKCM_SBS17b <- as.matrix(SKCM_SBS17b)


# Convert ensemble gene names to entrez


pathways <- fgsea::gmtPathways("~/c5.all.v7.1.symbols.gmt") # Download file containing pathways

# Create dataframe of gene names, including gene ID, gene symbol, and entrez ID
genes_SKCM_SBS17b <- ensembldb::select(EnsDb.Hsapiens.v86, columns = c("SYMBOL", "GENEID", "ENTREZID"), keytype = "GENEID", keys = rownames(SKCM_SBS17b))

# Remove duplicates
genes_SKCM_SBS17b2 <- genes_SKCM_SBS17b[!duplicated(genes_SKCM_SBS17b$GENEID),]

genes_SKCM_SBS17b2 <- as.data.frame(genes_SKCM_SBS17b2)

# Ensure only genes from original data frame are used
SKCM_SBS17b2 <- SKCM_SBS17b[which(rownames(SKCM_SBS17b) %in% (genes_SKCM_SBS17b2$GENEID)),]

# Ensure new dataframe contains all original gene names
all.equal(rownames(SKCM_SBS17b2), genes_SKCM_SBS17b$GENEID)


stats_SKCM_SBS17b <- SKCM_SBS17b2
names(stats_SKCM_SBS17b) <- genes_SKCM_SBS17b2$SYMBOL
#stats_SKCM_SBS17b <- stats_SKCM_SBS17b[!duplicated(stats_SKCM_SBS17b)]


stats_SKCM_SBS17b <- stats_SKCM_SBS17b[!duplicated(names(stats_SKCM_SBS17b))]

#length(stats_SKCM_SBS17b)

# fgsea requires the fold change be sorted in decreasing order
stats_SKCM_SBS17b <- sort.int(stats_SKCM_SBS17b, decreasing = TRUE)

stats_SKCM_SBS17b <- na.omit(stats_SKCM_SBS17b)
save(stats_SKCM_SBS17b, file = "stats_SKCM_SBS17b.Rdata")
stats_SKCM_SBS17b2 <- as.data.frame(stats_SKCM_SBS17b)
save(stats_SKCM_SBS17b2, file = "stats_SKCM_SBS17b2.Rdata")
stats_SKCM_SBS17b <- as.vector(stats_SKCM_SBS17b)

# Gene set enrichment analysis using fgsea

SKCM_SBS17b_gsea <- fgsea(pathways = pathways, stats = stats_SKCM_SBS17b, minSize = 1, maxSize = Inf, gseaParam = 1, eps = 0)
save(SKCM_SBS17b_gsea, file = "~/SKCM_SBS17b_gsea.Rdata")

SKCM_SBS17b_gsea_ordered <- (SKCM_SBS17b_gsea[order(padj), ])

fwrite(SKCM_SBS17b_gsea, file = "SKCM_SBS17b_gsea.txt", sep = "\t", sep2 = c("", " ", ""))


# Plotting a GSEA table using a built in function from fgsea package
topPathwaysUp <- SKCM_SBS17b_gsea[ES > 0][head(order(pval), n=10), SKCM_SBS17b_gsea$pathway]
topPathwaysDown <- SKCM_SBS17b_gsea[ES < 0][head(order(pval), n=10), SKCM_SBS17b_gsea$pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGSEA_SKCMSBS17b <- plotGseaTable(pathways[topPathways], stats_SKCM_SBS17b, SKCM_SBS17b_gsea, gseaParam=1)
save(plotGSEA_SKCMSBS17b, file = "~/plotGSEA_SKCMSBS17b.Rdata")

# Collapse pathways. This function collapses the list of enriched pathways to independent ones.
collapsedPathwaysSKCM <- collapsePathways(SKCM_SBS17b_gsea[order(pval)][padj < 0.01], pathways, stats_SKCM_SBS17b)

#mainPathwaysSKCM <- SKCM_SBS17b_gsea[pathway %in% collapsedPathwaysSKCM$mainPathwaysSKCM][order(-NES), pathway]

SKCMSBS17b_plotmain <- plotGseaTable(pathways[mainPathwaysSKCM], stats_SKCM_SBS17b, SKCM_SBS17b_gsea, gseaParam = 0.5)

save(SKCMSBS17b_plotmain, file = "~/SKCMSBS17b_plotmain.Rdata")

mainpathways_SKCMSBS17b <- collapsedPathwaysSKCM$mainPathways
mainpathways_SKCMSBS17b <- as.data.frame(mainpathways_SKCMSBS17b)
try <- formattable::formattable(mainpathways_SKCMSBS17b, format = "markdown")

# Create a table of the collapsed pathways in the same format as the main pathway table output by fgsea
mainpathways_SKCMSBS17b <- SKCM_SBS17b_gsea[SKCM_SBS17b_gsea$pathway %in% mainpathways_SKCMSBS17b$mainpathways_SKCMSBS17b]
mainpathways_SKCMSBS17b <- as.data.frame(mainpathways_SKCMSBS17b)
mainpathways_SKCMSBS17b <- mainpathways_SKCMSBS17b[order(NES)]
save(mainpathways_SKCMSBS17b, file = "~/SKCM_mainpath_df.Rdata")

# Saving data in two tables of upregulated and downregulated pathways
SKCM_upreg <- mainpathways_SKCMSBS17b[mainpathways_SKCMSBS17b$NES > 0]
SKCM_upreg <- SKCM_upreg[order(NES)]
save(SKCM_upreg, file = "~/SKCM_upreg.Rdata")
SKCM_downreg <- mainpathways_SKCMSBS17b[mainpathways_SKCMSBS17b$NES < 0]
SKCM_downreg <- SKCM_downreg[order(NES)]
save(SKCM_downreg, file = "~/SKCM_downreg.Rdata")

# Converting data frame to xtable, which allows the table to be printed into latex

print(xtable(mainpathways_SKCMSBS17b))
png("SKCM_table_GSEA.png", height = 5*nrow(mainpathways_SKCMSBS17b), width = 20*ncol(mainpathways_SKCMSBS17b))
S <- tableGrob(mainpathways_SKCMSBS17b)
grid.arrange(S)
dev.off()

# Subset dataframe of main collapsed pathways, only keeping columns I feel are vital for writeup
SKCM_x <- mainpathways_SKCMSBS17b[,c(1, 4:7)]

# Rounding pvals to four significant digits so they may fit into a table in latex
pval <- format.pval(mainpathways_SKCMSBS17b[,2], digits = 4)
p.adj <- format.pval(mainpathways_SKCMSBS17b[,3], digits = 4)
pval <- as.data.frame(pval)
p.adj <- as.data.frame(p.adj)


SKCM_xx <- cbind(SKCM_x, pval, p.adj)

# Create xtable
SKCM_xtable <- xtable(SKCM_xx)
print(xtable(SKCM_xtable), type = "latex")

```








```{r}

# Linear regression to find correlation between genes and signatures extracted from LUAD


LUAD_Norm <- read.delim("~/LUAD_Primary_Tumor_log2_TPM.txt", header = TRUE, sep = " ")


LUAD_exposures <- read.delim("~/LUAD.snvs.exposures.tsv", header = TRUE, sep = "\t")

LUAD_exposures <- LUAD_exposures[,colnames(LUAD_exposures) %in% c(colnames(LUAD_Norm))]
LUAD_Norm <- LUAD_Norm[,colnames(LUAD_Norm) %in% c(colnames(LUAD_exposures))]

rownames(LUAD_exposures) <-  c("SBS6", "SBS18", "SBS2", "SBS49", "SBS17b", "SBS26", "SBS4(1)", "SBS4(2)")

colnames(LUAD_exposures) <- gsub('.', '-', colnames(LUAD_exposures), fixed = TRUE)
colnames(LUAD_Norm) <- gsub('.', '-', colnames(LUAD_exposures), fixed = TRUE)


LUAD_Norm <- as.matrix(LUAD_Norm)
LUAD_exposures <- as.matrix(LUAD_exposures)

dim(LUAD_Norm)
dim(LUAD_exposures)
head(LUAD_exposures)
head(LUAD_Norm)

# Linear regression

reg_LUAD <- apply(LUAD_Norm, 1, function(x) lm(x ~ scale(LUAD_exposures[1,], scale = 508) + scale(LUAD_exposures[2,], scale = 508) + scale(LUAD_exposures[3,], scale = 508) + scale(LUAD_exposures[4,], scale = 508) + scale(LUAD_exposures[5,], scale = 508) + scale(LUAD_exposures[6,], scale = 508) + scale(LUAD_exposures[7,], scale = 508) + scale(LUAD_exposures[8,], scale = 508) - 1)$coefficients)


# Row five contains SBS17b
LUAD_SBS17b <- reg_LUAD[5,]
LUAD_SBS17b2 <- as.data.frame(LUAD_SBS17b)
LUAD_SBS17b2 <- as.matrix(LUAD_SBS17b)

# Checking to see if there is a normal distribution of log2 fold change
plot(density(LUAD_SBS17b))

# Convert ensemble gene names to entrez

pathways <- fgsea::gmtPathways("~/c5.all.v7.1.symbols.gmt")

# Create dataframe of gene names, including gene ID, gene symbol, and entrez ID
genes_LUAD_SBS17b <- ensembldb::select(EnsDb.Hsapiens.v86, columns = c("SYMBOL", "GENEID", "ENTREZID"), keytype = "GENEID", keys = rownames(LUAD_SBS17b2))

LUAD_SBS17b2 <- LUAD_SBS17b2[rownames(LUAD_SBS17b2) %in% (genes_LUAD_SBS17b$GENEID),]
LUAD_SBS17b2 <- as.data.frame(LUAD_SBS17b2)

genes_LUAD_SBS17b <- genes_LUAD_SBS17b[!duplicated(genes_LUAD_SBS17b$GENEID),]

#dim(genes_LUAD_SBS17b)
#dim(LUAD_SBS17b2)

# Ensure new dataframe contains all original gene names
all.equal(rownames(LUAD_SBS17b2), genes_LUAD_SBS17b$GENEID)


stats_LUAD_SBS17b <- LUAD_SBS17b2$LUAD_SBS17b
names(stats_LUAD_SBS17b) <- genes_LUAD_SBS17b$SYMBOL
#stats_LUAD_SBS17b <- stats_LUAD_SBS17b[!duplicated(stats_LUAD_SBS17b)]


stats_LUAD_SBS17b <- stats_LUAD_SBS17b[!duplicated(names(stats_LUAD_SBS17b))]

length(stats_LUAD_SBS17b)

# fgsea requires the fold change be sorted in decreasing order
stats_LUAD_SBS17b <- sort.int(stats_LUAD_SBS17b, decreasing = TRUE)

stats_LUAD_SBS17b <- na.omit(stats_LUAD_SBS17b)
save(stats_LUAD_SBS17b, file = "stats_LUAD_SBS17b.Rdata")
stats_LUAD_SBS17b2 <- as.data.frame(stats_LUAD_SBS17b)
save(stats_LUAD_SBS17b2, file = "stats_LUAD_SBS17b2.Rdata")
stats_LUAD_SBS17b <- as.vector(stats_LUAD_SBS17b)

# Gene set enrichment analysis using fgsea

LUAD_SBS17b_gsea <- fgsea(pathways = pathways, stats = stats_LUAD_SBS17b, minSize = 1, maxSize = Inf, gseaParam = 1, eps = 0)

head(LUAD_SBS17b_gsea[order(pval), ])

fwrite(LUAD_SBS17b_gsea, file = "LUAD_SBS17b_gsea.txt", sep = "\t", sep2 = c("", " ", ""))


# Plotting a GSEA table of the top under and overexpressed pathways using a built in function from fgsea package
topPathwaysUp <- LUAD_SBS17b_gsea[ES > 0][head(order(pval), n=10), LUAD_SBS17b_gsea$pathway]
topPathwaysDown <- LUAD_SBS17b_gsea[ES < 0][head(order(pval), n=10), LUAD_SBS17b_gsea$pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGSEA_LUADSBS17b <- plotGseaTable(pathways[topPathways], stats_LUAD_SBS17b, LUAD_SBS17b_gsea, gseaParam=1)
save(plotGSEA_LUADSBS17b, file = "~/plotGSEA_LUADSBS17b.Rdata")

# Collapse pathways

collapsedPathwaysLUAD <- collapsePathways(LUAD_SBS17b_gsea[order(pval)][padj < 0.01], pathways, stats_LUAD_SBS17b)

# Order pathways based on Normalized enrichment score
mainPathwaysLUAD <- LUAD_SBS17b_gsea[pathway %in% collapsedPathwaysLUAD$mainPathwaysLUAD][order(-NES), pathway]

LUADSBS17b_plotmain <- plotGseaTable(pathways[mainPathwaysLUAD], stats_LUAD_SBS17b, LUAD_SBS17b_gsea, gseaParam = 0.5)
save(LUADSBS17b_plotmain, file = "~/LUADSBS17b_plotmain.Rdata")

LUADSBS17b_mainpath <- collapsedPathwaysLUAD$mainPathways
LUADSBS17b_mainpath <- as.data.frame(LUADSBS17b_mainpath)

collapsedPathwaysLUAD$parentPathways
length(collapsedPathwaysLUAD$parentPathways)

# Create a table of the collapsed pathways in the same format as the main pathway table output by fgsea
mainpathways_LUADSBS17b <- LUAD_SBS17b_gsea[LUAD_SBS17b_gsea$pathway %in% LUADSBS17b_mainpath$LUADSBS17b_mainpath]
mainpathways_LUADSBS17b <- mainpathways_LUADSBS17b[order(NES)]
save(mainpathways_LUADSBS17b, file = "~/LUADSBS17b_pathtable.Rdata")
write.csv(mainpathways_LUADSBS17b, file = "~/LUAD_paths.csv")

# Saving data in two tables of upregulated and downregulated pathways
LUAD_upreg <- mainpathways_LUADSBS17b[mainpathways_LUADSBS17b$NES > 0]
LUAD_upreg <- LUAD_upreg[order(padj)]
save(LUAD_upreg, file = "~/LUAD_upreg.Rdata")
LUAD_downreg <- mainpathways_LUADSBS17b[mainpathways_LUADSBS17b$NES < 0]
LUAD_downreg <- LUAD_downreg[order(padj)]
save(LUAD_downreg, file = "~/LUAD_downreg.Rdata")


# Subset dataframe of main collapsed pathways, only keeping columns I feel are vital for writeup
luad_x <- mainpathways_LUADSBS17b[,c(1, 4:7)]

# Rounding pvals to four significant digits

pval <- format.pval(mainpathways_LUADSBS17b[,2], digits = 4)
p.adj <- format.pval(mainpathways_LUADSBS17b[,3], digits = 4)
pval <- as.data.frame(pval)
p.adj <- as.data.frame(p.adj)

luad_xx <- cbind(luad_x, pval, p.adj)

# Create xtable
library(xtable)
LUAD_xtable <- xtable(luad_xx)
print(xtable(LUAD_xtable), type = "latex")



```




```{r}


# Linear regression to find correlation between genes and signatures extracted from LUSC


LUAD_Norm <- read.delim("~/LUAD_Primary_Tumor_log2_TPM.txt", header = TRUE, sep = " ")


LUAD_exposures <- read.delim("~/LUAD.snvs.exposures.tsv", header = TRUE, sep = "\t")

LUAD_exposures <- LUAD_exposures[,colnames(LUAD_exposures) %in% c(colnames(LUAD_Norm))]
LUAD_Norm <- LUAD_Norm[,colnames(LUAD_Norm) %in% c(colnames(LUAD_exposures))]

rownames(LUAD_exposures) <-  c("SBS6", "SBS18", "SBS2", "SBS49", "SBS17b", "SBS26", "SBS4(1)", "SBS4(2)")

colnames(LUAD_exposures) <- gsub('.', '-', colnames(LUAD_exposures), fixed = TRUE)
colnames(LUAD_Norm) <- gsub('.', '-', colnames(LUAD_exposures), fixed = TRUE)


LUAD_Norm <- as.matrix(LUAD_Norm)
LUAD_exposures <- as.matrix(LUAD_exposures)

dim(LUAD_Norm)
dim(LUAD_exposures)
head(LUAD_exposures)
head(LUAD_Norm)

# Linear regression

reg_LUAD <- apply(LUAD_Norm, 1, function(x) lm(x ~ scale(LUAD_exposures[1,], scale = 508) + scale(LUAD_exposures[2,], scale = 508) + scale(LUAD_exposures[3,], scale = 508) + scale(LUAD_exposures[4,], scale = 508) + scale(LUAD_exposures[5,], scale = 508) + scale(LUAD_exposures[6,], scale = 508) + scale(LUAD_exposures[7,], scale = 508) + scale(LUAD_exposures[8,], scale = 508) - 1)$coefficients)


# Row five contains SBS17b
LUAD_SBS17b <- reg_LUAD[5,]
LUAD_SBS17b2 <- as.data.frame(LUAD_SBS17b)
LUAD_SBS17b2 <- as.matrix(LUAD_SBS17b)

# Checking to see if there is a normal distribution of log2 fold change
plot(density(LUAD_SBS17b))

# Convert ensemble gene names to entrez

pathways <- fgsea::gmtPathways("~/c5.all.v7.1.symbols.gmt")

# Create dataframe of gene names, including gene ID, gene symbol, and entrez ID
genes_LUAD_SBS17b <- ensembldb::select(EnsDb.Hsapiens.v86, columns = c("SYMBOL", "GENEID", "ENTREZID"), keytype = "GENEID", keys = rownames(LUAD_SBS17b2))

LUAD_SBS17b2 <- LUAD_SBS17b2[rownames(LUAD_SBS17b2) %in% (genes_LUAD_SBS17b$GENEID),]
LUAD_SBS17b2 <- as.data.frame(LUAD_SBS17b2)

genes_LUAD_SBS17b <- genes_LUAD_SBS17b[!duplicated(genes_LUAD_SBS17b$GENEID),]

#dim(genes_LUAD_SBS17b)
#dim(LUAD_SBS17b2)

# Ensure new dataframe contains all original gene names
all.equal(rownames(LUAD_SBS17b2), genes_LUAD_SBS17b$GENEID)


stats_LUAD_SBS17b <- LUAD_SBS17b2$LUAD_SBS17b
names(stats_LUAD_SBS17b) <- genes_LUAD_SBS17b$SYMBOL
#stats_LUAD_SBS17b <- stats_LUAD_SBS17b[!duplicated(stats_LUAD_SBS17b)]


stats_LUAD_SBS17b <- stats_LUAD_SBS17b[!duplicated(names(stats_LUAD_SBS17b))]

length(stats_LUAD_SBS17b)

# fgsea requires the fold change be sorted in decreasing order
stats_LUAD_SBS17b <- sort.int(stats_LUAD_SBS17b, decreasing = TRUE)

stats_LUAD_SBS17b <- na.omit(stats_LUAD_SBS17b)
save(stats_LUAD_SBS17b, file = "stats_LUAD_SBS17b.Rdata")
stats_LUAD_SBS17b2 <- as.data.frame(stats_LUAD_SBS17b)
save(stats_LUAD_SBS17b2, file = "stats_LUAD_SBS17b2.Rdata")
stats_LUAD_SBS17b <- as.vector(stats_LUAD_SBS17b)

# Gene set enrichment analysis using fgsea

LUAD_SBS17b_gsea <- fgsea(pathways = pathways, stats = stats_LUAD_SBS17b, minSize = 1, maxSize = Inf, gseaParam = 1, eps = 0)

head(LUAD_SBS17b_gsea[order(pval), ])

fwrite(LUAD_SBS17b_gsea, file = "LUAD_SBS17b_gsea.txt", sep = "\t", sep2 = c("", " ", ""))


# Plotting a GSEA table of the top under and overexpressed pathways using a built in function from fgsea package
topPathwaysUp <- LUAD_SBS17b_gsea[ES > 0][head(order(pval), n=10), LUAD_SBS17b_gsea$pathway]
topPathwaysDown <- LUAD_SBS17b_gsea[ES < 0][head(order(pval), n=10), LUAD_SBS17b_gsea$pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGSEA_LUADSBS17b <- plotGseaTable(pathways[topPathways], stats_LUAD_SBS17b, LUAD_SBS17b_gsea, gseaParam=1)
save(plotGSEA_LUADSBS17b, file = "~/plotGSEA_LUADSBS17b.Rdata")

# Collapse pathways

collapsedPathwaysLUAD <- collapsePathways(LUAD_SBS17b_gsea[order(pval)][padj < 0.01], pathways, stats_LUAD_SBS17b)

# Order pathways based on Normalized enrichment score
mainPathwaysLUAD <- LUAD_SBS17b_gsea[pathway %in% collapsedPathwaysLUAD$mainPathwaysLUAD][order(-NES), pathway]

LUADSBS17b_plotmain <- plotGseaTable(pathways[mainPathwaysLUAD], stats_LUAD_SBS17b, LUAD_SBS17b_gsea, gseaParam = 0.5)
save(LUADSBS17b_plotmain, file = "~/LUADSBS17b_plotmain.Rdata")

LUADSBS17b_mainpath <- collapsedPathwaysLUAD$mainPathways
LUADSBS17b_mainpath <- as.data.frame(LUADSBS17b_mainpath)

collapsedPathwaysLUAD$parentPathways
length(collapsedPathwaysLUAD$parentPathways)

# Create a table of the collapsed pathways in the same format as the main pathway table output by fgsea
mainpathways_LUADSBS17b <- LUAD_SBS17b_gsea[LUAD_SBS17b_gsea$pathway %in% LUADSBS17b_mainpath$LUADSBS17b_mainpath]
mainpathways_LUADSBS17b <- mainpathways_LUADSBS17b[order(NES)]
save(mainpathways_LUADSBS17b, file = "~/LUADSBS17b_pathtable.Rdata")
write.csv(mainpathways_LUADSBS17b, file = "~/LUAD_paths.csv")

# Saving data in two tables of upregulated and downregulated pathways
LUAD_upreg <- mainpathways_LUADSBS17b[mainpathways_LUADSBS17b$NES > 0]
LUAD_upreg <- LUAD_upreg[order(padj)]
save(LUAD_upreg, file = "~/LUAD_upreg.Rdata")
LUAD_downreg <- mainpathways_LUADSBS17b[mainpathways_LUADSBS17b$NES < 0]
LUAD_downreg <- LUAD_downreg[order(padj)]
save(LUAD_downreg, file = "~/LUAD_downreg.Rdata")


# Subset dataframe of main collapsed pathways, only keeping columns I feel are vital for writeup
luad_x <- mainpathways_LUADSBS17b[,c(1, 4:7)]

# Rounding pvals to four significant digits

pval <- format.pval(mainpathways_LUADSBS17b[,2], digits = 4)
p.adj <- format.pval(mainpathways_LUADSBS17b[,3], digits = 4)
pval <- as.data.frame(pval)
p.adj <- as.data.frame(p.adj)

luad_xx <- cbind(luad_x, pval, p.adj)

# Create xtable
library(xtable)
LUAD_xtable <- xtable(luad_xx)
print(xtable(LUAD_xtable), type = "latex")



```
