Principal component analysis (PCA) plots show clusters of samples based on their similarity. For these PCA plots we can see the similarity between the mutational weights of the signatures. The mutational weight is the number of mutations that a sample has, which are attributable to a particular mutational signature.

Also included is an attempt at a TSNE plot. These are used to visualize high-dimensional datasets in a low-dimentional space.

```{r}
library(devtools)
library(ggbiplot)

# It should be noted that the mutational weights matrix is also known as the 'exposures' matrix
SKCM_exposures <- read.csv("/path/to/file", header = TRUE, sep = "\t")

SKCM_exposures <- data.frame(SKCM_exposures, row.names = 1)

# The signatures extracted
rownames(SKCM_exposures) <- c("SBS49", "SBS7a(1)", "SBS6", "SBS7b", "SBS7a(2)", "SBS17b")
pcaSKCM <- prcomp(SKCM_exposures, scale. = TRUE)
summary(pcaSKCM)
str(pcaSKCM)
SKCM_pca <- ggbiplot(pcaSKCM, labels = rownames(SKCM_exposures), var.axes = FALSE) +
  theme_minimal() +
  ggtitle("PCA of signature weights: SKCM")

ggsave(SKCM_pca, file = "~/SKCM_pca2.png")
  
  
```



```{r}


LUAD_exposures <- read.csv("/path/to/file", header = TRUE, sep = "\t")

LUAD_exposures <- data.frame(LUAD_exposures, row.names = 1)

# Extracted signatures
rownames(LUAD_exposures) <- c("SBS6", "SBS18", "SBS2", "SBS49", "SBS17b", "SBS26", "SBS4(1)", "SBS4(2)")
pcaLUAD <- prcomp(LUAD_exposures, scale. = TRUE)
summary(pcaLUAD)
str(pcaLUAD)
pca_LUAD <- ggbiplot(pcaLUAD, labels = rownames(LUAD_exposures), var.axes = FALSE) +
  theme_minimal() +
  ggtitle("PCA of signature weights: LUAD")

ggsave(pca_LUAD, file = "~/LUAD_pca.png")


```


```{r}

LUSC_exposures <- read.csv("/path/to/file", header = TRUE, sep = "\t")

LUSC_exposures <- data.frame(LUSC_exposures, row.names = 1)

# Extracted signatures
rownames(LUSC_exposures) <- c("SBS4", "SBS5", "SBS15", "SBS7b", "SBS13", "SBS49")
pcaLUSC <- prcomp(LUSC_exposures, scale. = TRUE)
summary(pcaLUSC)
str(pcaLUSC)
LUSC_pca <- ggbiplot(pcaLUSC, labels = rownames(LUSC_exposures), var.axes = FALSE) +
  theme_minimal() +
  ggtitle("PCA's of signature weights: LUSC")

ggsave(LUSC_pca, file = "LUSC_pca.png")

# Arrange all PCA's into the same frame
allPCA <- SKCM_pca + pca_LUAD + LUSC_pca
ggsave(allPCA, file = "~/allPCA.png", width = 12, height = 7)
```



```{r}
library(Rtsne)

SKCM_SBS17b$genenames <- rownames(SKCM_SBS17b)
LUAD_SBS17b$genenames <- rownames(LUAD_SBS17b)
LUSC_SBS5$genenames <- rownames(LUSC_SBS5)

# Selecting genes enriched for signatures of unknown aetiology from each subtype
# Selected based on fold change of genes

SKCM_SBS17b2 <- SKCM_SBS17b[which(abs(SKCM_SBS17b$SKCM_SBS17b )>= 1 ),]

LUAD_SBS17b2 <- LUAD_SBS17b[which(abs(LUAD_SBS17b$LUAD_SBS17b)>= 1 ),]

LUSC_SBS52 <- LUSC_SBS5[which(abs(LUSC_SBS5$LUSC_SBS5)>= 1),]

SKCM_SBS17b2 <- SKCM_SBS17b2[which(SKCM_SBS17b2$SKCM_SBS17b) %in% LUAD_SBS17b2$LUAD_SBS17b,]


SKCM_Norm <- read.delim("/path/to/file", header = TRUE, sep = " ")
LUAD_Norm <- read.delim("/path/to/file", header = TRUE, sep = " ")
LUSC_Norm <- read.delim("/path/to/file", header = TRUE, sep = " ")

#SKCM_Norm <- as.data.frame(SKCM_Norm)
#LUAD_Norm <- as.data.frame(LUAD_Norm)
#LUSC_Norm <- as.data.frame(LUSC_Norm)

SKCM_Norm$genenames <- rownames(SKCM_Norm)
LUAD_Norm$genenames <- rownames(LUAD_Norm)
LUSC_Norm$genenames <- rownames(LUSC_Norm)

#dim(SKCM_Norm)
#dim(LUAD_Norm)


related_genes_SKCM <- SKCM_Norm[which(SKCM_Norm$genenames %in% LUAD_Norm$genenames),]
related_genes_SKCM <- related_genes_SKCM[which(related_genes_SKCM$genenames %in% LUSC_Norm$genenames),]

related_genes_LUAD <- LUAD_Norm[which(LUAD_Norm$genenames %in% related_genes_SKCM$genenames),]

related_genes_LUSC <- LUSC_Norm[which(LUSC_Norm$genenames %in% related_genes_LUAD$genenames),]


all.equal(rownames(related_genes_LUAD), rownames(related_genes_LUSC))
all.equal(rownames(related_genes_LUAD), rownames(related_genes_SKCM))

related_genes_SKCM2 <- related_genes_SKCM[which(related_genes_SKCM$genenames %in% SKCM_SBS17b2$genenames),]
related_genes_LUAD2 <- related_genes_LUAD[which(related_genes_LUAD$genenames %in% LUAD_SBS17b2$genenames),]
related_genes_LUSC2 <- related_genes_LUSC[which(related_genes_LUSC$genenames %in% LUSC_SBS52$genenames),]

related_genes_SKCM2 <- related_genes_SKCM2[which(related_genes_SKCM2$genenames %in% related_genes_LUAD2$genenames),]
related_genes_SKCM2 <- related_genes_SKCM2[which(related_genes_SKCM2$genenames %in% related_genes_LUSC2$genenames),]
related_genes_LUAD2 <- related_genes_LUAD2[which(related_genes_LUAD2$genenames %in% related_genes_SKCM2$genenames),]
related_genes_LUSC2 <- related_genes_LUSC2[which(related_genes_LUSC2$genenames %in% related_genes_SKCM2$genenames),]

all.equal(rownames(related_genes_SKCM2), rownames(related_genes_LUAD2))
all.equal(rownames(related_genes_LUAD2), rownames(related_genes_LUSC2))

related_genes_LUAD2$genenames <- NULL
related_genes_SKCM2$genenames <- NULL
related_genes_LUSC2$genenames <- NULL

dim(related_genes_SKCM2) # 103
dim(related_genes_LUAD2) # 513
dim(related_genes_LUSC2) # 501

SKCM_col <- rep("SKCM", 103)
LUAD_col <- rep("LUAD", 513)
LUSC_col <- rep("LUSC", 501)

colnames(related_genes_SKCM2) <- SKCM_col
colnames(related_genes_LUAD2) <- LUAD_col
colnames(related_genes_LUSC2) <- LUSC_col

cbind_norms <- cbind(related_genes_SKCM2, related_genes_LUAD2, related_genes_LUSC2)

tsne_plot <- Rtsne(cbind_norms, perplexity = 80, check_duplicates = FALSE, theta = 0.0, normalize = F, max_iter = 1000)

#related_genes <- SKCM_Norm[which(SKCM_SBS17b2$genenames %in% LUAD_SBS17b2$genenames),]
#related_genes <- related_genes[which(related_genes$genenames %in% LUSC_SBS52$genenames),]

#related_genes <- as.data.frame(related_genes)
#related_genes <- related_genes$related_genes[which(related_genes$related_genes) %in% rownames(LUSC_SBS5)]
```


