Principal component analysis (PCA) plots show clusters of samples based on their similarity. For these PCA plots we can see the similarity between the mutational weights of the signatures. The mutational weight is the number of mutations that a sample has, which are attributable to a particular mutational signature.

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
