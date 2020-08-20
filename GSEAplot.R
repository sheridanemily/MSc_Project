Code to create a dumbbell plot in R of the top ten upregulated and top ten downregulated pathways determined by the fgsea package


```{r}

# GMT file containing pathways
pathways <- fgsea::gmtPathways("~/c5.all.v7.1.symbols.gmt")


gene_list <- stats_SKCM_SBS17b # Vector of gene names and their respective log2 fold changes


GSEA = function(gene_list, pathways, pval) {
  set.seed(5432)
  library(dplyr)
  library(gage)
  library(fgsea)
  library(ggplot2)
 
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
   
  #myGO = fgsea::gmtPathways(pathways)
  
  
  fgRes <- fgsea(pathways = pathways,
                           stats = gene_list)
  #print(fgRes)
  fgRes <- fgRes[order(padj), ]

 
# Filter FGSEA by using gage results. Must be significant and in same direction to keep
  gaRes = gage::gage(stats, gsets=pathways, same.dir=TRUE, set.size = c(16508))
  #head(gaRes)
 
  # Upregulated pathways
  ups = as.data.frame(gaRes$greater) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
 
  # Downregulated pathways
  downs = as.data.frame(gaRes$less) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
 
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]

  # Collapse redundant pathways

  Up = collapsePathways(keepups[order(keepups$pval)][keepups$padj < 0.05], pathways = pathways, gene_list, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns[order(keepdowns$pval)][keepdowns$padj < 0.05], pathways, gene_list,  nperm = 500, pval.threshold = 0.05)
  fgRes = fgRes[ !is.na(match(fgRes$pathway, c( Up$mainPathways, Down$mainPathways))), ] %>% arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "pathways" , "")
 
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                      "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes SKCM SBS17b") +
    theme_minimal()
  g
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}
library(dplyr)
library(readxl)


res = GSEA(stats, pathways, pval = 0.05)
res$Results
res$Plot
```
