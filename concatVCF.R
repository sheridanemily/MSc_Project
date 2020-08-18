Objective of this function: concatenate all files for a given cancer subtype so that it may be processed by SigProfilerJulia. SigProfilerJulia requires that the input file contains the following columns- CHROM, POS, REF, ALT, SAMPLE. 

Read in all vcf files ending in .gz. 
Convert the file to vcr2tidy data frame. 
Set to correct data frame ($fix). 
Find the nrow of the file.
Find the sample name (the sample name is the name of the file).
Create a vector in the dataframe that lists the sample name. 
Write it into a table. 
Eliminate all but the five necessary rows (CHROM, POS, REF, ALT, SAMPLE).
Repeat for remaining file that are .gz. 
Concatenate the files.


```{r}
library(vcfR)


dataframe <- list()
folder <- list.files("/path/to/file", pattern = ".gz")
f <- file.path("/path/to/file", folder)
rows = numeric(length = length(f))
for (i in 1:length(f)) {
 
  vcfs <- read.vcfR(f[i], convertNA = TRUE)
  rows[i] <- nrow(vcfs)
  vcfs <- vcfR2tidy(vcfs)
  dataframe[[i]] <- vcfs$fix
  samples <- as.character(f[i])
samples
 
  SAMPLE <- substr(samples, 52, 63)
  # You will need to change the sample substring beginning/end points depending on the path to file length
  dataframe[[i]]$SAMPLE <- rep(SAMPLE, rows[i])
  dataframe[[i]]
  dataframe[[i]] <- dataframe[[i]][, c('CHROM', 'POS', 'REF', 'ALT', 'SAMPLE')]
  # Optional
  dataframe[[i]]
  write.table(dataframe[[i]], sep = "\t")
  

}

dataframe_new <- dataframe[[1]]
for (i in 2:length(f)) {
    dataframe_new <- rbind(dataframe_new, dataframe[[i]]) 
      }
write.table(dataframe_new, "file_name.txt", sep = "\t")

```
