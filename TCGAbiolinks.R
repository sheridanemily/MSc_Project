TCGAbiolinks is a package used for integrative analysis with GDC data, such as the TCGA data used for this project. In the following chunk the objective is to create a data frame of information pertaining to clinical data of patients, such as if the patient was male/female, smoker/non-smoker, etc. Also requires the package SummarizedExperiment, the default data structure used in TCGAbiolinks.


```{r}
library(TCGAbiolinks)
library(SummarizedExperiment)

# Read in normalized data, in this case SKCM
SKCM_Norm <- read.delim("/path/to/file", header = TRUE, sep = " ")

SKCM_patients <- colnames(SKCM_Norm)

# Search GDC data using GDCquery function to filter data using predetermined arguments in the function 
SKCM_query <- GDCquery(project = "TCGA-SKCM", data.category = "Gene expression", platform = "Illumina HiSeq", file.type = "normalized_results" , legacy = TRUE, barcode = SKCM_patients)
barcodes<- getResults(SKCM_query, cols="cases")
length(barcodes)
primarytumour <- TCGAquery_SampleTypes(barcode = barcodes, typesample = "TP")
length(primarytumour)
SKCM_querylater <- GDCquery(project = "TCGA-SKCM", data.category = "Gene expression", data.type = "Gene expression quantification", platform = "Illumina HiSeq", file.type = "normalized_results" , legacy = TRUE, barcode = primarytumour)
barcodesl<- getResults(SKCM_querylater, cols="cases")
length(barcodesl)
t <- GDCdownload(SKCM_querylater, files.per.chunk = 10)

# GDCprepare function download and prepares the data from GDC for analysis
clinicaldata <- GDCprepare(SKCM_querylater, save.filename = "clinicaldataSKCM.rda")

# Create the dataframe
clinicaldataframe <- as.data.frame(colData(clinicaldata))
clinicaldataframe
load("clinicaldataSKCM.rda", skcmclindata <- new.env())
clinicaldataframe <- skcmclindata[["data"]]
clindata <- colData(clinicaldataframe)
```




```{r}


LUAD_Norm <- read.delim("/path/to/file", header = TRUE, sep = " ")


LUAD_patients <- colnames(LUAD_Norm)

LUAD_patients <- gsub('.', '-', LUAD_patients, fixed = TRUE)

# Search GDC data using GDCquery function to filter data using predetermined arguments in the function 
LUAD_query <- GDCquery(project = "TCGA-LUAD", data.category = "Gene expression", platform = "Illumina HiSeq", file.type = "normalized_results" , legacy = TRUE, barcode = LUAD_patients)
barcodes <- getResults(LUAD_query, cols="cases")
primarytumour <- TCGAquery_SampleTypes(barcode = barcodes, typesample = "TP")
LUAD_querylater <- GDCquery(project = "TCGA-LUAD", data.category = "Gene expression", data.type = "Gene expression quantification", platform = "Illumina HiSeq", file.type = "normalized_results" , legacy = TRUE, barcode = primarytumour)
barcodesl<- getResults(LUAD_querylater, cols="cases")
t <- GDCdownload(LUAD_querylater, files.per.chunk = 10)

# GDCprepare function download and prepares the data from GDC for analysis
clinicaldata <- GDCprepare(LUAD_querylater, save.filename = "clinicaldataLUAD.rda")

# Create the dataframe
clinicaldataframe <- as.data.frame(colData(clinicaldata))

load("clinicaldataLUAD.rda", luaddata <- new.env())
clinicaldataframe <- luaddata[["data"]]
clindata <- colData(clinicaldataframe)
```
