################################################################################################################################################
# Importing transcript abundance data and creating DESeqDataSet
# Code modelled after: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# The above link is recommended here: https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
################################################################################################################################################




# Install packages --------------------------------------------------------

#Uncomment if packages haven't been installed
#EBS below have been installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximportData")
# biocLite("TxDb.Celegans.UCSC.ce11.refGene")
# biocLite("tximport")




# Load packages into environment ------------------------------------------

library(tximportData)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(tximport)


# sample files ------------------------------------------------------------

#sample 1: N2_S1_L005
s1 <- "emilysalmonrun2/quants/N2_S1_L005_quant/quant.sf"
file.exists(s1)

#sample 2: N2_S1_L005
s2 <- "emilysalmonrun2/quants/UNC130_S2_L005_quant/quant.sf"
file.exists(s2)

#sample 3: N2_S1_L005
s3 <- "emilysalmonrun2/quants/UNC130_Sep7-2_S3_L005_quant/quant.sf"
file.exists(s3)

files <- c(s1,s2,s3)


# tx2gene dataframe -------------------------------------------------------

#The reference for the below library is here: http://bioconductor.org/packages/3.7/data/annotation/html/TxDb.Celegans.UCSC.ce11.refGene.html
#There are other libraries listed here: http://bioconductor.org/packages/3.7/data/annotation/

library(TxDb.Celegans.UCSC.ce11.refGene)
# note there is also TxDb.Celegans.UCSC.ce6.ensGene

txdb <- TxDb.Celegans.UCSC.ce11.refGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")





# Import Salmon Quants ----------------------------------------------------

library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)






































##########################################################
#EBS import transcript estimates
#EBS The next few steps are all related to this main step
#########################################################
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

#EBS creation of named vector pointing to quantification files
#EBS This is done by first reading in a table that contains sample IDs and then combining with dir and "quant.sf.gz"
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))

#EBS Transcript IDs need to be associated with gene IDs. Salmon only provides the transcript IDs.
#EBS This is done using the 'AnnotationDbi' package. 
#EBS Create tx2gene data frame

########################################################################
#EBS NOT SURE WHERE TO SOURCE THE GENE IDS - WAS THIS ALREADY CREATED?
#EBS I'm using an ensembl genome, so there may be another option.
#EBS I've put what I think are two options of doing the same thing. 

#### OPTION 1 #####
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- select(txdb, k, "GENEID", "TXNAME")

#### OPTION 2 #####
#EBS This uses the 'ensembldb' package
#EBS I'm not quite sure what the annotation is for C.elegans, but following their patten it should be something like: Ensdb.celegans.vXX
#EBS I don't know what version it would be, or if it's available for C.elegans. We used release 38, so maybe it's version 38?
#EBS Is there a way to checl the ensembldb package for which transcriptomes it has?

# library(readr)
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# head(tx2gene)

#EBS Not sure what this is, but it's the next block of code
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)


#EBS Give transcript level abundances (this may be important when I want to look at alternative splicing)
#EBS I don't think it's important right now, but if it doesn't affect downstream steps I'd like to include it
# txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

#EBS Summarize matrices
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)

######################################
# Specific to Salmon
######################################

#EBS I'm not sure where in the above steps this is supposed to occur
#EBS This is meant to import the files, so it probably comes first????

files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

#EBS If using a different transcriptome, the txgene2 needs to be read in
tx2knownGene <- read_csv(file.path(dir, "tx2gene.csv"))
files <- file.path(dir, "sailfish", samples$run, "quant.sf")
names(files) <- paste0("sample", 1:6)
txi.sailfish <- tximport(files, type = "sailfish", tx2gene = tx2knownGene)
head(txi.sailfish$counts)

#EBS This pertains to previous version of salmon. Since we just downloaded it this weekend I don't think it's relevant
#EBS I left it here as a comment just in case
# txi <- tximport("quant.sf", type = "none", txOut = TRUE, txIdCol = "Name", abundanceCol = "TPM", 
#                countsCol = "NumReads", lengthCol = "Length", importer = function(x) read_tsv(x, 
#                                                                                              skip = 8))


#EBS There is a section addressing inferential replicates in Salmon. This doesn't apply to my experiment.


###########################################
# USE WITH DOWNSTREAM BIOCONDUCTOR PACKAGES
###########################################
#EBS There are three suggestions for differential gene expression analysis packages. I'll include instruction for both


###############
# EdgeR
###############

#EBS Creating a DGEList for use with EdgeR

library(edgeR)


cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)

#################
# DESeq2
#################

#EBS Creating a DESeqDataSet for use with DESeq2
#EBS I think this package is discussed in the other workflow

library(DESeq2)

#EBS Make sure row names of sampleTable align with column names of txi$counts
sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

###################
# limma-voom
###################

#EBS creating a date object for use with limma-voom

files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)
y <- calcNormFactors(y)
design <- model.matrix(~condition, data = sampleTable)
v <- voom(y, design)




#########################
# END                   #
#########################




