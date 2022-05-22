#install.packages("devtools")
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#install.packages("vcfR")
library(vcfR)
library(GAPIT3)

setwd("~/MA_Bioinformatics/Thesis_Project/05_GWAS")

#Step 1: Set data directory and import files
vcf <- read.vcfR( "manta_overlap.vcf", verbose = FALSE )

##numeric format
gt <- extract.gt(vcf)
myGD <- gsub('0/0', '0', gsub('1/1', '2', gsub('0/1', '1', gt)))
myGD <- t(myGD) #transpose
class(myGD) <- "numeric"
rownames(myGD) <- gsub("\\..*","",rownames(myGD))
myGD <- cbind(rownames(myGD), data.frame(myGD, row.names=NULL))
names(myGD)[1] <- "taxa"

##genetic map
myGM <- as.data.frame(vcf@fix[,c('ID','CHROM','POS')])

##phenotypic data
pheno <- read.table('days2FwH.txt',head=T)
pheno <- pheno[ , -which(names(pheno) %in% c("Country","GenePool"))]
#pheno <- pheno[pheno$Accession %in% myGD$taxa,]

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=pheno,
  GD=myGD,
  GM=myGM,
  PCA.total=4,
  SNP.MAF=0.05,
  Multiple_analysis=TRUE,
  model=c("MLM","MLMM","SUPER","FarmCPU","Blink")
)

