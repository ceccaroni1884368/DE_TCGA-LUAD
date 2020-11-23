#######################################
# Differential expression (DE) analysis
#######################################

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#
#BiocManager::install('EnhancedVolcano')
#install.packages("GeoTcgaData")
library(EnhancedVolcano)
library("GeoTcgaData")
library(R.matlab)

# Settings ----------------------------------------------------------------

BASE_FOLD = '~/GitHub/DE_TCGA-LUAD/'
DATA_FOLD = paste0(BASE_FOLD, 'TCGA-LUAD/data/')


# Load data ---------------------------------------------------------------

load(paste0(DATA_FOLD,'dataCN.RData'))


# Log2 transformation -----------------------------------------------------

data.N.log2 <- log2(data.N)
data.C.log2 <- log2(data.C)

## Compute the means of the samples of each condition
data.N.log2.mean <- apply(data.N.log2, 1, mean)
data.C.log2.mean <- apply(data.C.log2, 1, mean)

## Just get the maximum of all the means
limit <- max(data.N.log2.mean, data.C.log2.mean)

## Scatter plot
plot(data.C.log2.mean ~ data.N.log2.mean, xlab = "Normal", ylab = "Cancer",
     main = "Scatter", xlim = c(0, limit), ylim = c(0, limit))
## Diagonal line
abline(0, 1, col = "red")


# Compute fold-change -----------------------------------------------------

## Compute biological significance
## Difference between the means of the conditions
FC <- data.C.log2.mean - data.N.log2.mean

## Histogram of the fold differences
hist(FC, col = "gray")



# Compute statistical significance (using t-test) -------------------------

pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(data.C)) { # For each gene : 
  x <- data.N.log2[i,] # N of gene number i
  y <- data.C.log2[i,] # C of gene number i
  
  # Compute t-test between the two conditions
  t <- t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] <- t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] <- t$statistic
}


# Correction for multiple comparison --------------------------------------
## Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

## FDR
FDR <-p.adjust(pvalue, method = "fdr")

## Histogram of p-values (-log10)
hist(-log10(FDR), col = "gray")


# Thresholds selection ----------------------------------------------------

FC.threshold <- 2.5
FDR.threshold <- 0.05

## Volcano: put the biological significance (fold-change)
## and statistical significance (p-value) in one plot
plot(FC, -log10(FDR), main = "Volcano")

abline(v = FC.threshold, col = "blue", lwd = 3)
abline(v = -FC.threshold, col = "red", lwd = 3)
abline(h = -log10(FDR.threshold), col = "green", lwd = 3)



# Screen for the genes that satisfy the filtering criteria ----------------

## Fold-change filter for "biological" significance
filter.by.FC <- abs(FC) >= FC.threshold

## P-value filter for "statistical" significance
filter.by.FDR <- FDR <= FDR.threshold

## Combined filter (both biological and statistical)
filter.combined <- filter.by.FC & filter.by.FDR

data.N.DE <- data.N[filter.combined,]
data.C.DE <- data.C[filter.combined,]

dim(data.N.DE)
dim(data.C.DE)


# Volcano Plot ------------------------------------------------------------

res <- as.data.frame(cbind(FC, FDR))
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'FC',
                y = 'FDR',
                title = 'Cancer versus Normal',
                subtitle = '',
                caption = '',
                pCutoff = FDR.threshold,
                FCcutoff = FC.threshold,
                pointSize = 3.0,
                labSize = 0,
                colAlpha = 1)


# Save the DE genes -------------------------------------------------------

gene.ID.DE <- rownames(data.C.DE)
mean.FC.DE <- FC[filter.by.FC]
rm(list=setdiff(ls(), c("data.N.DE", "data.C.DE",
                        "FC.threshold", "FDR.threshold",
                        "gene.ID.DE", "mean.FC.DE",
                        "DATA_FOLD")))
gene.sym.DE <- id_conversion_vector("Ensembl_ID","symbol", gene.ID.DE)

## data.N.DE: Normal tissue data of Differentially expressed genes (DEGs)
## data.C.DE: Cancer tissue data of DEGs
## gene.ID.DE: gene IDs of DEGs
## gene.sym.DE: gene symbols of DEGs
## FDR.threshold: FDR max (parameter used to define DEGs)
## FC.threshold: minimun fold change (parameter used to define DEGs)
## mean.FC.DE: average Fold Chance of DEGs among the patients


write.csv(data.N.DE,
          paste0(DATA_FOLD,"dataN_DE.csv"))
write.csv(data.C.DE,
          paste0(DATA_FOLD,"dataC_DE.csv"))

save.image(file= paste0(DATA_FOLD,"dataCN_DE.RData"))

## Matlab export
filename <- paste0(DATA_FOLD,"data_DE.mat")
writeMat(filename,
         dataNDe = as.matrix(data.N.DE),
         dataCDe = as.matrix(data.C.DE),
         gene_IDDe = gene.ID.DE,
         gene_symDe = gene.sym.DE,
         fdrt = FDR.threshold,
         fct = FC.threshold,
         mFCDe = mean.FC.DE)
