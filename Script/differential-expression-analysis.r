
#######################################
# Differential expression (DE) analysis
#######################################


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
plot(data.C.log2.mean ~ data.N.log2.mean, xlab = "N", ylab = "C",
     main = "Scatter", xlim = c(0, limit), ylim = c(0, limit))
## Diagonal line
abline(0, 1, col = "red")


# Compute fold-change -----------------------------------------------------

## Compute biological significance
## Difference between the means of the conditions
fold <- data.N.log2.mean - data.C.log2.mean

## Histogram of the fold differences
hist(fold, col = "gray")



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


## Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

## Volcano: put the biological significance (fold-change)
## and statistical significance (p-value) in one plot
plot(fold, -log10(pvalue), main = "Volcano")

fold_cutoff <- 2.5
pvalue_cutoff <- 0.05
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)



# Screen for the genes that satisfy the filtering criteria ----------------

## Fold-change filter for "biological" significance
filter.by.fold <- abs(fold) >= fold_cutoff

## P-value filter for "statistical" significance
filter.by.pvalue <- pvalue <= pvalue_cutoff

## Combined filter (both biological and statistical)
filter.combined <- filter.by.fold & filter.by.pvalue

data.N.DE <- data.N[filter.combined,]
data.C.DE <- data.C[filter.combined,]

dim(data.N.DE)
dim(data.C.DE)


## Let's generate the volcano plot again,
## highlighting the significantly differential expressed genes
plot(fold, -log10(pvalue), main = "Volcano #2")
points (fold[filter.combined], -log10(pvalue[filter.combined]),
        pch = 16, col = "red")

## Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "Volcano #3")
points (fold[filter.combined & fold < 0],
        -log10(pvalue[filter.combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter.combined & fold > 0],
        -log10(pvalue[filter.combined & fold > 0]),
        pch = 16, col = "blue")


# Save the DE genes -------------------------------------------------------
write.csv(data.N.DE,
          paste0(DATA_FOLD,"dataN_DE.csv"))
write.csv(data.C.DE,
          paste0(DATA_FOLD,"dataC_DE.csv"))

