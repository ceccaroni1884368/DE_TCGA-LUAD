
# Settings ----------------------------------------------------------------

BASE_FOLD = '~/GitHub/DE_TCGA-LUAD/'
DATA_FOLD = paste0(BASE_FOLD, 'TCGA-LUAD/')


# Import Data -------------------------------------------------------------

rna.expr.data.C <- read.csv(paste0(DATA_FOLD, "/TCGA-LUAD_rna_expr_data_C.txt"),
                            sep="")
rna.expr.data.N <- read.csv(paste0(DATA_FOLD, "/TCGA-LUAD_rna_expr_data_N.txt"),
                            sep="")
rna.genes.C <- scan(paste0(DATA_FOLD, "/TCGA-LUAD_rna_genes_C.txt"), 
                    what="", sep="\n")
rna.genes.N <- scan(paste0(DATA_FOLD, "/TCGA-LUAD_rna_genes_N.txt"), 
                    what="", sep="\n")
rna.patients.C <- scan(paste0(DATA_FOLD, "/TCGA-LUAD_rna_patients_C.txt"), 
                       what="", sep="\n")
rna.patients.N <- scan(paste0(DATA_FOLD, "/TCGA-LUAD_rna_patients_N.txt"), 
                       what="", sep="\n")
clinical.data <- read.csv(paste0(DATA_FOLD, "/TCGA-LUAD_clinical_data.txt"),
                          header=FALSE)


# Clean DataFrame ------------------------------------------------

## Check if there are duplicated
check.duplicated.N <- duplicated(rna.patients.N)
cat(paste('Are there duplicates in rna.patients.N?', 
            any(check.duplicated.N)))

check.duplicated.C <- duplicated(rna.patients.C)
cat(paste('Are there duplicates in rna.patients.C?', 
            any(check.duplicated.C)))

## Select Cancer-Normal patients
rna.patients.CN <- intersect(rna.patients.C, rna.patients.N)

check.N.isin.CN <- sapply(rna.patients.N, 
                       function(elem) elem %in% rna.patients.CN)
check.C.isin.CN <- sapply(rna.patients.C, 
                          function(elem) elem %in% rna.patients.CN)

rna.expr.data.N.clean <- rna.expr.data.N[!check.duplicated.N & check.N.isin.CN]
rna.expr.data.C.clean <- rna.expr.data.C[!check.duplicated.C & check.C.isin.CN]

colnames(rna.expr.data.N.clean) <- lapply(colnames(rna.expr.data.N.clean), 
                                          function(elem) 
                                          gsub("[.]","-",substring(elem,1,12)))
colnames(rna.expr.data.C.clean) <- lapply(colnames(rna.expr.data.C.clean), 
                                          function(elem) 
                                          gsub("[.]","-",substring(elem,1,12)))

rna.expr.data.N.clean <- rna.expr.data.N.clean[rna.patients.CN]
rna.expr.data.C.clean <- rna.expr.data.C.clean[rna.patients.CN]


## Erase the genes with at least one zero values
row.sub.N <- apply(rna.expr.data.N.clean, 1, function(row) all(row !=0 ))
row.sub.C <- apply(rna.expr.data.C.clean, 1, function(row) all(row !=0 ))

rna.expr.data.N.clean <- rna.expr.data.N.clean[row.sub.C & row.sub.N,]
rna.expr.data.C.clean <- rna.expr.data.C.clean[row.sub.C & row.sub.N,]

## Save clean data sets
data.N <- rna.expr.data.N.clean
data.C <- rna.expr.data.C.clean

write.csv(data.N,
          paste0(DATA_FOLD,"data/dataN.csv"))
write.csv(data.C,
          paste0(DATA_FOLD,"data/dataC.csv"))

rm(list=setdiff(ls(), c("data.N", "data.C", "DATA_FOLD")))
save.image(file= paste0(DATA_FOLD,"data/dataCN.RData"))
