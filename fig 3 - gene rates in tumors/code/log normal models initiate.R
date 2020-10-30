library(logNormReg)
library(MCPcounter)
library(here)

## retrieve commonly shared genes 
filesList <- list.files(here("fig 3 - gene rates in tumors", "data"), ".rdata")
uniqueGene <- c()

for(i in seq_along(filesList)) { 
  fileName <- strsplit(filesList[i], "_res.rdata")[[1]]
  load(here("fig 3 - gene rates in tumors", "data", filesList[i]))
  if(length(uniqueGene)==0){
    uniqueGene <- c(uniqueGene, as.character(res$gene))
  }
  uniqueGene <- intersect(uniqueGene, as.character(res$gene))
}


## combine expression data from TCGA cancer types 
### assign each saved data into res
for(i in seq_along(filesList)) { 
  fileName <- strsplit(filesList[i], "_res.rdata")[[1]]
  load(here("output", outDir, filesList[i]))
  rownames(res) <- res$gene
  newres <- res[uniqueGene,c(2,3)]
  colnames(newres) <- c("gene", fileName)
  newres$gene <- as.character(newres$gene)
  assign(paste("res", i, sep = ""), newres)
}


### combine each res into res_combined
res_combined <- res1
for (i in setdiff(seq_along(filesList), 1)){
  res_combined <- inner_join(res_combined, get(paste("res", i, sep = "")), by = "gene") 
}

rownames(res_combined) <- res_combined[,1]
res_combined <- res_combined[,-1]




scores <- read.csv(here("fig 3 - gene rates in tumors", "data", "scores", "40425_2017_215_MOESM4_ESM.csv"))
rownames(scores) <- scores[,1]
scores <- 2^(scores[,-1]) # exponentiate the scores 

load(here("fig 3 - gene rates in tumors", "data", "scores", "MCPcounter-scores.rdata")) #scoresls
scoresls <- 2^scoresls

cancer <- "ACC" #
load(here("fig 3 - gene rates in tumors", "data", paste0(cancer, ".RData")))
dat <- get(paste0(cancer,".dat.subset"))$e
# standarize it to be mean of 100
dat <- dat[, -which(colSums(dat)==0)]
datnew <- apply(dat, 2, function(x) x*100/mean(x))

if(file.exists(here("output", paste0(cancer,"_res.rdata")))){
  load(here("output", paste0(cancer,"_res.rdata")))
  table(res$lgfit.convergence)
  unvgeneList <- res$gene[which(res$lgfit.convergence!=0)]
}

for (i in seq_len(length(unvgeneList))){
  gene <- unvgeneList[i]
  single.gene.expression <- datnew[, which(colnames(datnew)==gene)]
  single.gene.expression[which(single.gene.expression==0)] <-1
  X <- cbind(scores[rownames(datnew),], 
             scoresls[rownames(datnew), c("Fibroblasts", "Endothelial cells")])
  
  if(length(which(duplicated(colnames(datnew))))>0){
    datnew <- datnew[, -which(duplicated(colnames(datnew)))]
  }
  
  lgfit <- lognlm(single.gene.expression ~ as.matrix(X),
                  lik = FALSE,
                  method = "L-BFGS-B", 
                  lower = c(rep(0, ncol(X) + 1)), 
                  opt = "optim",
                  control=list(maxit=1000))
  
  if(lgfit$convergence == 52){
    X <- scale(X, center = FALSE)
    lgfit <- lognlm(single.gene.expression ~ as.matrix(X),
                    lik = FALSE,
                    method = "L-BFGS-B", 
                    lower = c(rep(0, ncol(X) + 1)), 
                    opt = "optim",
                    control=list(maxit=20000))    
  }
  
  res[which(res$gene==gene), ]  <- data.frame(i, gene, coef(lgfit)[1],
                                              lgfit$convergence,
                                              sbeta0 = coef(lgfit)[1]/sum(coef(lgfit)*c(1, colMeans(X))))
}
table(res$lgfit.convergence)

outDir2 = "set2"
if(length(which(res$lgfit.convergence!=0))>0){
  res <- res[-which(res$lgfit.convergence!=0), ]
}
print(cancer)
save(res, file = here("output", outDir2, paste0(cancer,"_res.rdata")))
