# Step 1: obtain coefficients from the log normal models 
library(logNormReg)
library(MCPcounter)
library(here)
library(dplyr)

# generate Fibroblasts and Endothelial cells scores from MCPcounter 
filesList <- list.files(here("fig 3 - gene rates in tumors", "data"), pattern = ".RData")
scoresls <- c()
count = 0

for(i in seq_along(filesList)) { 
  load(here("fig 3 - gene rates in tumors", "data", filesList[i]))
  cancer <- gsub(".RData", "", filesList[i])
  dat <- get(paste0(cancer,".dat.subset"))$e
  # standarize it to be mean of 100
  dat <- dat[, -which(colSums(dat)==0)]
  dat <- log2(dat + 1)
  
  fit <- MCPcounter.estimate(t(dat), featuresType="HUGO_symbols")
  if(nrow(fit)==10){
    scoresls <- rbind(scoresls, t(fit))
  }
  count = count + nrow(dat)
}

save(scoresls, file = here("fig 3 - gene rates in tumors", "data", "scores", "MCPcounter-scores.rdata"))

# retrive the scores for different cells from two files under the data/scores folder
scores <- read.csv(here("fig 3 - gene rates in tumors", "data", "scores", "40425_2017_215_MOESM4_ESM.csv"))
rownames(scores) <- scores[,1]
## exponentiate the scores
scores <- 2^(scores[,-1])  

load(here("fig 3 - gene rates in tumors", "data", "scores", "MCPcounter-scores.rdata")) #scoresls
## exponentiate the scores
scoresls <- 2^scoresls

## run the model for each cancer type 
for (cancer in gsub(".RData", "", filesList)){
  if(cancer %in% gsub("_res.rdata","", list.files(here("fig 3 - gene rates in tumors", "output")))){
    next
  }
  
  if (cancer == "UCEC"){
    next
  }
  
  load(here("fig 3 - gene rates in tumors", "data", paste0(cancer, ".RData")))
  dat <- get(paste0(cancer,".dat.subset"))$e
  # standarize it to be mean of 100
  dat <- dat[, -which(colSums(dat)==0)]
  datnew <- apply(dat, 2, function(x) x*100/mean(x))
  
  if(length(which(duplicated(colnames(datnew))))>0){
    datnew <- datnew[, -which(duplicated(colnames(datnew)))]
  }
  
  res <- c()
  # loop over all genes in the cancer data 
  for(i in seq_len(ncol(datnew))){
    gene <- colnames(datnew)[i]
    # extract the expression data for a give gene 
    single.gene.expression <- datnew[, which(colnames(datnew)==gene)]
    # shift the expression level by 1 if it is zero
    single.gene.expression[which(single.gene.expression==0)] <-1
    
    # find the common sample IDs between the score matrice and the TCGA dataset
    ind <- intersect(rownames(scores), rownames(datnew))
    single.gene.expression <- single.gene.expression[ind]
    
    # X is the combined cell scores from NanoString and MCPcounter
    X <- cbind(scores[ind,], 
               scoresls[ind, c("Fibroblasts", "Endothelial cells")])
    X <- scale(X, center = FALSE)
    
    # apply the log normal model 
    lgfit <- lognlm(single.gene.expression ~ as.matrix(X),
                    lik = FALSE,
                    method = "L-BFGS-B", 
                    lower = c(rep(0, ncol(X) + 1)), 
                    opt = "optim",
                    control=list(maxit=200)) 
    
    # collect the coeffecients, convergence status, and intercepts. 
    res <- rbind(res, data.frame(i, gene, coef(lgfit)[1], lgfit$convergence,
                                 sbeta0 = coef(lgfit)[1]/sum(coef(lgfit)*c(1, colMeans(X)))))
  }
  
  res_original <- res
  
  print(cancer)
  save(res, file = here("fig 3 - gene rates in tumors", "output", paste0(cancer,"_res.rdata")))
}

