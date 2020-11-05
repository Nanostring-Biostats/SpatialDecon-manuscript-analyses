# Step 2: continue converging the log normal models from Step 1
library(logNormReg)
library(MCPcounter)
library(here)

# retrive the scores for different cells from two files under the data/scores folder
scores <- read.csv(here("fig 3 - gene rates in tumors", "data", "scores", "40425_2017_215_MOESM4_ESM.csv"))
rownames(scores) <- scores[,1]
## exponentiate the scores 
scores <- 2^(scores[,-1]) 

load(here("fig 3 - gene rates in tumors", "data", "scores", "MCPcounter-scores.rdata")) 
## exponentiate the scores
scoresls <- 2^scoresls

## specify the cancer to continue converging the log normal model
cancer <- "ACC" #
load(here("fig 3 - gene rates in tumors", "data", paste0(cancer, ".RData")))
dat <- get(paste0(cancer,".dat.subset"))$e
## standarize it to be mean of 100
dat <- dat[, -which(colSums(dat)==0)]
datnew <- apply(dat, 2, function(x) x*100/mean(x))
load(here("fig 3 - gene rates in tumors", "output", paste0(cancer,"_res.rdata")))
## show number of unconverged genes
table(res$lgfit.convergence)
## extract the unconverged genes
unvgeneList <- res$gene[which(res$lgfit.convergence!=0)]
if(file.exists(here("output", paste0(cancer,"_res.rdata")))){
  load(here("output", paste0(cancer,"_res.rdata")))
  table(res$lgfit.convergence)
  unvgeneList <- res$gene[which(res$lgfit.convergence!=0)]
}

## refitting the log normal models on these unconverged genes
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
## re-examine whether there are genes remained unconverged
table(res$lgfit.convergence)

# update results for these unconverged genes
if(length(which(res$lgfit.convergence!=0))>0){
  res <- res[-which(res$lgfit.convergence!=0), ]
}

print(cancer)
save(res, file = here("output", paste0(cancer,"_res.rdata")))
