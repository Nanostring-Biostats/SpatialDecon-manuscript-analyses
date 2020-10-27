
rm(list = ls())
library(logNormReg)
library(ComplexHeatmap)
library(pheatmap)
library(SpatialDecon)
library(scales)
library(circlize)
library(umap)

#### load grid data from earlier script:

load("191 grid decon results.RData")
aoicols = c("chartreuse2", "darkcyan")  
names(aoicols) = c("Tumor", "TME")


# infer polygon:
source("spaceplot utils.R")
bound = getBoundary(annot$x, annot$y, marg = 0.1)

#### model expression ~ cells -------------------------

beta = res$beta

#use.cells = rownames(res$beta)[rowSums(res$beta) > 25]
use.cells = rownames(beta)#[rowSums(res$beta) > 25]

ests = array(NA, dim = c(nrow(snr), length(use.cells) + 1, 2), dimnames = list(rownames(snr), c("b0", use.cells), c("Tumor", 'TME')))

for (atype in c("Tumor", "TME")) {
  # define data subsets:
  if (atype == "TME") {
    genemat = snr[, roiindices$TME]
    cellmat = beta[use.cells, roiindices$TME]
  }
  if (atype == "Tumor") {
    genemat = snr[, roiindices$Tumor[!is.na(roiindices$Tumor)]]
    cellmat = beta[use.cells, roiindices$TME[!is.na(roiindices$Tumor)]]
  }
  # run models:
  for (gene in dimnames(ests)[[1]]) {
    # run lognlm:
    y = pmax(genemat[gene, ], 1)
    fit = lognlm(y ~ t(cellmat),
                 lik = FALSE,
                 #start = c(1, init),
                 method = "L-BFGS-B", 
                 lower = rep(0, ncol(cellmat) + 1), 
                 upper = rep(Inf, ncol(cellmat) + 1),
                 opt = "optim",
                 control = list(maxit = 1000))
    ests[gene, , atype] = fit$coefficients                  
  }
}

save(ests, file = "gene vs cells nlm results.RData")

# now get yhats:
yhat = snr * NA
for (ind in roiindices$TME) {
  yhat[, ind] = ests[, c("b0", use.cells), "TME"] %*% c(1, beta[use.cells, ind])
}
for (ind in setdiff(roiindices$Tumor, NA)) {
  yhat[, ind] = ests[, c("b0", use.cells), "Tumor"] %*% c(1, beta[use.cells, ind])
}

# and resids:
resids = log2(pmax(snr, 1)) - log2(pmax(yhat, 1))

# sanity check:
plot(yhat[1:100, ], snr[1:100, ])


#### panel a: cartooning of the concept ---------------------------------

# cartoon of predicted vs. observed, colored by residuals

# function to color by residuals:
rcols = colorRampPalette(c("darkblue", "grey70", "darkred"))(101)
rthresh = 1.5
rbreaks = seq(-rthresh, rthresh, length.out = 102)
colorresids = function(x) {
  x2 = pmax(pmin(x, rthresh), -rthresh)
  x3 = (x2 + rthresh) / (2 * rthresh)
  col = rcols[round(x3 * 100) + 1]
  return(col)
}

# 
gene = "COL1A1"
svg("results/cartoon of residuals.svg", width = 6, height = 6)
par(mar = c(5,5,2,1))
plot(yhat[gene, roiindices$TME],
     snr[gene, roiindices$TME],
     col = alpha(colorresids(resids[gene, roiindices$TME]), 0.9),
     pch = 16, cex = 2,
     xlab = "Expression predicted from cell abundances",
     ylab = "Observed expression", cex.lab = 1.5)
abline(0,1)
dev.off()


# get cluster order:
o100 = order(apply(resids[, roiindices$TME], 1, sd), decreasing = T)[1:100]
p0 = pheatmap(t(scale(t(log2(snr[, roiindices$TME]))))[o100, 1:100],
         #colorRampPalette(c("blue", "white", "orange"))(100),
         breaks = seq(-3,3,length.out = 101),
         show_rownames = F, show_colnames = F)


# snr:
svg("results/fig 6a - cartoon of expression.svg", width = 8, height = 6)
pheatmap(t(scale(t(log2(snr[, roiindices$TME]))))[o100, 1:100][p0$tree_row$order, p0$tree_col$order],
         #colorRampPalette(c("blue", "white", "orange"))(100),
         breaks = seq(-3,3,length.out = 101),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F, legend = F,
         border_color = NA)
dev.off()

# yhat:
svg("results/fig 6a - cartoon of yhat.svg", width = 8, height = 6)
pheatmap(t(scale(t(log2(pmax(yhat, 1)[, roiindices$TME]))))[o100, 1:100][p0$tree_row$order, p0$tree_col$order],
         #colorRampPalette(c("blue", "white", "orange"))(100),
         breaks = seq(-3,3,length.out = 101),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F, legend = F,
         border_color = NA)
dev.off()


# resids:
svg("results/fig 6a - cartoon of resids.svg", width = 8, height = 6)
pheatmap(resids[, roiindices$TME][o100, 1:100][p0$tree_row$order, p0$tree_col$order],
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
         breaks = seq(-2,2,length.out = 101),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F, legend = F,
         border_color = NA)
dev.off()

pheatmap(resids[, roiindices$TME],
         #col = rcols, 
         col = colorRampPalette(c("darkred", "white", "darkblue"))(100),
         breaks = rbreaks, #seq(-,3,length.out = 101),
         show_rownames = F, show_colnames = F)

# beta:
p1 = pheatmap((pmin(res$cell.counts[, roiindices$TME], 100)),
              col = viridis_pal(option = "B")(100),
              show_colnames = F, show_rownames = F,
              cluster_rows = T, cluster_cols = T,
              legend = F)
dev.off()
svg("results/fig 6a - cartoon beta heatmap.svg", width = 1.5, height = .5)
pheatmap((pmin(res$cell.counts[, roiindices$TME], 100))[p1$tree_row$order, p1$tree_col$order],
         col = viridis_pal(option = "B")(100),
         show_colnames = F, show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         legend = F)
dev.off()


barplot(snr[,1])

#### quantify model goodness-of-fit within TME: -------------------------

cors = c()
for (gene in rownames(yhat)) {
  # TME results:
  cors[gene] = cor(pmax(yhat[gene, roiindices$TME], 1), pmax(snr[gene, roiindices$TME], 1))
}

# pick exemplar genes:
plot(cors, apply(resids[, roiindices$TME], 1, sd), col = 0)
text(cors, apply(resids[, roiindices$TME], 1, sd), names(cors), cex = 0.5)
showgenes = c("MT1M", "CCL19", "ARG1", "PDCD1")



#### fig 6b: plot cor vs resid sd --------------------------
svg("results/genes - cor w cells vs. resid SD.svg")
par(mar = c(5,5,1,1))
plot(cors, apply(resids[, roiindices$TME], 1, sd), 
     xlim = c(0.2, 1),
     col = alpha("#000066", 0.2), pch = 16, cex = 1.5,
     xlab = "Correlation between predicted and observed expression",
     ylab = "SD of residuals from predicted expression",
     cex.lab = 1.5)
xbound = 0.7
ybound = 0.32
abline(v = xbound, lty = 2)
abline(h = ybound, lty = 2)
text(0.92, 0.025, "Redundant with\ncell abundance")
text(0.92, 1.35, "Confounded by\ncell abundance")
text(0.275, 0.025, "Low variability")
text(0.275, 1.35, "High variability,\nindependent of\ncell abundance")
text(cors[showgenes], apply(resids[, roiindices$TME], 1, sd)[showgenes], showgenes,
     col = "red", cex = 1.25)
dev.off()




#### fig 6c: plot exemplar genes -----------------------
svg("results/resids - example genes.svg")
par(mar = c(2,2,1,1))
par(mfrow = c(2,2))
yrange = 100
for (gene in showgenes) {
  plot(yhat[gene, roiindices$TME] / median(yhat[gene, roiindices$TME]),
       snr[gene, roiindices$TME] / median(snr[gene, roiindices$TME]),
       col = colorresids(resids[gene, roiindices$TME]),
       pch = 16,
       log = "xy", 
       ylim = c(.1, 10), xlim = c(0.2, 6),
       xlab = "", ylab = "") 
  #, xaxt = "n", yaxt = "n")
  #legend("top", legend = paste0(gene, "   "))
  text(1.2, 10, paste0(gene, "   "), cex = 1.5)
}
dev.off()




#### fig 6d-f: plot selected genes and their resids in space -----------------------------

finalgenes = c("CXCL13", "LYZ",  "CCL17")
use = annot$AOI.name == "TME"
svg("results/resids in space - selected genes.svg", height = 7, width = 4.5)
par(mfrow = c(3, 2))
for (gene in finalgenes) {
  
  # plot expected vs observed:
  par(mar = c(4, 4, 2, 0))
  plot(snr[gene, use] ~ yhat[gene, use],
       pch = 16, cex = 1.5, cex.lab = 1,
       col = colorresids(resids[gene, use]),
       xlab = paste0("Fitted ", gene),
       ylab = paste0("Observed ", gene))
  abline(0,1)
  
  # resid spaceplot:
  tempcols = c()
  for (i in 1:ncol(resids)) {
    tempcols[i] = rcols[1 + sum(rbreaks <= resids[gene, i] )]
  }
  names(tempcols) = colnames(resids)
  par(mar = c(4,2,2,2))
  plot(x = annot[roiindices$TME, "x"],
       y = annot[roiindices$TME, "y"],
       cex = snr.sub[gene, roiindices$TME] / max(snr.sub[gene, roiindices$TME]) * 5,
       pch = 16,
       #main = gene,
       col = alpha(tempcols[roiindices$TME], 0.7),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n",  bty = "n",
       xlim = c(0, 9000), ylim = c(0, 9000))
  polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
  text(4500, 8800, gene, cex = 0.85)
}
dev.off()



#### co-expression analysis -----------------------------------------------

# take a gene, identify to 20 most correlated genes with it, and show cor heatmaps in snr and resid space:

cormat = cor(log2(t(snr[, roiindices$TME])))
cormat.r = cor((t(resids[, roiindices$TME])))

o = hclust(dist(cormat.r))


# identify tightly-linked clusters:
highcor = rownames(cormat.r)[apply(cormat.r > 0.7, 1 , sum) > 2]


# correlation heatmaps - snr and resids:
showgenes = highcor  #c("HLA-DRA", "HLA-DRB", "CD74", "HLA-DQB1", "HLA-DPA1")
ha = rowAnnotation(gene = anno_mark(at = match(showgenes, rownames(cormat[o$order, o$order])), labels = showgenes))

#### fig 6g: correlation matrices in normalized and resid space  --------------------------------
pdf("results/correlation heatmaps - all genes.pdf", width = 12)
p1 = Heatmap(cormat[o$order, o$order],
             cluster_rows = F, cluster_columns = F,
             col = colorRamp2(breaks = seq(-1, 1, length.out = 100),
                              colors = colorRampPalette(c("darkblue", "white", "darkred"))(100)),
             show_row_names = F, show_column_names = F,
             show_heatmap_legend = FALSE)
p2 = Heatmap(cormat.r[o$order, o$order],
             cluster_rows = F, cluster_columns = F,
             col = colorRamp2(breaks = seq(-1, 1, length.out = 100),
                              colors = colorRampPalette(c("darkblue", "white", "darkred"))(100)),
             show_row_names = F, show_column_names = F,
             show_heatmap_legend = TRUE, 
             name = "Correlation",
             right_annotation = ha)
print(p1 + p2)
dev.off()


# cuttree to find co-expressed clusters:
cut = cutree(o, h = 4.25)
table(table(cut))
cutids = names(table(cut))[table(cut) > 5]

pheatmap(cormat.r[names(cut)[is.element(cut, cutids)], names(cut)[is.element(cut, cutids)]])

if (FALSE) {
  for (cid in cutids) {
    tempgenes = names(cut)[cut == cid]
    #frame(); legend("center", legend = tempgenes)
    pheatmap(sweep(ests[tempgenes,,"TME"], 2, rowMeans(rbind(1, res$beta.granular[, use])), "*"),
             cluster_cols = F)
  }
  
  for (cid in cutids) {
    tempgenes = names(cut)[cut == cid]
    #frame(); legend("center", legend = tempgenes)
    if(length(intersect(tempgenes, rownames(training.matrices$monaco))) > 1) {
      pheatmap(training.matrices$monaco[intersect(tempgenes, rownames(training.matrices$monaco)), , drop = F],
               cluster_cols = F)
    }
  }
}

#pdf("results/resid clusters in space2.pdf", width = 12)
#for (cid in cutids) {
for (cid in c(240, 6)) {
  
  genemat = snr.sub[, roiindices$TME]
  #genemat = snr.sub
  tempgenes = names(cut)[cut == cid]
  tempmat = genemat[tempgenes, , drop = F]
  tempmat = tempmat[rowSums(tempmat) > 1, , drop = F]
  tempmat = sweep(tempmat, 1, apply(tempmat, 1, mean), "/")
  
  
  svg(paste0("results/resid florets - ", cid, ".svg"))
  #graphics::layout(mat = matrix(1:2, 1), widths = c(3, 8))
  
  #par(mar = c(0,0,0,0))
  #frame()
  #legend("center", legend = tempgenes, cex = 1)
  
  
  # florets:
  par(mar = c(0,0,0,0))
  #for (i in roiindices$TME) {
  florets(#x = annot[i, "x"],
    #y = annot[i, "y"],
    #b = tempmat[, i], 
    x = annot[roiindices$TME, "x"],
    y = annot[roiindices$TME, "y"],
    b = tempmat, 
    border = NA,
    #col = rep(alpha(aoicols["TME"], 0.5), nrow(tempmat)),
    col = matrix(colorresids(resids[rownames(tempmat), colnames(tempmat)]), nrow = nrow(tempmat), ncol = ncol(tempmat)),
    #col = "red",
    cex = 3,
    xlim = c(0, 9000), ylim = c(0, 9000)) 
  polygon(bound$x, bound$y, col = rgb(0,0,0,0.07), border = NA)
  #text(4500, 8800, 
  #     paste0(tempgenes, collapse = ", "), 
  #     cex = 0.85)
  dev.off()
  
  svg(paste0("results/resid florets legend - ", cid, ".svg"), width = 2.5, height = 2.5)
  par(mar = c(0,0,0,0))
  plot(0,0,col=0,xlim = c(-2,2), ylim = c(-2,2),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  angles = seq(0, 2 * pi, length.out = length(tempgenes) + 1)
  for (j in 1:length(tempgenes)) {
    tempangles = seq(angles[j], angles[j + 1], length.out = 20)
    xt = cos(tempangles)
    yt = sin(tempangles)
    polygon(c(0, xt), c(0, yt), col = NA, border = "black", lwd = 0.5)
    text(median(xt) * 1.2, median(yt) * 1.2, tempgenes[j])
  }
  dev.off()
}



### how much of expression is attributed to each cell type?
tempgenes = c( "ACP5","APOC1","ATP6V0D2", "LIPA","TREM2","CYP27A1" )
sum.yhat = c()
for (gene in tempgenes) {
  sum.yhat = rbind(sum.yhat, rowMeans(sweep(rbind(1, res$beta.granular[, use]), 1, ests[gene, , "TME"], "*")))
}
rownames(sum.yhat) = tempgenes
round(sweep(sum.yhat, 1, rowSums(sum.yhat), "/"), 2)


tempgenes = c( "CD74","HLA-DPA1", "HLA-DQA1/2", "HLA-DQB1", "HLA-DRA", "HLA-DRB" )
sum.yhat = c()
for (gene in tempgenes) {
  sum.yhat = rbind(sum.yhat, rowMeans(sweep(rbind(1, res$beta.granular[, use]), 1, ests[gene, , "TME"], "*")))
}
rownames(sum.yhat) = tempgenes
round(sweep(sum.yhat, 1, rowSums(sum.yhat), "/"), 2)


