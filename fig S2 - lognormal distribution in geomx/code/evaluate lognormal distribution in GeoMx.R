##### benchmark lognormal vs. least-squares deconvolution methods on CPA data

rm(list = ls())
library(NormqPCR)
library(scales)
library(ggplot2)
library(ggthemes)
library(e1071)
library(here)

#### data loading ----------------------------------------
name = "ICP20th"
# load segment/AOI annotation:
annot = read.table(here("fig 5 - NSCLC tumor analysis/", paste0("data/", name, "_SegmentProperties.txt")), 
                   header = T, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "")
rownames(annot) = make.names(annot$Sample_ID)
head(annot)


# load normalized data:
negnorm = as.matrix(read.table(here("fig 5 - NSCLC tumor analysis/", paste0("data/", name, "_NegNorm_TargetCountMatrix.txt")),
                               header = T, row.names = 1))
negnormall = negnorm



# use just TME AOIs:
use = rownames(annot)[annot$AOI.name == "TME"]
negnorm = negnorm[, use]
annot = annot[use, ]

# linear and logscale data:
lindat = t(negnorm)
logdat = t(log2(pmax(negnorm, 1)))


# calculate skew for all genes, in both logscale and linear scale:
skew.lin = apply(lindat, 2, skewness)
skew.log = apply(logdat, 2, skewness)
denslin = density(skew.lin[!is.na(skew.lin)])
denslog = density(skew.log[!is.na(skew.log)])

svg("../results/skewness in lin vs. log - geomx.svg", width = 6, height = 4)
par(mar = c(4,4,.1,.1))
plot(denslin, col = 0, xlab = "Skewness of genes in\nmicroenvironment regions", ylab = "",
     cex.lab = 1.2, cex.axis = 0.7, main = "", ylim = c(0, max(denslog$y)))
polygon(denslin, col = alpha("grey50", 0.5), border = NA)
polygon(denslog, col = alpha("orange", 0.5), border = NA)
legend("topright", fill = alpha(c("grey50", "orange"), 0.5), 
       legend = c("Linear-scale data", "Log-transformed data"), cex = 1)
dev.off()

# summary stats:
mean(skew.lin > 2, na.rm = T)
mean(abs(skew.log) > 2, na.rm = T)
mean(skew.lin, na.rm = T)
mean(skew.log, na.rm = T)


# show mean-variance relationship on logscale, linear scale:
means.lin = apply(lindat[, !grepl("NegProbe", colnames(lindat))], 2, mean)
means.log = apply(logdat[, !grepl("NegProbe", colnames(logdat))], 2, mean)
sds.lin = apply(lindat[, !grepl("NegProbe", colnames(lindat))], 2, sd)
sds.log = apply(logdat[, !grepl("NegProbe", colnames(logdat))], 2, sd)

# summary stats:
range(sds.lin)
range(sds.lin)[2] / range(sds.lin)[1]
range(sds.log, na.rm = T)
range(sds.log, na.rm = T)[2] / range(sds.log, na.rm = T)[1]

svg("../results/mean vs. sd.svg", height = 3, width = 6.6)
par(mar = c(4,4,.1,1))
par(mfrow = c(1,2))
plot(log2(means.lin), log2(sds.lin), log = "", xlab = "Mean", ylab = "SD", 
     col = alpha("dodgerblue4", 0.2), pch = 16, xaxt = "n", yaxt = "n")
axis(1, at = log2(10^seq(-3, 5, 2)), labels = 10^(seq(-3, 5, 2)), cex.axis = 0.65)
axis(2, at = log2(10^seq(-15, 15, 2)), labels = 10^(seq(-15, 15, 2)), cex.axis = 0.75)
legend("top", legend = "linear-scale data    ", bty = "n")
lines(lowess(log2(means.lin[sds.lin > 0]), log2(sds.lin[sds.lin > 0])), col = "orange", lwd = 2)
plot(means.log, sds.log, xlab = "Mean", ylab = "SD",
     col = alpha("dodgerblue4", 0.2), pch = 16, cex.axis = 0.65)
legend("top", legend = "log2-scale data    ", bty = "n")
lines(lowess(means.log[!is.na(sds.log)], sds.log[!is.na(sds.log)]), col = "orange", lwd = 2)
dev.off()


# histograms of CD274
svg("../results/gene histograms.svg", height = 3.1, width = 6)
gene = "CD274"
par(mfrow = c(1, 2))
par(mar = c(4,1,1,.1))
hist(lindat[, gene], breaks = 30, col = alpha("grey50", 0.5), border = NA, main = "",
     xlab = paste0("Linear-scale ", gene, " expression"), ylab = "", 
     yaxt = "n", cex.axis = 0.65)
legend("center", legend = paste0("skewness = ", round(skewness(lindat[, gene]), 1)), bty = "n")
hist(logdat[, gene], breaks = 30, col = alpha("orange", 0.5), border = NA, main = "",
     xlab = paste0("Log-scale ", gene, " expression"), ylab = "", 
     yaxt = "n", cex.axis = 0.65)
legend("center", legend = paste0("skewness = ", round(skewness(logdat[, gene]), 1)), bty = "n")
dev.off()




### aside for reviewer response: show that using TME focused the problem:

# look at heteroscedasticity:
svg("../results/heteroscedasticity.svg")
istme = grepl("TME", colnames(negnormall))
par(mfrow = c(2,2))
plot(apply(negnormall, 1, mean), apply(negnormall, 1, sd), 
     xlab = "Mean linear-scale expression", ylab = "SD of linear-scale expression",
     main = "All regions", log = "xy")
plot(apply(log2(negnormall), 1, mean), apply(log2(negnormall), 1, sd), 
     xlab = "Mean log2-scale expression", ylab = "SD of log2-scale expression",
     main = "All regions") #, xlim = c(2, 5.5), ylim = c(0,2.5))
plot(apply(negnormall[, istme], 1, mean), apply(negnormall[, istme], 1, sd),
     xlab = "Mean linear-scale expression", ylab = "SD of linear-scale expression",
     main = "Microenvironment regions only", log = "xy") #, xlim = c(0,70), ylim = c(0,260))
plot(apply(log2(negnormall)[, istme], 1, mean), apply(log2(negnormall)[, istme], 1, sd),
     xlab = "Mean log2-scale expression", ylab = "SD of log2-scale expression",
     main = "Microenvironment regions only") #, xlim = c(2, 5.5), ylim = c(0,2.5))
dev.off()
range(apply(negnormall, 1, sd))
range(apply(negnormall[, istme], 1, sd))

svg("../results/skewness.svg", height = 4)
# look at skweness:
par(mfrow = c(1,2))
# all regions:
skew.lin = apply(negnormall, 1, skewness)
skew.log = apply(log2(pmax(negnormall, 1)), 1, skewness)
denslin = density(skew.lin[!is.na(skew.lin)])
denslog = density(skew.log[!is.na(skew.log)])
plot(denslin, main = "All regions",
     col = 0, xlab = "Skewness", ylab = "",
     cex.lab = 1.2, cex.axis = 0.7, xlim = c(-2, 10), ylim = c(0, max(denslog$y)))
polygon(denslin, col = alpha("grey50", 0.5), border = NA)
polygon(denslog, col = alpha("orange", 0.5), border = NA)
mean(skew.lin, na.rm = T)

# TME only:
skew.lin = apply(negnorm, 1, skewness)
skew.log = apply(log2(pmax(negnorm, 1)), 1, skewness)
denslin = density(skew.lin[!is.na(skew.lin)])
denslog = density(skew.log[!is.na(skew.log)])
plot(denslin,  main = "Microenvironment regions only", 
     col = 0, xlab = "Skewness", ylab = "",
     cex.lab = 1.2, cex.axis = 0.7 , xlim = c(-2, 10))
polygon(denslin, col = alpha("grey50", 0.5), border = NA)
polygon(denslog, col = alpha("orange", 0.5), border = NA)
legend("topright", fill = alpha(c("grey50", "orange"), 0.5), 
       legend = c("Linear-scale data", "Log-transformed data"), cex = .75)
dev.off()
mean(skew.lin, na.rm = T)



