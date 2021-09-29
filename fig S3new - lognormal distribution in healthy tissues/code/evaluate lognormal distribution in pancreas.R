##### benchmark lognormal vs. least-squares deconvolution methods on CPA data

rm(list = ls())
library(NormqPCR)
library(scales)
library(ggplot2)
library(ggthemes)
library(e1071)
library(readxl)

#### data loading ----------------------------------------

# load data:
annot = as.data.frame(read_xlsx("../data/pancreas Q3 Norm.xlsx", sheet = "SegmentProperties"))
raw0 = read_xlsx("../data/pancreas Q3 Norm.xlsx", sheet = "TargetCountMatrix")
raw = as.matrix(raw0[, -1])
rownames(raw) = raw0$TargetName
rm(raw0)

# normalize:
norm = sweep(raw, 2, apply(raw, 2, quantile, 0.9), "/") *  mean(apply(raw, 2, quantile, 0.9))

# use geometric AOIs for a well-controlled comparison:
use = annot$SegmentLabel == "Geometric Segment" 
annot = annot[use, ]
norm = norm[, use]

# linear and logscale data:
lindat = (norm)
logdat = (log2(pmax(norm, 1)))


# calculate skew for all genes, in both logscale and linear scale:
skew.lin = apply(lindat, 2, skewness)
skew.log = apply(logdat, 2, skewness)
denslin = density(skew.lin[!is.na(skew.lin)])
denslog = density(skew.log[!is.na(skew.log)])

svg("../results/skewness in lin vs. log - pancreas.svg", width = 6, height = 4)
par(mar = c(4,4,.1,.1))
plot(denslin, col = 0, xlab = "Skewness of genes in\nmicroenvironment regions", ylab = "",
     cex.lab = 1.2, cex.axis = 0.7, main = "", xlim = c(0, max(denslin$x)), ylim = c(0, max(denslog$y)))
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
means.lin = apply(lindat[, !grepl("NegProbe", colnames(lindat))], 1, mean)
means.log = apply(logdat[, !grepl("NegProbe", colnames(logdat))], 1, mean)
sds.lin = apply(lindat[, !grepl("NegProbe", colnames(lindat))], 1, sd)
sds.log = apply(logdat[, !grepl("NegProbe", colnames(logdat))], 1, sd)

# summary stats:
range(sds.lin)
range(sds.lin)[2] / range(sds.lin)[1]
range(sds.log, na.rm = T)
range(sds.log, na.rm = T)[2] / range(sds.log, na.rm = T)[1]

png("../results/mean vs. sd - pancreas.png", height = 3, width = 6.6, units = "in", res = 300)
par(mar = c(4,4,.1,1))
par(mfrow = c(1,2))
plot(log2(means.lin), log2(sds.lin), log = "", xlab = "Mean", ylab = "SD", 
     col = alpha("dodgerblue4", 0.2), pch = 16, xaxt = "n", yaxt = "n")
axis(1, at = log2(10^seq(-3, 5, 1)), labels = 10^(seq(-3, 5, 1)), cex.axis = 0.65)
axis(2, at = log2(10^seq(-15, 15, 1)), labels = 10^(seq(-15, 15, 1)), cex.axis = 0.75)
legend("top", legend = "linear-scale data    ", bty = "n")
lines(lowess(log2(means.lin[sds.lin > 0]), log2(sds.lin[sds.lin > 0])), col = "orange", lwd = 2)
plot(means.log, sds.log, xlab = "Mean", ylab = "SD",
     col = alpha("dodgerblue4", 0.2), pch = 16, cex.axis = 0.65)
legend("top", legend = "log2-scale data    ", bty = "n")
lines(lowess(means.log[!is.na(sds.log)], sds.log[!is.na(sds.log)]), col = "orange", lwd = 2)
dev.off()

