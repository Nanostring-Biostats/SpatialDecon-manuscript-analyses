##### figures to demonstrate log-scale nature of variability and importance of log-linear method

rm(list = ls())
library(e1071)
library(scales)

#### TCGA: demonstrate logscale nature of bio variability:

# load TCGA LUAD:
load("LUSC.RData")
# linear and logscale data:
lindat = LUSC.dat.subset$e
logdat = log2(LUSC.dat.subset$e)
rm(LUSC.dat.subset)


# calculate skew for all genes, in both logscale and linear scale:
skew.lin = apply(lindat, 2, skewness)
skew.log = apply(logdat, 2, skewness)
denslin = density(skew.lin[!is.na(skew.lin)])
denslog = density(skew.log[!is.na(skew.log)])

svg("skewness in lin vs. log.svg", width = 6, height = 4)
par(mar = c(4,4,.1,.1))
plot(denslin, col = 0, xlab = "Skewness of genes in TCGA LUAD", ylab = "",
     cex.lab = 1.2, cex.axis = 0.7, main = "", ylim = c(0, max(denslog$y)))
polygon(denslin, col = alpha("grey50", 0.5), border = NA)
polygon(denslog, col = alpha("orange", 0.5), border = NA)
legend("topright", fill = alpha(c("grey50", "orange"), 0.5), 
       legend = c("Linear-scale data", "Log-transformed data"), cex = 1)
dev.off()

# summary stats:
mean(skew.lin > 2, na.rm = T)
mean(abs(skew.log) > 2, na.rm = T)

hist(skew.lin, col = alpha("grey50", 0.5), border = NA, ylab = "Number of genes", xlab = "Skewness",
     breaks = seq(-5,25, length.out = 100), ylim = c(0, 4000))
hist(skew.log, add = TRUE, col = alpha("orange", 0.5), border = NA,
     breaks = seq(-5,25, length.out = 100))



# show mean-variance relationship on logscale, linear scale:
means.lin = apply(lindat, 2, mean)
means.log = apply(logdat, 2, mean)
sds.lin = apply(lindat, 2, sd)
sds.log = apply(logdat, 2, sd)

# summary stats:
range(sds.lin[means.lin > 0])
range(sds.lin[means.lin > 0])[2] / range(sds.lin[means.lin > 0])[1]
range(sds.log, na.rm = T)
range(sds.log, na.rm = T)[2] / range(sds.log, na.rm = T)[1]

svg("mean vs. sd.svg", height = 3, width = 6.6)
par(mar = c(4,4,.1,1))
par(mfrow = c(1,2))
plot(log2(means.lin), log2(sds.lin), log = "", xlab = "Mean", ylab = "SD", 
     col = alpha("dodgerblue4", 0.05), pch = 16, xaxt = "n", yaxt = "n")
axis(1, at = log2(10^seq(-3, 5, 2)), labels = 10^(seq(-3, 5, 2)), cex.axis = 0.65)
axis(2, at = log2(10^seq(-15, 15, 2)), labels = 10^(seq(-15, 15, 2)), cex.axis = 0.75)
legend("top", legend = "linear-scale data    ", bty = "n")
lines(lowess(log(means.lin[sds.lin > 0]), log(sds.lin[sds.lin > 0])), col = "orange", lwd = 2)
plot(means.log, sds.log, xlab = "Mean", ylab = "SD",
     col = alpha("dodgerblue4", 0.05), pch = 16, cex.axis = 0.65)
legend("top", legend = "log2-scale data    ", bty = "n")
lines(lowess(means.log[!is.na(sds.log)], sds.log[!is.na(sds.log)]), col = "orange", lwd = 2)
dev.off()


# histograms of CD274
svg("gene histograms.svg", height = 3.1, width = 6)
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

#### in a TCGA sample, show influence of every point in LM22
# (how do you do that without running one time per gene, leaving a single gene out?)
# (should prob just do that. also for svm.)






