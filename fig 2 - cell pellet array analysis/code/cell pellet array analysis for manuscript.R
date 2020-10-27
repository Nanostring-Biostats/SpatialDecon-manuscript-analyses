##### benchmark lognormal vs. least-squares deconvolution methods on CPA data

rm(list = ls())
library(logNormReg)
library(e1071)
library(NormqPCR)
library(scales)
library(ggplot2)
library(ggthemes)
source("deconvolution functions - NNLS and vSVR.R")

#### load CPA data: (HEK293 mixture) ---------------------

raw = t(as.matrix(read.csv("../data/HEK393-CCRF mixture raw data.csv", row.names = 1, header = T, stringsAsFactors = F)))
annot = read.csv("../data/HEK393-CCRF mixture AOI annotations.csv", row.names = 1, header = T, stringsAsFactors = F)
pnotes = read.csv("../data/HEK393-CCRF mixture probe notes.csv", stringsAsFactors = F, row.names = 1)


#### normalize raw data: -----------------------

# derive housekeepers:
topexpressing = rownames(raw)[order(rowMeans(raw), decreasing = T)[1:50]]
genorm = selectHKs(t(log2(raw)[topexpressing, ]), 
                   method = "geNorm", 
                   Symbols = topexpressing)
plot(genorm$variation)
hks = genorm$ranking[1:27]

# normalize:
annot$hk.factor = exp(colMeans(log(raw[hks, ])))
norm = sweep(raw, 2, annot$hk.factor, "/") * mean(annot$hk.factor)

write.csv(norm, file = "../results/HK normalized data.csv")

# estimate background in the normalized data: (separately for each probe pool):
probes.b = make.names(rownames(pnotes)[pnotes$Module == "bkp"])
probes.l = make.names(rownames(pnotes)[pnotes$Module == "lkp"])
bg = raw * 0
bg[probes.b, ] = sweep(bg[probes.b, ], 2, norm["NegProbe.bkp", ], "+")
bg[probes.l, ] = sweep(bg[probes.l, ], 2, norm["NegProbe.lkp", ], "+")


#### get expression profiles of pure cell lines: -------------------------
X = cbind(apply(norm[, annot$mixprop == 1], 1, median),
          apply(norm[, annot$mixprop == 0], 1, median))
colnames(X) = c("HEK293T", "CCRF-CEM")
X[, 2] = X[, 2] * exp(median(log(X[, 1]))) / exp(median(log(X[, 2])))

#### define gene subsets: --------------------------------------
plot(pmax(X, 1), log = "xy")
abline(0, 1)

genesets = list()
genesets[["all"]] = rownames(X)
genesets[["informative"]] = names(which(abs(log(X[, 1]) - log(X[, 2])) > 1.5))
genesets[["uninformative"]] = names(which(abs(log(X[, 1]) - log(X[, 2])) <= 1.5))
genesets[["nothigh"]] = names(which((apply(X, 1, max) < 1000)))

points(X[genesets[["informative"]], ], col = 2)
points(X[genesets[["uninformative"]], ], col = 3)
points(X[genesets[["nothigh"]], ], col = 4)

par(mfrow = c(1,4))
for (gsname in names(genesets)) {
  plot(pmax(X, 1), log = "xy", pch = 16, main = gsname,
       col = c(alpha("grey50", 0.3), alpha("orange", 0.5))[1 + is.element(rownames(X), genesets[[gsname]])])
  abline(0, 1)
}

# ggplot of gene sets:
hek = ccrf = gset = inset = c()
for (gname in names(genesets)) {
  hek = c(hek, X[, "HEK293T"])
  ccrf = c(ccrf, X[, "CCRF-CEM"])
  gset = c(gset, rep(gname, nrow(X)))
  inset = c(inset, c("out", "in")[1 + is.element(rownames(X), genesets[[gname]])])
}

setcols = c(alpha("grey60", 0.3), alpha("orange", 0.5))
names(setcols) = c("out", "in")

gsetcols = c("grey50", "#000066", "#006633",  "#FF0000", "#CC9900") #"#660033"(purpleish),
names(gsetcols) = c("out", "informative", "uninformative", "all", "nothigh")

plota = data.frame(hek = hek, ccrf = ccrf, gset = gset, inset = inset, stringsAsFactors = F)
plota$inset[plota$inset == "in"] = as.character(plota$gset)[plota$inset == "in"]
plota$gset[plota$gset == "all"] = "all genes"
plota$gset[plota$gset == "informative"] = "most informative genes"
plota$gset[plota$gset == "uninformative"] = "least informative genes"
plota$gset[plota$gset == "nothigh"] = "low-medium expressors"
plota$gset = factor(plota$gset, levels = c("most informative genes", "least informative genes", "all genes", "low-medium expressors"))
g = ggplot(plota, aes(x = hek, y = ccrf, col = inset)) + 
  geom_point(alpha = 0.5, size = 0.75) + 
  scale_color_manual(values = gsetcols) + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  facet_wrap(~gset, nrow = 1) +
  theme_few() +
  labs(x = "Expression in HEK293T", y = "Expression in CCRF-CEM") + 
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 13), 
        strip.text.y = element_text(size = 15)) 
  
svg("../results/fig 2a.svg", width = 9.8, height = 3.2)
print(g)
dev.off()

#### run deconvolution under loglinear, NNLS, v-SVR: ---------------------------------

res = list()
for (gsname in names(genesets)) {
  res[[gsname]] = list()
  tempY = norm[genesets[[gsname]], ]
  tempX = X[genesets[[gsname]], ]
  tempbg = bg[genesets[[gsname]], ]
  res[[gsname]][["NNLS"]] = deconUsingConstrainedLM(Y = tempY, X = tempX, bg = tempbg)
  res[[gsname]][["vSVR"]] = deconSVR(Y = tempY, X = tempX, bg = tempbg)
  res[[gsname]][["lognorm"]] = deconLNR(Y = tempY, X = tempX, bg = tempbg)
  save(res, file = "competing decon results.RData")
}


#### compare decon performance:


# assemble data frame for ggplot:
est = obs = meth = gset = c()
par(mfrow = c(length(res[[1]]), length(res)))
for (mname in names(res[[1]])) {
  for (gsname in names(res)) {
    tempbeta = res[[gsname]][[mname]]$beta
    est = c(est, tempbeta[1, ] / colSums(tempbeta))
    obs = c(obs, annot$mixprop)
    meth = c(meth, rep(mname, ncol(tempbeta)))
    gset = c(gset, rep(gsname, ncol(tempbeta)))
  }
}

plotdf = data.frame(est = est, obs = obs, meth = meth, gset = gset, stringsAsFactors = F)

# rename methods for plottig:
plotdf$meth[plotdf$meth == "lognorm"] = "log-normal"
plotdf$meth[plotdf$meth == "vSVR"] = "v-SVR"
plotdf$meth = factor(plotdf$meth, levels = c("NNLS", "v-SVR", "log-normal"))

# rename gene sets for plottig:
plotdf$gset0 = plotdf$gset
plotdf$gset[plotdf$gset == "all"] = "all genes"
plotdf$gset[plotdf$gset == "informative"] = "most informative genes"
plotdf$gset[plotdf$gset == "uninformative"] = "least informative genes"
plotdf$gset[plotdf$gset == "nothigh"] = "low-medium expressors"
plotdf$gset = factor(plotdf$gset, levels = c("most informative genes", "least informative genes", "all genes", "low-medium expressors"))

# plot true vs observed proportions:
g = ggplot(plotdf, aes(x = obs, y = est, col = gset0)) +
  scale_color_manual(values = gsetcols) + 
  #geom_point(col = I(alpha("dodgerblue4", 0.5)), size = I(2)) +
  geom_point(size = 2, alpha = 0.5) +
  facet_grid(meth ~ gset) + 
  theme_few() + 
  labs(x = "True proportion HEK293T", y = "Estimated proportion HEK293T") + 
  theme(axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text.x = element_text(size = 13), 
        strip.text.y = element_text(size = 15), 
        legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, col = rgb(0,0,0,0.3))
svg("../results/fig 2b.svg", width = 10, height = 7.5)
print(g)
dev.off()

## summary stats:
maxdelta = mads = mses = cors = matrix(NA, length(res[[1]]), length(res),
              dimnames = list(names(res[[1]]), names(res)))
par(mfrow = c(length(res[[1]]), length(res)))
for (mname in names(res[[1]])) {
  for (gsname in names(res)) {
    tempbeta = res[[gsname]][[mname]]$beta
    tempprop = tempbeta[1, ] / colSums(tempbeta)
    cors[mname, gsname] = cor(annot$mixprop, tempprop)
    mads[mname, gsname] = median(abs(annot$mixprop - tempprop))
    maxdelta[mname, gsname] = max(abs(annot$mixprop - tempprop))
    mses[mname, gsname] = sqrt(mean(annot$mixprop - tempprop)^2)
  }
}

round(cors,3)
round(mads,2)
round(maxdelta,2)
round(mses,3)
write.csv(round(mses, 3), file = "../results/rMSEs of decon fits.csv")


#### Assess influence of genes on decon:

# pick an even-ish mix:
oneaoi = which(annot$mixprop == 0.5)[1]

# create dataset with genes left out:
y.loo = norm[, oneaoi, drop = F]
for (i in 1:nrow(y.loo)) {
  y.loo = cbind(y.loo, replace(y.loo[, 1], i, NA))
}
colnames(y.loo) = c("all", rownames(y.loo))
bg.loo = sweep(y.loo * 0, 1, bg[, oneaoi], "+")

### now run decon, for all methods, with each gene removed once:
lres = list()
lres[["NNLS"]] = deconUsingConstrainedLM(Y = y.loo, X = X, bg = bg.loo)
lres[["lognorm"]] = deconLNR(Y = y.loo, X = X, bg = bg.loo)
lres[["vSVR"]] = deconSVR(Y = y.loo, X = X, bg = bg.loo)
colnames(lres$vSVR$beta) = colnames(lres$NNLS$beta)
for (name in names(lres)) {
  lres[[name]]$prop = lres[[name]]$beta[1, ] / colSums(lres[[name]]$beta)
}
# for each omitted gene, get change in cell proportions from complete gene set, and get squared delta in beta
propdeltas = betasse = c()
for (mname in names(lres)) {
  tempbeta = lres[[mname]]$beta
  tempprop = tempbeta[1, ] / colSums(tempbeta)
  
  propdeltas = cbind(propdeltas, tempprop - tempprop[1])
  betasse = cbind(betasse, sqrt(colSums(sweep(tempbeta, 1, tempbeta[, 1], "-")^2)))
}
colnames(propdeltas) = colnames(betasse) = names(lres)


# identify the gene with the biggest impact:
tempg = rownames(propdeltas)[order(abs(propdeltas[, "NNLS"]), decreasing = T)[1]]
print(c(lres$NNLS$prop[1], lres$NNLS$prop[tempg]))
tempg = rownames(propdeltas)[order(abs(propdeltas[, "vSVR"]), decreasing = T)[1]]
print(c(lres$vSVR$prop[1], lres$vSVR$prop[tempg]))
tempg = rownames(propdeltas)[order(abs(propdeltas[, "lognorm"]), decreasing = T)[1]]
print(c(lres$lognorm$prop[1], lres$lognorm$prop[tempg]))


# what are the properties of the worst vSVR fits?
hist(abs(propdeltas[, "vSVR"]))
bad = which(abs(propdeltas[, "vSVR"]) > 0.1)
order(abs(propdeltas[, "lognorm"]), decreasing = T)[1]
y.loo[names(bad),1]


# show influence with point size:
svg("../results/fig 2c - influence of genes when removed.svg", width = 13, height = 4)
layout(mat = matrix(1:4, 1), width = c(3.5, 3.5, 3.5, 2.5))
for (mname in c("NNLS", "vSVR", "lognorm")) {
  if (mname == "NNLS") {main = "NNLS"}
  if (mname == "lognorm") {main = "log-normal"}
  if (mname == "vSVR") {main = "v-SVR"}
  plot(X, pch = 16, cex = (abs(propdeltas[-1, mname]) / max(abs(propdeltas[-1, ]))) * 20 + .3, 
       log = "xy", col = alpha(gsetcols["all"], 0.5), main = main,
       xlab = paste0("Expression in ", colnames(X)[1]),
       ylab = paste0("Expression in ", colnames(X)[2]),
       cex.lab = 1.5, cex.main = 1.5)
  abline(0,1,col = "grey50")
}

par(mar = c(5,.1,2,.1))
plot(c(0,0), xlim = c(0, 1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = 0, bty = "n")
points(rep(.4, 5), seq(0, .75, length.out = 5), 
       cex = seq(0, 1, length.out = 5) * 20 + .3,
       pch = 16, col = alpha(gsetcols["all"], 0.5))
text(rep(.7, 5), seq(0, .75, length.out = 5), 
     labels = seq(0, signif(max(abs(propdeltas)), 0.2), length.out = 5))
text(.5, .9, labels = "Change in estimated cell\nproportion when removed", cex = 1.3)
dev.off()


