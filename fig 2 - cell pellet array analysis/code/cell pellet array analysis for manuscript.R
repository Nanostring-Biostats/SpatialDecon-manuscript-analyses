##### benchmark lognormal vs. least-squares deconvolution methods on CPA data

rm(list = ls())
library(logNormReg)
library(e1071)
library(NormqPCR)
library(scales)
library(ggplot2)
library(ggthemes)
source("deconvolution functions - NNLS and vSVR and DWLS.R")
source("DWLS_functions.R")

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

#write.csv(norm, file = "../results/HK normalized data.csv")

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


#### run deconvolution under loglinear, NNLS, v-SVR: ---------------------------------

if (FALSE) {
  res = list()
    res[["all"]] = list()
    tempY = norm[genesets[["all"]], ]
    tempX = X[genesets[["all"]], ]
    tempbg = bg[genesets[["all"]], ]
    res[["all"]][["NNLS"]] = deconUsingConstrainedLM(Y = tempY, X = tempX, bg = tempbg)
    res[["all"]][["vSVR"]] = deconSVR(Y = tempY, X = tempX, bg = tempbg)
    res[["all"]][["lognorm"]] = deconLNR(Y = tempY, X = tempX, bg = tempbg)
    res[["all"]][["dwls"]] = deconUsingDWLS(Y = tempY, X = tempX, bg = tempbg)
    save(res, file = "competing decon results with dwls.RData")
}
load("competing decon results with dwls.RData")


#### compare decon performance:

## plots:

# assemble data frame for ggplot:
est = obs = meth = gset = c()
par(mfrow = c(length(res[[1]]), length(res)))
for (mname in names(res[[1]])) {
  #for (gsname in names(res)) {
  for (gsname in c("all")) {
    tempbeta = res[[gsname]][[mname]]$beta
    est = c(est, tempbeta[1, ] / colSums(tempbeta))
    obs = c(obs, annot$mixprop)
    meth = c(meth, rep(mname, ncol(tempbeta)))
    gset = c(gset, rep(gsname, ncol(tempbeta)))
  }
}

plotdf = data.frame(est = est, obs = obs, meth = meth, gset = gset, stringsAsFactors = F)

# rename methods for plotting:
plotdf$meth[plotdf$meth == "lognorm"] = "log-normal"
plotdf$meth[plotdf$meth == "vSVR"] = "v-SVR"
plotdf$meth[plotdf$meth == "dwls"] = "DWLS"
plotdf$meth = factor(plotdf$meth, levels = c("NNLS", "v-SVR", "DWLS", "log-normal"))

# plot true vs observed proportions:
g = ggplot(plotdf, aes(x = obs, y = est)) +
  #scale_color_manual(values = gsetcols) + 
  #geom_point(col = I(alpha("dodgerblue4", 0.5)), size = I(2)) +
  geom_point(col = I("darkblue"), size = 2, alpha = 0.5) +
  facet_wrap(~meth, nrow = 1) + 
  theme_few() + 
  labs(x = "True proportion HEK293T", y = "Estimated proportion HEK293T") + 
  theme(axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 13), 
        strip.text.y = element_text(size = 15), 
        legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, col = rgb(0,0,0,0.3))
svg("../results/fig 2a.svg", width = 10, height = 3)
print(g)
dev.off()

## summary stats:
maxdelta = mads = mses = cors = matrix(NA, length(res[[1]]), length(res),
              dimnames = list(names(res[[1]]), names(res)))
par(mfrow = c(length(res[[1]]), length(res)))
for (mname in names(res[[1]])) {
#  for (gsname in names(res)) {
  for (gsname in c("all")) {
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
save(lres, file = "decon results in leave one out study.RData")
#save(list = ls(), file = "timesnap.RData")
lres[["vSVR"]] = deconSVR(Y = y.loo, X = X, bg = bg.loo)
colnames(lres$vSVR$beta) = colnames(lres$NNLS$beta)
save(lres, file = "decon results in leave one out study - all 3.RData")
lres[["dwls"]] = deconUsingDWLS(Y = y.loo, X = X, bg = bg.loo)
save(lres, file = "decon results in leave one out study - all 3.RData")

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
#rownames(propdeltas) = rownames(betasse) = c("all", genesubset)


# how big was the delta for each method when that gene was removed?


# identify the gene with the biggest impact:
tempg = rownames(propdeltas)[order(abs(propdeltas[, "NNLS"]), decreasing = T)[1]]
print(c(lres$NNLS$prop[1], lres$NNLS$prop[tempg]))
tempg = rownames(propdeltas)[order(abs(propdeltas[, "vSVR"]), decreasing = T)[1]]
print(c(lres$vSVR$prop[1], lres$vSVR$prop[tempg]))
tempg = rownames(propdeltas)[order(abs(propdeltas[, "dwls"]), decreasing = T)[1]]
print(c(lres$lognorm$prop[1], lres$lognorm$prop[tempg]))
tempg = rownames(propdeltas)[order(abs(propdeltas[, "lognorm"]), decreasing = T)[1]]
print(c(lres$lognorm$prop[1], lres$lognorm$prop[tempg]))


# what are the properties of the worst vSVR fits?
hist(abs(propdeltas[, "vSVR"]))
bad = which(abs(propdeltas[, "vSVR"]) > 0.1)
order(abs(propdeltas[, "lognorm"]), decreasing = T)[1]
y.loo[names(bad),1]


# show influence with point size:
#svg("../results/fig 2b - influence of genes when removed.svg", width = 13, height = 4)
#par(mfrow = c(1, length(lres) + 1))
layout(mat = matrix(1:5, 1), width = c(3.5, 3.5, 3.5, 3.5, 2.5))
for (mname in c("NNLS", "vSVR", "dwls", "lognorm")) {
  if (mname == "NNLS") {main = "NNLS"}
  if (mname == "lognorm") {main = "log-normal"}
  if (mname == "vSVR") {main = "v-SVR"}
  if (mname == "dwls") {main = "DWLS"}
  plot(X, pch = 16, cex = (abs(propdeltas[-1, mname]) / max(abs(propdeltas[-1, ]))) * 20 + .3, 
       log = "xy", col = alpha("darkblue", 0.5), main = main,
       xlab = paste0("Expression in ", colnames(X)[1]),
       ylab = paste0("Expression in ", colnames(X)[2]),
       cex.lab = 1.5, cex.main = 1.5)
  abline(0,1,col = "grey50")
}

par(mar = c(5,.1,2,.1))
plot(c(0,0), xlim = c(0, 1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = 0, bty = "n")
points(rep(.4, 5), seq(0, .75, length.out = 5), 
       cex = seq(0, 1, length.out = 5) * 20 + .3,
       pch = 16, col = alpha("darkblue", 0.5))
text(rep(.7, 5), seq(0, .75, length.out = 5), 
     labels = seq(0, signif(max(abs(propdeltas)), 0.2), length.out = 5))
text(.5, .9, labels = "Change in estimated cell\nproportion when removed", cex = 1.3)
#dev.off()

#save(list = ls(), file = "analysis snapshot.RData")
## ggplot version of the above for uniformity:


hek = ccrf = propdelta = model = c()
for (name in colnames(propdeltas)) {
  hek = c(hek, X[, "HEK293T"])
  ccrf = c(ccrf, X[, "CCRF-CEM"])
  propdelta = c(propdelta, propdeltas[rownames(X), name])
  model = c(model, rep(name, nrow(X)))
}
df = data.frame(hek = hek, ccrf = ccrf, meth = model, propdelta = abs(propdelta) + 0.01)

# rename methods for plotting:
df$meth[df$meth == "lognorm"] = "log-normal"
df$meth[df$meth == "vSVR"] = "v-SVR"
df$meth[df$meth == "dwls"] = "DWLS"
df$meth = factor(df$meth, levels = c("NNLS", "v-SVR", "DWLS", "log-normal"))


g = ggplot(df, aes(x = hek, y = ccrf)) + 
  geom_point(col = I("darkblue"), alpha = 0.35, size = I((abs(propdelta) / max(abs(propdelta))) * 20 + .1)) +
  facet_wrap(~meth, nrow = 1) + 
  theme_few() + 
  labs(x = "Expression in HEK293T", y = "Expression in CCRF-CEM") + 
  theme(axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 13), 
        strip.text.y = element_text(size = 15), 
        legend.position = "bottom") +
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') 

svg("../results/fig 2b - influence of genes when removed.svg", width = 10, height = 3.1)
print(g)
dev.off()

svg("../results/fig 2c legend.svg", width = 10, height = 1.25)
par(mar = c(0,0,0,0))
plot(c(0,0), xlim = c(0, 1), ylim = c(0.2,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = 0, bty = "n")
points(.1 + seq(.2, .8, length.out = 5), rep(.8, 5),
       cex = seq(0, 1, length.out = 5) * 7 + .3,
       pch = 16, col = alpha("darkblue", 0.35))
points(.1 + seq(.2, .8, length.out = 5), rep(.8, 5),
       cex = seq(0, 1, length.out = 5) * 7 + .3,
       pch = 1, col = alpha("darkblue", 0.45))
text(.1 + seq(.2, .8, length.out = 5), rep(.55, 5), 
     labels = seq(0, signif(max(abs(propdeltas)), 0.2), length.out = 5))
#text(.5, .1, labels = "Change in estimated cell proportion when removed", cex = 1.3)
text(.1, .65, labels = "Change in estimated\ncell proportion\nwhen removed", cex = 1)
dev.off()






# how does influence vary with 
par(mfrow = c(1, length(lres)))
for (mname in colnames(propdeltas)) {
  plot(abs(propdeltas[-1, mname]) ~ y.loo[, 1], log = "xy", ylim = c(1e-9, max(abs(propdeltas))), main = mname,
       xlab = "Expression level", ylab = "Change in estimated proportion when gene is removed", cex.lab = 1.5)
}
for (mname in colnames(propdeltas)) {
  plot(abs(propdeltas[-1, mname]) ~ y.loo[, 1], log = "x", ylim = c(1e-9, max(abs(propdeltas))), main = mname,
       xlab = "Expression level", ylab = "Change in estimated proportion when gene is removed", cex.lab = 1.5)
}

