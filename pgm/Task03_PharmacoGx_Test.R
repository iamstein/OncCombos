#This code comes directly from ??PharmacoGx R source

library(Biobase)
library(PharmacoGx)
library(xtable)
data("GDSCsmall")
data("CCLEsmall")
commonGenes <- intersect(fNames(GDSCsmall, "rna"),
                         fNames(CCLEsmall,"rna"))
common <- intersectPSet(list('CCLE'=CCLEsmall,
                             'GDSC'=GDSCsmall),
                        intersectOn=c("cell.lines", "drugs"), strictIntersect=TRUE)


GDSC.auc <- summarizeSensitivityProfiles(
  pSet=common$GDSC,
  sensitivity.measure='auc_published', 
  summary.stat="median",
  verbose=FALSE)
CCLE.auc <- summarizeSensitivityProfiles(
  pSet=common$CCLE,
  sensitivity.measure='auc_published', 
  summary.stat="median",
  verbose=FALSE)

GDSC.ic50 <- summarizeSensitivityProfiles(
  pSet=common$GDSC, 
  sensitivity.measure='ic50_published', 
  summary.stat="median",
  verbose=FALSE)
CCLE.ic50 <- summarizeSensitivityProfiles(
  pSet=common$CCLE, 
  sensitivity.measure='ic50_published', 
  summary.stat="median",
  verbose=FALSE)

GDSCexpression <- summarizeMolecularProfiles(common$GDSC, 
                                             cellNames(common$GDSC),
                                             mDataType="rna",
                                             features=commonGenes,
                                             verbose=FALSE)
CCLEexpression <- summarizeMolecularProfiles(common$CCLE, 
                                             cellNames(common$CCLE),
                                             mDataType="rna",
                                             features=commonGenes,
                                             verbose=FALSE)
gg <- fNames(common[[1]], 'rna')
cc <- cellNames(common[[1]])

ge.cor <- sapply(cc, function (x, d1, d2) {
  return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                     use="pairwise.complete.obs"))
}, d1=exprs(GDSCexpression), d2=exprs(CCLEexpression))
ic50.cor <- sapply(cc, function (x, d1, d2) {
  return (stats::cor(d1[, x], d2[ , x], method="spearman",
                     use="pairwise.complete.obs"))
}, d1=GDSC.ic50, d2=CCLE.ic50)
auc.cor <- sapply(cc, function (x, d1, d2) {
  return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                     use="pairwise.complete.obs"))
}, d1=GDSC.auc, d2=CCLE.auc)

w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor,
                         conf.int=TRUE, exact=FALSE)
w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor,
                         conf.int=TRUE, exact=FALSE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
              w1$p.value, w2$p.value)
boxplot(list("GE"=ge.cor,
             "AUC"=auc.cor,
             "IC50"=ic50.cor),
        main="Concordance between cell lines",
        ylab=expression(R[s]),
        sub=ss,
        ylim=yylim,
        col="lightgrey",
        pch=20,
        border="black")





###################################################
### code chunk number 8: load-cmap
###################################################
library(PharmacoGx)
require(xtable)
data(CMAPsmall)
drug.perturbation <- drugPerturbationSig(CMAPsmall, 
                                         mDataType="rna",
                                         verbose=FALSE)
data(HDAC_genes)

res <- apply(drug.perturbation[,,c("tstat", "fdr")],
             2, function(x, HDAC){ 
               return(connectivityScore(x=x, 
                                        y=HDAC[,2,drop=FALSE], 
                                        method="fgsea", nperm=100))
             }, HDAC=HDAC_genes)

rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=TRUE),]
xtable(res, 
       caption='Connectivity Score results for HDAC inhibitor gene signature.')


###################################################
### code chunk number 9: biomarkers
###################################################

data(CCLEsmall)
features <- fNames(CCLEsmall, "rna")[
  which(featureInfo(CCLEsmall,
                    "rna")$Symbol == "NQO1")]
sig.rna <- drugSensitivitySig(pSet=CCLEsmall, 
                              mDataType="rna", 
                              drugs=c("17-AAG"), 
                              features=features, 
                              sensitivity.measure="auc_published", 
                              molecular.summary.stat="median", 
                              sensitivity.summary.stat="median",
                              verbose=FALSE)
sig.mut <- drugSensitivitySig(pSet=CCLEsmall, 
                              mDataType="mutation", 
                              drugs=c("PD-0325901"), 
                              features="BRAF", 
                              sensitivity.measure="auc_published", 
                              molecular.summary.stat="and", 
                              sensitivity.summary.stat="median",
                              verbose=FALSE)
sig <- rbind(sig.rna, sig.mut)
rownames(sig) <- c("17-AAG + NQO1","PD-0325901 + BRAF")
colnames(sig) <- dimnames(sig.mut)[[3]]
xtable(sig, caption='P Value of Gene-Drug Association')


###################################################
### code chunk number 10: sessionInfo
###################################################
toLatex(sessionInfo())