library('MCMC.qpcr')

data = read.csv('loaP comp test 98+99+orig.csv', header=TRUE)
eff = read.csv('loaP comp test 98+99+orig eff.csv', header=TRUE)

qs = cq2counts(data=data, effic=eff, genecols=c(5:10), condcols=c(1:4))
qs$genotype = relevel(qs$genotype, ref='WT')
qs$xylose = relevel(qs$xylose, ref='N')
qs$culture = relevel(qs$culture, ref='WT')

cds = cq2log(
  data=data[1:36,],
  genecols=c(5:10), 
  condcols=c(1:4),
  effic=eff
)
cds$genotype = relevel(cds$genotype, ref='WT')
cds$xylose = relevel(cds$xylose, ref='N')
cds$culture = relevel(cds$culture, ref='WT')

naive = mcmc.qpcr(data=qs, fixed="culture", pr=T, pl=T, nitt=100000, thin=10, burnin=15000)
soft_norm = mcmc.qpcr(data=qs, fixed="culture", pr=T, pl=T, controls=c("rpoB", "gyrB"), normalize=TRUE, nitt=100000, thin=10, burnin=15000)
informed = mcmc.qpcr(data=qs, fixed="culture", pr=T, pl=T, controls=c("rpoB", "gyrB"), nitt=100000, thin=10, burnin=15000)
classic = mcmc.qpcr.classic(data=cds, fixed="culture", pr=T, pl=T, controls=c("rpoB", "gyrB"), nitt=100000, thin=10, burnin=15000)

sn_sum = HPDsummary(model=soft_norm,data=qs, genes=c("dfnA", "dfnG", "dfnM"))
sn_sum$ggPlot$layers <- sn_sum$ggPlot$layers[c(1,3)]
sn_sum$ggPlot+ylim(c(2.5, 16))+coord_flip()


doplot <- function (model, data, genes, xFactor, xFactorelem, groupFactor, lineWidth = 0.4, 
                           whiskerWidth = 0.2, pointSize = 2.5, 
                           ylab = "log(abundance)", legendPos = "bottom", posDodge = 0.3) 
{
  modelSummary = HPDsummary(model=model,data=data, xgroup=groupFactor, summ.plot=F)
  lower = NULL
  upper = NULL
  pd = position_dodge(posDodge)
  xFactor_ss=modelSummary$summary[modelSummary$summary[,xFactor] %in% xFactorelem,]
  gpl = ggplot(subset(xFactor_ss, gene %in% genes), aes_string(x = xFactor, 
                                                                         group = groupFactor, colour = groupFactor, y = "mean")) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), lwd = lineWidth, 
                  width = whiskerWidth, position = pd) +
    geom_point(position = pd, size = pointSize, aes_string(shape=xFactor), stroke=1) + theme_bw() + 
    #theme(axis.text.x=element_text(angle=90))+
    #theme(axis.text.y=element_text(angle=90))+
    coord_flip()+
    theme(strip.text.x = element_text(face = "italic", size = 12, 
                                      family = "serif")) + ylab(ylab) + theme(legend.position = legendPos)
  return(gpl)
}

trel=doplot(soft_norm, qs, c("dfnA", "dfnG", "dfnM"), xFactor = "culture", xFactorelem=c("dL", "WT"), groupFactor="gene", lineWidth = 0.8, whiskerWidth = 0.4)
trel+scale_shape_manual(values=c(1,19))

trel=trellisByGene(soft_norm, qs, c("dfnA", "dfnG", "dfnM"), xFactor = "culture", xFactorelem=c("dLcompNon", "dLcompInd"), groupFactor="gene", lineWidth = 0.8, whiskerWidth = 0.4)
trel+scale_shape_manual(values=c(19,1))



trellisByGene <- function (model, data, genes, xFactor, groupFactor, subplots, nrow = 1, lineWidth = 0.4, 
                           whiskerWidth = 0.2, pointSize = 2.5, 
                           ylab = "log(abundance)", legendPos = "bottom", posDodge = 0.3) 
{
  modelSummary = HPDsummary(model=model,data=data, xgroup=groupFactor, summ.plot=F)
  modelSummary$summary$Expt <- as.vector(subplots[modelSummary$summary$culture])
  print(as.vector(subplots[modelSummary$summary$culture]))
  print(modelSummary$summary$culture)
  lower = NULL
  upper = NULL
  pd = position_dodge(posDodge)
  print(modelSummary$summary)
  gpl = ggplot(subset(modelSummary$summary[!is.na(modelSummary$summary$Expt),], gene %in% genes), aes_string(x = xFactor, 
                                                                         group = groupFactor, 
                                                                         colour = "gene", y = "mean")) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), lwd = lineWidth, 
                  width = whiskerWidth, position = pd) +
    geom_point(position = pd, size = pointSize, aes_string(shape=xFactor), stroke=1) + theme_bw() + 
    facet_wrap(~Expt, nrow = nrow, scale = "free_x") + 
    theme(axis.text.x=element_text(angle=90))+
    theme(axis.text.y=element_text(angle=90))+
    theme(strip.text.x = element_text(face = "italic", size = 12, 
                                      family = "serif")) + ylab(ylab) + theme(legend.position = legendPos)
  return(gpl)
}
subplots = c(1,2,2,1)
names(subplots)<-c('dL', "dLcompInd", "dLcompNon", "WT")
trel = trellisByGene(soft_norm, qs, c("dfnA", "dfnG", "dfnM") ,xFactor="culture",groupFactor="culture", 
                     subplots=subplots, pointSize=3.5, whiskerWidth = 0.2, lineWidth = 0.8,
                     posDodge = 1)
trel+scale_shape_manual(values=c(1, 19, 1, 19), guide=F)


s0=HPDsummary(model=naive,data=qs,relative=TRUE, summ.plot=FALSE)
s1=HPDsummary(model=informed,data=qs,relative=TRUE, summ.plot=FALSE)
s2=HPDsummary(model=classic,data=cds,relative=TRUE, summ.plot=FALSE, genes=c('dfnA', 'dfnG', 'dfnM', 'loaP'))
s3=HPDsummary(model=soft_norm,data=qs,relative=TRUE, summ.plot=FALSE)

diagnostic.mcmc(
  model=naive,
  col="grey50",
  cex=0.8,
  main="Naive model diagnostics"
)

diagnostic.mcmc(
  model=informed,
  col="grey50",
  cex=0.8,
  main="Informed model diagnostics"
)

diagnostic.mcmc(
  model=soft_norm,
  col="grey50",
  cex=0.8,
  main="Soft-Normaliation model diagnostics"
)

diagnostic.mcmc(
  model=classic,
  col="grey50",
  cex=0.8,
  main="Classic model diagnostics"
)

st0=HPDsummary(model=naive,data=qs,relative=TRUE, genes=c('dfnA', 'dfnG', 'dfnM'))
st1=HPDsummary(model=informed,data=qs,relative=TRUE, genes=c('dfnA', 'dfnG', 'dfnM'))
st2=HPDsummary(model=classic,data=cds,relative=TRUE, genes=c('dfnA', 'dfnG', 'dfnM'))
st3=HPDsummary(model=soft_norm,data=cds,relative=TRUE, genes=c('dfnA', 'dfnG', 'dfnM'))
st0$ggPlot+labs(title='Naive Model - Fold-change relative to WT cells', x='Culture Type')+theme_grey()
st1$ggPlot+labs(title='Informed Model - Fold-change relative to WT cells', x='Culture Type')+theme_grey()
st2$ggPlot+labs(title='Classic Model - Fold-change relative to WT cells', x='Culture Type')+theme_grey()
st3$ggPlot+labs(title='Soft-Normalization Model - Fold-change relative to WT cells', x='Culture Type')+theme_grey()




#below is from previous WT/dLoaP analysis

plot.new()
HPDplot(
  model=soft_norm,
  factors="genotype.xylosedL.N",
  main="delta-loaP vs WT (comparison of models)",
  ylimits=c(-7.5,1),
)
HPDpoints(
  model=naive,
  factors="genotypedL",
  jitter=-0.15,
  col="orange2",
  pch=17
)
HPDpoints(
  model=informed,
  factors="genotypedL",
  jitter=0.15,
  col="blue2",
  pch=18
)
HPDpoints(
  model=classic,
  factors="genotypedL",
  jitter=0.30,
  col="green4",
  pch=16
)
legend("topright",c("Naive", "Soft-Norm", "Informed", "Classic"),lty=1,pch=c(17,1,18,16),bty="n", col=c("orange2", "black", "blue2", "green4"))


plot.new()
HPDplot(
  model=informed,
  factors="genotypedL",
  main="delta-loaP vs WT (comparison of models)",
  ylimits=c(-7.5,4),
)
HPDpoints(
  model=naive,
  factors="genotypedL",
  jitter=-0.15,
  col="orange2",
  pch=17
)
HPDpoints(
  model=soft_norm,
  factors="genotypedL",
  jitter=0.15,
  col="blue2",
  pch=18
)
HPDpoints(
  model=classic,
  factors="genotypedL",
  jitter=0.30,
  col="green4",
  pch=16
)
legend("bottomright",c("Naive", "Informed", "Soft-Norm", "Classic"),lty=1,pch=c(17,1,18,16),bty="n", col=c("orange2", "black", "blue2", "green4"))


HPDplot(
  model=informed,
  factors="plate2",
  main="plate effect (2v1)"
)  


