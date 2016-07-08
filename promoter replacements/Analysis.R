library('MCMC.qpcr')

setwd("~/Dropbox/Goodson et al PepsA/LoaP/2016 Goodson LoaP/data/qPCR/Pconst replacement noear/")

data = read.csv('Cqs.csv', header=TRUE)
eff = read.csv('Eff.csv', header=TRUE)

qs = cq2counts(data=data, effic=eff, genecols=c(4:10), condcols=c(1:3))
qs$strain = relevel(qs$strain, ref='99')
qs$xylose = relevel(qs$xylose, ref='m')

soft_norm = mcmc.qpcr(data=qs, fixed="strain+xylose+xylose:strain", controls=c("rpoB", "gyrB", "dnaG"), normalize=TRUE, nitt=15000, burnin=1500)

sn_sum = HPDsummary(model=soft_norm,data=qs, genes=c("loaP", "dfnA", "dfnG", "dfnM"))


trellisByGene(sn_sum,xFactor="xylose",groupFactor="strain",nrow=1)

trellisByGene2 <- function (model, data, genes, xFactor, groupFactor, subplots, nrow = 1, lineWidth = 0.4, 
                            whiskerWidth = 0.2, pointSize = 2.5, 
                            ylab = "log(abundance)", legendPos = "bottom", posDodge = 0.4) 
{
  modelSummary = HPDsummary(model=model,data=data, xgroup=groupFactor, summ.plot=F)
  modelSummary$summary$Expt <- as.vector(subplots[modelSummary$summary$strain])
  lower = NULL
  upper = NULL
  pd = position_dodge(posDodge)
  gpl = ggplot(subset(modelSummary$summary, gene %in% genes), aes_string(x = xFactor, 
                                                                         group = groupFactor, colour = "gene", y = "mean")) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), lwd = lineWidth, 
                  width = whiskerWidth, position = pd) +
    geom_point(position = pd, size = pointSize, aes_string(shape="xylose"), stroke=1) + theme_bw() + 
    #facet_wrap(~Expt, nrow = nrow, scale = "free_x") + 
    theme(axis.text.x=element_text(angle=90))+
    theme(axis.text.y=element_text(angle=90))+
    theme(strip.text.x = element_text(face = "italic", size = 12, 
                                      family = "serif")) + ylab(ylab) + theme(legend.position = legendPos)
  return(gpl)
}
subplots = c("Pconst", "Pconst-Leader", "PdfnA-Leader")

trel = trellisByGene2(soft_norm, qs, c('dfnA', 'dfnG', 'dfnM') ,xFactor="strain",groupFactor="xylose", subplots=subplots, pointSize=3.5, whiskerWidth = 0.3, lineWidth = 0.8)
trel+scale_shape_manual(values=c(1,19))+scale_color_manual(values=c("dfnA"="#48b1a7", "dfnG"="#c55d88", "dfnM"="#b69340"))
pdf("plot.pdf")