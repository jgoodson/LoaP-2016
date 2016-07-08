library('MCMC.qpcr')

#Set the working directory where files are
setwd("~/Dropbox/Goodson et al PepsA/LoaP/2016 Goodson LoaP/data/qPCR/YFP All Reporter Check/")

#Import csv files containing Cq and efficiency values
#Data file should contain both Cq values for all amplicons and also experimental variables as columns
#Data file should contain independent samples (technical or biological) as rows
#Technical replicates should have all experimental variables as well as sample ID identical
#Biological replicates should have all experimental variables identical with different sample IDs
#Required columns: "sample" (sample ID), at least one experimental variable, Cq columns for each amplicon
data = read.csv('YFP reporter Cqs.csv', header=TRUE, skip=1)
eff = read.csv('YFP reporter eff.csv', header=TRUE)

#Import cq data. Specify sample data as condcols and amplicons as genecols 
#Relevel all experimental variables so the reference level is correct
qs = cq2counts(data=data, effic=eff, condcols=c(1:5), genecols=c(6:11))
qs$Xylose = relevel(qs$Xylose, ref='Uninduced')
qs$Genotype = relevel(qs$Genotype, ref = 'P-L-Y')
soft_norm=mcmc.qpcr(data=qs, fixed="Genotype+Xylose+Genotype:Xylose", pr=T, pl=T, controls=c("rpoB", "gyrB", "dnaG"), normalize=TRUE, nitt=100000, burnin=10000)
sum_sn = HPDsummary(model=soft_norm,data=qs,relative=TRUE, genes=c('loaP', 'YFP'), xgroup="Xylose", x.order=c("Uninduced", "Induced"))
sum_sn$ggPlot+scale_fill_manual(values=c("#56B4E9", "#F6BF00"))
 

trellisByGene <- function (model, data, genes, xFactor, groupFactor, subplots, nrow = 1, lineWidth = 0.4, 
          whiskerWidth = 0.2, pointSize = 2.5, 
          ylab = "log(abundance)", legendPos = "bottom", posDodge = 0.3) 
{
  modelSummary = HPDsummary(model=model,data=data, xgroup=groupFactor, summ.plot=F)
  modelSummary$summary$Expt <- as.vector(subplots[modelSummary$summary$Genotype])
  lower = NULL
  upper = NULL
  pd = position_dodge(posDodge)
  gpl = ggplot(subset(modelSummary$summary, gene %in% genes), aes_string(x = xFactor, 
                                                group = groupFactor, colour = xFactor, y = "mean")) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), lwd = lineWidth, 
                  width = whiskerWidth, position = pd) +
    geom_point(position = pd, size = pointSize, aes_string(shape="Xylose"), stroke=1) + theme_bw() + 
    facet_wrap(~Expt, nrow = nrow, scale = "free_x") + 
    theme(axis.text.x=element_text(angle=90))+
    theme(axis.text.y=element_text(angle=90))+
    theme(strip.text.x = element_text(face = "italic", size = 12, 
                                      family = "serif")) + ylab(ylab) + theme(legend.position = legendPos)
  return(gpl)
}
subplots = c(1,1,2,2,3,3,3)
names(subplots)<-c('P-L-Y', 'P-LNT-Y', 'P-L-TT-Y', 'P-LNT-TT-Y', 'P-TL-Y', 'P-TLNT-Y', 'P-TL-TT-Y')
trel = trellisByGene(soft_norm, qs, c("YFP") ,xFactor="Genotype",groupFactor="Xylose", subplots=subplots, pointSize=3.5, whiskerWidth = 0.4, lineWidth = 0.8)
trel+scale_shape_manual(values=c(19,1))

