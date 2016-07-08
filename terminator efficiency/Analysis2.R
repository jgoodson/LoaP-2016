library('MCMC.qpcr')

setwd("~/Dropbox/Goodson et al PepsA/LoaP/2016 Goodson LoaP/data/qPCR/Both terminator readthrough/")

data = read.csv('Cqs2.csv', header=TRUE)
eff = read.csv('Eff.csv', header=TRUE)

qs = cq2counts(data=data[1:16,], effic=eff, genecols=c(4:10), condcols=c(1:3))
qs$Xylose = relevel(qs$Xylose, ref='Uninduced')
qs$sample = qs$Sample

ls = cq2log(data=data[1:16,], effic=eff, genecols=c(4:10), condcols=c(1:3))
ls$Xylose = relevel(ls$Xylose, ref='Uninduced')
ls$sample = ls$Sample

soft_norm = mcmc.qpcr(data=qs, fixed="Xylose", globalFixed=c("Batch"), controls=c("rpoB", "gyrB", "dnaG"), random="sample", 
                      normalize=T, geneSpecRes=F, 
                      pr=T, nitt=220000, burnin=20000)

sn_sum = HPDsummary(model=soft_norm, data=qs)

sn=soft_norm$Sol
colnames(sn)

termEff <- function(sol, gene, cond) {
  pA = sol[,paste("gene",gene,"A", sep='')]+sol[,paste("gene",gene,"A",":",cond, sep='')]+sol[,cond]
  mA = sol[,paste("gene",gene,"A", sep='')]
  
  pB = sol[,paste("gene",gene,"B", sep='')]+sol[,paste("gene",gene,"B",":",cond, sep='')]+sol[,cond]
  mB = sol[,paste("gene",gene,"B", sep='')]

  p_Eff = pA-pB
  p_Eff_mean = exp(1)**mean(p_Eff)*100
  p_Eff_ci = exp(1)**HPDinterval(p_Eff)*100
  m_Eff = mA-mB
  m_Eff_mean = exp(1)**mean(m_Eff)*100
  m_Eff_ci = exp(1)**HPDinterval(m_Eff)*100
  
  print(c(p_Eff_mean, p_Eff_ci))
  print(c(m_Eff_mean, m_Eff_ci))
  print(mcmc.pval(p_Eff, m_Eff, ptype='mcmc', sided=1), digits=8)
}

termEff(sn, "E", "XyloseInduced")

termEff(sn, "leader", "XyloseInduced")



