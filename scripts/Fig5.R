#!/usr/bin/env Rscript
############################
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
#}
library(data.table)
library(dplyr)
library(ggplot2)
library(EnvStats)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
args<-c(6000)
options(scipen=999)
W<-as.numeric(args[1])#window size

############################
############################
############################
############################
tag='phys_100000_0.0005'
readRDS('~/height_prediction/unweighted_prs/output/results.WHI.Rds')[[tag]]-> results.WHI
readRDS('~/height_prediction/unweighted_prs/output/results.JHS.Rds')[[tag]]-> results.JHS
readRDS('~/height_prediction/unweighted_prs/output/results.UKB_afr.Rds')[[tag]]-> results.UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/results.HRS_eur.Rds')[[tag]]-> results.HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/results.HRS_afr.Rds')[[tag]]-> results.HRS_afr
readRDS('~/height_prediction/unweighted_prs/output/B_WHI.Rds')[[tag]]->B_WHI
readRDS('~/height_prediction/unweighted_prs/output/B_JHS.Rds')[[tag]]->B_JHS
readRDS('~/height_prediction/unweighted_prs/output/B_UKB_afr.Rds')[[tag]]-> B_UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/B_HRS_eur.Rds')[[tag]]-> B_HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/B_HRS_afr.Rds')[[tag]]-> B_HRS_afr
readRDS('~/height_prediction/unweighted_prs/output/PGS3_UKB_afr.Rds')[[tag]]-> PGS3_UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/PGS3_WHI.Rds')[[tag]]-> PGS3_WHI
readRDS('~/height_prediction/unweighted_prs/output/PGS3_JHS.Rds')[[tag]]-> PGS3_JHS
readRDS('~/height_prediction/unweighted_prs/output/PGS3_HRS_eur.Rds')[[tag]]-> PGS3_HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/PGS3_HRS_afr.Rds')[[tag]]-> PGS3_HRS_afr

ALL2<-rbind(B_JHS[1:2,][, Dataset:='JHS_afr'], B_WHI[1:4,][, Dataset:='WHI_afr'], B_UKB_afr[1:4,][,Dataset:='UKB_afr'],B_HRS_afr[1:2,][, Dataset:='HRS_afr'], B_HRS_eur)
tmp<-1/c(var(results.JHS[[1]]$t), var(results.JHS[[2]]$t), var(results.WHI[[1]]$t), var(results.WHI[[2]]$t), var(results.WHI[[3]]$t), var(results.WHI[[4]]$t), var(results.UKB_afr[[1]]$t),var(results.UKB_afr[[2]]$t), var(results.UKB_afr[[3]]$t), var(results.UKB_afr[[4]]$t), var(results.HRS_afr[[1]]$t), var(results.HRS_afr[[2]]$t), var(results.HRS_eur$t))  #weighing lm by boostrap replicates.
my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
ALL2$Dataset<-factor(ALL2$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
cbind(ALL2, W=tmp)-> ALL2

shapes<-c(16,17,18,6,15)
plotA<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
geom_point(aes(shape=Dataset), size=5, fill="white", alpha=0.8) + stat_smooth(data=ALL2,method = "lm", mapping = aes(weight = sqrt(W)), col='darkgray', lty=2) +
geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
#geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +
scale_color_manual(values=my_colrs) +
scale_shape_manual(values=shapes)+
ylab(expression(paste("Partial-R"^"2"))) + xlab("European Ancestry Proportion") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 17), axis.title.x=element_text(size=17),axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=12),legend.position = c(0.13,0.86))
############################
############################
############################
############################
############################
############################
#read in frequencies
readRDS('~/height_prediction/epistasis/output/table_HRS_eur.Rds')->hrs_eur
readRDS('~/height_prediction/epistasis/output/table_HRS_afr.Rds')->hrs_afr
#calculate frequency difference per SNP

setkey(hrs_eur, CHR, POS)
setkey(hrs_afr, CHR, POS)
comb<-hrs_eur[hrs_afr, nomatch=0]
comb$Freq_effect_allele<-as.numeric(comb$Freq_effect_allele)
comb$i.Freq_effect_allele<-as.numeric(comb$i.Freq_effect_allele)
comb[,Freq_Diff:=Freq_effect_allele-i.Freq_effect_allele]
comb[, POS1:=POS]
comb[, POS2:=POS]
#betas_ukb<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')
final_plink<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
select(do.call(rbind, readRDS('~/height_prediction/gwas/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds')), CHR, POS, b, SE, p, N)-> prs_snps_hrs
coords_hrs<-prs_snps_hrs[,.(CHR,POS)][, POS1:=POS-(W/2)][,POS2:=POS+(W/2)]
coords_hrs<-split(coords_hrs, by='CHR')
lapply(coords_hrs, function(X) setkey(X, CHR, POS1, POS2))

RES<-vector('list', 22)

for(chr in 22:1){
        tp<-comb[CHR==chr]
        setkey(tp, CHR, POS1, POS2)
        coords_hrs[[chr]]$CHR<-as.integer(coords_hrs[[chr]]$CHR)
        setkey(coords_hrs[[chr]], CHR,POS1, POS2)
        foverlaps(tp, coords_hrs[[chr]], type='within', nomatch=0)[, Win:=paste0(POS1, "|",POS2)]-> res
        res[,MeanFreqDiff:=mean(Freq_Diff^2, na.rm=T), by="Win"]
        RES[[chr]]<-res
}
do.call(rbind, RES)-> RES2
RES2[, Data:='HRS']

RES2<-unique(RES2, by=c('CHR','Win'))

#calculate beta diff
merge(final_plink, prs_snps_hrs, by=c('CHR', 'POS'))-> betas_plink_hrs
betas_plink_hrs[, Beta_Diff:=b-PLINK]
betas_plink_hrs[, Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
betas_plink_hrs$CHR<-as.integer(betas_plink_hrs$CHR)
setkey(betas_plink_hrs,CHR, POS)
setkey(RES2, CHR, POS)
dt2<-betas_plink_hrs[RES2, nomatch=0]
dt2[,Quantile:=cut(MeanFreqDiff, breaks=quantile(MeanFreqDiff, probs=seq(0,1, by=0.2), na.rm=T), include.lowest=T)]
dt2[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
dt2[, MedianMeanFreqDiff:=median(MeanFreqDiff, na.rm=T), by=Quantile]
saveRDS(dt2, paste0('~/height_prediction/epistasis/output/dt2_', W, '.Rds'))

nrow(dt2)
my_dt<-dt2[Beta_Diff_Chisq<=15]
nrow(my_dt)
model_hrs<-lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2)
require(broom)
glance(model_hrs)
pval<-glance(model_hrs)$p.value
plotB<-ggplot(my_dt, aes(x=MeanFreqDiff, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col="black", alpha=0.6) +
geom_smooth(method='lm', se=T, col='darkgray', lty=2, lwd=1) +
geom_point(aes(x=MedianMeanFreqDiff, y=MeanBetaDiffChisq), col='red', cex=0.5) +
labs(x="Mean Squared Frequency Difference", y=expression(chi[diff]^2), cex=17) +
annotate("text", x=0.35, y=4, label=paste("p = ", round(pval,4))) +
theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 17), axis.title.x=element_text(size=17))


plot_grid(plotA, plotB,labels = c("A", "B"), nrow=1)
ggsave('figs/Fig5.png', width=12, height=4.5)
