#!/usr/bin/env Rscript
######################################################################################################################
#Name:		Plot_v2
#Args:		<summary statistics: gwas/sib_betas> <recombination map:AA/CEU> <win size:3000_100000> <test_set(WHI, HRS_afr, HRS_eur, JHS> 
#
#Last modified:	July 13-2020
#Usage:		Rscript --vanilla scripts/Fig3.r gwas CEU 20000 WHI
######################################################################################################################
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
#}
#Load packages #########################
suppressPackageStartupMessages({library("optparse")
library(ggplot2)
library(reshape2)
library(data.table)
library(reshape)
library(RColorBrewer)
library(MASS)
library(cowplot)
library(dplyr)
require(broom)
})
options(scipen=999)
source('scripts/fancy_scientific.R')
dtset<-c()
ci<-c()

cat('Analysing the ', args[1], ' summary statistics with ', args[2], ' recombination map and testing on ', args[3], '\n')
dtset<-vector('list', 5)
ci<-vector('list', 5)
names(dtset)<-c("ukb_afr", "WHI", "JHS", "HRS_afr",  "HRS_eur")
names(ci)<-c("ukb_afr", "WHI", "JHS", "HRS_afr",  "HRS_eur")
for(I in c("WHI","JHS","ukb_afr","HRS_eur", "HRS_afr")){
	dtset[[I]]<-vector('list', 6)
	ci[[I]]<-vector('list', 6)
	dtset[[I]]<-readRDS(paste0('strat_prs/output/part_R2_', I,"_",args[1], "_", args[3],"_", args[2], "_v2.Rds"))
	ci[[I]]<-readRDS(paste0('strat_prs/output/results_', I,"_",args[1], "_", args[3], "_", args[2], "_v2.Rds"))
}

cat('Done reading in stratified PRS results for ', args[2], ' map\n')
if(args[1]=='gwas'){ ##UKB_afr, WHI_afr, JHS_afr, HRS_afr,HRS_eur
	r2_vec<-c(0.04084170,0.03578924,0.0382532, 0.03103462, 0.1563672)
	r2_l_vec<-c(0.03255115,0.027783862,0.02252898,0.01906344,0.1438352)
	r2_u_vec<-c(0.04955051,0.04474530,0.05733168, 0.04574303,0.1680824)
} else if (args[1]=='sib_betas'){
#	r2_vec<-c(0.07511196,0.00967814,0.01993696,0.02717107,0.01060191)
}
cat('Read and edit quantiles\n')

my_lev<-levels(readRDS(paste0("strat_prs/output/rec_quant_", args[4], "_", args[1], "_", args[3], "_", args[2], "_v2.Rds"))[,Q:=cut(diff, breaks=quantile(diff), na.rm=T, include.lowest=T, dig.lab=1)]$Q)

df1<-rbind(data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='HRS_eur',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['HRS_eur']])/r2_vec[5],
Perc_L=c(unlist(lapply(ci[['HRS_eur']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['HRS_eur']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='JHS_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['JHS']])/r2_vec[3],
Perc_L=c(unlist(lapply(ci[['JHS']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['JHS']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='WHI_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['WHI']])/r2_vec[2],
Perc_L=c(unlist(lapply(ci[['WHI']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['WHI']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='UKB_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['ukb_afr']])/r2_vec[1],
Perc_L=c(unlist(lapply(ci[['ukb_afr']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['ukb_afr']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='HRS_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['HRS_afr']])/r2_vec[4],
Perc_L=c(unlist(lapply(ci[['HRS_afr']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['HRS_afr']], function(X) X$percent[5])))))
####
df2<-rbind(data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='HRS_eur',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['HRS_eur']]),
Perc_L=c(unlist(lapply(ci[['HRS_eur']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['HRS_eur']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='JHS_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['JHS']]),
Perc_L=c(unlist(lapply(ci[['JHS']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['JHS']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='WHI_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['WHI']]),
Perc_L=c(unlist(lapply(ci[['WHI']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['WHI']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='UKB_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['ukb_afr']]),
Perc_L=c(unlist(lapply(ci[['ukb_afr']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['ukb_afr']], function(X) X$percent[5])))),
data.table(Quantile=c("q1","q2","q3","q4"),
Dataset='HRS_afr',Win=c(rep(args[3], 4)),
R2=unlist(dtset[['HRS_afr']]),
Perc_L=c(unlist(lapply(ci[['HRS_afr']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['HRS_afr']], function(X) X$percent[5])))))
#
df1[, Set:='Relative'] #order ##UKB_afr, WHI_afr, JHS_afr, HRS_afr,HRS_eur
df1[, Perc_U:=Perc_U/c(rep(r2_u_vec[5],4), rep(r2_u_vec[3], 4), rep(r2_u_vec[2], 4), rep(r2_u_vec[1], 4), rep(r2_u_vec[4], 4))]
df1[, Perc_L:=Perc_L/c(rep(r2_l_vec[5],4), rep(r2_l_vec[3], 4), rep(r2_l_vec[2], 4), rep(r2_l_vec[1], 4), rep(r2_l_vec[4], 4))]
df2[, Set:='Absolute']
df2[, Win:=args[3]]
df1[, Win:=args[3]]
df4<-rbind(df1, df2)
#melt(df4, id=c("Quantile","Win", "Map", "Perc_L", "Perc_U", "Dataset"))-> df5
df4$Dataset<-factor(df4$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
df1$Dataset<-factor(df1$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
df4$Set<-factor(df4$Set, levels=c("Absolute", "Relative"))
pd <- position_dodge(0.5)
plot1<-ggplot(df4, aes(x=Quantile, y=R2, colour=Dataset)) + 
geom_point(position=pd) +
geom_errorbar(aes(ymin=Perc_L, ymax=Perc_U), position = pd) +
facet_wrap(. ~Set, scales='free_y') + 
labs(y=expression(paste("Partial-R"^"2")), x=paste0("Recombination Rate (cm/",as.character(as.integer(args[3])/10000),"Kb)"))+  
scale_colour_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=11), axis.text.y=element_text(size=11),legend.text=element_text(size=11),legend.position = "bottom", strip.text.x = element_text(size = 13))+
scale_x_discrete(labels=my_lev)
#
df2$Dataset<-factor(df2$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))

pd <- position_dodge(0.5)
plot2<-ggplot(df2, aes(x=Quantile, y=R2, colour=Dataset)) +
geom_point(position=pd) +
geom_errorbar(aes(ymin=Perc_L, ymax=Perc_U), position = pd) +
labs(y=expression(paste("Partial-R"^"2")), x=paste("Recombination Rate (cm/",as.character(as.integer(args[3])/10000),"Kb)"))+
scale_colour_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right", legend.direction='vertical', legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10))+
scale_x_discrete(labels=my_lev)
########################
########################
########################

if(args[1]=='gwas'){
	beta<-dplyr::select(do.call(rbind,readRDS(paste0('gwas/', args[4],'/output/hei_phys_100000_0.0005_v2.Rds'))), CHR, POS, b, MarkerName, Allele1, SE) %>% as.data.table
	beta1<-readRDS('loc_anc_analysis/output/final_plink.Rds')
} else {
#TO DO
}

merge(beta, beta1, by=c('CHR', 'POS'))-> beta2
cat('Read in betas\n')
BETA<-vector('list', 22)
prun<-"phys_100000_0.0005"
hei<-readRDS(paste0('gwas/', args[4], '/output/hei_', prun, '_v2.Rds'))
gc()
if(args[2]=='AA'){
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps/hm2/hm2/genetic_map_GRCh37_chr', chr,'.txt.gz'))[,CHR:=gsub("chr","",Chromosome)][, Chromosome:=NULL])-> rec #need to fix this path
for(chr in 1:22){colnames(rec[[chr]])<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')}
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps/maps_b37/maps_chr.', chr, '.gz')))-> maps #need to fix this path
for(chr in 1:22){colnames(maps[[chr]])[1]<-"POS"}
lapply(1:22, function(chr) setkey(rec[[chr]],POS))
lapply(1:22, function(chr) setkey(maps[[chr]],POS))
lapply(1:22, function(chr) maps[[chr]][rec[[chr]], nomatch=0])-> map
remove(maps)
} else if(args[2]=='CEU'){
lapply(1:22, function(chr) fread(paste0('/project/mathilab/data/maps/1000_genomes_maps/hg19/CEU/CEU_recombination_map_hapmap_format_hg19_chr_', chr, '.txt')))-> map
for(chr in 1:22){map[[chr]][,CHR:=gsub("chr", "", Chromosome)][,Chromosome:=NULL]}
for(chr in 1:22){colnames(map[[chr]])<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')}
}
lapply(1:22, function(chr) map[[chr]][,POS2:=POS])
lapply(1:22, function(chr) arrange(map[[chr]], CHR, POS))
lapply(1:22, function(chr) setDT(map[[chr]]))
cat('loading done\n')
for(chr in 22:1){
	cat('start chr ')
	cat(chr)
	cat('\n')
	rate.dist<-as.integer(args[3])
	betas<-beta2[CHR==chr]
	snps <- read.table(paste0("input/", args[4], "/", args[4], "_b37_strand_include_kgCY_chr", chr, ".phsnp"), as.is=TRUE)
        colnames(snps) <- c("ID", "CHR", "Map", "POS", "REF", "ALT")
	betas<-merge(snps, betas, by=c('CHR', 'POS'))
	if(args[2]=='AA'){
	REC.rate <- approxfun(map[[chr]]$POS, map[[chr]]$AA_Map, rule=2)
	} else if(args[2]=='CEU'){
	REC.rate <- approxfun(map[[chr]]$POS, map[[chr]]$MAP_cM, rule=2)
	cat('checkpoint \n')
	}
  	setDT(betas)
	betas$rec.rate <- 0
	for(i in 1:NROW(betas)){
    		cat(i, '\r')
    		rec.diff <- REC.rate(betas$POS[i]+(rate.dist/2))-REC.rate(betas$POS[i]-(rate.dist/2))
    		betas$rec.rate[i] <- rec.diff
	}
	betas[,Beta_Diff:=b-PLINK]
	betas[,Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
	#At this point you will restrict the betas to the SNPS that you are using in the PRS.
	BETA[[chr]]<-betas
	cat('chr ')
	cat(chr)
	cat(' done\n')
}
cat('sleep\n')
cat('Checkpoint after rec rate\n')
Sys.sleep(15)
do.call(rbind, BETA)-> BETA
setDT(BETA)
BETA[order(CHR, POS)]-> BETA
as.factor(BETA$CHR)-> BETA$CHR

saveRDS(BETA, file=paste0('gwas/', args[4], '/output/plink_', toupper(args[4]),"_", args[1], "_", args[2], "_", args[3],".Rds"))
#########################
########################
########################
#########################
########################
########################
BETA[,Quantile:=cut(rec.rate, breaks=quantile(rec.rate, probs=seq(0,1, by=0.05), na.rm=T), include.lowest=T)]
BETA[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
BETA[, MedianRecRate:=median(rec.rate, na.rm=T), by=Quantile]
lm_test0<-lm(Beta_Diff_Chisq~rec.rate, data=BETA)
glance(lm_test0)
pval<-round(glance(lm_test0)$p.value,3)
nrow(BETA)-sum(BETA$Beta_Diff_Chisq<=15) #10
BETA<-BETA[Beta_Diff_Chisq<=15]
nrow(BETA)
tp<-as.character(as.integer(args[3])/10000)
plot3<-ggplot(BETA, aes(x=rec.rate, y=Beta_Diff_Chisq)) + 
geom_point(cex=0.5, col='black', alpha=0.6) + 
geom_smooth(method='lm', se=T, lwd=1, col="darkgray", lty=2) + 
labs(y=expression(chi[diff]^2), x=paste0("Recombination Rate (cM/", tp, "Kb)")) + 
geom_point(aes(x=MedianRecRate, y=MeanBetaDiffChisq, col="red"), cex=0.5) + 
annotate("text", x=0.4, y=4, label=paste("p=", round(pval,4)), size=4) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none", legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))
plot1a<-plot1 + guides(fill=FALSE)
cat('CHECKPOINT 2\n')
c_plot<-readRDS(paste0('output/c_plot.Rds'))
my_list<-list(plot1, plot2, c_plot, plot3)
saveRDS(my_list, paste0('strat_prs/output/my_list_', args[2], '_',args[3], '.Rds'))
ld<-do.call(rbind, lapply(1:22, function(X) fread(paste0('eur_w_ld_chr/', X,'.l2.ldscore.gz'))))
#beta<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')
colnames(ld)[3]<-'POS'
colnames(ld)[2]<-'MarkerName.y'
setkey(ld, CHR, POS)
#beta$CHR<-as.integer(beta$CHR)
as.integer(BETA$CHR)-> BETA$CHR
setkey(BETA, CHR, POS)
#beta[ld,nomatch=0]-> test
test<-merge(BETA,ld, by=c('CHR', 'POS','MarkerName.y'))

test[,Quantile:=cut(L2, breaks=quantile(L2, probs=seq(0,1, by=0.05), na.rm=T), include.lowest=T)]
test[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
test[, MedianL2:=median(L2, na.rm=T), by=Quantile]
cat('CHECKPOINT 3\n')
#plot4<-ggplot(test, aes(x=MedianL2, y=MeanBetaDiffChisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm', se=F) + labs(y=expression(Mean~chi^2), x="LD Score" )
lm_test<-lm(Beta_Diff_Chisq~L2, data=test)
glance(lm_test)
pval<-glance(lm_test)$p.value
plot4<-ggplot(test, aes(x=L2, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col='black', alpha=0.6) + geom_smooth(method='lm', se=T, col='darkgray') + 
labs(y=expression(chi[diff]^2), x="LD Score" ) + 
geom_point(aes(x=MedianL2, y=MeanBetaDiffChisq, col="red"), cex=0.5) + 
annotate("text", x=270, y=5, label=paste("p=", round(pval,4)), size=4) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none", legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))

cat('CHECKPOINT 4\n')
if(args[2]=='AA'){
png('figs/Fig3.png',  res=600, unit="in", height=9, width=11)
plot_grid(plot1,plot_grid(c_plot,plot3,labels = c("C", "D"), nrow=1, align='h', rel_widths=c(1,1), axis="t"), nrow=2, labels="A", align="h", axis="l")
cat('sleep...\n')
dev.off()
cat('CHECKPOINT 5\n')
png(paste0('figs/FigS11_', args[2], '.png'), res=600, unit="in", height=9, width=7)
print(plot4)
dev.off()
cat('CHECKPOINT 6\n')
} else {
png(paste0('figs/FigS14_', args[2], '.png'),  res=600, width=7, height=9, units="in")
plot_grid(plot1,plot_grid(c_plot,plot3,labels = c("C", "D"), nrow=1, align='h', rel_widths=c(1,1), axis="t"), nrow=2, labels="A", align="h", axis="l")
cat('CHECKPOINT 7\n')
dev.off()
}
#The End
