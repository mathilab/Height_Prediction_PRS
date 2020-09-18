#!/usr/bin/env Rscript
library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scales)
source('scripts/Rsq_R2.R')
##panel figure

#R2 vs ancestry for HRS_eur, HRS_afr, UKB_afr imputed data

#barplot eur_afr, imp, non-imp, prunings
comp<-readRDS('imputed/output/comparison.Rds')

dt<-comp[Name %in% c("phys_500000_0.000005", "phys_100000_0.0005")][, .(Name, HRS_afr_imp, HRS_afr, HRS_eur_imp, HRS_eur)]
melt(dt)-> dt2
dt2$Name<-rep(c("100 Kb p<0.005","500 Kb p<0.000005"),4)
dt2[, Dataset:=c(rep('HRS_afr',4), rep('HRS_eur',4))]
dt2$Dataset<-factor(dt2$Dataset, levels=c("HRS_afr", "HRS_eur"))
dt2[,Method:=rep(c(rep("Imputed",2), rep("Genotyped",2)),2)]

ggplot(dt2, aes(x=Name, y=value, fill=Method)) + 
facet_wrap(~Dataset) + geom_bar(stat="identity",position=position_dodge(), alpha=0.8) + 
#scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + 
scale_fill_manual(values=c(brewer.pal(4, 'Set1')[4],"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position='bottom', legend.title=element_blank()) + labs(y=expression(Partial~R^2)) + scale_alpha_manual(values=c(0.1, 0.5))
ggsave('figs/test.pdf')


comp_ukb<-readRDS('imputed/output/comparison_ukb.Rds')
dt_ukb<-comp_ukb[Name %in% c("phys_500000_0.000005", "phys_100000_0.0005")][, .(Name, UKB_afr_imp, UKB_afr)]
melt(dt_ukb)-> dt2_ukb

dt2_ukb$Name<-rep(c("100 Kb p<0.005","500 Kb p<0.000005"),2)
dt2_ukb[, Dataset:=rep('UKB_afr',4)]
dt2_ukb[,Method:=c(rep("Imputed",2), rep("Genotyped",2))]


rbind(dt2, dt2_ukb)

rbind(dt2_ukb, dt2)-> dt_final
dt_final$Dataset<-factor(dt_final$Dataset, levels=c("UKB_afr", "HRS_afr", "HRS_eur"))

ggplot(dt_final, aes(x=Name, y=value, fill=Method)) + facet_wrap(~Dataset) + 
geom_bar(stat="identity",position=position_dodge(), alpha=0.8) + 
scale_fill_manual(values=c(brewer.pal(4, 'Set1')[c(1,4)],"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank(), axis.text.x = element_text(angle = 45)) + labs(y=expression(Partial~R^2)) + scale_alpha_manual(values=c(0.1, 0.5))
#ggsave('figs/Fig4_B.pdf')


#now do R2 vs anc
##############################################################
#combine all

##
readRDS(paste0('imputed/output/results.UKB_afr.Rds'))-> results.UKB_afr
readRDS(paste0('imputed//output/results.HRS_eur.Rds'))-> results.HRS_eur
readRDS(paste0('imputed/output/results.HRS_afr.Rds'))-> results.HRS_afr
readRDS(paste0('imputed/output/B_UKB_afr.Rds'))-> B_UKB_afr
readRDS(paste0('imputed/output/B_HRS_eur.Rds'))-> B_HRS_eur
readRDS(paste0('imputed/output/B_HRS_afr.Rds'))-> B_HRS_afr
readRDS(paste0('imputed/output/PGS3_UKB_afr.Rds'))-> PGS3_UKB_afr
readRDS(paste0('imputed/output/PGS3_HRS_eur.Rds'))-> PGS3_HRS_eur
readRDS(paste0('imputed/output/PGS3_HRS_afr.Rds'))-> PGS3_HRS_afr

I<-which(names(B_UKB_afr)=="100000_0.0005")
ALL2<-rbind(B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
tmp<-1/c(var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t))
cbind(ALL2, W=tmp)-> ALL2
ALL2$Dataset<-factor(ALL2$Dataset, levels=c("UKB_afr", "HRS_afr", "HRS_eur"))
shapes<-c(16,6,15)
my_plotA<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=5, fill="white", cex=0.8) + stat_smooth(data=ALL2,method = "lm", mapping = aes(weight = sqrt(W)), col='darkgray', lty=2) +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
 #       geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  
	scale_color_manual(values=c(brewer.pal(4, 'Set1')[c(1,4)],"#101010")) +
	scale_shape_manual(values=shapes)+
        ylab(expression(paste("Partial-R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=12))

dt_final[, label:=rep(c("100 Kb \n p<0.0005", "500 Kb \n p<0.0005"),6)]
my_plotB<-ggplot(dt_final, aes(x=label, y=value, fill=Method, label=label)) + 
facet_wrap(~Dataset, ncol = 1) + geom_bar(stat="identity",position=position_dodge(), alpha=0.8) + 
scale_fill_manual(values=c("#96a8b2", "#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank(), axis.title.y = element_text(size = 15), axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), legend.text=element_text(size=12)) + 
labs(y=expression(Partial-R^2)) + scale_alpha_manual(values=c(0.1, 0.5))

#additive genetic variance for imputed data and non-imputer data

gen_phen_var<-data.table(Method=c(rep("Genotyped", 2), rep("Imputed", 2)), Dataset=rep(c( "UKB", "HRS"), 2), Method2=c(rep("Genotyped", 2), rep("Imputed", 2)))
gen_phen_var[,value:=c(0.78,0.92,0.8039084, 0.86)]
gen_phen_var[,y_min := value*0.5, by = Dataset]
gen_phen_var[,y_max := value*1.5, by = Dataset]
my_plotC<-ggplot(gen_phen_var, aes(x=Method2, y=value, fill=Method))+
geom_bar(stat="identity",position=position_dodge(3), alpha=0.8) +
facet_wrap(~Dataset, ncol = 1, scales="free_y")+
scale_fill_manual(values=c("#96a8b2", "#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank(), axis.title.y = element_text(size = 15), axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0, 1))+
labs(y=expression(G[PRS]~ratio)) + scale_alpha_manual(values=c(0.1, 0.5))

pdf('figs/Fig4.pdf',height =9, width=7)
plot_grid(my_plotA, plot_grid(my_plotB, my_plotC, nrow=1, labels=c("B", "C")),nrow=2, labels="A")
dev.off()
