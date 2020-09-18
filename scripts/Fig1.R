library("optparse")
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(cowplot)
library(psychometric)
library(boot)
library(RColorBrewer)
library(cowplot)
options(scipen=999)

i<-'phys_100000_0.0005'
args<-'gwas'
readRDS(paste0(args[1], '/WHI/output/results.WHI.Rds'))[[i]]-> results.WHI
readRDS(paste0(args[1], '/JHS/output/results.JHS.Rds'))[[i]]-> results.JHS
readRDS(paste0(args[1], '/ukb_afr/output/results.UKB_afr.Rds'))[[i]]-> results.UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/results.HRS_eur.Rds'))[[i]]-> results.HRS_eur
readRDS(paste0('~/height_prediction/',args[1], '/HRS_afr/output/results.HRS_afr.Rds'))[[i]]-> results.HRS_afr
readRDS(paste0('~/height_prediction/',args[1], '/WHI/output/B_WHI.Rds'))[[i]]->B_WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/B_JHS.Rds'))[[i]]->B_JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/B_UKB_afr.Rds'))[[i]]-> B_UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/B_HRS_eur.Rds'))[[i]]-> B_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/B_HRS_afr.Rds'))[[i]]-> B_HRS_afr
readRDS(paste0('~/height_prediction/',args[1],'/ukb_afr/output/PGS3_UKB_afr.Rds'))[[i]]-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/',args[1],'/WHI/output/PGS3_WHI.Rds'))[[i]]-> PGS3_WHI
readRDS(paste0('~/height_prediction/',args[1],'/JHS/output/PGS3_JHS.Rds'))[[i]]-> PGS3_JHS
readRDS(paste0('~/height_prediction/',args[1],'/HRS_eur/output/PGS3_HRS_eur.Rds'))[[i]]-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/PGS3_HRS_afr.Rds'))[[i]]-> PGS3_HRS_afr

ALL<-rbind(B_JHS[1:2,][, Dataset:='JHS_afr'], B_WHI[1:4,][, Dataset:='WHI_afr'], B_UKB_afr[1:4,][,Dataset:='UKB_afr'],B_HRS_afr[1:2,][, Dataset:='HRS_afr'], B_HRS_eur)

tmp<-1/c(var(results.JHS[[1]]$t), var(results.JHS[[2]]$t), var(results.WHI[[1]]$t), var(results.WHI[[2]]$t), var(results.WHI[[3]]$t), var(results.WHI[[4]]$t), var(results.UKB_afr[[1]]$t),var(results.UKB_afr[[2]]$t), var(results.UKB_afr[[3]]$t), var(results.UKB_afr[[4]]$t), var(results.HRS_afr[[1]]$t), var(results.HRS_afr[[2]]$t), var(results.HRS_eur$t))  #weighing lm by boostrap replicates.

ALL[,W:=tmp]
ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
shapes<-c(16,17,18,6,15)

A_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
	geom_point(aes(shape=Dataset), size=5, fill="white", alpha=0.8) + stat_smooth(data=ALL,method = "lm", mapping = aes(weight = sqrt(W)), col='darkgray', lty=2, lwd=1) +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
#        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +
        scale_color_manual(values=my_colrs) + coord_cartesian(ylim = c(0, 0.16)) +
	scale_shape_manual(values=shapes)+
        ylab(expression(paste("Partial-R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), 
	axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15), legend.position=c(0.10,0.85))


#####################3
#Second plot
ALL2<-readRDS('~/height_prediction/loc_anc_analysis/output/eur_only.Rds')

ALL2[Dataset!= 'HRS_eur']-> ALL3
ALL2[, y_k1:=ALL2[Dataset=='HRS_eur']$R_sq*(Med_Eur_Anc^1)]
ALL2[, y_k2:=ALL2[Dataset=='HRS_eur']$R_sq*(Med_Eur_Anc^2)]
ALL2[, y_k1.5:=ALL2[Dataset=='HRS_eur']$R_sq*(Med_Eur_Anc^1.5)]
ALL2[, y_k1.8:=ALL2[Dataset=='HRS_eur']$R_sq*(Med_Eur_Anc^1.8)]
a<-ALL2[Dataset=='HRS_eur']$R_sq
dts<-data.table(X=rep(seq(0,1, 0.01),3), K=c(rep(1,101), rep(1.5,101), rep(2, 101)))[, Y:=a*(X^K)]

B_plot<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=5, fill="white", alpha=0.8) +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
	geom_path(data=dts, aes(x=X,y=Y), col="orange")+
	#geom_line(aes(x=Med_Eur_Anc,y=y_k1.5), lty=2, col="darkorange", pwd=2)+
	#geom_line(aes(x=Med_Eur_Anc,y=y_k1.8), lty=2, col="darkorange2", pwd=3)+
	#geom_line(aes(x=Med_Eur_Anc,y=y_k2), lty=2, col="darkorange4", pwd=4)+
	#stat_smooth(data=ALL3,method = "lm", mapping = aes(weight = sqrt(ALL3$W)), col='darkgray', lty=2, lwd=1) +
        scale_color_manual(values=my_colrs) + coord_cartesian(ylim = c(0, 0.16)) + 
	scale_shape_manual(values=shapes)+
        ylab(expression(paste("Partial-R"^"2"))) + xlab("European Ancestry Proportion") +
	annotate("text", x=0.39, y=0.07, label=paste("k=1"), size=4.5) +
	annotate("text", x=0.52, y=0.07, label=paste("k=1.5"), size=4.5) +
	annotate("text", x=0.72, y=0.07, label=paste("k=2"), size=4.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15),legend.position = "none")

        #print(B_plot)

#png('figs/Fig1_v1.png', res=300, unit="in", height=8, width=12)
#plot_grid(A_plot, B_plot, nrow=1, labels=c('A', 'B'))
#dev.off()

pdf('figs/Fig1.pdf',  width=7, height=9)
plot_grid(A_plot, B_plot, nrow=2, labels=c('A', 'B'))
dev.off()

