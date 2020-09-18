#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)
########################
PRS1_all_res<-readRDS(file='loc_anc_analysis/output/PRS1_all_res.Rds')
PRS2_all_res<-readRDS(file='loc_anc_analysis/output/PRS2_all_res.Rds')
PRS1_2_all_res<-readRDS(file='loc_anc_analysis/output/PRS1_2_all_res.Rds')

dt<-melt(PRS1_2_all_res,id=c('Alpha', 'Dataset', 'Test'))

ggplot(dt, aes(x=Alpha, y=value,colour=Test)) + facet_wrap(~Dataset) +
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + scale_color_manual(values=c("darkseagreen4", "darkslateblue", "darkorange3")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2))

ggsave('figs/PRS1_2.pdf')

#just WHI
PRS1_all_res<-readRDS(file='loc_anc_analysis/output/PRS1_all_res.Rds')
PRS2_all_res<-readRDS(file='loc_anc_analysis/output/PRS2_all_res.Rds')
PRS1_2_all_res<-readRDS(file='loc_anc_analysis/output/PRS1_2_all_res.Rds')

PRS1_2_all_res<-PRS1_2_all_res[Dataset=='WHI_afr']
dt<-melt(PRS1_2_all_res,id=c('Alpha', 'Dataset', 'Test'))

ggplot(dt, aes(x=Alpha, y=value,colour=Test)) + 
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + scale_color_manual(values=c("darkseagreen4", "darkslateblue", "darkorange3")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2))

ggsave('figs/PRS1_2_WHI.pdf')

##add PRS3

PRS3_WHI<-readRDS('loc_anc_analysis/output/part_R2_WHI_unscaled.Rds')
dt_PRS3_WHI<-data.table(Part_R2=unlist(PRS3_WHI[1:101]), Alpha=seq(0,1, 0.01), Test="PRS3", Dataset="WHI_afr")
PRS3_JHS<-readRDS('loc_anc_analysis/output/part_R2_JHS_unscaled.Rds')
dt_PRS3_JHS<-data.table(Part_R2=unlist(PRS3_JHS[1:101]), Alpha=seq(0,1, 0.01), Test="PRS3", Dataset="JHS_afr")
PRS3_HRS<-readRDS('loc_anc_analysis/output/part_R2_HRS_unscaled.Rds')
dt_PRS3_HRS<-data.table(Part_R2=unlist(PRS3_HRS[1:101]), Alpha=seq(0,1, 0.01), Test="PRS3", Dataset="HRS_afr")

PRS1_2_all_res<-readRDS(file='loc_anc_analysis/output/PRS1_2_all_res.Rds')
#test5[, Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]
PRS123<-rbind(PRS1_2_all_res, dt_PRS3_WHI, dt_PRS3_JHS, dt_PRS3_HRS)

melt(PRS123,id=c('Dataset', 'Alpha', 'Test'))->dt2

png("figs/multi_prs.png", unit='in', res=300, height=5, width=9)
ggplot(dt2, aes(x=Alpha, y=value,colour=Test))  + facet_wrap(~Dataset)+
#facet_wrap(~Dataset2) + 
geom_line(size=1.2) + 
coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + 
scale_color_manual(values=c("darkseagreen4", "darkslateblue", "deeppink4", "gray7")) +
theme(strip.text.x = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=14, angle=45), axis.text.y=element_text(size=14), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial-R^2))
dev.off()

dt2_WHI<-dt2[Dataset=='WHI_afr']
dt2_WHI[,value2:=ifelse(Test=='PRS3', value/PRS3_WHI[[1]], value/PRS3_WHI[[102]])]
dt2_JHS<-dt2[Dataset=='JHS_afr']
dt2_JHS[,value2:=ifelse(Test=='PRS3', value/PRS3_JHS[[1]], value/PRS3_JHS[[102]])]
dt2_HRS<-dt2[Dataset=='HRS_afr']
dt2_HRS[,value2:=ifelse(Test=='PRS3', value/PRS3_HRS[[1]], value/PRS3_HRS[[102]])]
dt3<-rbind(dt2_WHI, dt2_JHS, dt2_HRS)

pdf("figs/Fig2.pdf",width=7, height=4.5)
ggplot(dt3, aes(x=Alpha, y=value2,colour=Test))  + facet_wrap(~Dataset) +
geom_line(size=1.2) +
geom_hline(yintercept=1, linetype="dashed", color = "darkgray")+
coord_cartesian(xlim=c(0,0.5), ylim=c(0.5, 1.5)) +
scale_color_manual(values=c("darkseagreen4", "darkslateblue", "deeppink4", "gray7"), labels = c(expression("PRS"[c]^1), expression("PRS"[c]^2), expression("PRS"[c]^3))) +
#geom_hline(yintercept=1) +
theme(strip.text.x = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=14, angle=45), axis.text.y=element_text(size=14), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=16)) + 
#labs(x=expression(alpha), y=expression(Partial-~R^2~relative~change))
labs(x=expression(alpha), y=expression(paste("Partial-R"^"2"~"Relative Change")))
dev.off()

