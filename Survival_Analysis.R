
#--------------survival plots------------------------#

setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/")

set.seed(123)

library(dplyr)
library(survival)
library(survminer)
library(survival)


cln.order<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep = ";")

cbind(cln.order$metastasis,
      cln.order$relasp_status,
      cln.order$metastasis_date,
      cln.order$relaps_date)

cln.order[cln.order$FU_metastasis_latestControl>cln.order$FU_tod_latestControl,]

cln.order$FU_tod_latestControl_new<-ifelse(cln.order$FU_metastasis_latestControl>cln.order$FU_tod_latestControl,
                                           cln.order$FU_metastasis_latestControl,
                                           cln.order$FU_tod_latestControl)


#cln.order<-cln.order[!cln.order$newID %>% is.na(),]

death_type<-list(c(1:4),1,c(1:4),1,c(1:4))
wrt<-c( "Cluster", "Cluster","Group","Group")

survival.formula<-sapply(wrt, function(x) as.formula(paste('Surv(Time, Status)~', x)))

model.srv<-list()
p.val<-c()

col.cluster<-c("blue", "red3", "green4", "black", "orange")
col.group<-color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3")


theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 


par(mar = c(5, 4, 4, 1))

for(k in 1:4){
  
  cli.data<-cln.order[,c("oldID", "Status", "FU_tod_latestControl_new",  "Group", "Cluster")]
  names(cli.data)<-c("ID", "Status", "Time", "Group", "Cluster")
  
  bc_death<-cln.order$Status
  bc_death<-ifelse(bc_death %in% death_type[[k]], 1, ifelse(bc_death != 0, NA, 0))
  
  cli.data$Status<-bc_death
  
  #cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
  #cli.data<-na.omit(cli.data)
  cli.data$Time<-as.numeric(cli.data$Time)
  
  model.srv[[k]]<-survfit(survival.formula[[k]], data=cli.data)
  p.val<-c(p.val, surv_pvalue(survfit(survival.formula[[k]], data=cli.data))$pval)
  
  if(wrt[k]=="Cluster"){my.col=col.cluster; 
  legend=paste("MOC", 1:3, sep=""); 
  color="Cluster"}else{
    my.col=col.group; 
    legend=cln.order$Group %>% sort() %>% unique; 
    color="Group"}
  
  ggsurv<-ggsurvplot(model.srv[[k]], palette= my.col, 
                     size.lab=1.5, 
                     pval = paste("p =", round(p.val[k] %>% as.numeric(),digits = 3)),#formatC(p.val[k] %>% as.numeric(), format = "d", digits = 2) %>% round(digits = 3)), 
                     font.size=1.5,
                     surv.scale = "percent",
                     font.legend = list(size = 15, color = "black", face = "bold"),
                     font.tickslab = c(20, "plain", "black"),
                     risk.table = T, 
                     pval.method = T,
                     legend.title="",
                     font.x = c(20, "bold", "black"),
                     font.y = c(20, "bold", "black"),
                     pval.size = 7, pval.coord = c(60, .1),
                     conf.int = T,
                     legend.labs = legend,
                     ggtheme = theme,
                     risk.table.fontsize = 6,
                     conf.int.alpha=0.05,
                     xlab="Months"
  )
  
  ggsurv$table <- ggrisktable(model.srv[[k]], 
                              data = cli.data, 
                              color = color, 
                              palette=my.col,
                              y.text = T,   ylab = "",  xlab = "",
                              tables.theme = theme,
                              font.tickslab = c(15, "bold"),
                              legend.labs = legend,
                              #ggtheme = theme_survminer(),
                              fontsize = 6,
                              clip="off"
                              
  ) + coord_cartesian(xlim = c(0, 200), 
                      clip = 'on')
  
  ggsurv %>% print()
  
}


wald.p<-c()
comb<-list(c("MOC1","MOC2"),c("MOC1","MOC3"),c("MOC2","MOC3"))

for(i in 1:3){
  
  cli.data<-cln.order[,c("oldID", "Status", "FU_tod_latestControl_new",  "Group", "Cluster")] %>% as.data.frame()
  
  cli.data$MOC<-ifelse(cli.data$Cluster==1,"MOC1",
                       ifelse(cli.data$Cluster==2,"MOC2",
                              ifelse(cli.data$Cluster==3,"MOC3",NA)))
  
  names(cli.data)<-c("ID", "Status", "Time", "Group","Cluster", "MOC")
  
  bc_death<-cln.order$Status
  bc_death<-ifelse(bc_death %in% death_type[[2]], 1, ifelse(bc_death != 0, NA, 0))
  
  cli.data$Status<-bc_death
  
  #cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
  #cli.data<-na.omit(cli.data)
  cli.data$Time<-as.numeric(cli.data$Time)
  
  #cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
  cli.data<-cli.data %>% dplyr::filter(MOC%in%comb[[i]])
  #cli.data<-na.omit(cli.data)
  cli.data$Time<-as.numeric(cli.data$Time)
  
  cli.data$MOC<-cli.data$MOC %>% as.factor()
  
  model.srv1<-survfit(Surv(Time, Status) ~ MOC, data=cli.data)
  p.val1<-surv_pvalue(model.srv1, data=cli.data)
  #summary(model.srv1) %>% print()
  wald.p<-c(wald.p, round(p.val1$pval, digits = 3))
}
wald.p


cli.data<-cln.order[,c("oldID", "Status", "FU_tod_latestControl_new",  "Group", "Cluster")] %>% as.data.frame()

cli.data<-data.frame(cli.data, 
                     MOC= ifelse(cln.order$Cluster==1,"MOC1",
                                 ifelse(cln.order$Cluster==2,"MOC2",
                                        ifelse(cln.order$Cluster==3,"MOC3",NA))))


names(cli.data)<-c("ID", "Status", "Time", "Group", "Cluster", "MOC")

#bc_death<-cln.order$`Cause of death (1=breast cancer, 2=different disease, 3=unknown, 4=cause of death not updated (dead after des2017)`
bc_death<-ifelse(cli.data$Status %in% c(1:4), 1, ifelse(cli.data$Status != 0, NA, 0))

cli.data$Status<-bc_death
cli.data$MOC<-cli.data$MOC %>% as.factor()
cli.data$MOC<-relevel(cli.data$MOC, "MOC2")

#cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
#cli.data<-na.omit(cli.data)
cli.data$Time<-as.numeric(cli.data$Time)

res.cox <- coxph(Surv(Time, Status) ~ MOC, data =  cli.data, ties = "efron")
res.cox %>% summary()
print(forestmodel::forest_model(res.cox))






#--------------------Metastasis-------------------#

cln.order<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep = ";")

cln.order$HR<-ifelse(cln.order$ER=="pos"|cln.order$PR=="pos", "+",
                     ifelse(cln.order$ER %>% is.na, NA, "-"))

cln.order$TNBC<-ifelse(cln.order$HR=="-" & cln.order$HER2_status=="neg", "Yes",
                       ifelse(cln.order$HR %>% is.na | cln.order$HER2_status %>% is.na(),
                              NA, "No"))

cln.order$Age[cln.order$Cluster==3] %>% as.numeric() %>% mean(na.rm=T)
cln.order$Age[cln.order$Cluster==3] %>% as.numeric() %>% sd(na.rm=T)


write.csv(cln.order, "/mnt/work/Oslo2/data/Clinic_data_oslo2_HR&TNBC.csv", row.names=F, quote = F)



cln.order$FU_tod_latestControl_new<-ifelse(cln.order$FU_metastasis_latestControl>cln.order$FU_tod_latestControl,
                                           cln.order$FU_metastasis_latestControl,
                                           cln.order$FU_tod_latestControl)

cln.order$metastasis<-ifelse(cln.order$metastasis%in%c(NA, "0"),0, 1)
cln.order$metastasis[cln.order$newID %>% is.na()]<-NA

cln.order$BC_Meta<-ifelse(cln.order$metastasis==0, cln.order$Status, 1)

table(cln.order$relasp_status, cln.order$Cluster, useNA = "always")

c(22,124,2,7,3,134,4)*100/148

table(cln.order$HER2_status, cln.order$Cluster, useNA="always")



cbind(cln.order$TNBC, cln.order$HR)
table(cln.order$HR, useNA="always")*100/335






#cln.order<-cln.order[!cln.order$newID %>% is.na(),]

cli.data<-cln.order[,c("oldID", "BC_Meta", "FU_metastasis_latestControl",  "Group", "Cluster")]
names(cli.data)<-c("ID", "Status", "Time", "Group", "Cluster")

cli.data$Status<-ifelse(cli.data$Status%in%c(1), 1, ifelse(cli.data$Status != 0, NA, 0))

model.cluster<-survfit(Surv(Time, Status)~Cluster, data=cli.data)
model.sbgrp<-survfit(Surv(Time, Status)~Group, data=cli.data)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 



ggsurv<-ggsurvplot(model.cluster, palette= col.cluster, 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = paste("MOC", 1:3, sep=""),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.cluster, 
                            data = cli.data, 
                            color = "Cluster", 
                            palette=col.cluster,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = paste("MOC", 1:3, sep=""),
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()



ggsurv<-ggsurvplot(model.sbgrp, palette= col.group, 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cln.order$Group %>% sort() %>% unique,
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.sbgrp, 
                            data = cli.data, 
                            color = "Group", 
                            palette=col.group,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cln.order$Group %>% sort() %>% unique,
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()



#-----------------Some Table--------------------#

cln.order<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep = ";")
cln.order$FU_tod_latestControl_new<-ifelse(cln.order$FU_metastasis_latestControl>cln.order$FU_tod_latestControl,
                                           cln.order$FU_metastasis_latestControl,
                                           cln.order$FU_tod_latestControl)

cln.order$metastasis<-ifelse(cln.order$metastasis%in%c(NA, "0"),0, 1)
cln.order$metastasis[cln.order$newID %>% is.na()]<-NA

table(cln.order$Status, cln.order$Cluster, useNA="always")

6+23+15+1+7+15+9+1+5

table(cln.order$metastasis, cln.order$Status, useNA = "always")
cln.order$metastasis %>% table(useNA = "always")
1/60

table(cln.order$relasp_status, cln.order$Status, useNA = "always")
cln.order$relasp_status %>% table(useNA = "always")
table(cln.order$relasp_status, cln.order$Status, useNA = "always")[2,]*100/23

table(cln.order$relasp_status,cln.order$metastasis,  useNA = "always")




#---------------------DFS-------------------------#

cln.order<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep = ";")
cln.order$FU_tod_latestControl_new<-ifelse(cln.order$FU_metastasis_latestControl>cln.order$FU_tod_latestControl,
                                           cln.order$FU_metastasis_latestControl,
                                           cln.order$FU_tod_latestControl)

cln.order$metastasis<-ifelse(cln.order$metastasis%in%c(NA, "0"),0, 1)
cln.order$metastasis[cln.order$newID %>% is.na()]<-NA

cln.order$BC_Meta<-ifelse(cln.order$metastasis==0, cln.order$Status, 1)

cln.order$DFS_status<-ifelse(cln.order$relasp_status!=1, cln.order$BC_Meta , 1)
cln.order$DFS_status<-ifelse(cln.order$DFS_status%in%c(2:4), NA, cln.order$DFS_status)

cln.order$FU_DFS<-cbind(#cln.order$FU_metastasis_latestControl,
  cln.order$FU_tod_latestControl_new,
  cln.order$FU_relapse_latestControl, 
) %>% apply(MARGIN = 1, FUN = min)

#cln.order<-cln.order[!cln.order$newID %>% is.na(),]

cli.data<-cln.order[,c("oldID", "relasp_status", "FU_relapse_latestControl",  "Group", "Cluster")]
names(cli.data)<-c("ID", "Status", "Time", "Group", "Cluster")

cli.data$Status<-ifelse(cli.data$Status%in%c(1), 1, ifelse(cli.data$Status != 0, NA, 0))

model.cluster<-survfit(Surv(Time, Status)~Cluster, data=cli.data)
model.sbgrp<-survfit(Surv(Time, Status)~Group, data=cli.data)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 



ggsurv<-ggsurvplot(model.cluster, palette= col.cluster, 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = paste("MOC", 1:3, sep=""),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.cluster, 
                            data = cli.data, 
                            color = "Cluster", 
                            palette=col.cluster,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = paste("MOC", 1:3, sep=""),
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()



ggsurv<-ggsurvplot(model.sbgrp, palette= col.group, 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cln.order$Group %>% sort() %>% unique,
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.sbgrp, 
                            data = cli.data, 
                            color = "Group", 
                            palette=col.group,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cln.order$Group %>% sort() %>% unique,
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()

cli.data1<-cli.data[cli.data$Group%in%c("LumA", "LumB"),]
cli.data1$Status<-ifelse(cli.data1$Status%in%c(1), 1, ifelse(cli.data$Status != 0, NA, 0))

model.cluster<-survfit(Surv(Time, Status)~Cluster, data=cli.data1)
model.sbgrp<-survfit(Surv(Time, Status)~Group, data=cli.data1)



theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 



ggsurv<-ggsurvplot(model.cluster, palette= col.cluster[c(2,3)], 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cli.data1$Cluster %>% sort() %>% unique,
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.cluster, 
                            data = cli.data, 
                            color = "Cluster", 
                            palette=col.cluster[c(2,3)],
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cli.data1$Cluster %>% sort() %>% unique,
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()



ggsurv<-ggsurvplot(model.sbgrp, palette= col.group[c(3,4)], 
                   size.lab=1.5, 
                   pval = T, 
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cli.data1$Group %>% sort() %>% unique(),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.sbgrp, 
                            data = cli.data, 
                            color = "Group", 
                            palette=col.group[c(3,4)],
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cli.data1$Group %>% sort() %>% unique,
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            clip="off",
                            xlim = c(0, 200)
                            
) 
ggsurv %>% print()




























clinic<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/data_clinical_patient.txt", header = T, sep="\t")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val_cluster.txt", header = T)

merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  merge.clinic$SEX=="Female" 
  & merge.clinic$HISTORY_NEOADJUVANT_TRTYN=="No"
),]

#table----#

clinic1<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/data_clinical_sample.txt", header = T, sep="\t")

table(merge.clinic$PATH_M_STAGE, merge.clinic$cluster, useNA="always")

c(7,255,62)*100/324


merge.clinic %>% nrow()
merge.clinic$RADIATION_THERAPY %>% table()

which(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$PERSON_NEOPLASM_CANCER_STATUS=="With Tumor")
merge.clinic$OS_STATUS %>% table()
merge.clinic$Status<-ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="1:DEAD WITH TUMOR", 1,
                            ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="0:ALIVE OR DEAD TUMOR FREE", 2,
                                   ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="", 3,0)))
merge.clinic$OS_STATUS<-ifelse(merge.clinic$OS_STATUS == "1:DECEASED", 1, 0)
merge.clinic$SUBTYPE<-gsub("BRCA_", "", merge.clinic$SUBTYPE)
merge.clinic$SUBTYPE<-gsub("Her2", "HER2", merge.clinic$SUBTYPE)
merge.clinic$SUBTYPE<-ifelse(merge.clinic$SUBTYPE=="", NA, merge.clinic$SUBTYPE)

merge.clinic$CANCER_TYPE_ACRONYM
#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="White"),]

merge.clinic$DFS_STATUS<-ifelse(merge.clinic$DFS_STATUS == "0:DiseaseFree", 0,
                                ifelse(merge.clinic$DFS_STATUS == "1:Recurred/Progressed", 1, NA))
merge.clinic$Status %>% table()
#----------survival-------------
col.cluster<-c("blue", "red3", "green4", "black", "orange")
col.group<-color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3")


cli.data<-merge.clinic[,c("PATIENT_ID", "OS_STATUS",  "OS_MONTHS", "cluster", "SUBTYPE")] #DSS status for only BC, OS status for all cause

#cli.data$Status<-ifelse(cli.data$Status%in%c(2:3), NA,cli.data$Status)
names(cli.data)<-c("ID", "Status", "Time", "Cluster", "Subtype")

#cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()
cli.data$Time<-as.numeric(cli.data$Time)


end.point<-200
cli.data1 <- cli.data %>% 
  mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
         Status2 = ifelse(Time >= end.point, Status, Status))

model.srv<- survfit(Surv(Time2, Status2)~Cluster, data=cli.data1)
model.srv.sub<- survfit(Surv(Time2, Status2)~Subtype, data=cli.data1)

surv_pvalue(model.srv, data=cli.data1)$pval.txt
#model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 

p.value <- surv_pvalue(fit = model.srv, data = cli.data1)$pval
ggsurv<-ggsurvplot(model.srv, palette= col.cluster, 
                   size.lab=1.5, 
                   xlim=c(0,150),
                   #pval = T,
                   pval = paste("p =", round(surv_pvalue(fit = model.srv, data = cli.data1)$pval, digits = 3)),
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = paste("MOC", 1:3, sep=""),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.srv,
                            xlim=c(0,150),
                            data = cli.data, 
                            color = "Cluster", 
                            palette=c("blue", "red3", "green4"),
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = paste("MOC", 1:3, sep=""),
                            #ggtheme = theme_survminer(),
                            fontsize = 6
                            
)
ggsurv %>% print()


ggsurv<-ggsurvplot(model.srv.sub, palette= col.group, 
                   size.lab=1.5, 
                   xlim=c(0,150),
                   #pval = T,
                   pval = paste("p =", round(surv_pvalue(fit = model.srv.sub, data = cli.data1)$pval, digits = 3)),
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cli.data1$Subtype %>% unique() %>% sort(),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
                   
)

ggsurv$table <- ggrisktable(model.srv.sub, 
                            data = cli.data1, 
                            color = "Subtype", 
                            palette=col.group,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cli.data1$Subtype %>% unique() %>% sort(),
                            #ggtheme = theme_survminer(),
                            fontsize = 6,
                            xlim=c(0,150),
                            
)
ggsurv %>% print()








clinic<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/data_clinical_patient.txt", header = T, sep="\t")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/val_cluster.txt", header = T)
clinic1<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/data_clinical_sample.txt", header = T, sep="\t")


merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")
merge.clinic$OS_MONTHS<-gsub(",", "\\.", merge.clinic$OS_MONTHS) %>% as.numeric()

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  merge.clinic$SEX=="Female"
  &merge.clinic$CHEMOTHERAPY=="NO"###
  &merge.clinic$RADIO_THERAPY=="NO"###
),]

a<-merge(merge.clinic, clinic1, "PATIENT_ID", all.x=T)
a<-a[a$cluster==3,]
dim(a)
ifelse(c(a$ER_STATUS=="Positive"|a$PR_STATUS=="Positive")&a$HER2_STATUS=="Negative","+","-") %>% table()

271/311

merge.clinic.new<-merge(merge.clinic, clinic1, by="PATIENT_ID", all.x = T)

merge.clinic.new$HR<-ifelse(c(merge.clinic.new$ER_STATUS=="Positive"|merge.clinic.new$PR_STATUS=="Positive"),"+","-")

merge.clinic.new$TNBC<-ifelse(merge.clinic.new$HER2_STATUS=="Negative"&merge.clinic.new$HR=="-", "-", "+")

table(merge.clinic.new$TNBC, merge.clinic.new$cluster, useNA="always")

table(merge.clinic.new$HR, useNA="always")

c(600, 94, 0, 63, 631, 0, 64, 630, 0)*100/694


merge.clinic$OS_STATUS<-ifelse(merge.clinic$OS_STATUS == "1:DECEASED", 1, 0)
merge.clinic$DSS_STATUS <- ifelse(merge.clinic$VITAL_STATUS == "Died of Disease", 1, 
                                  ifelse(merge.clinic$VITAL_STATUS == "Died of Other Causes", NA, 0))
merge.clinic$RFS_STATUS <- ifelse(merge.clinic$RFS_STATUS == "0:Not Recurred", 0, 
                                  ifelse(merge.clinic$RFS_STATUS == "1:Recurred", 1, NA))

merge.clinic$CLAUDIN_SUBTYPE<-ifelse(merge.clinic$CLAUDIN_SUBTYPE%in%c("claudin-low","NC"), NA, merge.clinic$CLAUDIN_SUBTYPE)
merge.clinic$CLAUDIN_SUBTYPE<-gsub("Her2", "HER2", merge.clinic$CLAUDIN_SUBTYPE)

#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="White"),]
merge.clinic$RFS_MONTHS<-gsub(",", "\\.", merge.clinic$RFS_MONTHS) %>% as.numeric()
#----------survival-------------

cli.data<-merge.clinic[,c("PATIENT_ID", "OS_STATUS",  "OS_MONTHS", "cluster", "CLAUDIN_SUBTYPE")] 

names(cli.data)<-c("ID", "Status", "Time", "Cluster", "Subtype")

#cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()
cli.data$Time<-as.numeric(cli.data$Time)

end.point<-200
cli.data1 <- cli.data %>% 
  mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
         Status2 = ifelse(Time >= end.point, Status, Status))

model.srv<- survfit(Surv(Time2, Status2)~Cluster, data=cli.data1)
model.srv.sub<-survfit(Surv(Time2, Status2)~Subtype, data=cli.data1)

surv_pvalue(model.srv, data=cli.data1)$pval.txt
#model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "grey90"),
               panel.grid.minor = element_line(colour = "grey90"),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.text = element_text(size = 20, color = "black", face = "bold"),
               legend.key.size = unit(1.5, 'cm'),
               axis.text.y = element_text(size=20)) 

ggsurv<-ggsurvplot(model.srv, palette= col.cluster, 
                   size.lab=1.5, 
                   #pval = T,
                   pval = paste("p =", round(surv_pvalue(fit = model.srv, data = cli.data1)$pval, digits = 3)),
                   xlim=c(0,200),
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = paste("MOC", 1:3, sep=""),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
                   
                   
)

ggsurv$table <- ggrisktable(model.srv, 
                            data = cli.data, 
                            color = "Cluster", 
                            palette=c("blue", "red3", "green4"),
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = paste("MOC", 1:3, sep=""),
                            #ggtheme = theme_survminer(),
                            fontsize = 6
                            
)
ggsurv %>% print()

ggsurv<-ggsurvplot(model.srv.sub, palette= col.group, 
                   size.lab=1.5, 
                   #pval = T,
                   pval = paste("p =", round(surv_pvalue(fit = model.srv.sub, data = cli.data1)$pval, digits = 3)),
                   xlim=c(0,200),
                   font.size=1.5,
                   surv.scale = "percent",
                   font.legend = list(size = 15, color = "black", face = "bold"),
                   font.tickslab = c(20, "plain", "black"),
                   risk.table = T, 
                   pval.method = F,
                   legend.title="",
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   pval.size = 7, pval.coord = c(60, .1),
                   conf.int = T,
                   legend.labs = cli.data1$Subtype %>% unique() %>% sort(),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.05,
                   xlab="Months"
)

ggsurv$table <- ggrisktable(model.srv.sub, 
                            data = cli.data1, 
                            color = "Subtype", 
                            palette=col.group,
                            y.text = T,   ylab = "",  xlab = "",
                            tables.theme = theme,
                            font.tickslab = c(15, "bold"),
                            legend.labs = cli.data1$Subtype %>% unique() %>% sort(),
                            #ggtheme = theme_survminer(),
                            fontsize = 6
                            
)
ggsurv %>% print()



#---------------table Update------------------------#

cln.order<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep = ";")


table(cln.order$Status, cln.order$Cluster, useNA = "always")

cbind(table(cln.order$Status, cln.order$Cluster, useNA = "always")[,1]*100/50,
      table(cln.order$Status, cln.order$Cluster, useNA = "always")[,2]*100/137,
      table(cln.order$Status, cln.order$Cluster, useNA = "always")[,3]*100/148)


cln.order$Status %>% table(useNA = "always")*100/335



clinic<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/data_clinical_patient.txt", header = T, sep="\t")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val_cluster.txt", header = T)

merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  merge.clinic$SEX=="Female" 
  & merge.clinic$HISTORY_NEOADJUVANT_TRTYN=="No"
  #merge.clinic$RADIATION_THERAPY=="No"
),]


merge.clinic$OS_STATUS %>% table(useNA = "always")

merge.clinic$DSS_STATUS %>% table(useNA = "always")

paste(merge.clinic$OS_STATUS, merge.clinic$DSS_STATUS) %>% table()

merge.clinic$Status<-ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="1:DEAD WITH TUMOR", 1,
                            ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="0:ALIVE OR DEAD TUMOR FREE", 2,
                                   ifelse(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$DSS_STATUS=="", 3,0)))


table(merge.clinic$Status, useNA="always")
table(merge.clinic$Status, merge.clinic$cluster, useNA="always")


