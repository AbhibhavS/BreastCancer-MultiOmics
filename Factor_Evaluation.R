set.seed(123)

library(MOFA2)
library(ggplot2)
library(gplots)
library(tidyverse)
library(corrplot)
library(faraway)
library(survival)
library(survminer)


model <- load_model("./MOFA_ALL/MOFA_allgene_20F.hdf5")
factor<-read.csv("./Factor_evaluation/FA_factor20.csv", header = T)
clinic1<-read.csv("./Factor_evaluation/FA_clinic.csv", header = T)
clinic<-read.csv("/mnt/work/Oslo2/data/Clinic_data_oslo2.csv", header = T, sep=";")

clinic$FU_tod_latestControl_new<-ifelse(clinic$FU_metastasis_latestControl>clinic$FU_tod_latestControl,
                                        clinic$FU_metastasis_latestControl,
                                        clinic$FU_tod_latestControl)

clinic$metastasis<-ifelse(clinic$metastasis%in%c(NA, "0"),0, 1)
clinic$metastasis[clinic$newID %>% is.na()]<-NA

clinic$BC_Meta<-ifelse(clinic$metastasis==0, clinic$Status, 1)

#cln.order<-cln.order[!cln.order$newID %>% is.na(),]

factor.clinic<-data.frame(oldID= clinic1$newID,
                          factor)

minmax<-function(x){
  (x-min(x, na.rm = T))/(max(x, na.rm=T)-min(x, na.rm = T))
}

#factor.clinic[,-1]<-apply(factor.clinic[,-1], 2, minmax)

clinic.for.surv<-data.frame(oldID = clinic$oldID,
                            FU = clinic$FU_metastasis_latestControl,
                            status = clinic$BC_Meta,
                            age = clinic$Age,
                            BMI = clinic$BMI,
                            Tumorsize = clinic$TumorSize,
                            Grade = clinic$Grade,
                            N_status = clinic$N_status
)

clinic.new<-merge(factor.clinic, clinic.for.surv, by.x = "oldID", all = T)


my.data<-cbind(#clinic$status..0.alive..1.dead., 
  ifelse(clinic.new$status%in%c(2:4), NA, clinic.new$status),
  #ifelse(clinic.new$status%in%c(2:4), NA, clinic.new$status),
  clinic.new$FU %>% as.numeric(),
  #clinic$Age %>% as.numeric(),
  #clinic$BMI,
  #clinic$TumorSize %>% as.numeric(),
  #clinic$Grade,
  #clinic$N_status,
  clinic.new[,2:21])

colnames(my.data)<-c("status", "time", 
                     #"Age", 
                     #"BMI",
                     #"Tumoursize",
                     #"Grade",
                     #"Node_Status",
                     paste("Factor", 1:ncol(factor), sep = ""))

#my.data$status<-ifelse(clinic$Cause.of.death..1.breast.cancer..2.different.disease..3.unknown..4.cause.of.death.not.updated..dead.after.des2017.==1,1,0)

#-------------------Multi Variate Cox----------------------------#
res.cox.crude <- coxph(Surv(time, status) ~ ., data =  my.data, ties = "efron")
res.cox.crude %>% summary()
ggforest(res.cox.crude, fontsize = 1.5) %>% print()

#test.ph <- cox.zph(res.cox.crude)
#ggcoxzph(test.ph)

#-----------Age Adjusted

my.data<-cbind(#clinic$status..0.alive..1.dead., 
  ifelse(clinic.new$status%in%c(2:4), NA, clinic.new$status),
  #ifelse(clinic.new$status%in%c(2:4), NA, clinic.new$status),
  clinic.new$FU %>% as.numeric(),
  clinic.new$age %>% as.numeric(),
  #clinic.new$BMI,
  #clinic$TumorSize %>% as.numeric(),
  #clinic$Grade,
  clinic.new[,2:21])

colnames(my.data)<-c("status", "time", 
                     "Age", 
                     #"BMI",
                     #"grade",
                     #"Tumoursize",
                     paste("Factor", 1:ncol(factor), sep = ""))

#my.data$status<-ifelse(clinic$Cause.of.death..1.breast.cancer..2.different.disease..3.unknown..4.cause.of.death.not.updated..dead.after.des2017.==1,1,0)

#-------------------Multi Variate Cox----------------------------#
res.cox.age <- coxph(Surv(time, status) ~ ., data =  my.data, ties = "efron")
summary(res.cox.age)

ggforest(res.cox.age, fontsize = 1.5) %>% print()


#-----------All Adjusted
my.data<-cbind(#clinic$status..0.alive..1.dead., 
  #ifelse(clinic.new$status%in%c(NA), NA, 1),
  ifelse(clinic.new$status%in%c(2:4), NA, clinic.new$status),
  clinic.new$FU %>% as.numeric(),
  clinic.new$age %>% as.numeric(),
  clinic.new$BMI,
  clinic.new$Grade %>% as.numeric(),
  ifelse(clinic.new$N_status%in%c("pN0"),0,
         ifelse(clinic.new$N_status%in%c("pN2"),2,
                ifelse(clinic.new$N_status%in%c("pN3"),3,
                       ifelse(clinic.new$N_status%in%c(NA),NA,1)))),
  clinic.new$Tumorsize %>% as.numeric(),
  
  clinic.new[,2:21])

colnames(my.data)<-c("status", "time", 
                     "Age", 
                     "BMI",
                     "Grade",
                     "Node",
                     "Tumor Size",
                     paste("Factor", 1:ncol(factor), sep = ""))
#my.data$status<-ifelse(clinic$Cause.of.death..1.breast.cancer..2.different.disease..3.unknown..4.cause.of.death.not.updated..dead.after.des2017.==1,1,0)

#-------------------Multi Variate Cox----------------------------#
res.cox.age <- coxph(Surv(time, status) ~ ., data =  my.data, ties = "efron")
summary(res.cox.age)

ggforest(res.cox.age, fontsize = 1) %>% print()


#test.ph <- cox.zph(res.cox.age)
#ggcoxzph(test.ph)
#-------------Univariate Survival------------------#

#----------------------------------------------------#

clinic<-read.csv("./Factor_evaluation/FA_clinic.csv", header = T)
my.clinic<-clinic[,c(10,11,13,14,16, 17:19, 23,9)] 
colnames(my.clinic)<-c("Histology",
                       "Tumorsize_mm", 
                       "grade","n_status",
                       "M_status",
                       "HER","ER", "PR",
                       "BMI",
                       "Age") 

c(clinic$Follow.up.time..time.to.death.or.latest.control. %>% as.numeric() /12) %>% summary()
clinic$status..0.alive..1.dead. %>% table()

my.clinic$Histology<-ifelse(my.clinic$Histology == "Ductal", 1, 2)
my.clinic$Tumorsize_mm<-my.clinic$Tumorsize_mm %>% as.numeric()
my.clinic$grade<-my.clinic$grade %>% as.numeric()
my.clinic$n_status<-my.clinic$n_status %>% as.factor() %>% as.numeric()
my.clinic$M_status<-my.clinic$M_status %>% as.factor() %>% as.numeric()
my.clinic$HER<-my.clinic$HER %>% as.factor() %>% as.numeric()
my.clinic$ER<-my.clinic$ER %>% as.factor() %>% as.numeric()
my.clinic$PR<-my.clinic$PR %>% as.factor() %>% as.numeric()
my.clinic$BMI<-my.clinic$BMI %>% as.numeric()

corr<-NULL; pval<-NULL
for(i in 1:ncol(my.clinic)){
  for(j in 1:20){
    a<-cor.test(factor[,j], my.clinic[,i], method = "spearman", use = "complete.obs")
    corr<-c(corr, a$estimate)
    pval<-c(pval, a$p.value)
  }
}


corr<-matrix(corr, byrow = F, nrow=20)
pval<-matrix(pval, byrow = F, nrow=20)

row.names(corr)<-paste("Factor", 1:20, sep="")
colnames(corr)<-gsub("_","-",colnames(my.clinic))

row.names(pval)<-paste("Factor", 1:20, sep="")
colnames(pval)<-gsub("_","-",colnames(my.clinic))

corrplot(corr=t(corr), p.mat = t(pval), col= colorRampPalette(c("pink","orange","red","black"))(n = 4)
         ,bg = 'lemonchiffon1',tl.col = 'black',tl.cex =1.3)




clinic<-read.csv("./Factor_evaluation/FA_clinic.csv", header = T)
my.clinic<-clinic[,c(10,11,13,14,16, 17:19, 23)] 
colnames(my.clinic)<-c("Histology",
                       "Tumorsize_mm", 
                       "grade","n_status",
                       "M_status",
                       "HER","ER", "PR",
                       "BMI") 

#---------Histology
pval.Histology<-NULL

histo<-ifelse(my.clinic$Histology == "Ductal", 1, 2)

for(i in 1:ncol(factor)){
  w<-wilcox.test(factor[,i][which(my.clinic$Histology == "Ductal")], 
                 factor[,i][which(my.clinic$Histology != "Ductal")])
  pval.Histology<-c(pval.Histology, w$p.value)
}
which(pval.Histology <= 0.001)

#--------Tumor Size
my.clinic$Tumorsize_mm <- my.clinic$Tumorsize_mm %>% as.numeric()

Tumorsize_mm<-ifelse(my.clinic$Tumorsize_mm %in% c(1:ncol(factor)), 1, 
                     ifelse(my.clinic$Tumorsize_mm %in% c(21:50), 2,
                            ifelse(my.clinic$Tumorsize_mm > 50, 2, NA)))
pval.Tumor_Size<-NULL
for(i in 1:ncol(factor)){
  w<-wilcox.test(factor[,i][which(Tumorsize_mm == 1)], 
                 factor[,i][which(Tumorsize_mm != 1)])
  pval.Tumor_Size<-c(pval.Tumor_Size, w$p.value)
}
which(pval.Tumor_Size <= 0.01)

#------Grade
pval.Grade<-NULL
for(i in 1:ncol(factor)){
  k<-kruskal.test(factor[,i], 
                  my.clinic$grade %>% as.factor())
  pval.Grade<-c(pval.Grade, k$p.value)
}
which(pval.Grade <= 0.05)


#------Node Status
node_status<-my.clinic$n_status
node_status<-gsub("\\(|a|mi|)", "", node_status)
node_status<-gsub("3", "2", node_status)

pval.NS<-NULL
for(i in 1:ncol(factor)){
  k<-kruskal.test(factor[,i], 
                  node_status %>% as.factor())
  pval.NS<-c(pval.NS, k$p.value)
}
which(pval.NS <= 0.01)


#------ER
pvalER<-NULL
for(i in 1:ncol(factor)){
  w<-wilcox.test(factor[,i][which(my.clinic$ER == "pos")], 
                 factor[,i][which(my.clinic$ER == "neg")])
  pvalER<-c(pvalER, w$p.value)
}
which(pvalER <= 0.05)


#------PR
pvalPR<-NULL
for(i in 1:ncol(factor)){
  w<-wilcox.test(factor[,i][which(my.clinic$PR == "pos")], 
                 factor[,i][which(my.clinic$PR == "neg")])
  pvalPR<-c(pvalPR, w$p.value)
}
which(pvalPR <= 0.05)

#------HER
pvalHER<-NULL
for(i in 1:ncol(factor)){
  w<-wilcox.test(factor[,i][which(my.clinic$HER == "pos")], 
                 factor[,i][which(my.clinic$HER == "neg")])
  pvalHER<-c(pvalHER, w$p.value)
}
which(pvalHER <= 0.05)


#------BMI
BMI_cat <- ifelse(my.clinic$BMI >= 30, 1, 
                  ifelse( my.clinic$BMI < 30 & my.clinic$BMI >= 25, 2, 
                          ifelse(my.clinic$BMI < 25 & my.clinic$BMI >= 18.5, 3,
                                 ifelse(my.clinic$BMI < 18, 4, NA)))) %>% as.factor()

pval.BMI<-NULL
for(i in 1:ncol(factor)){
  k<-kruskal.test(factor[,i], 
                  BMI_cat %>% as.factor())
  pval.BMI<-c(pval.BMI, k$p.value)
}
which(pval.BMI <= 0.05)


ptable<-cbind(pval.Histology,
              pval.Tumor_Size,
              pval.Grade,
              pval.NS,
              pvalER, 
              pvalPR, 
              pvalHER)

sample_metadata <- data.frame(sample = samples_names(model)[[1]],
                              Histology = histo %>% as.factor(),
                              Tumorsize = Tumorsize_mm,
                              grade = my.clinic$grade %>% as.numeric(),
                              Nodestatus = node_status,
                              ER = my.clinic$ER,
                              PR = my.clinic$PR,
                              HER2 = my.clinic$HER,
                              Fact = factor
)


samples_metadata(model) <- sample_metadata

integer0_test <- function(data) {
  
  if(identical(data, integer(0))) {
    return(T)
  }
  else {return(F)}
}

for(i in 1:7){
  
  if(which(ptable[,i] < 0.05) %>% integer0_test){next}
  
  p <- plot_factor(model, 
                   factors = c(which(ptable[,i] < 0.05)),
                   color_by = colnames(sample_metadata)[-1][i],
                   dot_size = 3,        # change dot size
                   dodge = T,           # dodge points with different colors
                   legend = T,          # remove legend
                   add_violin = T,      # add violin plots,
                   violin_alpha = 0.25, # transparency of violin plots
                   show_missing = F,
                   scale = F,
                   color_name = colnames(sample_metadata)[-1][i]
                   
  )
  print(p)
}


pval.mat<-cbind(pval.BMI,
                pval.Grade,
                pval.Histology,
                pval.NS,
                pval.Tumor_Size,
                pvalER,
                pvalHER,
                pvalPR)

pval.mat<-ifelse(pval.mat<=0.01,4, 
                 ifelse(pval.mat<=0.05 & pval.mat>0.01,3, 
                        ifelse(pval.mat<=0.1 & pval.mat>0.05,2,1)))

rownames(pval.mat)<-paste("Factor", 1:ncol(factor))
colnames(pval.mat)<-c("BMI", "Grade", "Histology", "Node Status",
                      "Tumor Size", "ER", "HER", "PR")


#write.table(pval.mat %>% as.data.frame(), "sup1_tab1_A&N.csv",row.names = FALSE, sep=",")


corrplot(pval.mat %>% t(), is.corr = F, 
         col= colorRampPalette(c("pink","orange","red","black"))(n = 4)
         ,bg = 'lemonchiffon1',tl.col = 'black',tl.cex =1.3)







