
setwd("********/Oslo2/data/")

set.seed(123)

library(MOFA2)
library(readxl)
library(ggplot2)
library(gplots)
library(cluster)
library(tidyverse)
library(survival)
library(dplyr)
library(corrplot)
library(siggenes)
library(RecordLinkage)
library(DMwR2)
library(e1071)
library(Matrix)
library(magrittr)
library(randomForest)
library(caret)
library(pROC)
library(factoextra)
library(faraway)
library(survminer)
library(dplyr)


#subgroup
data.subg<-as.data.frame(read_excel(dir()[grep("Intrinsic", dir())], sheet = 1))

#clinical data
data.cln<-as.data.frame(read_excel(dir()[grep("Clinic", dir())], sheet = 1))

#rppa
data.rppa<-as.data.frame(read_excel(dir()[grep("RPPAgr", dir())], sheet = 1))

#metabolomic
data.metab<-as.data.frame(read_excel(dir()[grep("metabo",dir())], sheet = 1))

#protein
data.prot<-as.data.frame(read_excel(dir()[grep("protein",dir())], sheet = 1))

#gene data
data.gene<-as.data.frame(read_excel(dir()[grep("Gene",dir())], sheet = 1, col_names =F)) #"C:\\Users\\abhibhas\\AppData\\Roaming\\MobaXterm\\slash\\RemoteFiles\\2231052_4_3\\GeneExpressionData.xlsx"
data.gene<-as.matrix(t(data.gene))%>%as.data.frame()
colnames(data.gene)<-data.gene[1,]
data.gene<-data.gene[-1,]

max.sample<-max(c(data.gene$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric(),
                  data.subg$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric(),
                  data.rppa$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric(),
                  data.metab$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric(),
                  data.prot$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric(),
                  data.cln$newID%>%gsub(pattern = "BC", replacement = "")%>%as.numeric()))

index<-cbind(paste("BC", 1:max.sample, sep = ""), c(1:max.sample))%>%as.data.frame()
colnames(index)<-c("newID", "index")

#put all data frames into list
df_list <- list(index, data.subg, data.cln, data.rppa, data.metab, data.prot, data.gene)

#merge all data frames in list
merge.data<-df_list %>% purrr::reduce(full_join, by='newID')

duplicate<-replace(merge.data$newID, merge.data$newID%>%duplicated(), NA)%>%is.na()%>%which() #duplicate entries
merge.data<-merge.data[-duplicate,]

#------------All data list----------------#
m.data.gene<-apply(merge.data[,colnames(data.gene)[-1]], 2, as.numeric)%>%t()
m.data.prot<-apply(merge.data[,colnames(data.prot)[-1]], 2, as.numeric)%>%t()
m.data.metab<-apply(merge.data[,colnames(data.metab)[-c(1:2)]], 2, as.numeric)%>%t()
m.data.cln<-merge.data[,c(colnames(data.subg), colnames(data.cln), colnames(data.rppa))]
metab.cluster<-merge.data[,colnames(data.metab)[2]]

#----------cleaning----------------
rm(data.subg)
rm(data.cln)
rm(data.rppa)
rm(data.metab)
rm(data.prot)
rm(merge.data)
rm(df_list)

#write.csv(m.data.gene, ".\\functional_analysis\\merge_gene.csv")
#write.csv(m.data.metab, ".\\functional_analysis\\merge_metab.csv")
#write.csv(m.data.prot, ".\\functional_analysis\\merge_prot.csv")

#write.csv(t(data[[3]]), ".\\ml_data/Gene.csv", row.names = F)


#------------------------Gene Filtering---------------------#

#data.gene<-as.data.frame(t(m.data.gene))
#data.gene$label<-m.data.cln$Grade%>%as.numeric()%>%as.factor()
#write.csv(data.gene, "geneforfeatureselection.csv", row.names = F)

#data.gene<-na.omit(data.gene)

#------------Functions----------

welch<-function(x){
  a<-t.test(x[which(data.gene$label %in% c(1,2))], x[which(data.gene$label==3)], var.equal = T)#welch t tailed test
  return(a$p.value)
}                                                                 

wilcox<-function(x){
  a<-wilcox.test(x[which(data.gene$label %in% c(1))], x[which(data.gene$label==3)], var.equal = T)#wilcoxn run tailed test
  return(a$p.value)
}

stnd<-function(x){
  (x-mean(x, na.rm=T))/sd(x, na.rm=T)
}
m.data.metab<-apply(m.data.metab, 2, stnd)

minmax<-function(x){
  (x-min(x))/(max(x)-min(x))
}
#-----------------naming the DATA list-----------------------
colnames(m.data.gene)<-colnames(m.data.metab)<-colnames(m.data.prot)<-m.data.cln$newID
data<-list(m.data.metab, m.data.prot, m.data.gene)

names(data)<-c("metabolite", "protein", "Gene")

#for(i in 1:3){
# write.csv(t(data[[i]]), paste(names(data)[i], ".csv", sep = ""), row.names = F)
#}

#-------------------MOFA Simulatio--------------------
condition<-F
if(condition){
  #--------------MOFA---------------------------------
  num_fac<-c(10, 15, 20, 30)
  
  for(i in seq_along(num_fac)){
    
    MOFAobject <- create_mofa(data)
    
    plot_data_overview(MOFAobject)
    
    data_opts <- get_default_data_options(MOFAobject)
    head(data_opts)
    
    model_opts <- get_default_model_options(MOFAobject)
    head(model_opts)
    
    model_opts$num_factors<-num_fac[i]
    
    train_opts <- get_default_training_options(MOFAobject)
    head(train_opts)
    
    MOFAobject <- prepare_mofa(
      object = MOFAobject,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    outfile = file.path(getwd(),paste("MOFA_allgene_", as.character(num_fac[i]),"F.hdf5", sep=""))
    MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
    
  }
}

