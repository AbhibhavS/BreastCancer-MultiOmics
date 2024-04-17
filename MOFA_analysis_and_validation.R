setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/")

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



#--------------------------------------------------------------#


#-----Load Trained Model-----------------

model <- load_model("MOFA_allgene_20F.hdf5")
factors <- get_factors(model, as.data.frame = T) # factor matrix
factor_mat<-matrix(factors$value, ncol=length(unique(factors$factor)), byrow = F)
plot_data_overview(model, colors = c("darkorange3","goldenrod","plum2"), )

plot_variance_explained(model, x="view", y="factor")
#write.csv(factor_mat, "FA_factor30.csv", row.names = F)


#-------------------MODEL Properties Plots-----------------
plot.condition<-T

if(plot.condition){
  
  head(model@cache$variance_explained$r2_total[[1]])
  p<-barplot(model@cache$variance_explained$r2_total[[1]],
             col = "skyblue", ylim = c(0,60),
             main = "Variance explained per view")
  text(x = p, 
       y = c(as.numeric(model@cache$variance_explained$r2_total[[1]]))+ 2, 
       labels = round(model@cache$variance_explained$r2_total[[1]],2))
  
  #---------------------#
  var.fac<-model@cache$variance_explained$r2_per_factor[[1]]
  model@cache$variance_explained$r2_total$group1 %>% sum()
  
  var_exp_per_view<-cbind(var.fac, rowSums(var.fac))
  sum(rowSums(var.fac)/sum(rowSums(var.fac)))
  
  
  
  #----------------------#
  plot(c(15.84, 21, 22.29, 23.14), xlab = "factor_number",
       xlim=c(1,4), ylim=c(13,55), yaxt="n", ylab="Total Variance explained",
       type = "o", col="red", lwd=2, xaxt="n")
  points(c(29.18, 34.39, 37.41, 40), type="o", col="blue", lwd=2)
  points(c(38.62, 44.4, 47.78, 52.55), type="o", col="black", lwd=2)
  
  axis(1, labels = c(10,15,20,30), at=1:4)
  axis(2, labels = seq(13,55,3), at=seq(13,55,3))
  grid()
  legend(legend = c("Metabolomic", "Protein", "Gene"),
         col = c("Red", "Blue", "black"),
         "topleft",
         lwd=5)
  
  #----------------------------#
  cor.p<-psych::corr.test(factor_mat %>% as.data.frame(), adjust="fdr")
  
  M<-cor(factor_mat, use = "complete.obs")
  
  corrplot(M,col = c('white', 'black'),
           bg = 'gold2', insig = 'pch', pch = "",
           p.mat = cor.p$p, sig.level = 0.05)
  
  #------------------------------#
  p<-barplot(t(as.matrix(var.fac)),
             beside = T,
             col = c("red","orange","blue"), ylim = c(0,15),
             main = "Variance explained per factor")
  
  legend("topright",
         x = c(40, 50), y = c(1, 12),
         legend =names(data), 
         fill = c("red","orange","blue"),
         text.width = 1,
         cex=2,
         text.font=9,
         bty="n",
         box.lwd=2,
         ncol = 1)
  
  vif(factor_mat)
}


#------Clustering----------------


silhouette_score <- function(k){
  km <- kmeans(factor_mat[,c(1,2,13)], centers = k, nstart=10)
  ss <- silhouette(km$cluster, dist(factor_mat[,c(2,13)]))
  mean(ss[, 3])
}
k <- 2:10
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters',
     ylab='Average Silhouette Scores', frame=FALSE)


gap_stat <- clusGap(factor_mat[,c(1,2,13)], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


kmc<-kmeans(factor_mat[,c(1,2,13)], 3, nstart = 25)
my.cluster<-kmc$cluster


fviz_cluster(kmc, data = factor_mat[,c(1,2)],
             palette = c("blue","red","green", "black", "orange"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw(),
             xlab = "", ylab = "", main = ""
)

#write.csv(my.cluster, "./ml_data/cluster.csv", row.names = F)
my.cluster<-read.csv("./ml_data/cluster.csv", header = T)[,1]

#----------survival-------------

death_type<-list(c(1:4),1,c(1:4),1,c(1:4))
wrt<-c("Group","Group", "Cluster", "Cluster", "Grade")

survival.formula<-sapply(wrt,
                         function(x) as.formula(paste('Surv(Time, Status)~', x)))
model.srv<-list()
p.val<-c()
for(k in 1:5){
  
  cli.data<-m.data.cln[,c(1,24,27,13,2)]
  cli.data<-cbind(cli.data, my.cluster)
  names(cli.data)<-c("ID", "Status", "Time","Grade", "Group", "Cluster")
  
  bc_death<-m.data.cln$`Cause of death (1=breast cancer, 2=different disease, 3=unknown, 4=cause of death not updated (dead after des2017)`
  bc_death<-ifelse(bc_death %in% death_type[[k]], 1, ifelse(bc_death != 0, NA, 0))
  
  status1<-as.numeric(as.factor(bc_death[order(my.cluster)]))
  cli.data$Status<-bc_death
  
  #cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
  #cli.data<-na.omit(cli.data)
  cli.data$Time<-as.numeric(cli.data$Time)
  
  
  model.srv[[k]]<-survfit(survival.formula[[k]], data=cli.data)
  p.val<-c(p.val, surv_pvalue(survfit(survival.formula[[k]], data=cli.data))$pval.txt)
  
}

for(k in 1:5){
  col.cluster<-c("blue","red","green", "black", "orange")
  col.group<-c("green","red", "Blue","Orange","Purple")
  
  ggsurvplot(model.srv[[k]], palette= col.cluster, 
             size.lab=1.5, 
             pval = p.val[k], 
             font.size=1.5,
             surv.scale = "percent",
             font.legend = list(size = 15, color = "black", face = "bold"),
             font.tickslab = c(20, "plain", "black"),
             risk.table = T, 
             pval.method = T,
             legend.title="",
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             pval.size = 7, pval.coord = c(60, .1),
             conf.int = T,
  ) %>% print()
}



#-------------------Factors propertirs---------------------------#

node_status<-m.data.cln$N.status
node_status<-gsub("\\(|a|mi|)", "", node_status)
node_status<-gsub("3", "2", node_status)

t_status<-m.data.cln$N.status

m.data.cln$M.status %>% table()

sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  status = m.data.cln$`status (0=alive, 1=dead)`,
  grade = m.data.cln$Grade %>% as.numeric(),
  pam50 = m.data.cln$Group,
  n.status = node_status
)

samples_metadata(model) <- sample_metadata

p <- plot_factor(model, 
                 factors = c(1,2,13),
                 color_by = "n.status",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25, # transparency of violin plots
                 show_missing = F
)

# The output of plot_factor is a ggplot2 object that we can edit

print(p)



#---------------------Heat Map---------------
ref<-m.data.cln$newID[order(my.cluster)]

sorted<-function(x){
  x<-t(x)[order(my.cluster),]
}

data.order<-lapply(data, sorted)
#saveRDS(data.order, file = "data_order_all.rds")
cln.order<-m.data.cln[order(my.cluster),]
#saveRDS(cln.order, file = "cln_order_all.rds")



my_palette <-colorRampPalette(c("blue","white",
                                "red"))(n = 5)
#my_palette <-colorspace::sequential_hcl

weights <- get_weights(model, views = "all", factors = "all")


heat.data <- log2(data.order[[1]]/colMeans(data.order[[1]], na.rm = T))

pval<-NULL
for(k in 1:ncol(heat.data)){
  a<-t.test(heat.data[which(my.cluster[order(my.cluster)]==1),k],
            heat.data[which(my.cluster[order(my.cluster)]==2),k],
            var.equal = F)
  pval[k]<-a$p.value
}
heat.data<- heat.data[,order(pval, decreasing = F)]

heat.data<-heat.data[,order(colMeans(heat.data, na.rm = T))]

par(bg="white")
heatmap.2( heat.data, 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col = colorspace::heat_hcl(3,c=c(150,200), l = c(30,90), power = c(1/5, 1.5)),
           
           #col= colorRampPalette(c("violet","blue","green",
           #                       "yellow", "orange","red"))(n = 50),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(heat.data),
           labCol = colnames(heat.data),
           main = "metabolite",
           na.color = "Black")



par(bg="white")
heatmap.2( cbind(sort(my.cluster), 
                 sort(my.cluster)), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue","red", "green"))(n = 3),
           #col=c("blue","red","green", "black", "orange"),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order[[1]]),
           labCol = colnames(data.order[[1]]),
           main = "metabolite")


#-------------------------------------------------------------------#

pam50<-as.numeric(as.factor(m.data.cln$Group[order(my.cluster)]))
heatmap.2( cbind(pam50, pam50), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue",
                                   "red", "skyblue", 
                                   "orange", "purple"))(n = 5),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order[[1]]),
           labCol = colnames(data.order[[1]]),
           main = "PAM50",
           na.color = "black")



rppa<-as.numeric(as.factor(m.data.cln$Named.Group[order(my.cluster)]))

heatmap.2( cbind(rppa, rppa), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue", "red",
                                   "yellow", 
                                   "purple", "skyblue"))(n = 5),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order[[1]]),
           labCol = colnames(data.order[[1]]),
           main = "RPPA",
           na.color = "black")




bc_death<-m.data.cln$`Cause of death (1=breast cancer, 2=different disease, 3=unknown, 4=cause of death not updated (dead after des2017)`
bc_all_death<-ifelse(bc_death %in% c(1:4), 1, ifelse(is.na(bc_death), NA, 0))
bc_1_death<-ifelse(bc_death %in% c(1), 1, ifelse(is.na(bc_death), NA, 0))

bcdeath<-ifelse(bc_death %in% c(1), 1, ifelse(is.na(bc_death), NA, ifelse(bc_death %in% c(2,3,4),2,0)))
bcdeath<-bcdeath[order(my.cluster)]

heatmap.2( cbind(bcdeath,bcdeath), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("black", "red", "green"))(n = 3),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order[[1]]),
           labCol = colnames(data.order),
           main = "death",
           na.color = "white")

status1to4<-as.numeric(as.factor(bc_all_death[order(my.cluster)]))-1
status1<-as.numeric(as.factor(bc_1_death[order(my.cluster)]))
status<-as.numeric(as.factor(bcdeath[order(my.cluster)]))

status1to4[which(my.cluster[order(my.cluster)]==1)]<-mean(status1to4[which(my.cluster[order(my.cluster)]==1)], na.rm=T)
status1to4[which(my.cluster[order(my.cluster)]==2)]<-mean(status1to4[which(my.cluster[order(my.cluster)]==2)], na.rm=T)
status1to4[which(my.cluster[order(my.cluster)]==3)]<-mean(status1to4[which(my.cluster[order(my.cluster)]==3)], na.rm=T)

heatmap.2( cbind(status1to4, status1to4), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("lightblue", "blue", "darkblue"))(n = 10),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order),
           labCol = colnames(data.order),
           main = "death",
           na.color = "white")



ER<-as.numeric(m.data.cln$`ER (1=negative (<1%), 2= 1-10%, 3=10-50%, 4=>=50%, 5=positive, 98 = not defined)`[order(my.cluster)])
ER[ER==98 | ER ==0]<-NA

heatmap.2( cbind(ER, ER), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           #col = colorRampPalette(c("blue","pink"))(n = 2),
           col= colorRampPalette(c("blue","lightblue","white", "pink", "red"))(n = 5),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order),
           labCol = colnames(data.order),
           main = "ER",
           na.color = "black")

PR<-as.numeric(m.data.cln$`PR (1=negative (<1%), 2= 1-10%, 3=10-50%, 4=>=50%, 98 = not defined)`[order(my.cluster)])
PR[PR==98]<-NA

heatmap.2( cbind(PR, PR), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue","lightblue", "pink", "red"))(n = 4),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order),
           labCol = colnames(data.order),
           main = "PR",
           na.color = "black")


her<-as.numeric(m.data.cln$Her2[order(my.cluster)])
heatmap.2( cbind(her, her), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue","white", "red"))(n = 20),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order[[1]]),
           labCol = colnames(data.order),
           main = "HER",
           na.color = "black")


condition.sam<-T

if(condition.sam){
  
  setwd("./functional_analysis")
  set.seed(123)
  
  data.prot<-read.csv("merge_prot.csv", header = T)
  data.gene<-read.csv("merge_gene.csv", header = T)
  data.metab<-read.csv("merge_metab.csv", header = T)
  prot_name<-read.table("prot_name.txt", header = F, sep = "\t", fill = T)
  
  name.prot<-prot_name[,2]
  name.gene<-data.gene[,1]
  name.metab<-data.metab[,1]
  
  data.prot<-data.prot[,-1]
  data.gene<-data.gene[,-1]
  data.metab<-data.metab[,-1]
  
  my.cluster<-read.csv("cluster.csv", header = T)[,1]
  
  my.data<-data.gene
  my.name<-name.gene
  
  na_index<-which(is.na(my.data[1,]))
  
  my.data<-my.data[,-na_index]
  my.cluster.sam<-my.cluster[-na_index]
  
  my.data<-apply(my.data, 2, as.numeric)
  
  sam.out.all <- sam(my.data, my.cluster.sam,rand = 123, gene.names = my.name)
  
  sum.sam.out.all <- summary(sam.out.all, 4)
  
  plot(sam.out.all, 4)
  
  name1<-my.name[sum.sam.out.all@mat.sig[,1]]
  #write.table(name1, "DEG_gene.txt", row.names = F)
}

condition.gene.plot<-F

if(condition.gene.plot){
  
  pval<-NULL
  for(k in 1:ncol(data.order[[3]])){
    a<-t.test(data.order[[3]][which(my.cluster[order(my.cluster)]==1),k],
              data.order[[3]][which(my.cluster[order(my.cluster)]==2),k],
              var.equal = F)
    pval[k]<-a$p.value
  }
  
  #gene.heat<- data.order[[3]][,order(apply(data.order[[3]], 2, mean, na.rm = T), decreasing = F)]
  gene.heat<- data.order[[3]][,order(pval, decreasing = F)[which(p.adjust(sort(pval, decreasing = F), method = "bonferroni")<=0.05)]]
  
  
  gene.heat<-rbind(gene.heat[my.cluster[order(my.cluster)]==1,][order(gene.heat[my.cluster[order(my.cluster)]==1,1], na.last = T),],
                   gene.heat[my.cluster[order(my.cluster)]==2,][order(gene.heat[my.cluster[order(my.cluster)]==2,1], na.last = T),],
                   gene.heat[my.cluster[order(my.cluster)]==3,][order(gene.heat[my.cluster[order(my.cluster)]==3,1], na.last = T),])
  
  heatmap.2( gene.heat, 
             trace="none", 
             dendrogram = "none",
             margins = c(5,9), 
             density.info=c("none"),
             key = T, Rowv = F, Colv = F,
             cexRow = 0.6, cexCol = 0.7, 
             key.xlab = "importance",
             key.ylab = "Density",
             col= colorRampPalette(c("blue3","white","red3"))(n = 50),
             font.lab=9,
             xlab = "Protein_Residue",
             ylab = "DNA_Residue",
             labRow = row.names(data.order[[1]]),
             labCol = colnames(data.order[[1]]),
             main = "GeneXpr1",
             na.color = "black")
  
  #write.csv(gene.heat, ".\\ml_data/Gene.csv", row.names = F)
  #write.csv(data.order[[2]], ".\\ml_data/Prot.csv", row.names = F)
}

plot_data_heatmap(model,
                  view = ,         # view of interest
                  factor = 1,             # factor of interest
                  features = 50,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)


#----------Weights-------------------

condition.weight<-F

if(condition.weight){
  
  plot_weights(model, abs = FALSE)
  
  weights$metabolite
  
  plot_top_weights(model,
                   view = c(3),
                   factor = 2,
                   nfeatures = 10
  )
  
  plot_top_weights(model,
                   view = c(3:1),
                   factor = 3,
                   nfeatures = 50
  )
  
  plot_top_weights(model,
                   view = "metabolite",
                   factor = 1,
                   nfeatures = 50
  )
  
  plot_top_weights(model,
                   view = "protein",
                   factor = 1,
                   nfeatures = 50
  )
  
  
  plot_top_weights(model,
                   view = "Gene",
                   factor = 1,
                   nfeatures = 45
  )
  
  
  model.tsne <- run_tsne(model)
  
  plot_dimred(model.tsne,
              method = "UMAP",
              scale="none")
  
  plot_weights(model,
               view = "protein",
               factor = 1,
               nfeatures = 10,     # Number of features to highlight
               scale = T,          # Scale weights from -1 to 1
               abs = F             # Take the absolute value?
  )
  
  plot_weights(model,
               view = "metabolite",
               factor = 1,
               nfeatures = 10,     # Number of features to highlight
               scale = T,          # Scale weights from -1 to 1
               abs = F             # Take the absolute value?
  )
  
  plot_weights(model,
               view = "Gene",
               factor = 13,
               nfeatures = 10,     # Number of features to highlight
               scale = T,          # Scale weights from -1 to 1
               abs = F             # Take the absolute value?
  )
  
}

##################################TCGA Validation (preprosessing)#############################################33
setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/")
my.cluster<-read.csv("./ml_data/cluster.csv", header = T)[,1]
data.order<-readRDS("data_order_all.rds")
cln.order<-readRDS("cln_order_all.rds")

name1<-read.table("DEG_gene.txt", header = T)
name1<-name1$x
setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/")


label<-read.csv("cluster.csv", header = T)[,1]%>%as.numeric()%>%as.factor()%>%sort()
#--------------------------Protein data preparation----------------------------#

data.prot<-as.data.frame(read.csv("proteinoslo.csv", header = T))
data.rppa<-read.table("data_rppa.txt", header = T, sep ="\t")

data.rppa.name<-strsplit(data.rppa$Composite.Element.REF, "\\|")%>%lapply(function(x){x[2]})%>%unlist()
data.prot.name<-strsplit(colnames(data.prot), "\\_\\.")%>%lapply(function(x){x[1]})%>%unlist()

a<-b<-c<-NULL
for(i in seq_along(data.prot.name)){
  
  a[i]<-data.prot.name[i]
  b[i]<-data.rppa.name[order(levenshteinSim(a[i], data.rppa.name), decreasing = T)[1]]
  c[i]<-data.rppa$Composite.Element.REF[order(levenshteinSim(a[i], data.rppa.name), decreasing = T)[1]]
  
}

cbind(a,b)
row.names(data.rppa)<-data.rppa$Composite.Element.REF

a<-gsub("_","-",a)
b<-gsub("_","-",b)

data.rppa.common<-data.rppa[c[which(levenshteinSim(a,b)>0.85)],]%>%as.data.frame()
data.rppa.common<-t(na.omit(data.rppa.common))
data.rppa.common<-data.rppa.common[-1,]
data.rppa.common<-apply(data.rppa.common,2,as.numeric)

data.prot.common<-data.prot[,colnames(data.prot)[which(levenshteinSim(a,b)>0.85)][-1]]

dim(data.prot.common)
dim(data.rppa.common)

row.names(data.rppa.common)<-colnames(data.rppa)[-1]

#write.csv(data.rppa.common, "val2018_prot_test.csv", row.names=T)
#write.csv(data.prot.common, "val2018_prot_train.csv", row.names=F)

#------------------------Gene data preparation---------------------------#

#data.oslo<-t(read.csv("geneoslo.csv", header = F))
data.oslo<-data.order[[3]][,name1] %>% t()
data.oslo<-cbind(name1, data.oslo)

data.val<-read.csv("gene_zscore.csv", header = T)%>%as.data.frame()
data.val<-data.val[,-2]

colnames(data.oslo)<-paste("V", 1:ncol(data.oslo), sep="")
colnames(data.val)[1]<-"V1"

data.oslo<-as.data.frame(data.oslo)
data.val<-as.data.frame(data.val)

#-----------------------------------------------------
data.val<-na.omit(data.val)


data<-merge(data.oslo, data.val, by="V1", all=F)


oslo.gene.common<-data[, colnames(data.oslo)]
val.gene.common<-data[, colnames(data.val)]

write.table(val.gene.common$V1, "common.txt", row.names = F, quote = F)
#row.names(val.gene.common)<-make.names(val.gene.common[,1], unique = TRUE)
#val.gene.common<-val.gene.common[,-c(1:2)]

val.gene.common[,-1]<-apply(val.gene.common[,-1], 2, as.numeric)
val.gene.common[,-1]<-preprocessCore::normalize.quantiles(as.matrix(val.gene.common[,-1]))

#row.names(data)<-row.names(val.gene.common)
#colnames(data)<-colnames(val.gene.common)

#write.csv(val.gene.common, "val2018_gene_test.csv", row.names=F)
#write.csv(oslo.gene.common, "val2018_gene_train.csv", row.names=F)

#-----------------------------------------------------------------------
val.gene.common[which(val.gene.common$V1%in%imp[1:20]),]

###################################TCGA Valiadtion##################################################

set.seed(1)

stnd<-function(x){
  (x-mean(x))/sd(x)
}

minmax<-function(x){
  (x-min(x))/(max(x)-min(x))
}
prot.test<-read.csv("val2018_prot_test.csv", header = T, row.names = 1)
prot.train<-read.csv("val2018_prot_train.csv", header = T)
#gene.test<-read.csv("val2018_gene_test.csv", header = T)
#gene.train<-read.csv("val2018_gene_train.csv", header = T)
gene.test<-val.gene.common
gene.train<-oslo.gene.common

#write.table(gene.test[,1] %>% as.data.frame(),
#           "/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/TCGA_gene.csv", 
#          row.names=F, col.names = F)
#write.table(prot.test %>% colnames() %>% as.data.frame(),
#           "/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/TCGA_prot.csv", 
#          row.names=F, col.names = F)



label<-read.csv("cluster.csv", header = T)[,1]%>%as.numeric()%>%as.factor()%>%sort()

#-------------------------------------------
prot.test %>% dim()
gene.test %>% dim()

cbind(rownames(gene.train),row.names(gene.test))
cbind(colnames(prot.train), colnames(prot.test))

row.names(gene.test)<-make.names(gene.test$V1, unique = TRUE) #gene.test$V1
gene.test<-gene.test[,-1]

row.names(gene.train)<-make.names(gene.train$V1, unique = TRUE) #gene.train$V1
gene.train<-gene.train[,-1]


mldata<-cbind(prot.train, t(gene.train))%>%apply(MARGIN = 2, as.numeric)%>%as.data.frame()
colnames(mldata)<-paste("V", 1:ncol(mldata), sep ="")
mldata$label<-label
mldata<-na.omit(mldata)


mldata[,1:(ncol(mldata)-1)]<-apply(mldata[,1:(ncol(mldata)-1)], 2, stnd)


#--------------------------------------------

class1<-mldata[which(mldata$label==1),]
class2<-mldata[which(mldata$label==2),]
class3<-mldata[which(mldata$label==3),]

class_list<-list(class1, class2, class3)

cv<-1:5
bin.svm.l<-list()
svm.cm<-list()
bin.rf<-list()
rf.cm<-list()
pls.cm<-list()
bin.pls<-list()
cv.label<-list()


for(k in 1:5){
  train_list<-test_list<-list(NULL)
  
  for(i in 1:3){
    set.seed(1)
    ss<-sample(cv, size = nrow(class_list[[i]]), prob = rep(0.2, 5), replace = T )
    train_list[[i]] <- class_list[[i]][which(ss!=cv[k]),]
    test_list[[i]] <- class_list[[i]][which(ss==cv[k]),]
  }
  
  
  train<-do.call(rbind, train_list)
  train<-train[sample(1:nrow(train)),]
  
  test<-do.call(rbind, test_list)
  test<-test[sample(1:nrow(test)),]
  
  #------------------------------------------------
  #detach(package:MOFA2, unload = TRUE)
  x<-T
  if(x){
    classifier.svm.l = svm(formula = label ~ .,
                           data = as.data.frame(train),
                           type = 'C-classification',
                           kernel = "linear",
                           probability=TRUE)
    
    pred.svm.l = predict(classifier.svm.l, newdata = test[,-ncol(test)])
    a<-predict(classifier.svm.l, newdata = test[,-ncol(test)],probability = T)
    
    bin.svm.l[[k]]<-predict(classifier.svm.l, newdata = test[,-ncol(test)],probability = T)
    
    bin.svm.l[[k]]<-attr(a, "probabilities")[,c(3,1,2)]
    svm.cm[[k]]<-confusionMatrix(pred.svm.l,
                                 as.factor(test[,ncol(test)]),
                                 mode = "everything")
    
    
    #----------------PLS-DA
    
    model.plsda<-plsda(train[,-ncol(train)], 
                       as.factor(train$label), ncomp = 3,
                       probMethod = "softmax",
                       type="proba")
    
    pred.plsda<- predict(model.plsda, 
                         test[,-ncol(test)])
    bin.pls[[k]]<-predict(model.plsda, 
                          test[,-ncol(test)])
    
    bin.pls[[k]]<-predict(model.plsda, 
                          test[,-ncol(test)], 
                          type="prob")
    
    prob.plsda<-matrix(as.matrix(a), ncol=3, byrow = F)[,2]
    
    pls.cm[[k]]<-confusionMatrix(pred.plsda, 
                                 as.factor(test$label), 
                                 mode = "everything")
    cv.label[[k]]<-test$label
  }
  #-------------randomforest
  
  model.rf<-randomForest(y= as.factor(as.numeric(as.matrix(train[,ncol(train)]))), 
                         x= as.matrix(train[,-ncol(train)]), 
                         probability=T)
  
  pred.rf <- predict(model.rf, test)
  bin.rf[[k]] <- predict(model.rf, test, type="prob")
  #bin.rf[[k]]<-predict(model.rf, test)
  
  rf.cm[[k]]<-confusionMatrix(pred.rf, as.factor(test$label), mode = "everything")
  mod.cv.rf<-model.rf
}

par(mfrow=c(1,3))
for(i in 1:3){
  
  svm_pred<-lapply(bin.svm.l, function(x){x[,as.character(i)]}) %>% unlist()
  lab1<-ifelse(unlist(cv.label)==i,1,0)
  roc_svm<-roc(lab1 %>% as.numeric(),svm_pred %>% as.numeric(), percent = T)
  plot(roc_svm, 
       grid=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       main= paste("ROC MOC",i," vs MOCs", sep=""),
       xaxt="n", yaxt="n",
       cex.lab=1.3,
       font.lab=9,
       print.auc.y=40,
       print.auc.x=24, 
       print.auc=F,
       legacy.axes=T,
       percent = T,
       col="red",
       lwd=5)
  
  
  rf_pred<-lapply(bin.rf, function(x){x[,as.character(i)]}) %>% unlist()
  roc_rf<-roc(lab1, rf_pred,percent=T)
  plot(roc_rf,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=30,
       print.auc.x=24, 
       print.auc=F,
       col="darkgoldenrod1", 
       lwd=5)
  
  pls_pred<-lapply(bin.pls, function(x){x %>% matrix(ncol=3, byrow=F)})
  pls_pred<-lapply(pls_pred, function(x){x[,i]}) %>% unlist()
  
  roc_pls<-roc(lab1, pls_pred,percent=T)
  plot(roc_pls,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=20,
       print.auc.x=24, 
       print.auc=F,
       col="deepskyblue4", 
       lwd=5)
  
  all_pred<-cbind(pls_pred, rf_pred, svm_pred) %>% rowMeans()
  roc_pls1<-roc(lab1, all_pred,percent=T)
  plot(roc_pls1,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=50,
       print.auc.x=24, 
       print.auc=F,
       col="black", 
       lwd=5)
  
  text<-c(paste("SVM","(AUC=", as.character(round(roc_svm$auc,2)),"%)"),
          paste("RF","(AUC=", as.character(round(roc_rf$auc,2)),"%)"),
          paste("PLS","(AUC=", as.character(round(roc_pls$auc,2)),"%)"),
          paste("Ensemble","(AUC=", as.character(round(roc_pls1$auc,2)),"%)"))
  
  legend("bottomright",
         legend = text,
         lty=1,
         lwd=4,
         cex=1,
         col = c("red","darkgoldenrod1", "deepskyblue4","black"),
         bty="n",
         text.font = 9)
  
  axis(1, at=seq(100,0,-10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
  axis(2, at=seq(0,100,10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
  
}

for(k in 1:5){
  capture.output(
    pls.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","TCGA_pls_",k,".csv", sep=""))
  capture.output(
    rf.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","TCGA_cm_",k,".csv", sep=""))
  capture.output(
    svm.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","TCGA_svm_",k,".csv", sep=""))
}
#-------------------------------------------------------------------

data<-merge(apply(prot.test, 2, stnd), t(gene.test), all = F, by = "row.names")
data.name<-data[,1]
data<-data[,-1]%>%apply(MARGIN = 2, as.numeric)
colnames(data)<-paste("V",1:ncol(data), sep="")

pred.rf <- predict(model.rf, data)
pred.svm.l <- predict(classifier.svm.l, data)
pred.plsda <- predict(model.plsda, data)

val.cluster<-round(rowMeans(cbind(pred.svm.l, pred.rf, pred.plsda)), 0)
data %>% dim()
#------------------------------------------------------------------

val.cluster<-as.data.frame(cbind(data.name, val.cluster))
colnames(val.cluster)<-c("PATIENT_ID", "cluster")
row.names(val.cluster)<-NULL

xx<-gsub( "\\.01", "", val.cluster$PATIENT_ID)
PATIENT_ID<-gsub("\\.", "-", xx)
val.cluster$PATIENT_ID<-PATIENT_ID
#write.table(PATIENT_ID, "pid.txt", row.names = F)
#write.table(val.cluster, "val_cluster.txt", row.names = F)
#-------------------------------------------------------------------

clinic<-read.csv("clinic2018.csv", header = T)
merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  #merge.clinic$RADIATION_THERAPY=="No" &
  merge.clinic$HISTORY_NEOADJUVANT_TRTYN=="No"),]

merge.clinic.white<-merge.clinic[which(merge.clinic$RACE=="White"&
                                         #merge.clinic$RADIATION_THERAPY=="No" &
                                         merge.clinic$HISTORY_NEOADJUVANT_TRTYN=="No"),]

#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="White"),]

#----------survival-------------
cli.data<-merge.clinic[,c(1,8,9,2)]

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()

cli.data$Time<-as.numeric(cli.data$Time)
model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

ggsurvplot(model.srv, palette= c("blue", "red", "green", 
                                 "black", "purple", "orange"), 
           size.lab=1.5, pval = T, font.size=1.5,
           surv.scale = "percent",
           font.legend = list(size = 15, color = "black", face = "bold"),
           font.tickslab = c(20, "plain", "black"),
           risk.table = T,
           legend.title="TCGA",
           font.x = c(20, "bold.italic", "black"),
           font.y = c(20, "bold.italic", "black"),
           pval.size = 7, pval.coord = c(60, .1)) %>% print()

survdiff(Surv(Time, Status)~Cluster, data=cli.data)



cli.data<-merge.clinic.white[,c(1,8,9,2)]

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()

cli.data$Time<-as.numeric(cli.data$Time)
model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

ggsurvplot(model.srv, palette= c("blue", "red", "green", 
                                 "black", "purple", "orange"), 
           size.lab=1.5, pval = T, font.size=1.5,
           surv.scale = "percent",
           font.legend = list(size = 15, color = "black", face = "bold"),
           font.tickslab = c(20, "plain", "black"),
           risk.table = T,
           legend.title="TCGA White Population",
           font.x = c(20, "bold.italic", "black"),
           font.y = c(20, "bold.italic", "black"),
           pval.size = 7, pval.coord = c(60, .1)) %>% print()

survdiff(Surv(Time, Status)~Cluster, data=cli.data)

#table(merge.clinic$OS_STATUS)
#--------------------------------------------------------------------------------

table(clinic$SEX)
nrow(clinic)

data.heat<-merge.clinic[order(merge.clinic$cluster),]


heatmap.2( cbind(data.heat$cluster%>%as.numeric(),
                 data.heat$cluster%>%as.numeric()), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue",
                                   "red", "green"))(n = 3),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.heat),
           labCol = colnames(data.heat),
           main = "cluster")

#-------------------------------------------------------------------#

pam50<-data.heat$SUBTYPE
pam50<-as.factor(pam50)%>%as.numeric()

heatmap.2( cbind(pam50, pam50), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("black","blue",
                                   "red", "skyblue", 
                                   "orange", "purple"))(n = 6),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.heat),
           labCol = colnames(data.heat),
           main = "PAM50",
           na.color = "black")



status<-data.heat$OS_STATUS
#status<-data.heat$VITAL_STATUS

status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==1)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==1)], na.rm = T)
status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==2)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==2)], na.rm = T)
status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==3)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==3)], na.rm = T)


heatmap.2( cbind(status, status), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("lightblue", "blue", "darkblue"))(n = 10),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           #labRow = row.names(data.order),
           #labCol = colnames(data.order),
           main = "death",
           na.color = "white")




############################METABRIC#########################################

setwd("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/")

set.seed(111)

#---------------Metabolomic Data----------------------
#data.oslo<-t(read.csv("geneoslo.csv", header = F))

data.val<-read.csv("gene_zscore.csv", header = T)
#write.csv(data.val, "gene_zscore.csv", row.names=F)

label<-read.csv("cluster.csv", header = T)[,1]%>%as.numeric()%>%as.factor()%>%sort()

colnames(data.val)[1]<-"V1"

data.val<-as.data.frame(data.val)

#-----------------------------------------------------

data<-full_join(data.oslo[,1:2], data.val, by="V1")
data<-na.omit(data)

data<-data[,-c(2:3)]
names<-data[,1]
data<-data[,-1]
data<-apply(data, 2, as.numeric)
data<-preprocessCore::normalize.quantiles(as.matrix(data))
row.names(data)<-names
colnames(data)<-colnames(data.val)[-c(1:2)]

#write.csv(data, "gene_val.csv", row.names = F)
#write.csv(names, "/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/meta_gene.csv", row.names = F)
#---------------------------------------------------------
row.names(data.oslo)<-data.oslo[,1]
data.oslo<-data.oslo[,-1]

#---------------------------------------------------------
mldata<-apply(t(data.oslo[names,]), 2, as.numeric)%>%as.data.frame()
mldata$label<-label
mldata<-na.omit(mldata)

mldata[,1:(ncol(mldata)-1)]<-apply(mldata[,1:(ncol(mldata)-1)], 2, stnd)

row.names(data)<-paste("V", c(1:ncol(t(data))), sep ="")
colnames(mldata)[-ncol(mldata)]<-paste("V", c(1:ncol(t(data))), sep ="")

#---------------data spliting-----------------

class1<-mldata[which(mldata$label==1),]
class2<-mldata[which(mldata$label==2),]
class3<-mldata[which(mldata$label==3),]

class_list<-list(class1, class2, class3)

cv<-1:5
bin.svm.l<-list()
svm.cm<-list()
bin.rf<-list()
rf.cm<-list()
pls.cm<-list()
bin.pls<-list()
cv.label<-list()


for(k in 1:5){
  train_list<-test_list<-list(NULL)
  
  for(i in 1:3){
    set.seed(1)
    ss<-sample(cv, size = nrow(class_list[[i]]), prob = rep(0.2, 5), replace = T )
    train_list[[i]] <- class_list[[i]][which(ss!=cv[k]),]
    test_list[[i]] <- class_list[[i]][which(ss==cv[k]),]
  }
  
  
  train<-do.call(rbind, train_list)
  train<-train[sample(1:nrow(train)),]
  
  test<-do.call(rbind, test_list)
  test<-test[sample(1:nrow(test)),]
  
  #------------------------------------------------
  #detach(package:MOFA2, unload = TRUE)
  x<-T
  if(x){
    classifier.svm.l = svm(formula = label ~ .,
                           data = as.data.frame(train),
                           type = 'C-classification',
                           kernel = "linear",
                           probability=TRUE)
    
    pred.svm.l = predict(classifier.svm.l, newdata = test[,-ncol(test)])
    a<-predict(classifier.svm.l, newdata = test[,-ncol(test)],probability = T)
    
    bin.svm.l[[k]]<-predict(classifier.svm.l, newdata = test[,-ncol(test)],probability = T)
    
    bin.svm.l[[k]]<-attr(a, "probabilities")[,c(3,1,2)]
    svm.cm[[k]]<-confusionMatrix(pred.svm.l,
                                 as.factor(test[,ncol(test)]),
                                 mode = "everything")
    
    
    #----------------PLS-DA
    
    model.plsda<-plsda(train[,-ncol(train)], 
                       as.factor(train$label), ncomp = 3,
                       probMethod = "softmax",
                       type="proba")
    
    pred.plsda<- predict(model.plsda, 
                         test[,-ncol(test)])
    bin.pls[[k]]<-predict(model.plsda, 
                          test[,-ncol(test)])
    
    bin.pls[[k]]<-predict(model.plsda, 
                          test[,-ncol(test)], 
                          type="prob")
    
    prob.plsda<-matrix(as.matrix(a), ncol=3, byrow = F)[,2]
    
    pls.cm[[k]]<-confusionMatrix(pred.plsda, 
                                 as.factor(test$label), 
                                 mode = "everything")
    cv.label[[k]]<-test$label
  }
  #-------------randomforest
  
  model.rf<-randomForest(y= as.factor(as.numeric(as.matrix(train[,ncol(train)]))), 
                         x= as.matrix(train[,-ncol(train)]), 
                         probability=T)
  
  pred.rf <- predict(model.rf, test)
  bin.rf[[k]] <- predict(model.rf, test, type="prob")
  #bin.rf[[k]]<-predict(model.rf, test)
  
  rf.cm[[k]]<-confusionMatrix(pred.rf, as.factor(test$label), mode = "everything")
  mod.cv.rf<-model.rf
}


for(k in 1:5){
  capture.output(
    pls.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","meta_pls_",k,".csv", sep=""))
  capture.output(
    rf.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","meta_cm_",k,".csv", sep=""))
  capture.output(
    svm.cm[[k]],
    file = paste("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/","meta_svm_",k,".csv", sep=""))
}

par(mfrow=c(1,3))
for(i in 1:3){
  
  svm_pred<-lapply(bin.svm.l, function(x){x[,as.character(i)]}) %>% unlist()
  lab1<-ifelse(unlist(cv.label)==i,1,0)
  roc_svm<-roc(lab1 %>% as.numeric(),svm_pred %>% as.numeric(), percent = T)
  plot(roc_svm, 
       grid=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       main= paste("ROC MOC",i," vs MOCs", sep=""),
       xaxt="n", yaxt="n",
       cex.lab=1.3,
       font.lab=9,
       print.auc.y=40,
       print.auc.x=24, 
       print.auc=F,
       legacy.axes=T,
       percent = T,
       col="red",
       lwd=5)
  
  
  rf_pred<-lapply(bin.rf, function(x){x[,as.character(i)]}) %>% unlist()
  roc_rf<-roc(lab1, rf_pred,percent=T)
  plot(roc_rf,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=30,
       print.auc.x=24, 
       print.auc=F,
       col="darkgoldenrod1", 
       lwd=5)
  
  pls_pred<-lapply(bin.pls, function(x){x %>% matrix(ncol=3, byrow=F)})
  pls_pred<-lapply(pls_pred, function(x){x[,i]}) %>% unlist()
  
  roc_pls<-roc(lab1, pls_pred,percent=T)
  plot(roc_pls,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=20,
       print.auc.x=24, 
       print.auc=F,
       col="deepskyblue4", 
       lwd=5)
  
  all_pred<-cbind(pls_pred, rf_pred, svm_pred) %>% rowMeans()
  roc_pls1<-roc(lab1, all_pred,percent=T)
  plot(roc_pls1,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=50,
       print.auc.x=24, 
       print.auc=F,
       col="black", 
       lwd=5)
  
  text<-c(paste("SVM","(AUC=", as.character(round(roc_svm$auc,2)),"%)"),
          paste("RF","(AUC=", as.character(round(roc_rf$auc,2)),"%)"),
          paste("PLS","(AUC=", as.character(round(roc_pls$auc,2)),"%)"),
          paste("Ensemble","(AUC=", as.character(round(roc_pls1$auc,2)),"%)"))
  
  legend("bottomright",
         legend = text,
         lty=1,
         lwd=4,
         cex=1,
         col = c("red","darkgoldenrod1", "deepskyblue4","black"),
         bty="n",
         text.font = 9)
  
  axis(1, at=seq(100,0,-10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
  axis(2, at=seq(0,100,10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
  
}
#-------------------------------------------------------------------

pred.rf <- predict(model.rf, t(data))
pred.svm.l <- predict(classifier.svm.l, t(data))
pred.plsda <- predict(model.plsda, t(data))

val.cluster<-round(rowMeans(cbind(pred.svm.l, pred.rf, pred.plsda)), 0)
#------------------------------------------------------------------

val.cluster<-as.data.frame(cbind(names(val.cluster), val.cluster))
colnames(val.cluster)<-c("PATIENT_ID", "cluster")
row.names(val.cluster)<-NULL

#xx<-gsub( "\\.01", "", val.cluster$PATIENT_ID)
PATIENT_ID<-gsub("\\.", "-", val.cluster$PATIENT_ID)
val.cluster$PATIENT_ID<-PATIENT_ID
#write.table(PATIENT_ID, "pid.txt", row.names = F)
#write.table(val.cluster, "val_cluster.txt", row.names = F)
#-------------------------------------------------------------------

clinic<-read.csv("clinic_data.csv", header = T)
merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")
merge.clinic[which(merge.clinic$VITAL_STATUS==2), "VITAL_STATUS"]<-1


#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="WHITE"),]
#----------survival-------------

cli.data<-merge.clinic[,c(1,12,6,2)]

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

cli.data<-na.omit(cli.data)

cli.data$Time<-as.numeric(cli.data$Time)
model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

ggsurvplot(model.srv, palette= c("blue", "red", "green", 
                                 "black", "purple", "pink", "orange"), 
           size.lab=1.5, pval = T, font.size=1.5,
           surv.scale = "percent",
           font.legend = list(size = 15, color = "black", face = "bold"),
           font.tickslab = c(20, "plain", "black"),
           risk.table = T,
           legend.title="All Death",
           font.x = c(20, "bold.italic", "black"),
           font.y = c(20, "bold.italic", "black"),
           pval.size = 7, pval.coord = c(60, .1))

survdiff(Surv(Time, Status)~Cluster, data=cli.data)



clinic<-read.csv("clinic_data.csv", header = T)
merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")
merge.clinic[which(merge.clinic$VITAL_STATUS==2), "VITAL_STATUS"]<-NA

merge.clinic<-merge.clinic[which(merge.clinic$RADIO_THERAPY=="NO"),]
#----------survival-------------

cli.data<-merge.clinic[,c(1,12,6,2)]

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

cli.data<-na.omit(cli.data)

cli.data$Time<-as.numeric(cli.data$Time)
model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)

ggsurvplot(model.srv, palette= c("blue", "red", "green", 
                                 "black", "purple", "pink", "orange"), 
           size.lab=1.5, pval = T, font.size=1.5,
           surv.scale = "percent",
           font.legend = list(size = 15, color = "black", face = "bold"),
           font.tickslab = c(20, "plain", "black"),
           risk.table = T,
           legend.title="BC Death",
           font.x = c(20, "bold.italic", "black"),
           font.y = c(20, "bold.italic", "black"),
           pval.size = 7, pval.coord = c(60, .1))

survdiff(Surv(Time, Status)~Cluster, data=cli.data)





#--------------------------------------------------------------------------------


library(gplots)

data.heat<-merge.clinic[order(merge.clinic$cluster),]


heatmap.2( cbind(data.heat$cluster%>%as.numeric(),
                 data.heat$cluster%>%as.numeric()), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue",
                                   "red", "green"))(n = 3),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.heat),
           labCol = colnames(data.heat),
           main = "cluster")

#-------------------------------------------------------------------#

pam50<-data.heat$CLAUDIN_SUBTYPE
pam50<-as.factor(pam50)%>%as.numeric()

heatmap.2( cbind(pam50, pam50), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("blue",
                                   "red", "skyblue", 
                                   "orange", "purple",
                                   "pink","yellow"))(n = 7),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.heat),
           labCol = colnames(data.heat),
           main = "PAM50",
           na.color = "black")



status<-data.heat$OS_STATUS
status<-data.heat$VITAL_STATUS

status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==1)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==1)], na.rm = T)
status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==2)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==2)], na.rm = T)
status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==3)]<-mean(status[which(merge.clinic$cluster[order(merge.clinic$cluster)]==3)], na.rm = T)


heatmap.2( cbind(status, status), 
           trace="none", 
           dendrogram = "none",
           margins = c(5,9), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 0.7, 
           key.xlab = "importance",
           key.ylab = "Density",
           col= colorRampPalette(c("lightblue", "blue", "darkblue"))(n = 10),
           font.lab=9,
           xlab = "Protein_Residue",
           ylab = "DNA_Residue",
           labRow = row.names(data.order),
           labCol = colnames(data.order),
           main = "death",
           na.color = "white")



clitation("siggenes") %>% print(bib = T)



