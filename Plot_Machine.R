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
library(reshape)
library(waffle)
library(hrbrthemes)
library(ggthemes)
library(forestmodel)

#-----Load Trained Model-----------------

model <- load_model("MOFA_allgene_20F.hdf5")

model@data <- model@data[c("metabolite","protein", "Gene")]


#----------partition

views_names(model)<-c("Zmetabolite","Yprotein", "XGene")
plotplot<-function (object, covariate = 1, colors = NULL, show_covariate = FALSE, 
                    show_dimensions = TRUE) {
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (sum(get_dimensions(object)[["N"]]) > 10000) 
    warning("This function is inefficient with large number of cells...")
  if (length(object@data) == 0) 
    stop("Data not found")
  M <- get_dimensions(object)[["M"]]
  G <- get_dimensions(object)[["G"]]
  if (M == 1 & G == 1) 
    warning("This function is not useful when there is just one view and one group")
  if (!.hasSlot(object, "covariates") || any(object@dimensions[["C"]] < 
                                             1, is.null(object@covariates))) 
    covariate <- NULL
  if (!is.null(covariate)) {
    if (is.numeric(covariate)) {
      if (covariate > object@dimensions[["C"]]) 
        stop("Covariate index out of range")
      covariate <- covariates_names(object)[covariate]
    }
    if (!is.character(covariate) | !covariate %in% covariates_names(object)) 
      stop("Covariate mispecified. Please read the documentation")
    covari <- .set_xax(object, covariate)
  }
  if (is.null(colors)) {
    palette <- c("#FF7F50", "#D95F02", "#377EB8", "#E6AB02", 
                 "#31A354", "#7570B3", "#E7298A", "#66A61E", "#A6761D", 
                 "#666666", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", 
                 "#FFFF33", "#A65628", "#F781BF", "#1B9E77")
    if (M < 18) 
      colors <- palette[seq_len(M)]
    else colors <- rainbow(M)
    names(colors) <- views_names(object)
  }
  else {
    if (length(colors) != M) 
      stop("Length of 'colors' does not match the number of views")
    if (is.null(names(colors))) {
      names(colors) <- views_names(object)
    }
    else {
      stopifnot(sort(names(colors)) == sort(views_names(object)))
    }
  }
  tmp <- lapply(object@data, function(m) sapply(m, function(g) apply(g, 
                                                                     2, function(x) !all(is.na(x)))))
  ovw <- do.call(cbind, lapply(seq_len(M), function(m) {
    do.call(rbind, lapply(tmp[[m]], as.data.frame))
  }))
  rownames(ovw) <- object@samples_metadata$sample
  colnames(ovw) <- views_names(object)
  ovw$sample <- object@samples_metadata$sample
  ovw$group <- object@samples_metadata$group
  to.plot <- reshape::melt(ovw, id.vars = c("sample", "group"), var = c("view"))
  if (!is.null(covariate)) {
    to.plot <- left_join(to.plot, covari, by = "sample")
    to.plot$sample <- factor(to.plot$sample, levels = unique(to.plot$sample[order(to.plot$covariate_value)]))
  }
  else {
    to.plot$sample <- factor(to.plot$sample, levels = rownames(ovw))
  }
  n <- length(unique(to.plot$sample))
  to.plot$combi <- ifelse(to.plot$value, as.character(to.plot$view), 
                          "missing")
  if (show_dimensions) {
    to.plot$ntotal <- paste("N=", sapply(object@data[[1]], 
                                         function(e) ncol(e))[as.character(to.plot$group)], 
                            sep = "")
    to.plot$ptotal <- paste("D=", sapply(object@data, function(e) nrow(e[[1]]))[as.character(to.plot$view)], 
                            sep = "")
    if (length(unique(to.plot$group)) == 1) {
      to.plot <- mutate(to.plot, view_label = paste(view, 
                                                    ptotal, sep = "\n"), group_label = ntotal)
    }
    else {
      to.plot <- mutate(to.plot, view_label = paste(view, 
                                                    ptotal, sep = "\n"), group_label = paste(group, 
                                                                                             ntotal, sep = "\n"))
    }
  }
  else {
    to.plot <- mutate(to.plot, view_label = view, group_label = group)
  }
  to.plot$group_label <- factor(to.plot$group_label, levels = unique(to.plot$group_label))
  p <- ggplot(to.plot, aes_string(x = "sample", y = "view_label", 
                                  fill = "combi")) + geom_tile() + scale_fill_manual(values = c(missing = "darkblue", 
                                                                                                colors)) + guides(fill = "none") + facet_wrap(~group_label, 
                                                                                                                                              scales = "free_x", nrow = length(unique(to.plot$view_label))) + 
    theme(panel.background = element_rect(fill = "white"), 
          text = element_text(size = 14), axis.line = element_blank(), 
          axis.ticks = element_blank(), axis.title = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_text(color = "black"), 
          strip.background = element_blank(), panel.grid = element_blank())
  if (show_covariate) {
    p2 <- ggplot(to.plot, aes_string(x = "sample", y = "covariate_value")) + 
      geom_point(size = 0.5) + theme_bw() + theme(text = element_text(size = 10), 
                                                  axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
                                                  axis.text.x = element_blank(), strip.background = element_blank(), 
                                                  strip.text = element_blank()) + ylab(covariate) + 
      facet_wrap(~group_label, ncol = 1, scales = "free_x")
    if (object@dimensions["G"] == 1) {
      p <- cowplot::plot_grid(p, p2, align = "v", ncol = 1, 
                              rel_heights = c(1, 0.2))
    }
    else {
      p <- cowplot::plot_grid(p, p2, align = "h", nrow = 1, 
                              rel_widths = c(1, 1))
    }
  }
  return(p)
}

plotplot(model, colors = c("darkorange3","goldenrod","plum2"))




vard<-model@cache$variance_explained$r2_total[[1]]
#vard<-vard[c(2,1,3)]

p<-barplot(vard, beside = T,
           xaxt="n", #yaxt="n",
           ylim=c(0,50), xlim=c(0,4),
           #space = c(0.01,3),
           border = "black",
           col = c("darkorange3","goldenrod","plum2"),
           cex.axis = 3, font=2)

text(x = p, cex=2.4, font=2,
     #y = c(as.numeric(vard))-11,
     y=c(11.6,22,22),
     labels = paste(c("Metabolite","Protein","Transcriptome"), paste(round(vard,2), "%", sep=""), sep = "   "), srt=90)



heatheat<-function (object, x = "view", y = "factor", split_by = NA, plot_total = FALSE, 
                    factors = "all", min_r2 = 0, max_r2 = NULL, legend = TRUE, 
                    use_cache = TRUE, ...) 
{
  if (length(unique(c(x, y, split_by))) != 3) {
    stop(paste0("Please ensure x, y, and split_by arguments are different.\n", 
                "  Possible values are `view`, `group`, and `factor`."))
  }
  if (is.na(split_by)) 
    split_by <- setdiff(c("view", "factor", "group"), c(x, 
                                                        y, split_by))
  if ((use_cache) & .hasSlot(object, "cache") && ("variance_explained" %in% 
                                                  names(object@cache))) {
    r2_list <- object@cache$variance_explained
  }
  else {
    r2_list <- calculate_variance_explained(object, factors = factors, 
                                            ...)
  }
  r2_mk <- r2_list$r2_per_factor
  r2_mk_df <- melt(lapply(r2_mk, function(x) melt(as.matrix(x), 
                                                  varnames = c("factor", "view"))), id.vars = c("factor", 
                                                                                                "view", "value"))
  colnames(r2_mk_df)[ncol(r2_mk_df)] <- "group"
  if ((length(factors) == 1) && (factors[1] == "all")) {
    factors <- factors_names(object)
  }
  else {
    if (is.numeric(factors)) {
      factors <- factors_names(object)[factors]
    }
    else {
      stopifnot(all(factors %in% factors_names(object)))
    }
    r2_mk_df <- r2_mk_df[r2_mk_df$factor %in% factors, ]
  }
  r2_mk_df$factor <- factor(r2_mk_df$factor, levels = factors)
  r2_mk_df$group <- factor(r2_mk_df$group, levels = groups_names(object))
  r2_mk_df$view <- factor(r2_mk_df$view, levels = views_names(object))
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])
  if (!is.null(min_r2)) 
    r2_mk_df$value[r2_mk_df$value < min_r2] <- 0.001
  min_r2 = 0
  if (!is.null(max_r2)) {
    r2_mk_df$value[r2_mk_df$value > max_r2] <- max_r2
  }
  else {
    max_r2 = max(r2_mk_df$value)
  }
  p1 <- ggplot(r2_mk_df, aes_string(x = x, y = y), lwd=3) + geom_tile(aes_string(fill = "value"), size=0.7,
                                                                      color = "black") + facet_wrap(as.formula(sprintf("~%s", 
                                                                                                                       split_by)), nrow = 1) + labs(x = "", y = "", title = "") + 
    scale_fill_gradientn(colors = c("white", "darkred"), 
                         guide = "colorbar", limits = c(min_r2, max_r2)) + 
    guides(fill = guide_colorbar("Var. (%)")) + theme(axis.text.x = element_text(size = rel(0), 
                                                                                 color = "black"), axis.text.y = element_text(size = rel(0), 
                                                                                                                              color = "black"), axis.line = element_blank(),
                                                      axis.ticks = element_blank(), 
                                                      panel.background = element_blank(), strip.background = element_blank(), 
                                                      strip.text = element_text(size = rel(2.5)),
                                                      #legend.position = "bottom",
                                                      legend.key.size = unit(1, 'cm'), #change legend key size
                                                      legend.key.height = unit(1, 'cm'), #change legend key height
                                                      legend.key.width = unit(1, 'cm'), #change legend key width
                                                      legend.title = element_text(size=24), #change legend title font size
                                                      legend.text = element_text(size=40))+
    scale_size_manual(values = c(1, 2))
  
  if (isFALSE(legend)) 
    p1 <- p1 + theme(legend.position = "none")
  if (length(unique(r2_mk_df[, split_by])) == 1) 
    p1 <- p1 + theme(strip.text = element_blank())
  if (plot_total) {
    r2_m_df <- melt(lapply(r2_list$r2_total, function(x) lapply(x, 
                                                                function(z) z)), varnames = c("view", "group"), value.name = "R2")
    colnames(r2_m_df)[(ncol(r2_m_df) - 1):ncol(r2_m_df)] <- c("view", 
                                                              "group")
    r2_m_df$group <- factor(r2_m_df$group, levels = MOFA2::groups_names(object))
    r2_m_df$view <- factor(r2_m_df$view, levels = views_names(object))
    min_lim_bplt <- min(0, r2_m_df$R2)
    max_lim_bplt <- max(r2_m_df$R2)
    p2 <- ggplot(r2_m_df, aes_string(x = x, y = "R2")) + 
      geom_bar(stat = "identity", fill = "deepskyblue4", 
               color = "black", width = 0.9) + facet_wrap(as.formula(sprintf("~%s", 
                                                                             split_by)), nrow = 1) + xlab("") + ylab("Variance explained (%)") + 
      scale_y_continuous(limits = c(min_lim_bplt, max_lim_bplt), 
                         expand = c(0.005, 0.005)) + theme(axis.ticks.x = element_blank(), 
                                                           axis.text.x = element_text(size = rel(2), color = "black"), 
                                                           axis.text.y = element_text(size = rel(2), color = "black"), 
                                                           axis.title.y = element_text(size = rel(2), color = "black"), 
                                                           axis.line = element_line(size = rel(2), color = "black"), 
                                                           panel.background = element_blank(), strip.background = element_blank(), 
                                                           strip.text = element_text(size = rel(2)))
    if (length(unique(r2_m_df[, split_by])) == 1) 
      p2 <- p2 + theme(strip.text = element_blank())
    plot_list <- list(p1, p2)
  }
  else {
    plot_list <- p1
  }
  return(plot_list)
}

getwd()
heatheat(model, x="view", y="factor", plot_total = F)

AN<-F
source("factor_evaluation.R")

tooth<-data.frame(get_factors(model, factors = c(1,2,13), 
                              as.data.frame = F)[[1]],
                  model@samples_metadata[,c("ER","PR","HER2","grade")])
tooth[,1:3]<-scale(tooth[,1:3])
tooth<-tooth[-which(tooth$ER %>% is.na()),]

margin<-theme(plot.margin = unit(c(-1,0.5,-3,-1), "cm"))

p1<-ggplot(tooth, aes(x=ER, y=Factor1, fill=ER)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor1") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5,
                     label.y = min(tooth$Factor1),size=8)


p2<-ggplot(tooth, aes(x=ER, y=Factor2, fill=ER)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor2") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5, 
                     label.y = min(tooth$Factor2), size=8)

p3<-ggplot(tooth, aes(x=ER, y=Factor13, fill=ER)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=50, face="bold"),
        axis.title=element_text(size=40,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor13") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 1.5, 
                     vjust = 0, size=8)

gridExtra::grid.arrange(p1, p2, ncol = 1)


#gridExtra::grid.arrange(p1, p2, nrow = 1)

#--------------------------------------------#

tooth<-data.frame(get_factors(model, factors = c(1,2,13), 
                              as.data.frame = F)[[1]],
                  model@samples_metadata[,c("ER","PR","HER2","grade")])
tooth[,1:3]<-scale(tooth[,1:3])
tooth<-tooth[-which(tooth$PR %>% is.na()),]

margin<-theme(plot.margin = unit(c(-1,0.5,-3,-1), "cm"))

p4<-ggplot(tooth, aes(x=PR, y=Factor1, fill=PR)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor1") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5,
                     label.y = min(tooth$Factor1),size=8)

p5<-ggplot(tooth, aes(x=PR, y=Factor2, fill=PR)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor2") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5, 
                     label.y = min(tooth$Factor2), size=8)

p6<-ggplot(tooth, aes(x=PR, y=Factor13, fill=PR)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=50, face="bold"),
        axis.title=element_text(size=40,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor13") +  margin +
  theme(plot.title = element_text(hjust = 0.5, vjust = -4, size=40, face="bold"))

gridExtra::grid.arrange(p4, p5, p6, ncol = 1)

#--------------------------------------------#

tooth<-data.frame(get_factors(model, factors = c(1,2,13), 
                              as.data.frame = F)[[1]],
                  model@samples_metadata[,c("ER","PR","HER2","grade")])
tooth[,1:3]<-scale(tooth[,1:3])
tooth<-tooth[-which(tooth$HER2 %>% is.na()),]

margin<-theme(plot.margin = unit(c(-1,0.5,-3,-1), "cm"))

p7<-ggplot(tooth, aes(x=HER2, y=Factor1, fill=HER2)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor1") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5,
                     label.y = min(tooth$Factor1),size=8)

p8<-ggplot(tooth, aes(x=HER2, y=Factor2, fill=HER2)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor2") + margin +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  stat_compare_means(aes(label = paste0("p: ", ..p.format..)),
                     label.x = 0.5, 
                     label.y = min(tooth$Factor2), size=8)

p9<-ggplot(tooth, aes(x=HER2, y=Factor13, fill=HER2)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2"),
                    labels= c("Negative","Positve"))+
  theme(axis.text=element_text(size=50, face="bold"),
        axis.title=element_text(size=40,face="bold")) +
  scale_x_discrete(labels= c("", "")) +
  xlab("")+ylab("")+
  ggtitle("Factor13") + margin +
  theme(plot.title = element_text(hjust = 0.5, vjust = -4, size=40, face="bold"))

gridExtra::grid.arrange(p7, p8, ncol = 1)

#gridExtra::grid.arrange(p1, p2, nrow = 1)
#-------------------------------------#

tooth<-data.frame(get_factors(model, factors = c(1,2,13), 
                              as.data.frame = F)[[1]],
                  model@samples_metadata[,c("ER","PR","HER2","grade")])
tooth[,1:3]<-scale(tooth[,1:3])
tooth<-tooth[-which(tooth$grade %>% is.na()),]
tooth$grade<-tooth$grade %>% factor()

margin<-theme(plot.margin = unit(c(-1,0.5,-3,-1), "cm"))

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3"))
p10<-ggplot(tooth, aes(x=grade, y=Factor1, fill=grade)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2", "yellowgreen"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2","yellowgreen"),
                    labels= c("Grade 1","Grade 2", "Grade 3"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "","")) +
  xlab("")+ylab("")+
  ggtitle("Factor1") +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  margin+
  stat_compare_means(comparisons = my_comparisons, size=6)



p11<-ggplot(tooth, aes(x=grade, y=Factor2, fill=grade)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2","yellowgreen"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2","yellowgreen"),
                    labels= c("Grade 1","Grade 2", "Grade 3"))+
  theme(axis.text=element_text(size=30, face="bold"),
        axis.title=element_text(size=30,face="bold")) +
  scale_x_discrete(labels= c("", "","")) +
  xlab("")+ylab("")+
  ggtitle("Factor2") +
  theme(plot.title = element_text(hjust = 1, vjust = -10, 
                                  size=20, face="bold"))+
  margin+
  stat_compare_means(comparisons = my_comparisons, size=6)


p12<-ggplot(tooth, aes(x=grade, y=Factor13, fill=grade)) +
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("lightsalmon", "skyblue2","yellowgreen"))+
  geom_boxplot(width=0.1, color="black", lwd=2, show.legend = F) +
  theme(legend.position="none",
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=40),
        legend.key = element_rect(fill = NA, color = NA)) +
  scale_fill_manual(name = "", 
                    values=c("lightsalmon", "skyblue2","yellowgreen"),
                    labels= c("Grade 1","Grade 2", "Grade 3"))+
  theme(axis.text=element_text(size=50, face="bold"),
        axis.title=element_text(size=40,face="bold")) +
  scale_x_discrete(labels= c("", "","")) +
  xlab("")+ylab("")+
  ggtitle("Factor13") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -4, size=40, face="bold"))+margin

gridExtra::grid.arrange(p10, p11,ncol = 1)


gridExtra::grid.arrange(p1,p2,p4,p5,
                        p7,p8,
                        p10,p11,ncol = 2)




gridExtra::grid.arrange(a,b,c,d, ncol=2)

pl = replicate(3, ggplot(), FALSE)

#--------------------------------------------------------#
minmax <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
for(f in c(1,2,13)){
  W1<-get_weights(model, factors = f, views = 1, 
                  as.data.frame = TRUE)
  W1$feature<-c("??-glucose", "Ascorbate", "Lactate", "Tyrosine", "Glycine", 
                "Myoinositol", "Taurine", "Scylloinositol", 
                "GlyPcho_phosphocholine", "Phosphocholine", "Choline", 
                "Creatine", "Glutathione", "Glutamine", "Succinate", 
                "Glutamate", "Acetate", "Alanine")
  data.frame(c("??-glucose", "Ascorbate", "Lactate", "Tyrosine", "Glycine", 
               "Myoinositol", "Taurine", "Scylloinositol", 
               "GlyPcho_phosphocholine", "Phosphocholine", "Choline", 
               "Creatine", "Glutathione", "Glutamine", "Succinate", 
               "Glutamate", "Acetate", "Alanine"))
  W1.new<-W1[order(abs(W1$value), decreasing = T),]
  
  W1.new$value<-W1.new$value %>% abs() %>% minmax()
  
  
  W2<-get_weights(model, factors = f, views = 2, 
                  as.data.frame = TRUE)
  W2.new<-W2[order(abs(W2$value), decreasing = T),]
  W2.new$value<-W2.new$value %>% abs() %>% minmax()
  
  
  W3<-get_weights(model, factors = f, views = 3, 
                  as.data.frame = TRUE)
  W3.new<-W3[order(abs(W3$value), decreasing = T),]
  W3.new$value<-W3.new$value %>% abs() %>% minmax()
  
  W<-c(W1.new$value[1:18], W2.new$value[1:25], W3.new$value[1:25])
  label<-c(W1.new$feature[1:18], W2.new$feature[1:25] %>% as.character(), W3.new$feature[1:25]%>% as.character())
  
  label<-lapply(strsplit(label %>% as.character(), split = "_\\("), function(x){x[[1]][1]}) %>% unlist()
  label<-lapply(strsplit(label %>% as.character(), split = "_p"), function(x){x[[1]][1]}) %>% unlist()
  names(W)<-label
  
  par(mar=c(3,16.8,0.01,2))
  barplot(rev(W), horiz = T, space=1.205,
          col = c(rep("darkorange3",18),
                  rep("goldenrod",25),
                  rep("plum2",25)) %>% rev(),
          yaxt="n", width = 1,
          xaxt="n",
          main=paste("Factor ",f))
  
  axis(2, at = seq(1,150, length.out=68) %>% rev, las=1, font = 2,
       labels = names(W), 
       col = "black", 
       tick=T, cex.axis=2.2)
  axis(1, at = seq(0,1,0.2), las=1, font = 2,
       labels = seq(0,1,0.2), 
       col = "black", 
       tick=T, cex.axis=2)
  abline(v=0.9, lwd=3, col="darkblue", lty=3)
}


#-------------------------------------------------#

aa<-summary(res.cox.crude)

#saveRDS(aa, "/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/data/aa.Rdata")

aa<-cbind(aa$conf.int,aa$coefficients)

# make some mock data
colo<-ifelse(aa[,ncol(aa)] %>% round(digits = 2) <=0.05, "red", "black")
colo1<-ifelse(aa[,ncol(aa)] %>% round(digits = 2) <=0.05, "black", "grey25")

aaa<-aa[,c(1,4,3,9)] %>% as.data.frame()
colnames(aaa)<-c("coeff", "upper", "lower", "pval")
aaa$variable<-paste("Fact.", sprintf("%02d", 1:20), sep="")
aaa$odd<-rep(1,20)

aaa %>%
  ggplot(aes(x = variable, y = coeff, group = 1)) +
  geom_polygon(aes(y = upper), fill = "plum1", alpha = 0.9) +
  geom_polygon(aes(y = lower), fill = "white", alpha = 0.8) +
  
  geom_polygon(fill = NA, colour = "purple", lwd=3) +
  coord_polar() +
  theme(panel.grid.minor = element_blank()) + 
  labs(x = "", y = "")+
  theme(
    #axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=15, face = "bold", colour = "black"),
    #legend.key = element_blank(),
    #legend.title = element_blank(),
    legend.background = element_rect(color="#ffffff", fill="transparent"), ### neu check !!!
    panel.background = element_rect(fill = "white", colour = "white", size = 3, linetype = "solid"),
    panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "grey25"),
    axis.text.x=element_text(size=18, face = "bold", colour = colo, angle = 0),
    #axis.text.x=element_text(size=25, face = "bold", colour = colo, angle = 10)
    #plot.margin = margin(0.3, 1, 0.3, 1, unit="inch"),
    #panel.grid.minor = element_blank()
  )+
  geom_polygon(aes(y = odd), fill = NA, color="red", lwd=1, lty=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

# make some mock data
bb<-summary(res.cox.age)
bb<-cbind(bb$conf.int,bb$coefficients)

bb<-bb[-c(1),]
#colo<-ifelse(bb[,ncol(bb)] %>% round(digits = 2) <=0.05, "red", "grey25")
colo1<-ifelse(bb[,ncol(bb)] %>% round(digits = 2) <=0.05, "black", "grey25")

bbb<-bb[,c(1,4,3,9)] %>% as.data.frame()
colnames(bbb)<-c("coeff", "upper", "lower", "pval")

bbb$variable<-paste("Fact.", sprintf("%02d", 1:20), sep="")
bbb$odd<-rep(1,20)

bbb %>%
  ggplot(aes(x = variable, y = coeff, group = 1)) +
  geom_polygon(aes(y = upper), fill = "lightgoldenrod", alpha = 0.9) +
  geom_polygon(aes(y = lower), fill = "white", alpha = 0.8) +
  
  geom_polygon(fill = NA, colour = "orange2", lwd=3) +
  coord_polar() +
  theme(panel.grid.minor = element_blank()) + 
  labs(x = "", y = "")+
  theme(
    #axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=15, face = "bold", colour = "black"),
    #legend.key = element_blank(),
    #legend.title = element_blank(),
    legend.background = element_rect(color="#ffffff", fill="transparent"), ### neu check !!!
    panel.background = element_rect(fill = "white", colour = "white", size = 3, linetype = "solid"),
    panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "grey25"),
    axis.text.x=element_text(size=18, face = "bold", colour = colo, angle = 0),
    #plot.margin = margin(0.3, 1, 0.3, 1, unit="inch"),
    #panel.grid.minor = element_blank()
  )+
  geom_polygon(aes(y = odd), fill = NA, color="red", lwd=1, lty=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))



my.cluster<-read.csv("./ml_data/cluster.csv", header = T)[,1]
factors <- get_factors(model, as.data.frame = T) # factor matrix
factor<-matrix(factors$value, ncol=length(unique(factors$factor)), byrow = F)
#
my.cluster.code<-ifelse(my.cluster==1, "MOC1", ifelse(my.cluster==2, "MOC2", "MOC3"))
factor.group<-data.frame(factor, cluster=my.cluster.code %>% as.factor())

col.code<-c("blue", "red3", "darkgreen")
col.code.br<-ifelse(my.cluster==1, "lightblue", ifelse(my.cluster==2, "red", "yellowgreen"))

factor.group %>%
  ggplot(aes(x = X1,
             y = X2))+
  ggforce::geom_mark_ellipse(aes(fill = cluster),
                             colour = NA,
                             alpha = 0.2,
                             expand = unit(0.5,"mm"),
                             size = 0,
                             show.legend=F)+
  geom_point(aes(x = X1, y = X2, color = cluster), size=4)+
  scale_fill_manual(breaks = unique(factor.group$cluster) %>% rev(),
                    values = col.code)+
  scale_colour_manual(breaks = unique(factor.group$cluster) %>% rev(),
                      values = col.code)+
  #geom_point(shape = 1,size = 0.1,colour = col.code.br, stroke = 1.5)+
  labs(title = "", 
       x = "MOFA+ Latent Factor1\n Metabolites (3.9%), Proteins (15.1%), Transcripts (9.8%)", 
       y = "MOFA+ Latent Factor2\n Metabolites (5.9%), Proteins (6.8%), Transcripts (5.4%)", 
       color = "Multi-Omics\nClusters")+
  #theme_bw() +
  theme(axis.text.x = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 20, face = "bold", color = "black"),
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.y = element_text(size = 20, face = "bold", color = "black"),
        #plot.title = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 18, face = "bold", color = "black"),
        legend.position=c(.88,0.12))


#------------------------------------------------------------#
my.cluster<-read.csv("./ml_data/cluster.csv", header = T)[,1]
data.order<-readRDS("data_order_all.rds")
cln.order<-readRDS("cln_order_all.rds")

clinic.oslo<-data.frame(cln.order, cluster=sort(my.cluster) %>% as.factor())
na.index<-clinic.oslo$Group %>% is.na() %>% which()
clinic.oslo<-clinic.oslo[-na.index,]

clinic.oslo$Group<-clinic.oslo$Group %>% as.factor()
clinic.oslo<-clinic.oslo[order(clinic.oslo$cluster, clinic.oslo$Group),]

color1<-c("palevioletred4","deepskyblue","darkorange3","palevioletred1","yellow4")
tab2<-clinic.oslo %>% group_by(cluster) %>% summarise(n = n())

gap<-10

new.clst<-c(rep(1,tab2$n[1]),rep(0,gap),rep(2,tab2$n[2]),rep(0,gap),rep(3,tab2$n[3]), rep(0,gap))
color.inside<-ifelse(new.clst==1,"blue", ifelse(new.clst==2, "red3", ifelse(new.clst==3, "darkgreen","white")))

color.outside<-new.clst
color.outside[which(new.clst!=0)]<-clinic.oslo$Group %>% as.character()
color.outside<-ifelse(color.outside=="Basal", color1[1], 
                      ifelse(color.outside=="Her2", color1[2],
                             ifelse(color.outside=="LumA", color1[3],
                                    ifelse(color.outside=="LumB", color1[4],
                                           ifelse(color.outside=="Normal", color1[5],"white")))))

pie(rep(1, length(color.inside)), border = NA, radius = 0.9, labels = NA,
    col=color.outside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.37, labels = NA, 
    col = "white")
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.35, labels = NA, 
    col = color.inside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.20, labels = NA, 
    col = "white")
legend(0.9, 1, 
       c( "Basal", "Her2", "Luminal A", "Luminal B", "Normal"), 
       cex = 4, fill = color1, box.lwd = 0,
       bty="n", text.font=2, border = color1,
       text.col = color1)




clinic.oslo<-data.frame(cln.order, cluster=sort(my.cluster) %>% as.factor())
na.index<-clinic.oslo$Named.Group %>% is.na() %>% which()
clinic.oslo<-clinic.oslo[-na.index,]

clinic.oslo$Group<-clinic.oslo$Named.Group %>% as.factor()
clinic.oslo<-clinic.oslo[order(clinic.oslo$cluster, clinic.oslo$Named.Group),]

color1<-c("palevioletred4","deepskyblue","darkorange3","palevioletred1","yellow4")
tab2<-clinic.oslo %>% group_by(cluster) %>% summarise(n = n())

gap<-10

new.clst<-c(rep(1,tab2$n[1]),rep(0,gap),rep(2,tab2$n[2]),rep(0,gap),rep(3,tab2$n[3]), rep(0,gap))
color.inside<-ifelse(new.clst==1,"blue", ifelse(new.clst==2, "red3", ifelse(new.clst==3, "darkgreen","white")))

color.outside<-new.clst
color.outside[which(new.clst!=0)]<-clinic.oslo$Named.Group %>% as.character()
color.outside<-ifelse(color.outside=="Basal", color1[1], 
                      ifelse(color.outside=="Her2", color1[2],
                             ifelse(color.outside=="Liminal", color1[3],
                                    ifelse(color.outside=="ReacI", color1[4],
                                           ifelse(color.outside=="ReacII", color1[5],"white")))))

pie(rep(1, length(color.inside)), border = NA, radius = 0.9, labels = NA,
    col=color.outside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.37, labels = NA, 
    col = "white")
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.35, labels = NA, 
    col = color.inside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.20, labels = NA, 
    col = "white")

legend(0.9, 1, 
       c( "Basal", "Her2", "Luminal", "Reactive I", "Reactive II"), 
       cex = 4, fill = color1, box.lwd = 0,
       bty="n", text.font=2, border = color1,
       text.col = color1)



#----------------------------------------------#


data.order<-readRDS("data_order_all.rds")
cln.order<-readRDS("cln_order_all.rds")

clinic.oslo<-data.frame(cln.order, cluster=sort(my.cluster) %>% as.factor())
na.index<-c(clinic.oslo$ER..1.negative...1....2..1.10...3.10.50...4...50...5.positive..98...not.defined. %>% is.na() %>% which(),
            c(clinic.oslo$ER..1.negative...1....2..1.10...3.10.50...4...50...5.positive..98...not.defined.==98) %>% which())
clinic.oslo$Group_er<-clinic.oslo$ER..1.negative...1....2..1.10...3.10.50...4...50...5.positive..98...not.defined.
clinic.oslo$Group_er[na.index]<-98
clinic.oslo$Group_er<-ifelse(clinic.oslo$Group_er==2,1,clinic.oslo$Group_er)

clinic.oslo$Group<-clinic.oslo$Group_er %>% as.factor()
clinic.oslo<-clinic.oslo[order(clinic.oslo$cluster, clinic.oslo$Group_er),]

color1<-c("pink1","hotpink2","hotpink3","hotpink4","grey25")
tab2<-clinic.oslo %>% group_by(cluster) %>% summarise(n = n())

gap<-10

new.clst<-c(rep(1,tab2$n[1]),rep(0,gap),rep(2,tab2$n[2]),rep(0,gap),rep(3,tab2$n[3]), rep(0,gap))
color.inside<-ifelse(new.clst==1,"blue", ifelse(new.clst==2, "red3", ifelse(new.clst==3, "darkgreen","white")))

color.outside<-new.clst
color.outside[which(new.clst!=0)]<-clinic.oslo$Group_er %>% as.character()
color.outside<-ifelse(color.outside=="1", color1[1], 
                      ifelse(color.outside=="3", color1[2],
                             ifelse(color.outside=="4", color1[3],
                                    ifelse(color.outside=="5", color1[4],
                                           ifelse(color.outside=="98", color1[5],"white")))))

pie(rep(1, length(color.inside)), border = NA, radius = 0.9, labels = NA,
    col=color.outside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.37, labels = NA, 
    col = "white")
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.35, labels = NA, 
    col = color.inside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.20, labels = NA, 
    col = "white")
legend(0.9, 1, title = "ER cell %", title.col = "Black",
       c( "Negative (<1%)", "1-10%", "10-50%", "Positive (???50%)", "NA"), 
       cex = 3.5, fill = color1, box.lwd = 0,
       bty="n", text.font=2, border = color1,
       text.col = color1)



clinic.oslo<-data.frame(cln.order, cluster=sort(my.cluster) %>% as.factor())
na.index<-c(clinic.oslo$PR..1.negative...1....2..1.10...3.10.50...4...50...98...not.defined. %>% is.na() %>% which(),
            c(clinic.oslo$PR..1.negative...1....2..1.10...3.10.50...4...50...98...not.defined.==98) %>% which())
clinic.oslo$Group_er<-clinic.oslo$PR..1.negative...1....2..1.10...3.10.50...4...50...98...not.defined.
clinic.oslo$Group_er[na.index]<-98

clinic.oslo$Group<-clinic.oslo$Group_er %>% as.factor()
clinic.oslo<-clinic.oslo[order(clinic.oslo$cluster, clinic.oslo$Group_er),]

color1<-c("pink1","hotpink2","hotpink3","hotpink4","grey25")
tab2<-clinic.oslo %>% group_by(cluster) %>% summarise(n = n())

gap<-10

new.clst<-c(rep(1,tab2$n[1]),rep(0,gap),rep(2,tab2$n[2]),rep(0,gap),rep(3,tab2$n[3]), rep(0,gap))
color.inside<-ifelse(new.clst==1,"blue", ifelse(new.clst==2, "red3", ifelse(new.clst==3, "darkgreen","white")))

color.outside<-new.clst
color.outside[which(new.clst!=0)]<-clinic.oslo$Group_er %>% as.character()
color.outside<-ifelse(color.outside=="1", color1[1], 
                      ifelse(color.outside=="2", color1[2],
                             ifelse(color.outside=="3", color1[3],
                                    ifelse(color.outside=="4", color1[4],
                                           ifelse(color.outside=="98", color1[5],"white")))))

pie(rep(1, length(color.inside)), border = NA, radius = 0.9, labels = NA,
    col=color.outside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.37, labels = NA, 
    col = "white")
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.35, labels = NA, 
    col = color.inside)
par(new=T)
pie(rep(1,length(color.inside)), border = NA, radius = 0.20, labels = NA, 
    col = "white")
legend(0.9, 1, title = "PR cell %", title.col = "Black",
       c( "Negative (<1%)", "1-10%", "10-50%", "Positive (???50%)", "NA"), 
       cex = 3.5, fill = color1, box.lwd = 0,
       bty="n", text.font=2, border = color1,
       text.col = color1)



#----------------------------------------------------#

data.order<-readRDS("data_order_all.rds")
cln.order<-readRDS("cln_order_all.rds")
model <- load_model("MOFA_allgene_20F.hdf5")

cln.order$Age %>% as.numeric() %>% sd(na.rm=T)
cln.order$Age %>% as.numeric() %>% mean(na.rm=T)

clinic.oslo<-data.frame(cln.order, cluster=sort(my.cluster) %>% as.factor())
clinic.oslo<-clinic.oslo[order(clinic.oslo$cluster, clinic.oslo$Group),]

df<-clinic.oslo %>% group_by(cluster, Group) %>% summarise(n = n()) %>%mutate(freq = n*100 / sum(n)) %>% as.data.frame()

newrow<-data.frame(cluster=c(1,1,1),
                   Group=c("LumA", "LumB", "Normal"),
                   n=c(0,0,0),
                   freq=c(0,0,0))

df <- rbind(df[1:2,], newrow, df[-(1:2),])

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3", "grey45")
w1 <- waffle(df %>% filter(cluster==1) %>% dplyr::select(Group,freq), rows = 4, colors = color, size=1.5)
w2 <- waffle(df %>% filter(cluster==2) %>% dplyr::select(Group,freq), rows = 4, colors = color, size=1.5)
w3 <- waffle(df %>% filter(cluster==3) %>% dplyr::select(Group,freq), rows = 4, colors = color, size=1.5)

# Combine the plots
iron(w1, w2, w3)


#--------------------------------------------------#

index.na<-c(data.order[[1]][,1] %>% is.na() %>% which(),
            data.order[[2]][,1] %>% is.na() %>% which(),
            data.order[[3]][,1] %>% is.na() %>% which()) %>% unique()


prot<-data.order$metabolite[-index.na, ]
#prot<-prot[,-ncol(prot)] %>% as.matrix() %>% log2()
colnames(prot)<-c("??-glucose", "Ascorbate", "Lactate", "Tyrosine", "Glycine", 
                  "Myoinositol", "Taurine", "Scylloinositol", 
                  "GlyPcho", "Pcho", "Choline", 
                  "Creatine", "Glutathione", "Glutamine", "Succinate", 
                  "Glutamate", "Acetate", "Alanine")

prot<-apply(prot,2,function(x){log2(x/mean(x, na.rm=T))})
prot<-prot[,order(colSums(prot), decreasing = F)]

par(oma=c(6,0,0,2))
labCol = lapply(strsplit(colnames(prot),
                         "_\\("), function(x){x[1]}) %>% unlist()
labCol<-as.expression(lapply(labCol, function(a) bquote(bold(.(a)))))
heatmap.2( prot %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 2, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(10,15),
           lwid = c(10,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(100, "gold","black", "magenta") ,
           col= colorRampPalette(c("black","purple4","mediumpurple4","yellow1",
                                   "orangered1", "darkred"))(n = 100),
           
           font.lab=3,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol ,
           main = "metabolite",
           na.color = "grey85")

pr<-cbind(prot, clust=sort(my.cluster)[-index.na] %>% as.factor()) %>% as.data.frame()
pval<-c()
for(i in 1:ncol(prot)){
  a<-aov(pr[,ncol(pr)]~pr[,i])
  b<-summary(a)
  pval<-c(pval, b[[1]][["Pr(>F)"]][1])
}
d1<-data.frame(metabolite=colnames(prot), pvalue=pval %>% p.adjust() %>% round(digits = 3))



W1<-get_weights(model, factors = 1, views = 2, 
                as.data.frame = TRUE)
W1.new<-W1[order(abs(W1$value), decreasing = T),]
W1.new$value<-W1.new$value %>% abs() %>% minmax()


W2<-get_weights(model, factors = 2, views = 2, 
                as.data.frame = TRUE)
W2.new<-W2[order(abs(W2$value), decreasing = T),]
W2.new$value<-W2.new$value %>% abs() %>% minmax()

imp<-c(W1.new[1:20,1], W2.new[1:20,1]) %>% unique()
imp<-c(W1.new[1:25,1])

prot<-data.frame(data.order[[2]],
                 cluster=sort(my.cluster))

colnames(prot)<- gsub("_\\.", "_(", prot %>% colnames())
colnames(prot)<- gsub("\\._", ")_", prot %>% colnames())

prot<-prot[-index.na, c(colnames(prot)%in%c("p53_(V)_R", c(imp[1:21]) %>% as.character()) %>% which())]


#prot<-prot[,-ncol(prot)] %>% as.matrix() %>% log2()

prot<-prot[,order(apply(prot,2,var), decreasing = T)]

pal <- colorpanel(15, "green", "black", "red") 

par(oma=c(6,0,0,2))
labCol = lapply(strsplit(colnames(prot),
                         "_\\("), function(x){x[1]}) %>% unlist()
labCol<-as.expression(lapply(labCol, function(a) bquote(bold(.(a)))))
heatmap.2( prot %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 2, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(10,15),
           lwid = c(10,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(100, "gold","black", "magenta") ,
           col= colorRampPalette(c("orangered","orange","yellow", "black",
                                   "turquoise4","turquoise1", "turquoise1") %>% rev())(n = 100),
           font.lab=3,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol ,
           main = "metabolite",
           na.color = "grey85")
pr<-cbind(prot, clust=sort(my.cluster)[-index.na]) %>% as.data.frame()
pval<-c()
for(i in 1:ncol(prot)){
  a<-aov(pr[,ncol(pr)]~pr[,i])
  b<-summary(a)
  pval<-c(pval, b[[1]][["Pr(>F)"]][1])
}
d2<-data.frame(protein=colnames(prot), pvalue=pval %>% p.adjust() %>% round(digits = 100))


W1<-get_weights(model, factors = 1, views = 3, 
                as.data.frame = TRUE)
W1.new<-W1[order(abs(W1$value), decreasing = T),]
W1.new$value<-W1.new$value %>% abs() %>% minmax()


W2<-get_weights(model, factors = 2, views = 3, 
                as.data.frame = TRUE)
W2.new<-W2[order(abs(W2$value), decreasing = T),]
W2.new$value<-W2.new$value %>% abs() %>% minmax()

imp<-c(W1.new[1:20,1], W2.new[1:20,1]) %>% unique()
imp<-c(W1.new[1:20,1])

prot<-data.frame(data.order[[3]],
                 cluster=sort(my.cluster))

#colnames(prot)<- gsub("_\\.", "_(", prot %>% colnames())
#colnames(prot)<- gsub("\\._", ")_", prot %>% colnames())

prot<-prot[-index.na, c(colnames(prot)%in%imp[1:20] %>% which())]
#prot<-prot[,-ncol(prot)] %>% as.matrix() %>% log2()

prot<-prot[,order(colSums(prot),apply(prot,2,var), decreasing = T)]

pal <- colorpanel(15, "green", "black", "red") 

par(oma=c(6,0,0,2))
labCol = lapply(strsplit(colnames(prot),
                         "_\\("), function(x){x[1]}) %>% unlist()
labCol<-as.expression(lapply(labCol, function(a) bquote(bold(.(a)))))
heatmap.2( prot %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(10,15),
           lwid = c(10,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col= colorRampPalette(c("gold","yellow","black", "magenta", "pink") %>% rev())(n = 100),
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")
pr<-cbind(prot, clust=sort(my.cluster)[-index.na]) %>% as.data.frame()
pval<-c()
for(i in 1:ncol(prot)){
  a<-aov(pr[,ncol(pr)]~pr[,i])
  b<-summary(a)
  pval<-c(pval, b[[1]][["Pr(>F)"]][1])
}
d3<-data.frame(gene=colnames(prot), pvalue=pval %>% p.adjust() %>% round(digits = 100))



par(oma=c(6,0,0,2))
clst<-cbind(sort(my.cluster)[-index.na],
            sort(my.cluster)[-index.na])

heatmap.2( clst %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = T, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(1,15),
           lwid = c(1,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col=colorpanel(3, "blue", "red3", "darkgreen") ,
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")




#-----------------------------------------------------------------------#

cln.order<-readRDS("cln_order_all.rds")
my.cluster<-read.csv("./ml_data/cluster.csv", header = T)[,1]
death_type<-list(c(1:4),1,c(1:4),1,c(1:4))
wrt<-c("Group","Group", "Cluster", "Cluster", "Grade")

survival.formula<-sapply(wrt,
                         function(x) as.formula(paste('Surv(Time, Status)~', x)))
model.srv<-list()
p.val<-c()
for(k in 1:5){
  
  cli.data<-cln.order[,c(1,24,27,13,2)]
  cli.data<-cbind(cli.data, my.cluster %>% sort())
  names(cli.data)<-c("ID", "Status", "Time","Grade", "Group", "Cluster")
  
  bc_death<-cln.order$`Cause of death (1=breast cancer, 2=different disease, 3=unknown, 4=cause of death not updated (dead after des2017)`
  bc_death<-ifelse(bc_death %in% death_type[[k]], 1, ifelse(bc_death != 0, NA, 0))
  
  cli.data$Status<-bc_death
  
  #cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
  #cli.data<-na.omit(cli.data)
  cli.data$Time<-as.numeric(cli.data$Time)
  
  
  model.srv[[k]]<-survfit(survival.formula[[k]], data=cli.data)
  p.val<-c(p.val, surv_pvalue(survfit(survival.formula[[k]], data=cli.data))$pval.txt)
  
}

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

for(i in 3:4){
  ggsurv<-ggsurvplot(model.srv[[i]], palette= col.cluster, 
                     size.lab=1.5, 
                     pval = p.val[i], 
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
                     legend.labs = paste("MOC", 1:3, sep=""),
                     ggtheme = theme,
                     risk.table.fontsize = 6,
                     conf.int.alpha=0.15
  )
  
  ggsurv$table <- ggrisktable(model.srv[[i]], 
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
}

#---------------------wald test-----------#
wald.p<-c()
comb<-list(c("MOC1","MOC2"),c("MOC1","MOC3"),c("MOC2","MOC3"))

for(i in 1:3){
  
  cli.data<-cln.order[,c(1,42,27,13,2,31)] %>% as.data.frame()
  
  cli.data$MOC<-ifelse(cli.data$Cluster==1,"MOC1",
                       ifelse(cli.data$Cluster==2,"MOC2",
                              ifelse(cli.data$Cluster==3,"MOC3",NA)))%>% sort()
  
  names(cli.data)<-c("ID", "Status", "Time","Grade", "Group","Cluster", "MOC")
  
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
  summary(model.srv1) %>% print()
  wald.p<-c(wald.p, p.val1$pval.txt)
}
wald.p


n<-3
ref.grp<-2
df<-expand.grid(1:n, 1:n)
#df<-df[-which(df$Var1==n),]

line.coord<-combinat::combn(1:n, 2)

#line.coord<-line.coord[,-c(which(line.coord[1,]!=ref.grp & line.coord[2,]!=ref.grp))]

line_df <- data.frame(
  x = c(1:(n)),
  y = line.coord[1,],
  xend = c(1:(n)),
  yend = line.coord[2,]
)

exp.grid<-expand.grid(1:(n),1:n) %>% as.matrix()
exp.grid<-paste(exp.grid[,1], exp.grid[,2])

my.grid<-cbind(c(1:(n),1:(n)) %>% sort(), c(line.coord))
my.grid<-paste(my.grid[,1], my.grid[,2])

colo<-rep("black", n*(n))
colo[which(exp.grid%in%my.grid)]<-"white"

ggplot(df, aes(Var1, Var2, fill = Var2 %>% as.factor())) +
  geom_tile(color = "white" , size=1) +
  #geom_text(aes(label = Var2), size = 10) +
  scale_y_reverse() +
  geom_point(size=45, color=colo)+
  #geom_line(color = "black", size=3)+
  coord_equal() +
  theme_void() +
  scale_fill_manual(values = c("blue", "red3", "green4","orange","purple"), guide = guide_none())+
  #geom_vline(yintercept=0, color = "grey", size=2)+
  geom_segment(
    data = line_df, 
    mapping = aes(x=x, y=y, xend=xend, yend=yend), 
    inherit.aes = FALSE,
    size=10, color="white"
  )
#+theme(axis.text.x = element_text(angle = 90, hjust = 1, lebel = c(wald.p %>% round(digits = 2))))


#------------------------------------------------#

cli.data$MOC<-cli.data$Cluster %>% as.factor()
levels(cli.data$MOC)<-c("MOC1", "MOC2", "MOC3")
cli.data$MOC<-relevel(cli.data$MOC, "MOC2")

res.cox <- coxph(Surv(Time, Status) ~ MOC, data =  cli.data, ties = "efron")
print(forest_model(res.cox))

#---------------------Figure 4 c----------------------------#

cli.data<-cln.order[,c(1,24,27,13,2)]

cli.data<-data.frame(cli.data, 
                     MOC= ifelse(my.cluster==1,"MOC1",
                                 ifelse(my.cluster==2,"MOC2",
                                        ifelse(my.cluster==3,"MOC3",NA)))%>% sort() %>% as.factor())


names(cli.data)<-c("ID", "Status", "Time","Grade", "Group", "MOC")

bc_death<-cln.order$`Cause of death (1=breast cancer, 2=different disease, 3=unknown, 4=cause of death not updated (dead after des2017)`
bc_death<-ifelse(bc_death %in% c(1:4), 1, ifelse(bc_death != 0, NA, 0))

cli.data$Status<-bc_death
cli.data$MOC<-relevel(cli.data$MOC, "MOC2")

#cli.data<-cli.data[which(m.data.cln$Grade %in% c("1", "2")),]
#cli.data<-na.omit(cli.data)
cli.data$Time<-as.numeric(cli.data$Time)

res.cox <- coxph(Surv(Time, Status) ~ MOC, data =  cli.data, ties = "efron")
res.cox %>% summary()
print(forest_model(res.cox))


#----------------------------------------------------------------#
for(i in 1:2){
  ggsurv<-ggsurvplot(model.srv[[i]], palette= col.group, 
                     size.lab=1.5, 
                     pval = p.val[i], 
                     font.size=1.5,
                     surv.scale = "percent",
                     font.legend = list(size = 16, color = "black", face = "bold"),
                     font.tickslab = c(20, "plain", "black"),
                     risk.table = T, 
                     pval.method = T,
                     legend.title="",
                     font.x = c(20, "bold", "black"),
                     font.y = c(20, "bold", "black"),
                     pval.size = 7, pval.coord = c(60, .1),
                     conf.int = T,
                     legend.labs = cli.data$Group %>% unique() %>% sort(),
                     ggtheme = theme,
                     risk.table.fontsize = 6,
                     alpha = 100,
                     conf.int.alpha=0.15,
                     xlab="Month"
  )
  
  ggsurv$table <- ggrisktable(model.srv[[i]], 
                              data = cli.data, 
                              color = "Group", 
                              palette=color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3"),
                              y.text = T,   ylab = "",  xlab = "",
                              tables.theme = theme,
                              font.tickslab = c(15, "bold"),
                              legend.labs = cli.data$Group %>% unique() %>% sort(),
                              #ggtheme = theme_survminer(),
                              fontsize = 6
                              
  )
  
  ggsurv %>% print()
}

ggsurvp


#---------------------------------------------------------------------------#

#validation

clinic<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/data_clinical_patient.txt", header = T, sep="\t")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val_cluster.txt", header = T)

merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  merge.clinic$SEX=="Female" 
  & merge.clinic$HISTORY_NEOADJUVANT_TRTYN=="No"
  #merge.clinic$RADIATION_THERAPY=="No"
),]
merge.clinic %>% nrow()
merge.clinic$RADIATION_THERAPY %>% table()

which(merge.clinic$OS_STATUS=="1:DECEASED" & merge.clinic$PERSON_NEOPLASM_CANCER_STATUS=="With Tumor")
merge.clinic$OS_STATUS %>% table()

merge.clinic$OS_STATUS<-ifelse(merge.clinic$OS_STATUS == "1:DECEASED", 1, 0)
merge.clinic$DSS_STATUS <- ifelse(merge.clinic$OS_STATUS == 1 & merge.clinic$DSS_STATUS == "1:DEAD WITH TUMOR", 1, 
                                  ifelse(merge.clinic$OS_STATUS == 1 & merge.clinic$DSS_STATUS == "0:ALIVE OR DEAD TUMOR FREE", NA, 0))

merge.clinic$SUBTYPE<-gsub("BRCA_", "", merge.clinic$SUBTYPE)
merge.clinic$SUBTYPE<-gsub("Her2", "HER2", merge.clinic$SUBTYPE)
merge.clinic$SUBTYPE<-ifelse(merge.clinic$SUBTYPE=="", NA, merge.clinic$SUBTYPE)

merge.clinic$CANCER_TYPE_ACRONYM
#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="White"),]

#----------survival-------------
col.cluster<-c("blue", "red3", "green4", "black", "orange")
col.group<-color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3")


cli.data<-merge.clinic[,c("PATIENT_ID", "DSS_STATUS",  "OS_MONTHS", "cluster", "SUBTYPE")] #DSS status for only BC, OS status for all cause

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

ggsurv<-ggsurvplot(model.srv, palette= col.cluster, 
                   size.lab=1.5, 
                   xlim=c(0,150),
                   #pval = T,
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
                   conf.int.alpha=0.15
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
                   xlim=c(0,150),
                   #pval = T,
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
                   legend.labs = cli.data1$Subtype %>% unique() %>% sort(),
                   ggtheme = theme,
                   risk.table.fontsize = 6,
                   conf.int.alpha=0.15
                   
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



df<-cli.data %>% group_by(Cluster, Subtype) %>% summarise(n = n()) %>%mutate(freq = n*100 / sum(n)) %>% as.data.frame()

newrow<-data.frame(Cluster=c(1,1,1,2),
                   Subtype=c("LumA", "LumB", "Normal", "Normal"),
                   n=c(0,0,0,0),
                   freq=c(0,0,0,0))

df <- rbind(df[1:2,], newrow, df[-(1:2),])

color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3", "grey45")
w1 <- waffle(df %>% filter(Cluster==1) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)
w2 <- waffle(df %>% filter(Cluster==2) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)
w3 <- waffle(df %>% filter(Cluster==3) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)

# Combine the plots
iron(w1, w2, w3)




cli.data<-merge.clinic[,c("PATIENT_ID", "DSS_STATUS",  "OS_MONTHS", "cluster")] 

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

#cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()

cli.data$Time<-as.numeric(cli.data$Time)

end.point<-200
cli.data1 <- cli.data %>% 
  mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
         Status2 = ifelse(Time >= end.point, Status, Status))

model.srv<- survfit(Surv(Time2, Status2)~Cluster, data=cli.data1)

ggsurv<-ggsurvplot(model.srv, palette= col.cluster, 
                   size.lab=1.5, 
                   xlim = c(0, 150),
                   #pval = T,
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
                   conf.int.alpha=0.15
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
ggsurv



#----------survival-------------

wald.p<-c()
comb<-list(c(1,2),c(1,3),c(2,3))

for(i in 1:3){
  
  cli.data<-merge.clinic[,c("PATIENT_ID", "DSS_STATUS",  "OS_MONTHS", "cluster", "SUBTYPE")] 
  #cli.data<-na.omit(cli.data)
  names(cli.data)<-c("ID", "Status", "Time", "Cluster", "Subgroup")
  
  cli.data<-cli.data %>% filter(Cluster%in%comb[[i]])
  cli.data<-na.omit(cli.data)
  cli.data$Cluster<-cli.data$Cluster%>%as.factor()
  cli.data$Time<-as.numeric(cli.data$Time)
  
  end.point<-175
  cli.data1 <- cli.data %>% 
    mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
           Status2 = ifelse(Time >= end.point, Status, Status))
  
  model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)
  
  p.val1<-surv_pvalue(model.srv, data=cli.data)
  summary(model.srv) %>% print()
  wald.p<-c(wald.p, p.val1$pval.txt)
}
wald.p


#----------------------------------------------------------------------------------------#

clinic<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/data_clinical_patient.txt", header = T, sep="\t")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/val_cluster.txt", header = T)

merge.clinic<-left_join(val.cluster, clinic, by = "PATIENT_ID")
merge.clinic$OS_MONTHS<-gsub(",", "\\.", merge.clinic$OS_MONTHS) %>% as.numeric()

#merge.clinic<-merge.clinic[which(merge.clinic$RACE %in% c("WHITE","[Not Available]")),]
merge.clinic<-merge.clinic[which(#merge.clinic$RACE=="White"&
  merge.clinic$SEX=="Female"
  #&merge.clinic$HORMONE_THERAPY=="NO"#
  &merge.clinic$CHEMOTHERAPY=="NO"###
  &merge.clinic$RADIO_THERAPY=="NO"###
  #&merge.clinic$BREAST_SURGERY==""
),]
merge.clinic %>% nrow()

merge.clinic$OS_STATUS<-ifelse(merge.clinic$OS_STATUS == "1:DECEASED", 1, 0)
merge.clinic$DSS_STATUS <- ifelse(merge.clinic$VITAL_STATUS == "Died of Disease", 1, 
                                  ifelse(merge.clinic$VITAL_STATUS == "Died of Other Causes", NA, 0))

merge.clinic$CLAUDIN_SUBTYPE<-ifelse(merge.clinic$CLAUDIN_SUBTYPE%in%c("claudin-low","NC"), NA, merge.clinic$CLAUDIN_SUBTYPE)
merge.clinic$CLAUDIN_SUBTYPE<-gsub("Her2", "HER2", merge.clinic$CLAUDIN_SUBTYPE)

#merge.clinic<-merge.clinic[which(merge.clinic$RACE=="White"),]

#----------survival-------------

cli.data<-merge.clinic[,c("PATIENT_ID", "DSS_STATUS",  "OS_MONTHS", "cluster", "CLAUDIN_SUBTYPE")] 

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
                   pval = T,
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
                   risk.table.fontsize = 6
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
                   pval = T,
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
                   conf.int.alpha=0.15
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


df<-cli.data %>% group_by(Cluster, Subtype) %>% summarise(n = n()) %>%mutate(freq = n*100 / sum(n)) %>% as.data.frame()

newrow<-data.frame(Cluster=c(1,1,1),
                   Subtype=c("LumA", "LumB", "Normal"),
                   n=c(0,0,0),
                   freq=c(0,0,0))

df <- rbind(df[1:2,], newrow, df[-(1:2),])

color<-c("violetred3", "mediumpurple4", "goldenrod", "darkgoldenrod4", "lightsalmon3", "grey45")
w1 <- waffle(df %>% filter(Cluster==1) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)
w2 <- waffle(df %>% filter(Cluster==2) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)
w3 <- waffle(df %>% filter(Cluster==3) %>% select(Subtype,freq), rows = 4, colors = color, size=1.5)

# Combine the plots
iron(w1, w2, w3)



wald.p<-c()
comb<-list(c(1,2),c(1,3),c(2,3))

for(i in 1:3){
  
  cli.data<-merge.clinic[,c("PATIENT_ID", "OS_STATUS",  "OS_MONTHS", "cluster", "CLAUDIN_SUBTYPE")] 
  
  names(cli.data)<-c("ID", "Status", "Time", "Cluster", "Subtype")
  
  cli.data<-cli.data %>% filter(Cluster%in%comb[[i]])
  cli.data<-na.omit(cli.data)
  cli.data$Cluster<-cli.data$Cluster%>%as.factor()
  cli.data$Time<-as.numeric(cli.data$Time)
  
  end.point<-175
  cli.data1 <- cli.data %>% 
    mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
           Status2 = ifelse(Time >= end.point, Status, Status))
  
  model.srv<- survfit(Surv(Time, Status)~Cluster, data=cli.data)
  
  p.val1<-surv_pvalue(model.srv, data=cli.data)
  summary(model.srv) %>% print()
  wald.p<-c(wald.p, p.val1$pval.txt)
}
wald.p









cli.data<-merge.clinic[,c("PATIENT_ID", "DSS_STATUS",  "OS_MONTHS", "cluster")] 

names(cli.data)<-c("ID", "Status", "Time", "Cluster")

cli.data<-na.omit(cli.data)

cli.data$Cluster<-cli.data$Cluster%>%as.factor()

cli.data$Time<-as.numeric(cli.data$Time)

end.point<-300
cli.data1 <- cli.data %>% 
  mutate(Time2 = ifelse(Time >= end.point, end.point, Time),
         Status2 = ifelse(Time >= end.point, Status, Status))

model.srv<- survfit(Surv(Time2, Status2)~Cluster, data=cli.data1)
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
                   xlim=c(0,200),
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
                   conf.int.alpha=0.15
                   
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
                            fontsize = 6,
                            xlim=c(0,200)
                            
)
ggsurv %>% print()
























prot<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val2018_gene_test.csv", header = T, sep=",")
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val_cluster.txt", header = T)
prot<-t(prot[which(prot$V1%in%imp[1:20]),])
colnames(prot)<-prot[1,]
prot<-prot[-1,]
prot<-data.frame(prot)
nm<-row.names(prot)
prot<-apply(prot, 2, as.numeric)
row.names(prot)<-nm
prot<-prot %>% as.data.frame()

rownames(val.cluster)<-paste(gsub("-", ".", val.cluster$PATIENT_ID), ".01", sep="")
prot$name<-rownames(prot)
val.cluster$name<-rownames(val.cluster)

prot<-left_join(val.cluster, prot, by="name")
prot<-prot[order(prot$cluster, decreasing = F),-c(1:3)]


nm1<-c("ROPN1", "C4orf7", "KRT6B", "CRYAB",
       "GABRP", "ANKRD30A", "AGR3", "TFF3",
       "MLPH", "ESR1", "MAPT", "AGR2", "NAT1",
       "SCUBE2",  "GATA3", "FOXA1", "CA12",  "AR") %>% rev()

prot<-prot[,nm1]

#prot<-prot[,order(colSums(prot),apply(prot,2,var), decreasing = F)] %>% as.data.frame()

par(oma=c(6,0,0,2))
labCol = lapply(strsplit(colnames(prot),
                         "_\\("), function(x){x[1]}) %>% unlist()
labCol<-as.expression(lapply(labCol, function(a) bquote(bold(.(a)))))
heatmap.2( prot %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = F, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(1,15),
           lwid = c(1,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col= colorRampPalette(c("gold","yellow","black", "magenta", "pink") %>% rev())(n = 100),
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")

par(oma=c(6,0,0,2))
heatmap.2( cbind(val.cluster$cluster[order(val.cluster$cluster, decreasing = F)],
                 val.cluster$cluster[order(val.cluster$cluster, decreasing = F)])
           %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = F, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(1,15),
           lwid = c(1,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col= c("blue","red", "green"),
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")




prot<-read.csv("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/gene_zscore.csv", header = T)
val.cluster<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/metabric.data/val_cluster.txt", header = T)
#prot<-read.table("/mnt/work/workbench/abhibhas/oslo2_as/project/Oslo2/val.data.2018/val2018_gene_test.csv", header = T, sep=",")
prot<-t(prot[which(prot$V1%in%imp[1:20]),])
colnames(prot)<-prot[1,]
prot<-prot[-c(1:2),]
prot<-data.frame(prot)
nm<-row.names(prot)
prot<-apply(prot, 2, as.numeric)
row.names(prot)<-nm
prot<-prot %>% as.data.frame()

rownames(val.cluster)<-gsub("-", ".", val.cluster$PATIENT_ID)
prot$name<-rownames(prot)
val.cluster$name<-rownames(val.cluster)

prot<-left_join(val.cluster, prot, by="name")
prot<-prot[order(prot$cluster, decreasing = F),-c(1:3)]
prot<-prot[,nm1]

#prot<-prot[,order(colSums(prot),apply(prot,2,var), decreasing = T)] %>% as.data.frame()

par(oma=c(6,0,0,2))
labCol = lapply(strsplit(colnames(prot),
                         "_\\("), function(x){x[1]}) %>% unlist()
labCol<-as.expression(lapply(labCol, function(a) bquote(bold(.(a)))))
heatmap.2( prot %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = F, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(1,15),
           lwid = c(1,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col= colorRampPalette(c("gold","yellow","black", "magenta", "pink") %>% rev())(n = 100),
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")


par(oma=c(6,0,0,2))
heatmap.2( cbind(val.cluster$cluster[order(val.cluster$cluster, decreasing = F)],
                 val.cluster$cluster[order(val.cluster$cluster, decreasing = F)])
           %>% as.matrix(), 
           trace="none", 
           dendrogram = "none",
           margins = c(9,3), 
           density.info=c("none"),
           key = F, Rowv = F, Colv = F,
           cexRow = 0.6, cexCol = 1.9, 
           key.xlab = "",
           key.ylab = "",
           lhei = c(1,15),
           lwid = c(1,15),
           #col = colorspace::diverge_hsv(5, power = 0.3),
           #col=colorpanel(20, "gold","black", "magenta") ,
           col= c("blue","red", "green"),
           font.lab=6,
           font.axis=2,
           xlab = "",
           ylab = "",
           labRow = "",
           labCol = labCol,
           main = "metabolite",
           na.color = "grey85")






cor.t














#-----------------------------------Supplymentry--------------------------#

my.vif<-c()
my.sum<-c()
models <- c("MOFA_allgene_10F.hdf5", "MOFA_allgene_15F.hdf5", 
            "MOFA_allgene_20F.hdf5", "MOFA_allgene_30F.hdf5")

for(i in 1:4){
  model<-load_model(models[i])
  factors <- get_factors(model, as.data.frame = T) # factor matrix
  factor_mat<-matrix(factors$value, ncol=length(unique(factors$factor)), byrow = F)
  my.vif[i]<-mean(abs(vif(factor_mat)))
  
  var.fac<-model@cache$variance_explained$r2_per_factor[[1]]
  my.sum<-rbind(my.sum, colSums(var.fac))
}


plot(1:4, my.vif, ylim=c(1.05,1.3), font.lab=2, cex.lab=1.5,
     pch=19, cex=6, lwd=2, lty=1, xaxt="n", font.axis=2,
     col="red", xlab="Factor Number", ylab="Average VIF")
lines(1:4, my.vif,lwd=10, lty=1)
grid()
axis(1, at=c(1:4), labels = c(10,15,20,30), font.axis=2, cex.axis=2)

plot(rep(1:4,3), c(my.sum), ylim=c(13,57), font.lab=2, cex.lab=1.5,
     pch=19, cex=6, lwd=2, lty=1, xaxt="n", font.axis=2,
     col="grey50", xlab="Factor Number", ylab="Total Variance explained")
lines(1:4, my.sum[,1], lwd=10, lty=1, col="darkorange3")
lines(1:4, my.sum[,2], lwd=10, lty=1, col="goldenrod")
lines(1:4, my.sum[,3], lwd=10, lty=1, col="plum2")
grid()
axis(1, at=c(1:4), labels = c(10,15,20,30), font.axis=2, cex.axis=2)

legend(legend = c("Metabolomic", "Protein", "Transcriptome"),
       col = c("darkorange3","goldenrod","plum2"),
       "topleft",
       lwd=5)


model<-load_model(models[3])
wgt<-get_weights(model)



write.table( apply( wgt[[1]], 2, minmax) %>% as.data.frame(), "sup2_tab1.csv",row.names = T, sep=",")
write.table( apply( wgt[[2]], 2, minmax) %>% as.data.frame(), "sup2_tab2.csv",row.names = T, sep=",")
write.table( apply( wgt[[3]], 2, minmax) %>% as.data.frame(), "sup2_tab3.csv",row.names = T, sep=",")











