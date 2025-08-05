library(ggplot2)
library("data.table") 
library(R.matlab)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)

library("dplyr")
library(rstatix)

folder = r"{D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\matlab\data}"
out_folder = r"{D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\R\figures}"  

# Variable names
nov_names <- c(
  `FALSE` = "Familiar" ,
  `TRUE` = "Shifted"
)

sal_names <- c(
  `FALSE` = "Sparse" ,
  `TRUE` = "Abundant"
)

##Fig 5A - number of place cells
dataset = readMat(file.path(folder,"mean_pc_perc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==1]

pc_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Place cells (%)', x=NULL, fill=NULL)
pc_plot

# stats
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

##Fig 5B - MI of place cells
dataset = readMat(file.path(folder,"allAveMI.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==1]

mi_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Mutual information', x=NULL, fill=NULL)
mi_plot

lmm = lmer(value ~ sal*nov +(1|mouse) +(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Fig 5C -Bayesian decoder
dataset = readMat(file.path(folder,"all_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==1]

error_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Mean error (cm)', x=NULL, fill=NULL)
error_plot

# stats
lmm = lmer(value ~ sal*nov +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Fig 5D - Error difference when excluding PC
dataset = readMat(file.path(folder,"pc_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==1]

pc_diff = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Error difference (cm)', x=NULL, fill=NULL)
pc_diff

# Stats
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

stat.test <- perf_data %>%
  group_by(env) %>%
  wilcox_test(value~1, mu=0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
write.csv(tt, "test.csv")


## Fig 5E -  Object-vector cells
dataset = readMat(file.path(folder,"mean_ovc_perc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==1]

ovc_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Object-vector cells (%)', x=NULL, fill=NULL)
ovc_plot

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Fig 5F - MI OVC
dataset = readMat(file.path(folder,"allAveMI_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==1]

mi_plot_ovc = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Mutual information', x=NULL, fill=NULL)
mi_plot_ovc

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Fig 5G - Overlap in PCs
dataset = readMat(file.path(folder,"pc_overlap_raw.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1]

raw_place_cell_no = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='% cells PCs in \nboth environments', x=NULL, fill=NULL)+guides(fill="none")
raw_place_cell_no

#stats
lmm = lmer(value ~factor(comp) +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Fig 5H Relative number of place cells  
dataset = readMat(file.path(folder,"pc_overlap_p_cond_diff.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1]

overlap_plot_diff = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Place cell overlap \n(% difference from chance)', x=NULL, fill=NULL)+guides(fill="none")
overlap_plot_diff

lmm = lmer(value ~factor(comp) +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Fig 5I -  Cross decoder results
dataset = readMat(file.path(folder,"all_cross_decoder_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==-1]

all_decoder = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Error (cm)', x=NULL, fill=NULL)+guides(fill="none")
all_decoder

#stats
lmm = lmer(value ~comp +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Fig 5J -  Cross decoder results - control
dataset = readMat(file.path(folder,"all_cross_decoder_diff.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==-1]

all_decoder_difference = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Error difference (cm)', x=NULL, fill=NULL)+guides(fill="none")
all_decoder_difference

#stats
lmm = lmer(value ~comp +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

stat.test <- perf_data %>%
  group_by(comp) %>%
  wilcox_test(value~1, mu=0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
write.csv(tt, "test.csv")


all_plots = 
  (pc_plot|mi_plot|error_plot)/
  (pc_diff|ovc_plot|mi_plot_ovc)/
(raw_place_cell_no|overlap_plot_diff)/
(all_decoder|all_decoder_difference)


plot_all = all_plots +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()

plot_all


ggsave(file.path(out_folder,'fig5_training_effect_v2.png'), plot_all, height = 10, width = 10)
ggsave(file.path(out_folder,'fig5_training_effect_v2.pdf'), plot_all, height = 10, width = 10)

