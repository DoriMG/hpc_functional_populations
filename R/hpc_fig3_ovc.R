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

loc_names <- c(
  `1` = "Before obj." ,
  `2` = "Around obj.",
  `3` = "After obj."
)

## Fig 3A-B example plots
dataset = readMat(file.path(folder,"ovc_example_pcm_ex1.mat"))
single_cell = dataset[[1]]
trav_data <- data.table("traversal" = single_cell[,2], "location" = single_cell[,3]*4, "value" = single_cell[,1])
trav_data$value = trav_data$value-min(trav_data$value )
trav_data$value = trav_data$value/max(trav_data$value )

pcm = ggplot(trav_data, aes(x=location, y=traversal, fill=value)) +
  geom_tile(aes(fill = value))+
  scale_fill_continuous(type = "viridis") +
  labs(y ='Traversal', x='Location (cm)',title="Place map", fill=expression("Norm. " ~ frac(Delta*F, F)))
pcm

dataset = readMat(file.path(folder,"ovc_example_ovcm_ex1.mat"))
single_cell = dataset[[1]]
trav_data <- data.table("traversal" = single_cell[,2], "location" = (single_cell[,3]*4)-100, "value" = single_cell[,1])
trav_data$value = trav_data$value-min(trav_data$value )
trav_data$value = trav_data$value/max(trav_data$value )

ovcm = ggplot(trav_data, aes(x=location, y=traversal, fill=value)) +
  geom_tile(aes(fill = value))+
  scale_fill_continuous(type = "viridis") +
  labs(y ='Traversal', x='Location from object (cm)',title="Object-vector map", fill=expression("Norm. " ~ frac(Delta*F, F)))
ovcm

ex_plot = (pcm|ovcm)+
  plot_layout(guides = "collect")


## No pcs
dataset = readMat(file.path(folder,"mean_ovc_perc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

pc_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Object-vector cells (%)', x=NULL, fill=NULL)
pc_plot

# Stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Figure 3D - Mutual information
dataset = readMat(file.path(folder,"allAveMI_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

mi_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Mutual information', x=NULL, fill=NULL)
mi_plot

# Stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

# Figure 3E - Correlation with 2 objects
dataset =readMat(file.path(folder,"corr_ovc_with_obj.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = factor(perf_data[,2]), "sal" = factor(perf_data[,3]), "shuffle" = perf_data[,4], "training" = perf_data[,5], "mouse" = factor(perf_data[,6]))
perf_data = perf_data[perf_data$training==2]

control_overlap = ggplot(data=perf_data, aes(x=sal, y=value, fill=factor(shuffle))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(1,2),
                   labels=c("Sparse","Abundant")) + 
  scale_fill_manual(labels = c('Data', 'Shuffled'), values=c("#6962AD", "#96E9C6")) +
  labs(y ='', x=NULL, fill=NULL)
control_overlap

lmm = lmer(value ~ sal*shuffle +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

em_res = emmeans(lmm,  ~sal*shuffle,adjust = "tukey")
write.csv(pairs(em_res, by='sal'), "test.csv")



# Figure 3F - Location of OVC centres
dataset = readMat(file.path(folder,"pf_centres_binned_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "loc" = factor(perf_data[,4]),"training" =factor(perf_data[,5]), "mouse" = perf_data[,6])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]


binned_centres = ggplot(data=perf_data, aes(x=loc, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 

  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  #scale_x_discrete(breaks=c(1,2,3,4),
  #labels=c("Control obj.","Cue obj.", "Reward zone", "Between objects")) + 
  scale_x_discrete(breaks=c(1,2,3), labels = loc_names)+
  labs(y ='Object-vector cell centres (%)', x=NULL, fill=NULL, colour=NULL)
binned_centres

lmm = lmer(value ~ sal*loc +(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


em_res = emmeans(lmm,  ~loc,adjust = "tukey")
write.csv(pairs(em_res), "test.csv")


# Figure 3G -Mean error OVC

dataset = readMat(file.path(folder,"all_obj_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

obj_error_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Mean error (cm)', x=NULL, fill=NULL)
obj_error_plot

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data[perf_data$nov])
anova(lmm)
write.csv(anova(lmm), "test.csv")


##  Fig 2H -Bayesian decoder contribution of OVCs
dataset = readMat(file.path(folder,"ovc_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

perf_data = perf_data[perf_data$nov]

ovc_diff =  ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Error difference (cm)', x=NULL, fill=NULL)
ovc_diff

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Combine the plots
plot_temp = (ex_plot) /
  (pc_plot|mi_plot|control_overlap)/
  (binned_centres)/
  (obj_error_plot|ovc_diff)

plot_all= plot_temp  +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()
plot_all

ggsave(file.path(out_folder,'fig3_OVCs.png'), plot_all, height = 11, width = 8)
ggsave(file.path(out_folder,'fig3_OVCs.pdf'), plot_all, height = 11, width = 8)

