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
  `1` = "Control obj." ,
  `2` = "Cue obj.",
  `3` = "Reward zone",
  `4` = "Between objects"
)

## Fig 2A - all place cells in a single environment example 
dataset = readMat(file.path(folder,"fig3_pc_example_all_cells.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "val" = perf_data[,1], "cell" = perf_data[,2], "loc" = perf_data[,3]*4)
# Normalize
perf_data$value = perf_data$value-min(perf_data$value )
perf_data$value = perf_data$value/max(perf_data$value )

pcm = ggplot(perf_data, aes(x=loc, y=cell, fill=val)) +
  geom_tile(aes(fill = val))+
  scale_fill_continuous(type = "viridis") +
  labs(y ='Cell', x='Location (cm)', fill=expression("Norm. " ~ frac(Delta*F, F)))
pcm

## Fig 2B - Number of place cells

dataset = readMat(file.path(folder,"mean_pc_perc.mat"))
perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
# Filter to only include trained animals
perf_data = perf_data[perf_data$training==2]

pc_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(values = c('#C59AD4', '#FFC19C'), labels = c("Sparse","Abundant"))+
  labs(y ='Place cells (%)', x=NULL, fill=NULL)
pc_plot

# Stats
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)



## Fig 2C - Place field location
dataset = readMat(file.path(folder,"pf_centres_binned_pc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "loc" = factor(perf_data[,4]),"training" =factor(perf_data[,5]), "mouse" = perf_data[,6])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

binned_centres = ggplot(data=perf_data, aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  facet_wrap(~ loc, labeller = labeller(loc = loc_names), ncol = 4) +
  scale_fill_manual(values = c('#C59AD4', '#FFC19C'), labels = c("Sparse","Abundant"))+
  scale_x_discrete(breaks=c(FALSE,TRUE), labels = nov_names)+
  labs(y ='Place cell centres (%)', x=NULL, fill=NULL, colour=NULL)
binned_centres

# stats
lmm = lmer(value ~ nov*sal*loc +(1|dataset) + (1|mouse), data =perf_data)
anova(lmm)
aggregate(value ~ nov*loc  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
write.csv(anova(lmm), "test.csv")


## Fig 2D - single cell across traversals
dataset = readMat(file.path(folder,"fig3_pc_example_one_cell.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "val" = perf_data[,1], "cell" = perf_data[,2], "loc" = perf_data[,3]*4)

perf_data$value = perf_data$value-min(perf_data$value )
perf_data$value = perf_data$value/max(perf_data$value )

pc_one_cell = ggplot(perf_data, aes(x=loc, y=cell, fill=val)) +
  geom_tile(aes(fill = val))+
  scale_fill_continuous(type = "viridis") +
  labs(y ='Traversal', x='Location (cm)', fill=expression("Norm. " ~ frac(Delta*F, F)))
pc_one_cell

##  Fig 2E - Mutual information in place cells
dataset = readMat(file.path(folder,"allAveMI.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

mi_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(values = c('#C59AD4', '#FFC19C'), labels = c("Sparse","Abundant"))+
  labs(y ='Mutual information', x=NULL, fill=NULL)
mi_plot


lmm = lmer(value ~ sal*nov +(1|mouse) +(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


##  Fig 2F -Bayesian decoder error
dataset = readMat(file.path(folder,"all_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

error_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(values = c('#C59AD4', '#FFC19C'), labels = c("Sparse","Abundant"))+
  labs(y ='Mean error (cm)', x=NULL, fill=NULL)
error_plot

lmm = lmer(value ~ sal*nov +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

aggregate(value ~ sal*nov  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
em_res = emmeans(lmm,  ~nov,adjust = "tukey")
contrast(em_res)
pairs(em_res, by='sal')
write.csv(pairs(em_res, by='sal'), "test.csv")


##  Fig 2F -Bayesian decoder contribution of PCs
dataset = readMat(file.path(folder,"pc_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

pc_diff = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(values = c('#C59AD4', '#FFC19C'), labels = c("Sparse","Abundant"))+
  labs(y ='Error (PCs vs random cells) (cm)', x=NULL, fill=NULL)
pc_diff

# stats
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
em_res = emmeans(lmm,  ~nov,adjust = "tukey")
contrast(em_res)
pairs(em_res)
write.csv(anova(lmm), "test.csv")

stat.test <- perf_data %>%
  group_by(env) %>%
  wilcox_test(value~1, mu=0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
write.csv(tt, "test.csv")

## Combine the plots
plot_temp = (pcm|pc_plot)/
  binned_centres/ 
(pc_one_cell|mi_plot)/
  (error_plot|pc_diff)
  
plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()

ggsave(file.path(out_folder,'fig2_PCs_v2.png'), plot_all, height = 10, width = 8)
ggsave(file.path(out_folder,'fig2_PCs_v2.pdf'), plot_all, height = 10, width = 8)

