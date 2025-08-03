library(ggplot2)
library("data.table") 
library(R.matlab)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)


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

stat.test <- perf_data %>%
  wilcox_test(value~sal, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test



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

stat.test <- perf_data %>%
  wilcox_test(value~sal, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


# Figure 3E - Correlation with 2 objects
dataset =readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\corr_ovc_with_obj.mat"))

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
aggregate(value ~ sal*shuffle , perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
write.csv(anova(lmm), "test.csv")

em_res = emmeans(lmm,  ~sal*shuffle,adjust = "tukey")
contrast(em_res)
pairs(em_res, by='sal')
write.csv(pairs(em_res, by='sal'), "test.csv")

stat.test <- perf_data %>%
  group_by(sal) %>%
  wilcox_test(value~shuffle, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


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
aggregate(value ~ sal*loc , perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
write.csv(anova(lmm), "test.csv")

### Decoder stuff

## Errors relative to obj
dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\all_error_obj.mat"))

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

lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data[perf_data$nov])
anova(lmm)
write.csv(anova(lmm), "test.csv")
aggregate(value ~ sal  + (1|mouse), perf_data[perf_data$nov], function(x) c(mean = mean(x), sd = sd(x),length = length(x)))

stat.test <- perf_data %>%
  wilcox_test(value~sal) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test


## Error diff for ovc

dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\error_diff_ovc.mat"))

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

lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")
aggregate(value ~ sal  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))

stat.test <- perf_data[perf_data$nov] %>%
  wilcox_test(value~sal) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

stat.test <- perf_data[perf_data$nov] %>%
  group_by(sal) %>%
  wilcox_test(value~1, mu=0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
## move   stab_plot/ to obc figure
write.csv(tt, "test.csv")

plot_temp = (tot_plot) /
  (pc_plot|mi_plot|control_overlap)/
  (binned_centres)/
  (obj_error_plot|ovc_diff)
  

plot_temp

plot_all= plot_temp  +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


ggsave(file.path(folder,'v3\\figures\\fig3_OVCs_v2.png'), plot_all, height = 11, width = 8)
ggsave(file.path(folder,'v3\\figures\\fig3_OVCs_v2.pdf'), plot_all, height = 11, width = 8)


# Supplementary ######################################################################################


## overalp

dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\ovc_pc_overlap.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "dataset" = perf_data[,2], "env" = perf_data[,3], "comp" = factor(perf_data[,4]),"training" =perf_data[,5],  "mouse" = perf_data[,6])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

pc_v_ovc_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=comp)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = "stack", width=0.8) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_discrete(labels = c('OVC only','OVC + PC'))+
  labs(y ='Cells (%)', x=NULL, fill=NULL)
pc_v_ovc_plot


lmm = lmer(value ~ sal*comp +(1|mouse)+(1|dataset), data =perf_data[perf_data$nov])
anova(lmm)

perf_data = perf_data[perf_data$comp==2]
pc_v_ovc_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='% of object-vector cells', x=NULL, fill=NULL)
pc_v_ovc_plot
ag <- aggregate(value ~ sal*training , perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
ag
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data[perf_data$nov])
anova(lmm)

## Characteristics
## PF width
dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\allAvePfs_ovc_v2.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]
perf_data = perf_data[perf_data$env<5]

pf_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Object-vector field width', x=NULL, fill=NULL)
pf_plot
aggregate(value ~ sal , perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))

data_temp = perf_data
lmm = lmer(value ~ sal +(1|mouse) +(1|dataset), data =data_temp)
anova(lmm)

## Stability
dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\allAveStab_ovc_v2.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], 
                         "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]
perf_data = perf_data[perf_data$env<5]

stab_plot =ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Stability', x=NULL, fill=NULL)
stab_plot

data_temp = perf_data
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =data_temp)
anova(lmm)

## IO
dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\allIO_ovc_v2.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5], "ntravs"= perf_data[,6])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]
perf_data = perf_data[perf_data$env<5]

io_plot =ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Out/in ratio', x=NULL, fill=NULL)
io_plot
aggregate(value ~ training , perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
data_temp = perf_data
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =data_temp)
anova(lmm)


## Bayesian vs random
random = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\all_error_obj_rand.mat"))

perf_random = random[[1]]
perf_random <- data.table( "value" = perf_random[,1], "dataset" = perf_random[,2], "env" = perf_random[,3])

perf_random$nov = perf_random$env >2
perf_random$sal = perf_random$env ==2 | perf_random$env ==4
perf_random$random = 1

dataset = readMat(file.path(folder,"batch_analysis\\2point0\\R\\data\\all_error_obj.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data$random = 0

all_data = rbind(perf_data, perf_random)
rand_plot_ovc = ggplot(data=all_data, aes(x=sal, y=value, fill=factor(random))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),alpha=0.5, stroke = 0,shape=16,size=1) +  
  facet_wrap(~ nov, labeller = labeller(nov = nov_names)) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant")) + 
  scale_fill_manual(labels = c('Data', 'Shuffled'), values=c("#6962AD", "#96E9C6")) +
  labs(y ='Prediction error (cm)', x=NULL, fill=NULL)
rand_plot_ovc

lmm = lmer(value ~ sal*nov + random +(1|dataset), data =all_data)
anova(lmm)

for(n in 0:1)
{
  for(s in 0:1)
  {
    print(shapiro.test(all_data[all_data$nov == n&all_data$sal == s& all_data$random == 0]$value))
    print(shapiro.test(all_data[all_data$nov == n&all_data$sal == s& all_data$random == 1]$value))
  }
}

stat.test <- all_data %>%
  group_by(nov, sal) %>%
  wilcox_test(value~random, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

plot_temp = (pc_v_ovc_plot|pf_plot|
               io_plot)/(
                 stab_plot|
                   rand_plot_ovc)


plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

plot_all


save_folder = r"{\\gpfs.corp.brain.mpg.de\bark\personal\grijseelsd\Papers\Sussex\hippocampal_pops\v3\figures}"
ggsave(file.path(save_folder,'sfig4_ovcs.png'), plot_all, height = 8, width = 12)
ggsave(file.path(save_folder,'sfig4_ovcs.pdf'), plot_all, height = 8, width = 12)



#################################Supplementary################################################################
source("E:\\Dropbox (Brain Energy Lab)\\Everything\\Dori\\Analysis\\batch_analysis\\2point0\\R\\ovc_example_plot.R")

## example plot
plot5 = ovc_example_plot(5)
plot_temp = ovc_example_plot(5)/ovc_example_plot(6)/
  ovc_example_plot(7)/ovc_example_plot(3)/
  ovc_example_plot(2)/ovc_example_plot(4)

plot_left =   ovc_example_plot(2)/ovc_example_plot(6)/
  ovc_example_plot(7)
plot_right = ovc_example_plot(3)/ovc_example_plot(5)/ovc_example_plot(4)
plot_temp = (plot_left)|(plot_right)
plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('E:\\Dropbox (Brain Energy Lab)\\Everything\\Dori\\Analysis\\batch_analysis\\2point0\\R\\fig6_ovc_supp.png', plot_all, height = 8, width = 12)

ggsave('E:\\Dropbox (Brain Energy Lab)\\Everything\\Dori\\Analysis\\batch_analysis\\2point0\\R\\fig6_ovc_supp.eps', plot_all, device="eps", height = 8, width = 12)





plot_temp = (pf_plot|io_plot)/
  (stab_plot|plot_spacer())/
  binned_centres+ plot_layout(heights = c(1,1,1))
plot_temp

plot_all= plot_temp  +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave('E:\\Dropbox (Brain Energy Lab)\\Everything\\Dori\\Analysis\\batch_analysis\\2point0\\R\\fig6_ovc_supp_2.png', plot_all, height = 8, width =8)

ggsave('E:\\Dropbox (Brain Energy Lab)\\Everything\\Dori\\Analysis\\batch_analysis\\2point0\\R\\fig6_ovc_supp_2.eps', plot_all, height = 8, width =8)


