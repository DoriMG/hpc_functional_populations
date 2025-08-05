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

## S Fig 4A - Place field width
dataset = readMat(file.path(folder,"allAvePfs_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

pf_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Object-vector field width', x=NULL, fill=NULL)
pf_plot

# stats
lmm = lmer(value ~ sal +(1|mouse) +(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## S Fig 4B- In/out ratio
dataset = readMat(file.path(folder,"allIO_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

io_plot =ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Out/in ratio', x=NULL, fill=NULL)
io_plot

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## S Fig 4C- Stability
dataset = readMat(file.path(folder,"allAveStab_ovc.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], 
                         "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data =  perf_data[perf_data$nov]
perf_data = perf_data[perf_data$training==2]

stab_plot =ggplot(data=perf_data, aes(x=sal, y=value, fill=sal)) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitter(width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +  
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant"))+
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Stability', x=NULL, fill=NULL)
stab_plot

# stats
lmm = lmer(value ~ sal +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Bayesian vs random
random = readMat(file.path(folder,"all_obj_error_rand.mat"))

perf_random = random[[1]]
perf_random <- data.table( "value" = perf_random[,1], "dataset" = perf_random[,2], "env" = perf_random[,3])

perf_random$nov = perf_random$env >2
perf_random$sal = perf_random$env ==2 | perf_random$env ==4
perf_random$random = 1

dataset = readMat(file.path(folder,"all_obj_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data$random = 0

all_data = rbind(perf_data, perf_random)
all_data = all_data[all_data$nov == 1]
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


# stats
lmm = lmer(value ~ sal + random +(1|dataset), data =all_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

stat.test <- all_data %>%
  group_by( sal) %>%
  wilcox_test(value~random, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
write.csv(tt, "test.csv")

## Combine the plots
plot_temp = (pf_plot|
               io_plot)/(
                 stab_plot|
                   rand_plot_ovc)


plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()

plot_all

ggsave(file.path(out_folder,'sfig4_ovcs.png'), plot_all, height = 8, width = 12)
ggsave(file.path(out_folder,'sfig4_ovcs.pdf'), plot_all, height = 8, width = 12)

