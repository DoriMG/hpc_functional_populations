library(ggplot2)
library("data.table") 
library(R.matlab)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)


folder = r"{D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\matlab\data}"
out_folder = r"{D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\R\figures}"  


nov_names <- c(
  `FALSE` = "Familiar" ,
  `TRUE` = "Novel"
)


## S Fig 2A - Place field width
dataset = readMat(file.path(folder,"allAvePfs.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

pf_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Place field width (cm)', x=NULL, fill=NULL)
pf_plot

# Stats
summary(perf_data)
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =perf_data)
anova(lmm)

write.csv(anova(lmm), "test.csv")


## S Fig 2B- Place field stability
dataset = readMat(file.path(folder,"allAveStab.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

stab_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Stability', x=NULL, fill=NULL)
stab_plot

# Stats
lmm = lmer(value ~ sal*nov +(1|mouse) +(1|dataset), data =perf_data)
anova(lmm)
aggregate(value ~ training  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
write.csv(anova(lmm), "test.csv")

## S Fig 2C- In/out activity ratio
dataset = readMat(file.path(folder,"allIO.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

io_plot = ggplot(data=perf_data[perf_data$env<5], aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=nov_names) + 
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  labs(y ='Out/in ratio', x=NULL, fill=NULL)
io_plot

# Stats
data_temp = perf_data[perf_data$env<5]
lmm = lmer(value ~ sal*nov +(1|mouse)+(1|dataset), data =data_temp)
anova(lmm)
aggregate(value ~ training  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))
write.csv(anova(lmm), "test.csv")

## Bayesian vs random
random = readMat(file.path(folder,"all_error_rand.mat"))

perf_random = random[[1]]
perf_random <- data.table( "value" = perf_random[,1], "dataset" = perf_random[,2], "env" = perf_random[,3])

perf_random$nov = perf_random$env >2
perf_random$sal = perf_random$env ==2 | perf_random$env ==4
perf_random$random = 1

dataset = readMat(file.path(folder,"all_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data$random = 0

all_data = rbind(perf_data, perf_random)
rand_plot = ggplot(data=all_data, aes(x=sal, y=value, fill=factor(random))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) +
  facet_wrap(~ nov, labeller = labeller(nov = nov_names)) +
  scale_x_discrete(breaks=c(FALSE,TRUE),
                   labels=c("Sparse","Abundant"))+ 
  scale_fill_manual(labels = c('Data', 'Shuffled'), values=c("#6962AD", "#96E9C6")) +
  labs(y ='Prediction error (cm)', x=NULL, fill=NULL)
rand_plot

lmm = lmer(value ~ sal*nov + random +(1|dataset), data =all_data)
anova(lmm)


perf_data$random = factor(perf_data$random)
stat.test <- all_data %>%
  group_by(env) %>%
  wilcox_test(value~random, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
tt = as.data.frame(stat.test)
## move   stab_plot/ to obc figure
write.csv(tt, "test.csv")


plot_temp = (pf_plot|
               io_plot)/(
                 stab_plot|
                   rand_plot)


plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()
plot_all


ggsave(file.path(out_folder,'sfig2_PCs_v2.png'), plot_all, height = 6, width = 8)
ggsave(file.path(out_folder,'sfig2_PCs_v2.pdf'), plot_all, height = 6, width = 8)

