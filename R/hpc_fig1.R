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

sal_names <- c(
  `FALSE` = "Sparse" ,
  `TRUE` = "Salient"
)

loc_names <- c(
  `1` = "Control obj." ,
  `2` = "Cue obj.",
  `3` = "Reward zone",
  `4` = "Between objects"
)

## Figure 1E: Steps completed over days
dataset = readMat(file.path(folder,"steps_over_days.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "mouse" = factor(perf_data[,2]), "time" = perf_data[,3])

day_plot = ggplot(perf_data, aes(x=time, y=value, colour=mouse)) +
  geom_line() +
  labs(y ='Step completed', x='Days')+theme(legend.position= "none") +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
day_plot


## Figure 1F: % correct choice per environment

dataset = readMat(file.path(folder,"perf_per_env.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, # Convert to %
                         "dataset" = factor(perf_data[,2]), 
                         "env" = perf_data[,3], 
                         "training" = factor(perf_data[,4]), 
                         "mouse"=factor(perf_data[,5]))
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$nov==TRUE]

perf_env = ggplot(data=perf_data, aes(x=sal, y=value, fill=factor(training))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(FALSE,TRUE),
                     labels=c("Sparse","Salient")) + 
  scale_fill_brewer(labels = c('Untrained', 'Trained'),palette = "Set2")+
  labs(y ='Correct choice (%)', x=NULL, fill=NULL)
perf_env

# Stats test for 1F
lmm = lmer(value ~ sal*training + (1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Figure 1G: % correct choice per displacement 
dataset = readMat(file.path(folder,"all_perc_p_disp.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, 
                         "dataset" = factor(perf_data[,2]), 
                         "disp" = factor(perf_data[,3]), 
                         "training" = factor(perf_data[,4]), 
                         "mouse"=factor(perf_data[,5]))

disp_all = ggplot(data=perf_data, aes(x=factor(disp), y=value, fill=factor(training))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2),
                   labels=c("< 60 cm",">60 cm")) + 
  scale_fill_brewer(labels = c('Untrained', 'Trained'),palette = "Set2")+
  labs(y ='Correct choice (%)', x=NULL, fill=NULL)
disp_all

# Stats test for 1G
lmm = lmer(value ~ training*disp + (1|mouse) +(1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Figure 1H: Velocity by environment
dataset =  readMat(file.path(folder,"speed_p_env.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "training" = perf_data[,4], "mouse" = perf_data[,5])
perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4

speed_plot = ggplot(data=perf_data, aes(x=sal, y=value, fill=factor(training))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  facet_wrap(~ nov, labeller = labeller(nov = nov_names)) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Sparse","Abundant")) + 
  labs(y ='Velocity (cm/s)', x=NULL, fill=NULL)+
  scale_fill_brewer(labels = c('Untrained', 'Trained'),palette = "Set2")
speed_plot

# Stats Fig 1H
lmm = lmer(value ~ sal*nov*training +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Figure 1H: Velocity in trained animals depending on position
dataset = readMat(file.path(folder,"vel_p_bin.mat"))
perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "loc" = factor(perf_data[,4]),"training" =perf_data[,5], "mouse" = perf_data[,6])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==2]

vel_bin = ggplot(data=perf_data, aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  facet_wrap(~ loc, labeller = labeller(loc = loc_names), ncol = 4) +
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  #scale_x_discrete(breaks=c(1,2,3,4),
  #labels=c("Control obj.","Cue obj.", "Reward zone", "Between objects")) + 
  scale_x_discrete(breaks=c(FALSE,TRUE), labels = nov_names)+
  labs(y ='Velocity (cm/s)', x=NULL, fill=NULL, colour=NULL)
vel_bin

# Stats tests
lmm = lmer(value ~ nov*sal*loc +(1|mouse), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")
aggregate(value ~ nov*loc  + (1|mouse), perf_data, function(x) c(mean = mean(x), sd = sd(x),length = length(x)))

em_res = emmeans(lmm,  ~loc*nov, adjust = "tukey")
contrast(em_res, adjust = "tukey")
tt = as.data.frame(pairs(em_res,adjust = "tukey"))
tt
write.csv(tt, "test.csv")

# Create final figure
plot_temp = (plot_spacer() |  plot_spacer() ) /
  (day_plot|(perf_env|disp_all))/(speed_plot|vel_bin)

plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()

plot_all
ggsave(file.path(out_folder,'fig1_behaviors_v3.png'), plot_all, height = 8, width = 8)
ggsave(file.path(out_folder,'fig1_behaviors_v3.pdf'), plot_all, height = 8, width = 8)

