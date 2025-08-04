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

ap_names <- c(
  `1` = "Passive" ,
  `2` = "Active"
)

loc_names <- c(
  `1` = "Control obj." ,
  `2` = "Cue obj.",
  `3` = "Reward zone",
  `4` = "Between objects"
)

## S Fig 1A Illustrate training steps
dataset = readMat(file.path(folder,"example_behaviour.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, "session" = factor(perf_data[,2]), "control" = factor(perf_data[,3]), "ap" = perf_data[,4])

behav_exp = ggplot(perf_data, aes(x=session, y=value, group=control, color=control)) +
  geom_line() +
  facet_wrap(~ ap, labeller = labeller(ap = ap_names)) +
  labs(y ='% rewards obtained', x='Session')+theme(legend.position= "none") +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
behav_exp


## S Fig 1B Velocity in untrained animals
dataset = readMat(file.path(folder,"vel_p_bin.mat"))
perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], "dataset" = perf_data[,2], "env" = perf_data[,3], "loc" = factor(perf_data[,4]),"training" =perf_data[,5], "mouse" = perf_data[,6])

perf_data$nov = perf_data$env >2
perf_data$sal = perf_data$env ==2 | perf_data$env ==4
perf_data = perf_data[perf_data$training==1]

vel_bin_untrained = ggplot(data=perf_data, aes(x=nov, y=value, fill=factor(sal))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  facet_wrap(~ loc, labeller = labeller(loc = loc_names), ncol = 4) +
  scale_fill_manual(labels = c('Sparse', 'Abundant'),values = c("#FFC19C", "#C59AD4"))+
  #scale_x_discrete(breaks=c(1,2,3,4),
  #labels=c("Control obj.","Cue obj.", "Reward zone", "Between objects")) + 
  scale_x_discrete(breaks=c(FALSE,TRUE), labels = nov_names)+
  labs(y ='Velocity (cm/s)', x=NULL, fill=NULL, colour=NULL)+ theme_classic()
vel_bin_untrained

# Stats for S Fig 1B
lmm = lmer(value ~ nov*sal*loc +(1|mouse), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


# Create final figure
plot_temp = (behav_exp)/
  (vel_bin_untrained)

plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()

plot_all
ggsave(file.path(folder,'sfig1.png'), plot_all, height = 8, width = 8)
ggsave(file.path(folder,'sfig1.pdf'), plot_all, height = 8, width = 8)

