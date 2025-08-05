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

## Fig 4B Raw place cell overlap
dataset = readMat(file.path(folder,"pc_overlap_raw.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==2] # only include trained animals

raw_place_cell_no = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='% cells PCs in \nboth environments', x=NULL, fill=NULL)+guides(fill="none")
raw_place_cell_no

# Stats
lmm = lmer(value ~factor(comp) +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Fig 4B Relative number of place cells
dataset = readMat(file.path(folder,"pc_overlap_p_cond_diff.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1]*100, 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==2]

overlap_plot_diff = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Place cell overlap \n(% difference from chance)', x=NULL, fill=NULL)+guides(fill="none")
overlap_plot_diff

# stats
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

## Fig 4D -  Both familiar example
dataset = readMat(file.path(folder,"fig4_cell_10_58.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "traversal" = factor(perf_data[,2]), 
                         "loc" = perf_data[,3]*4, 
                         "env" = factor(perf_data[,4]))

same_novel = ggplot(perf_data, aes(x=loc, y=value, fill=env)) +
  stat_summary(fun=mean, geom='line', alpha=1,aes(color=env),linewidth=1.5) +
  stat_summary(fun.data = mean_cl_normal, geom="ribbon", alpha=0.1) +
  labs(y ='Relative fluorescence change', x='Location (cm)')+theme(legend.position= "none") +
  scale_color_manual(values=c("#9ADE7B", "#508D69", "#FF8F8F","#995555"),labels=c("Sparse-Familiar","Abundant-Familiar", "Sparse-Novel", "Abundant-Novel"))+
  scale_fill_manual(values=c("#9ADE7B", "#508D69", "#FF8F8F","#995555"),labels=c("Sparse-Familiar","Abundant-Familiar", "Sparse-Novel", "Abundant-Novel"))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
same_novel

## Fig 4E - Both abundant example
dataset = readMat(file.path(folder,"fig4_cell_17_30.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "traversal" = factor(perf_data[,2]), 
                         "loc" = perf_data[,3]*4, 
                         "env" = factor(perf_data[,4]))

same_sal = ggplot(perf_data, aes(x=loc, y=value, fill=env)) +
  stat_summary(fun=mean, geom='line', alpha=1,aes(color=env),linewidth=1.5) +
  stat_summary(fun.data = mean_cl_normal, geom="ribbon", alpha=0.1) +
  labs(y ='Relative fluorescence change', x='Location (cm)')+
  scale_color_manual(values=c("#9ADE7B", "#508D69", "#FF8F8F","#995555"),labels=c("Sparse-Familiar","Abundant-Familiar", "Sparse-Novel", "Abundant-Novel"))+
  scale_fill_manual(values=c("#9ADE7B", "#508D69", "#FF8F8F","#995555"),labels=c("Sparse-Familiar","Abundant-Familiar", "Sparse-Novel", "Abundant-Novel"))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
same_sal



## Fig 3F -  Correlation all place maps
dataset = readMat(file.path(folder,"all_overlap_corr.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1]

all_overlap = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Correlation in place map \n(Fisher corrected)', x=NULL, fill=NULL)+guides(fill="none")
all_overlap

#stats
lmm = lmer(value ~comp +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)

## Fig 3G -  Correlation PC place maps
dataset = readMat(file.path(folder,"pc_overlap_corr.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1]

pc_overlap = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Correlation in place map \n(Fisher corrected)', x=NULL, fill=NULL)+guides(fill="none")
pc_overlap

# stats
lmm = lmer(value ~comp +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")

## Fig 3H -  Cross decoder results
dataset = readMat(file.path(folder,"all_cross_decoder_error.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1] # only trained animals (untrained is -1)

all_decoder = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Error (cm)', x=NULL, fill=NULL)+guides(fill="none")
all_decoder

# stats
lmm = lmer(value ~comp +(1|mouse) + (1|dataset), data =perf_data)
anova(lmm)
write.csv(anova(lmm), "test.csv")


## Fig 3I -  Cross decoder results - control
dataset = readMat(file.path(folder,"all_cross_decoder_diff.mat"))

perf_data = dataset[[1]]
perf_data <- data.table( "value" = perf_data[,1], 
                         "dataset" = perf_data[,2], 
                         "comp" = factor(perf_data[,3]),
                         "expected" = perf_data[,4], 
                         "training" = factor(perf_data[,4]), 
                         "mouse" = perf_data[,5])
perf_data = perf_data[perf_data$training==1]

all_decoder_difference = ggplot(data=perf_data, aes(x=comp, y=value, fill=factor(comp))) +
  stat_summary(fun=mean, geom='bar', alpha=1, position = position_dodge(width=0.8), width=0.8,colour=NA) +
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.3, position = position_dodge(width=0.8)) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha=0.5, stroke = 0,shape=16,size=1) + 
  scale_x_discrete(breaks=c(1,2,3),
                   labels=c("Same abundance","Same novelty", "Different")) +
  scale_fill_manual(values = c("#7372ae", "#81b882", "#cecdcd"))+
  labs(y ='Error difference (cm)', x=NULL, fill=NULL)+guides(fill="none")
all_decoder_difference

# stats
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


## Combine the plots
plot_temp =(plot_spacer()|raw_place_cell_no|overlap_plot_diff)/
              (same_sal|same_novel)/
  (all_overlap|pc_overlap)/
  (all_decoder|all_decoder_difference)
  

plot_temp  


plot_all = plot_temp +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme_classic()
plot_all

ggsave(file.path(out_folder,'figure4_overlap.png'), plot_all, height = 10, width = 10)
ggsave(file.path(out_folder,'figure4_overlap.pdf'), plot_all, height = 10, width = 10)




