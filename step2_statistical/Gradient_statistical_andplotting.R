setwd('D:/youyi_fucha/first_test/gradient_voxel/final_re')
library(bruceR)
data <- import('Gradient_use_subjects.xlsx')
data <- import('gradient_chayi.xlsx')
library(PupillometryR)
library(ggsignif)

names(data)

behave_ <- data %>% filter(group=="SUB") %>% select(c(7:8,11:18))
behave_$overeatting_time_perweek <- as.numeric(behave_$overeatting_time_perweek )
behave_$hate_food_time_perweek <- as.numeric(behave_$hate_food_time_perweek )

IFG_R <- data %>% ggplot(aes(x=group,y=g1_RIFG,fill=group,color=group)) + 
  geom_flat_violin(position=position_nudge(x=.1, y=0), alpha=.8, width=.7,trim = TRUE) +
  geom_boxplot( width=.1, alpha = 0.5)+ theme_classic()+geom_signif(
    comparisons = list(c("CON", "SUB")),test="t.test",
    map_signif_level = TRUE, textsize = 6,color="black"
  )+xlab("Group")+ylab("Gradient 1")+ggsci::scale_color_npg()+ggsci::scale_fill_npg()

ggsave(filename = "IFG_R.png",IFG_R,width = 180, height = 140, dpi = 300, units = "mm", device='png')

Precunus <- data %>% ggplot(aes(x=group,y=g2_Precunus,fill=group,color=group)) + 
  geom_flat_violin(position=position_nudge(x=.1, y=0), alpha=.8, width=.7,trim = TRUE) +
  geom_boxplot( width=.1, alpha = 0.5)+ theme_classic()+geom_signif(
    comparisons = list(c("CON", "SUB")),test="t.test",
    map_signif_level = TRUE, textsize = 6,color="black"
  )+xlab("Group")+ylab("Gradient 2")+ggsci::scale_color_npg()+ggsci::scale_fill_npg()

ggsave(filename = "Precunus.png",IFG_R,width = 180, height = 140, dpi = 300, units = "mm", device='png')

LIFG <- data %>% ggplot(aes(x=group,y=g2_LIFG,fill=group,color=group)) + 
  geom_flat_violin(position=position_nudge(x=.1, y=0), alpha=.8, width=.7,trim = TRUE) +
  geom_boxplot( width=.1, alpha = 0.5)+ theme_classic(base_size = 20)+geom_signif(
    comparisons = list(c("CON", "SUB")),test="t.test",
    map_signif_level = TRUE, textsize = 6,color="black"
  )+xlab("Group")+ylab("Gradient 2")+ggsci::scale_color_npg()+ggsci::scale_fill_npg()

ggsave(filename = "LIFG.png",LIFG,width = 180, height = 140, dpi = 300, units = "mm", device='png')

data$g1_RIFG <- -data$g1_RIFG
DEBQ_33_external <- lm(DEBQ_33_external~g1_RIFG*group,
                    data)
GLM_summary(DEBQ_33_external)


interactions::sim_slopes(DEBQ_33_external,pred = g1_RIFG,modx = group)
interactions::interact_plot(DEBQ_33_external,pred = g1_RIFG,modx = group,interval = T)
Precunus<- data %>% ggplot(aes(x=group,y=g2,fill=group,color=group)) + 
  geom_flat_violin(position=position_nudge(x=.1, y=0), alpha=.8, width=.7,trim = TRUE) +
  geom_boxplot( width=.1, alpha = 0.5)+ theme_classic()+geom_signif(
    comparisons = list(c("CON", "SUB")),test="t.test",
    map_signif_level = TRUE, textsize = 6,color="black"
  )+xlab("Group")+ylab("Gradient 2")+ggsci::scale_color_npg()+ggsci::scale_fill_npg()

ggsave(filename = "Precunus.png",IFG_R,width = 180, height = 140, dpi = 300, units = "mm", device='png')

DEBQ_33_total <- lm(DEBQ_33_total~g1+g2+g1_range+g2_range+g1:group+g2:group+
             g1_range:group+g2_range:group,
           data)
summary(DEBQ_33_total)
EDI_BN <- lm(EDI_BN~g1+g2+g1_range+g2_range+g1:group+g2:group+
               g1_range:group+g2_range:group,
             data)
summary(EDI_BN)
library(sjPlot)
library(sjmisc)
library(ggplot2)
plot_model(EDI_BN, type = "pred", terms = c("g2_range", "group"))+theme_classic()
plot_model(DEBQ_33_total, type = "pred", terms = c("g2_range", "group"))+theme_classic()

interactions::sim_slopes(DEBQ_33_total,pred = g2_range,modx = group)

