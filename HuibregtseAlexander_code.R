# Choking study - biomarker data analysis ----

library(datawizard)
library(tidyverse)
library(openxlsx)
library(Hmisc)
library(ggpubr)
library(car)
library(effectsize)

setwd("/Users/meganhuibregtse/Library/CloudStorage/OneDrive-IndianaUniversity/Sexual Behavior + Brain Health [PILOT PROJECT]/Choking & neurobio pilot project/Manuscripts/Biomarkers/2 J Neurotrauma (5 biomarkers)/supplemental")
# Gather the data ----

# 4plex data
simoa4 <- read.csv("Biomarker_data2.csv", header = TRUE)
# S100B data
s100b <- read.xlsx("Biomarker_data.xlsx",sheet = 2, colNames = TRUE)
# participant characteristics
ptdata <- read.csv("pt_characteristics_pilot.csv", header = TRUE)

## Merge data together ----
d1 <- Merge(s100b, simoa4, ptdata, id = ~SubID, all = FALSE, verbose = TRUE)

## Data cleaning ----
head(d1)
library(dplyr)
library(stringr)
d1 <- d1 %>%
  mutate_at(c('SubID','sex','hispanic','marital','student','covid','tbi_ever','tbi_recent','modsevtbi','sex_orient','sex_active','ever_choked'),as.factor)

d1 <- d1 %>%
  mutate(race2 = case_when(
    hispanic == 'yes' ~ 'hispanic',
    str_detect(race,"asian") ~ 'non-hispanic asian',
    str_detect(race,"alaska") ~ 'non-hispanic AI/AN',
    str_detect(race,"white") ~ 'non-hispanic white',
    str_detect(race,"black") ~'non-hispanic black',), .after = race)%>%
  mutate_at(c('race2'),as.factor) %>% 
  mutate(grp = case_when(chokedfreq_1 >= 4 ~ "Hx",
                         is.na(chokedfreq_1) ~ "NoHx"), .after = SubID) %>%
  mutate_at('grp', as.factor)

# Summarize for table 1, calculate group differences ----
library(psych)
library(rstatix)
library(janitor)
d1 %>%
  group_by(grp) %>%
  dplyr::summarize(n = n(), mean = mean(age), sd = sd(age), median = median(age), iqr = IQR(age), age_quants = quantile(age, c(.25,.75)), q = c(.25, .75)) 

shapiro.test(d1$age)
#not normal
d1 %>% wilcox_test(age ~ grp)
#age differs by group

d1 <- d1 %>%
  mutate(student_bin = case_when(
    student == "undergrad" ~ 0,
    TRUE ~ 1
  ))
library(ltm)
biserial.cor(d1$age, d1$student, use = c("all.obs"), level = 2)
d1 %>% cor_test(vars=c(age,student_bin))

d1 %>% 
  group_by(grp,race2) %>%
  dplyr::summarize(n = n())
tab <- janitor::tabyl(d1,grp,race2)
janitor::fisher.test(tab)

d1 %>% 
  group_by(grp,student) %>%
  dplyr::summarize(n = n())
tab <- tabyl(d1,grp,student)
fisher.test(tab)

d1 <- d1 %>%
  mutate(mtbi_bin = as.factor(mtbi_count))
d1 %>%
  group_by(grp,mtbi_bin) %>%
  dplyr::summarize(n=n())
tab <- tabyl(d1,mtbi_bin,grp)
fisher.test(tab, simulate.p.value=TRUE)

d1 %>%
  group_by(grp) %>%
  dplyr::summarize(mean = mean(phq9_score), sd = sd(phq9_score), median = median(phq9_score), iqr = IQR(phq9_score), phq_quants = quantile(phq9_score, c(.25,.75)), q = c(.25, .75))
shapiro.test(d1$phq9_score)
d1 %>%
  wilcox_test(phq9_score ~ grp)


d1 %>%
  group_by(grp) %>%
  dplyr::summarize(mean = mean(gad7_score), sd = sd(gad7_score), median = median(gad7_score), iqr = IQR(gad7_score), gad_quants = quantile(gad7_score, c(.25,.75)), q = c(.25, .75))
shapiro.test(d1$gad7_score)
d1 %>%
  group_by(grp) %>%
  shapiro_test(gad7_score)
d1 %>%
  wilcox_test(gad7_score ~ grp)

d1 %>%
  group_by(grp) %>%
  dplyr::summarize(mean = mean(audit_score), sd = sd(audit_score), median = median(audit_score), iqr = IQR(audit_score))
shapiro.test(d1$audit_score)
d1 %>%
  group_by(grp) %>%
  shapiro_test(audit_score)
d1 %>%
  t_test(audit_score~grp, var.equal = T)

# Choking history description (also for Table 1) ----
d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(chokedfreq_1), sd = sd(chokedfreq_1), median = median(chokedfreq_1), iqr = IQR(chokedfreq_1), cf1_quants = quantile(chokedfreq_1, c(.25,.75)), q = c(.25, .75))
d1 %>% shapiro_test(chokedfreq_1)

d1 %>% shapiro_test(chokedfreq_2)
d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(chokedfreq_2), sd = sd(chokedfreq_2), median = median(chokedfreq_2), iqr = IQR(chokedfreq_2), cf2_quants = quantile(chokedfreq_2, c(.25,.75)), q = c(.25, .75))

d1 %>% shapiro_test(chokedfreq_12)
d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(chokedfreq_12), sd = sd(chokedfreq_12), median = median(chokedfreq_12), iqr = IQR(chokedfreq_12), cf12_quants = quantile(chokedfreq_12, c(.25,.75)), q = c(.25, .75))


d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(intensity_1), sd = sd(intensity_1), median = median(intensity_1), iqr = IQR(intensity_1))

d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(intensity_2), sd = sd(intensity_2), median = median(intensity_2), iqr = IQR(intensity_2))

d1 %>%
  filter(grp == "Hx") %>%
  dplyr::summarize(mean = mean(physicalresponse_freqscore), sd = sd(physicalresponse_freqscore), median = median(physicalresponse_freqscore), iqr = IQR(physicalresponse_freqscore))

# Model Fitting ----
model_m <- lm(cbind(GFAP, NfL, Tau, UCHL1, S100B) ~ grp, data = d1)
summary(model_m)
Manova(model_m)

d1 %>%
  dplyr::select(c(grp, SubID, NfL)) %>%
  group_by(grp) %>%
  identify_outliers(NfL)
# two outliers

d1 %>%
  dplyr::select(c(grp, SubID, GFAP)) %>%
  group_by(grp) %>%
  identify_outliers(GFAP)
#one outlier

d1 %>%
  dplyr::select(c(grp, SubID, Tau)) %>%
  group_by(grp) %>%
  identify_outliers(Tau)
#three outliers, one is extreme 

d1 %>%
  dplyr::select(c(grp, SubID, UCHL1)) %>%
  group_by(grp) %>%
  identify_outliers(UCHL1)
#one outlier

d1 %>%
  dplyr::select(c(grp, SubID, S100B)) %>%
  group_by(grp) %>%
  identify_outliers(S100B)
# one extreme outlier 

# multivariate outliers
d1 %>%
  dplyr::select(c(grp, SubID, S100B, NfL, Tau, GFAP, UCHL1)) %>%
  #group_by(grp) %>%
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE) %>%
  as.data.frame()
# one outlier

#univariate normality?
d1 %>%
  group_by(grp) %>%
  shapiro_test(S100B, NfL, Tau, GFAP, UCHL1) %>%
  arrange(variable)
#tau is not normally distributed in either group... S100B is not normally distributed in the Hx group... NfL is not normally distributed in the NoHx group

d1 <- d1 %>%
  mutate(
    NfL_log = log(NfL),
    S100B_log = log(S100B),
    Tau_log = log(Tau)
  )

d1 %>%
  group_by(grp) %>%
  shapiro_test(S100B_log, NfL_log, Tau_log, GFAP, UCHL1) %>%
  arrange(variable)
#assumption met now

#model with transformed tau, NfL, and S100B
model_m <- lm(cbind(GFAP, NfL_log, Tau_log, UCHL1, S100B_log) ~ grp, data = d1)
summary(model_m)
Manova(model_m)

d1 %>%
  dplyr::select(c(grp, SubID, NfL_log)) %>%
  group_by(grp) %>%
  identify_outliers(NfL_log)

d1 %>%
  dplyr::select(c(grp, SubID, GFAP)) %>%
  group_by(grp) %>%
  identify_outliers(GFAP)

d1 %>%
  dplyr::select(c(grp, SubID, Tau_log)) %>%
  group_by(grp) %>%
  identify_outliers(Tau_log)

d1 %>%
  dplyr::select(c(grp, SubID, UCHL1)) %>%
  group_by(grp) %>%
  identify_outliers(UCHL1)

d1 %>%
  dplyr::select(c(grp, SubID, S100B_log)) %>%
  group_by(grp) %>%
  identify_outliers(S100B_log)

#no extreme univariate outliers

# multivariate outliers
d1 %>%
  dplyr::select(c(grp,S100B_log, NfL_log, Tau_log, GFAP, UCHL1)) %>%
  #group_by(grp) %>%
  doo(~mahalanobis_distance(.)) %>%
  filter(is.outlier == TRUE)
#no multivariate outliers

ggqqplot(d1, "S100B_log", facet.by = "grp", ylab = "Log S100B", ggtheme = theme_minimal())
ggqqplot(d1, "NfL_log", facet.by = "grp", ylab = "Log NfL", ggtheme = theme_minimal())
ggqqplot(d1, "Tau_log", facet.by = "grp", ylab = "Log Tau", ggtheme = theme_minimal())
ggqqplot(d1, "UCHL1", facet.by = "grp", ylab = "UCHL1", ggtheme = theme_minimal())
ggqqplot(d1, "GFAP", facet.by = "grp", ylab = "GFAP", ggtheme = theme_minimal())

#multivariate normality
library(mvnormalTest)
mardia(d1[, c("S100B_log", "NfL_log", "GFAP", "UCHL1", "Tau_log")])$mv.test
#assumption met

#homogeneity of covariances
box_m(d1[, c("S100B_log", "NfL_log", "GFAP", "UCHL1", "Tau_log")], d1$grp)
#assumption met

#homogeneity of variance
d1 %>% 
  gather(key = "variable", value = "value", GFAP, NfL_log, Tau_log, UCHL1, S100B_log) %>%
  group_by(variable) %>%
  levene_test(value ~ grp)
#assumption met


### Adjust for age and AUDIT----
model_m <- lm(cbind(GFAP, NfL_log, Tau_log, UCHL1, S100B_log) ~ grp, data = d1)
summary(model_m)
Manova(model_m)
model_m_adj <- lm(cbind(GFAP, NfL_log, Tau_log, UCHL1, S100B_log) ~ age + audit_score + grp, data = d1)
summary(model_m_adj)
Manova(model_m_adj)
summary.aov(model_m_adj, type = 2)
eta_squared(model_m_adj)

model_anova <- Anova(model_m_adj, type = 2)
eta_squared(model_anova)


# Figures ----
## Figure 2 ----
library(ggsignif)
library(ggsci)

p_S100B_log <- ggplot(d1,aes(x=grp, y = S100B_log)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("Log S100B concentration (pg/mL)") +
  geom_signif(comparisons=list(c("Hx", "NoHx")), annotations="**",
              y_position = 5.2, tip_length = 0, vjust=0.4, textsize = 5)+
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_S100B <- ggplot(d1,aes(x=grp, y = S100B)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("S100B concentration (pg/mL)") +
  geom_signif(comparisons=list(c("Hx", "NoHx")), annotations="**",
              y_position = 5.2, tip_length = 0, vjust=0.4, textsize = 5)+
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_NfL_log <- ggplot(d1,aes(x=grp, y = NfL_log)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("Log NfL concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_NfL <- ggplot(d1,aes(x=grp, y = NfL)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("NfL concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_Tau_log <- ggplot(d1,aes(x=grp, y = Tau_log)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("Log Tau concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_Tau <- ggplot(d1,aes(x=grp, y = Tau)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("Tau concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_GFAP <- ggplot(d1,aes(x=grp, y = GFAP)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("GFAP concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

p_UCHL1 <- ggplot(d1,aes(x=grp, y = UCHL1)) + 
  geom_boxplot(aes(color = grp, fill = grp), show.legend = F, outlier.shape = NA, alpha = 0.8, width = .5) +
  #geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = .25, show.legend = F)+
  geom_jitter(width = .08, aes(color=grp), size = 1, show.legend = F)+
  scale_color_manual(values = c("#374e55","#F0B356"))+
  scale_fill_manual(values = c("#B5C9CF","#F8DDB5"))+
  xlab("Group") +
  ylab("UCHL1 concentration (pg/mL)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=11), axis.title.y = element_text(size=11))

library(patchwork)
(p_S100B_log | p_NfL_log | p_GFAP) / 
  (p_Tau_log | p_UCHL1 | plot_spacer()) + plot_annotation(tag_levels = 'A')
ggsave("Fig2_v1.tiff", width = 7, height = 7, bg ="white", dpi = 400)

(p_S100B | p_NfL | p_GFAP) / 
  (p_Tau | p_UCHL1 | plot_spacer()) + plot_annotation(tag_levels = 'A')

## Figure 3 ----
### ROC analysis ----
library(pROC)
S100Broc <- roc(d1, grp, S100B_log, levels = c("NoHx", "Hx"), ci=TRUE, ci.alpha=0.95, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, print.auc=TRUE, axes=TRUE)
sens.ci<- ci.se(S100Broc)
plot(sens.ci,type="shape",col='lightblue')

ggroc(S100Broc, legacy.axes = T) + 
  geom_abline(slope =1, intercept = 0, color = "red",linetype = "dashed", size = .25)+ 
  coord_equal() +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  annotate(geom = "text", x = .65, y = .2, label="AUC: 0.824\n(95%CI 0.670-0.977)") +
  theme_minimal()

ggsave("Fig3.tiff", width = 4, height = 4, bg ="white", dpi = 400)