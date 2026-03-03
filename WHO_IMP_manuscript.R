library(tidyverse)
library(rstatix)
library(readxl)
library(Hmisc)
library(ggpubr)
library(clusrank)
library(ggrepel)
library(gtsummary)

setwd("S:/Group/AF_LC_RT_share/WHO IMP CEIRS infant study/Write up papers conferences/comm med revision")

# READ demographic data-----------
demo <- read_excel("WHOIMP_dataset.xlsx",sheet = "demographics_prior_exposures")
## Table 1 demographics------
table1 <- demo %>% 
  filter (StudyYear == 1,  Classification %in% c("prior vaccination", "prior A(H3N2)")) %>%
  select(Classification,Sex,`AgeAtVax(months)`,`AgePriorVax_Inf(months)`,Interval,StartYear) %>% 
  tbl_summary(by=Classification,
              missing_text="Missing") %>% 
  add_p(list(all_continuous() ~ "wilcox.test", 
             all_categorical() ~ "fisher.test"))
table1

table1 %>% 
  as_hux_xlsx("table1.xlsx")

stats_age <- demo %>%
  filter (StudyYear == 1,  Classification %in% c("prior vaccination", "prior A(H3N2)")) %>%
  select (`AgeAtVax(months)`, Classification) %>%
  wilcox_test(`AgeAtVax(months)`~ Classification)
stats_age

stats_age_prior <- demo %>%
  filter (StudyYear == 1,  Classification %in% c("prior vaccination", "prior A(H3N2)")) %>%
  select (`AgePriorVax_Inf(months)`, Classification) %>%
  wilcox_test(`AgePriorVax_Inf(months)`~ Classification)
stats_age_prior

stats_interval <- demo %>%
  filter (StudyYear == 1,  Classification %in% c("prior vaccination", "prior A(H3N2)")) %>%
  select (Interval, Classification) %>%
  wilcox_test(Interval ~ Classification)
stats_interval

# READ data HI antibodies against A(H1N1) and A(H3N2) viruses----
data <- read_excel("WHOIMP_dataset.xlsx",sheet = "H1_H3_HI_titre")

data$l2hi <- log(data$Titer)/log(2)
data_renamed <- data %>%
  mutate(
    pid = PID, timepoint = Visit2Y, visit = Visit, sex = factor(Sex, levels = c("0","1"), labels =c("female","male")),
    age = AgeAtVax, clade = Clade, virus = Virus,vax_strain = vax_strain,
    abbr=Abbr, cohort = factor(Cohort), class = factor(Classification, levels =c("no/unclear prior","prior influenza A", "prior A(H3N2)", "prior vaccination","prior infection & vaccination")),
    virus_year = VirusYear,study_year = StudyYear, sample_year = SampleYear, start = StartYear, 
    virus_order= virus_order, virus_year = VirusYear,cell = Egg_Cell, stype = factor(Subtype, levels=c("H3N2","H1N1")),
    titre = Titer,l2hi = l2hi,gmr = ratio, vaccine_day = vaccine_day
  )
data_renamed$l2gmr <- log(data_renamed$gmr)/log(2)

data_extra <- data_renamed %>%
  mutate(
    virus = fct_reorder(virus, virus_order),
    abbr = fct_reorder(abbr,virus_order),
    timepoint_lbl = factor(
      timepoint, 1:8, c("1-d0", "1~d7", "1~d28","2-d0","2~d7","2~d28","3-d0","3~d28")
    )
  ) 

# add titre category variable
data_extra <- data_extra %>%
  mutate(
    exceed20 = case_when(
      titre > 20 ~ 1,
      titre < 40 ~ 0,
    ))%>%
  mutate(
    exceed40 = case_when(
      titre > 40 ~ 1,
      titre < 80 ~ 0,
    ),
    virus = fct_reorder(virus, virus_year)
  )


#exclude 9-07 who had both prior infection and prior vaccination
data_extra <- data_extra %>%
  filter (pid != "9-07")

###restructure to wide look at fold rise and effect of prevax titre------
ratios <-subset(data_extra,select = c("pid","cohort","class","study_year","visit","cell","titre","stype","abbr", "vax_strain"))
ratios$cohortC <- as.factor(ratios$cohort)
ratios <- ratios %>% 
  group_by(abbr,study_year, visit,pid) %>% 
  pivot_wider(names_from = visit, values_from = titre)
ratios <- ratios %>% 
  unnest (cols = c(`1`, `2`,`3`))
ratios$FRd7<- ratios$`2`/ratios$`1`
ratios$FRd28<- ratios$`3`/ratios$`1`
ratios$conv7 [ratios$FRd7 >=4] <- "Yes"
ratios$conv7 [ratios$FRd7 <4] <- "No"
ratios$conv28 [ratios$FRd28 >=4] <- "Yes"
ratios$conv28 [ratios$FRd28 <4] <- "No"

summarise_logmean <- function(vec, round_to = 0) {
  vec <- na.omit(vec)
  total <- length(vec)
  log_vec <- log(vec)
  logmean <- mean(log_vec)
  logse <- sd(log_vec) / sqrt(total)
  logmargin <- 1.96 * logse
  loglow <- logmean - logmargin
  loghigh <- logmean + logmargin
  mean <- exp(logmean)
  low <- exp(loglow)
  high <- exp(loghigh)
  tibble(total, mean, low, high)
}

##FIGURE 1 all strains by prior cohort----
box_all_H1strains_V1 <- ggplot(subset(data_extra, stype=="H1N1" & timepoint == 1 & cohort %in% c("prior infection","prior vaccination", "both prior")
                                      & pid != "9-06"), aes(abbr, titre,color = class, shape = pid,label = "")) +
  facet_grid(~ cohort) +
  geom_jitter(size=1.8,width=0.2,height = 0.1, show.legend = FALSE) +
  geom_text (hjust=0, vjust=0)+
  scale_shape_manual(values = c(0,16,1,2,6,5,9,8,3,16,10,1,17,5,6,8,9,3,16,17,18)) +
  scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:15)) +
  coord_cartesian(ylim=c(4, 4000)) +
  geom_abline(intercept= log10(40), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  scale_color_manual(values = c("red2","#ff00ff", "#69ba4c","skyblue2"))+
  theme_classic() + 
  ggtitle("A/H1N1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 10, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")+
  stat_summary(fun.y=mean, geom="point", shape=95, size=7, color="black")

box_all_H1strains_V1

box_all_H3strains_V1 <- ggplot(subset(data_extra, stype=="H3N2" & timepoint == 1 & cohort %in% c("prior infection","prior vaccination","both prior")
                                      & pid != "9-06"& !abbr %in% c("Sw13","NC14","Ka17","Sy18","Ncas18")), aes(abbr, titre,color = class, shape = pid,label = "")) +
  facet_grid(~ cohort) +
  geom_jitter(size=1.8,width=0.2,height = 0, show.legend = FALSE) +
  geom_text (hjust=0, vjust=0)+
  scale_shape_manual(values = c(0,16,1,2,6,5,9,8,3,16,10,1,17,5,6,8,9,3,16,17,18)) +
  scale_y_log10("Baseline HI antibody titre", breaks = 5 * 2^(0:15)) +
  coord_cartesian(ylim=c(4, 4000)) +
  geom_abline(intercept= log10(40), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  scale_color_manual(values = c("red2","#ff00ff", "#69ba4c","skyblue2"))+
  theme_classic() + 
  ggtitle("A/H3N2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")+
  stat_summary(fun.y=mean, geom="point", shape=95, size=7, color="black")
box_all_H3strains_V1

baseline_titres <- ggarrange(box_all_H3strains_V1,box_all_H1strains_V1,heights=10, align = c("h"),ncol=2,nrow=1)
baseline_titres
ggsave("baseline_titres.pdf",baseline_titres, unit = "cm", width = 24, height = 11)


## FIGURE 2a,b HI titres against H3 H1 vax strains------
Figure2_titres_plot <- data_extra %>%
  filter(timepoint <4 & vax_strain==1) %>%
  ggplot(aes(timepoint_lbl, titre, colour = class, fill= class,shape = pid ,linetype = pid )) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("A(H1N1) reactive HI antibody titre", breaks = 5 * 2^(0:15)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(3,22000)) +
  facet_grid(~ stype +class) +
  geom_line(aes(group = pid), alpha = 2,lwd = 0.3) +
  geom_jitter(width=0.1,height=0.02, size=2, alpha = 0.5) +
  scale_color_manual(values = c("#FF9933","#ff00ff","#ff00ff", "#69ba4c"))+
  scale_shape_manual(values = c(0,1,2,5,6,8,9,3,10,0,1,2,5,6,8,9,3,16,17,18)) +
  scale_linetype_manual(values = c("dotted","dashed","solid","solid","solid","dashed",
                                   "dotted","dashed","solid","dotted","dashed","dotted",
                                   "dashed","solid","dotted","dashed","dotted","solid"))

Figure2_titres_plot
ggsave("Figure2_line__HI_titres_ex_907.pdf", Figure2_titres_plot, units="cm", width = 28, height = 10.6)


Figure2_2_titres_plot <- data_extra %>%
  filter(timepoint <4 & vax_strain==1) %>%
  ggplot(aes(timepoint_lbl, titre, colour = class, fill= class)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:15)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(3,22000)) +
  facet_grid(~ stype + class) +
  # geom_line(aes(group = pid), alpha = 2,lwd = 0.4) +
  # geom_jitter(width=0.1,height=0.02, size=2, shape = c(1), alpha = 0.5) +
  scale_color_manual(values = c("#FF9933","#ff00ff","#ff00ff","#69ba4c"))+
  geom_boxplot(aes(group = paste0(timepoint_lbl)),colour = "darkgrey",
               width = 2.5, alpha = 0.2, fill = NA, outlier.alpha = 0.5, lwd = 0.3,outlier.shape = NA)
Figure2_2_titres_plot

ggsave("Figure2_3_box_HI_titres.pdf", Figure2_2_titres_plot, units="cm", width = 28, height = 10.6)


## Figure 2 a,b stats gmts gmrs vax strains ------------------

clusWilcox.test (l2hi~ cohort + cluster(pid)+stratum(timepoint), data= subset(data_extra, vax_strain == 1 & cohort <3 & timepoint<4))

gmts_HI <- data_extra %>% 
  filter(vax_strain ==1 & study_year ==1)%>%
  group_by(stype, class, timepoint_lbl) %>% 
  summarise(.groups = "drop", summarise_logmean(titre))
gmts_HI

gmrs_HI_7 <- ratios %>% 
  filter(vax_strain ==1 & study_year ==1)%>%
  group_by(stype, class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRd7))
gmrs_HI_7

gmrs_HI_28 <- ratios %>% 
  filter(vax_strain ==1  & study_year ==1)%>%
  group_by(stype, class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRd28))
gmrs_HI_28

stats_vax_strain <- subset(data_extra,vax_strain==1 & class %in% c("prior A(H3N2)","prior vaccination") & timepoint<4)%>% 
  group_by(timepoint_lbl, stype) %>%
  select (class, l2hi) %>%
  wilcox_test(l2hi ~ class)
stats_vax_strain

stats_vax_strain_gmr_7 <- subset(ratios, vax_strain == 1 & study_year ==1 & class %in% c("prior A(H3N2)","prior vaccination")) %>% 
  group_by(stype) %>%
  select (class,FRd7) %>%
  wilcox_test(FRd7 ~ class)
stats_vax_strain_gmr_7

stats_vax_strain_gmr_28 <- subset(ratios, vax_strain == 1 & study_year ==1 & class %in% c("prior A(H3N2)","prior vaccination")) %>% 
  group_by(stype) %>%
  select (class,FRd28) %>%
  wilcox_test(FRd28 ~ class)
stats_vax_strain_gmr_28

##Table S4 restructure to wide to compare H1 and H3  ------
H1vH3 <-subset(data_extra,select = c("pid","class","cohort","titre","l2hi","stype", "visit","vax_strain","study_year"))
H1vH3 <- H1vH3 %>% 
  filter(vax_strain == 1 & study_year == 1 & visit >1)%>% 
  group_by(pid, visit, class) %>% 
  pivot_wider(names_from = stype, values_from = c(titre,l2hi))

H1vH3$logratio<- H1vH3$l2hi_H3N2/H1vH3$l2hi_H1N1
H1vH3$ratio<- H1vH3$titre_H3N2/H1vH3$titre_H1N1
write.csv(H1vH3,"H1vH3ratio.csv")

###Table S4 geom means and 95% CI ----
gmt_H1H3ratios <- H1vH3 %>% 
  group_by(visit, class) %>% 
  summarise(.groups = "drop", summarise_logmean(ratio))
gmt_H1H3ratios

stats_H1H3_logratio <- H1vH3%>% 
  filter(class %in% c("prior A(H3N2)", "prior vaccination"))%>% 
  group_by(visit) %>%
  select (class,logratio) %>%
  wilcox_test(logratio ~ class)
stats_H1H3_logratio 

stats_H1H3_ratio <- H1vH3%>% 
  filter(class %in% c("prior A(H3N2)", "prior vaccination"))%>% 
  group_by(visit) %>%
  select (class,ratio) %>%
  wilcox_test(ratio ~ class)
stats_H1H3_ratio 

H1vH3_pi <- subset (H1vH3, class == "prior A(H3N2)")
H1vH3_pv <- subset (H1vH3, class == "prior vaccination")

t.test(H1vH3_pi$l2hi_H1N1, H1vH3_pi$l2hi_H3N2, paired = T)
t.test(H1vH3_pv$l2hi_H1N1, H1vH3_pv$l2hi_H3N2, paired = T)

### restructure ratios wide to compare H1 and H3 (not included in mansucript)-----
H1vH3_FR <-subset(ratios,select = c("pid","class","FRd28","stype","vax_strain","study_year"))
H1vH3_FR <- H1vH3_FR %>% 
  filter(vax_strain == 1)%>% 
  group_by(pid,class, study_year) %>% 
  pivot_wider(names_from = stype, values_from = c(FRd28))
H1vH3_FR$stype_ratio<- H1vH3_FR$H3N2/H1vH3_FR$H1N1

mean_H1H3FRratios <- H1vH3_FR  %>% 
  group_by(class) %>% 
  summarise(.groups = "drop", summarise_logmean(stype_ratio))
mean_H1H3FRratios

##Figure 3 titre line graphs vax strains ---------------
line_vax_strain <- ggplot(subset(data_extra, vax_strain == 1  & class %in% c("prior A(H3N2)","prior vaccination") & timepoint<7 &
                                   pid %in%c("8-07","9-01","9-02","9-03","9-04","9-10","9-11","9-12")), 
                          aes(timepoint_lbl, titre, group = pid, color = pid,linetype = pid, shape = class,label = "")) +
  facet_wrap(~ stype+ class) +
  geom_line(alpha = 1, show.legend = FALSE, size = 0.6)+
  geom_point(size=3, show.legend = FALSE) +
  geom_text (hjust=0, vjust=0,show.legend = FALSE)+
  scale_shape_manual(values = c(0, 1,2,3)) +
  scale_color_manual(values = c("#000000","yellowgreen", "#56B4E9", "#0072B2","#D55E00","green3","blue2","firebrick2", "#E69F00","#F0E442", "#CC79A7",
                                "magenta","aquamarine2", "green3","orange","grey33","violetred2" ))+
  scale_linetype_manual(values = c(1,2,1,6,1,4,1,2,1,4,1,2,6,1,2,1,6,1,4,1,2,1,4,1,2,6)) +
  scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:16), limits = c(1,6000)) +
  geom_abline(intercept=40, slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")
line_vax_strain
ggsave("Fig3 line_vax_strain_HI.pdf",line_vax_strain, unit = "cm", width = 12, height = 12.5)

## Figure S4 egg cell titre ratios------
gmrs_e_c <- data_extra %>% 
  filter(clade == "3c2a1b + 131K"  & class %in% c("prior A(H3N2)","prior vaccination") & 
           ! is.na (gmr) & study_year ==1)%>%
  group_by(cell, class, visit) %>% 
  summarise(.groups = "drop", summarise_logmean(ratio))
gmrs_e_c

line_gmr_ec <- ggplot(subset(data_extra,
                             clade == "3c2a1b + 131K"  & class %in% c("prior A(H3N2)","prior vaccination") & 
                               ! is.na (gmr) & study_year ==1), 
                      aes(cell, gmr, color = class,label = "")) +
  facet_wrap(~  class+visit ) +
  geom_line(aes(group = pid),alpha = 1, show.legend = FALSE, size = 0.6,position = position_dodge(width = 0.2),)+
  geom_point(aes(group = pid),position = position_dodge(width = 0.2),size=3) +
  geom_text (hjust=0, vjust=0,show.legend = FALSE)+
  scale_shape_manual(values = c(1,2,6,8,9,3,10,0,1,2,5,6,9,3,16,17,18)) +
  scale_color_manual(values = c("#ff00ff", "#69ba4c"))+
  scale_linetype_manual(values = c(1,2,1,6,1,4,1,2,1,4,1,2,6,1,2,1,6,1,4,1,2,1,4,1,2,6)) +
  scale_y_log10("HI antibody titre fold-rise", breaks = 1 * 2^(-1:7)) +
  coord_cartesian(ylim=c(0.25, 64)) +
  scale_x_continuous (breaks = c(0,1), labels =c("egg", "cell")) +
  geom_abline(intercept= log10(4), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=9, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=10),
        strip.text = element_text(size=12),
        legend.position="right")+
  stat_summary(aes(label=round(..y..,1)),fun.y="mean", geom="text", 
               position=position_nudge(x = 0, y = 0), fontface="bold", color = "black",size =3.5)
line_gmr_ec

line_gmr_ec <- ggplot(subset(data_extra,
                             clade == "3c2a1b + 131K"  & class %in% c("prior A(H3N2)","prior vaccination") & 
                               ! is.na (gmr) & study_year ==1), 
                      aes(cell, gmr, color = class,shape= pid, linetype=pid,label = "")) +
  facet_wrap(~  class+visit ) +
  geom_line(aes(group = pid),alpha = 1, show.legend = FALSE, size = 0.6,position = position_dodge(width = 0.2),)+
  geom_point(aes(group = pid),position = position_dodge(width = 0.2),size=3) +
  geom_text (hjust=0, vjust=0,show.legend = FALSE)+
  scale_shape_manual(values = c(1,2,6,8,9,3,10,0,1,2,5,6,9,3,16,17,18)) +
  scale_color_manual(values = c("#ff00ff", "#69ba4c"))+
  scale_linetype_manual(values = c(1,2,1,6,1,4,1,2,1,4,1,2,6,1,2,1,6,1,4,1,2,1,4,1,2,6)) +
  scale_y_log10("HI antibody titre fold-rise", breaks = 1 * 2^(-1:7)) +
  coord_cartesian(ylim=c(0.25, 64)) +
  scale_x_continuous (breaks = c(0,1), labels =c("egg", "cell")) +
  geom_abline(intercept= log10(4), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=9, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=10),
        strip.text = element_text(size=12),
        legend.position="right")
line_gmr_ec

ggsave ("FigS4_ec_gmr.pdf", line_gmr_ec, units = "cm", height = 18, width = 16)

stats_gmrs_ec <- subset(data_extra, clade == "3c2a1b + 131K"   & class %in% c("prior A(H3N2)","prior vaccination") & !is.na (gmr) & study_year ==1 & stype=="H3N2")%>% 
  group_by(visit, class) %>%
  select (cell, gmr) %>%
  wilcox_test(gmr ~ cell)
stats_gmrs_ec

## Average titres across all viruses within a subtype for Table S2-------------
summarise_mean <- function(arr) {
  arr <- na.omit(arr)
  mn <- mean(arr)
  se <- sd(arr) / sqrt(length(arr))
  err_margin <- qnorm(0.975) * se
  tibble(
    mean = mn,
    low = mn - err_margin,
    high = mn + err_margin
  )
}

V1_mean_titre <- data_extra %>%
  filter (timepoint==1, !abbr %in% c("Sw13","NC14","Ka17","Sy18","Ncas18"))%>%
  group_by(stype,pid, class)%>%
  summarise(summarise_mean(titre), .groups = "drop")
V1_mean_titre
write.csv(V1_mean_titre, "V1_mean_titre.csv", row.names=FALSE)

##Breadth--------
###Table 2 exceed20 (~ 40+ ~ seropositive)---------------
prop_ex20 <- data_extra %>%
  group_by(stype,study_year, pid,class,timepoint) %>%
  summarise(percent = sum(exceed20)/sum(titre!="NA"))
prop_ex20

meanbreadth_ex20 <- prop_ex20 %>%
  filter (study_year ==1) %>%
  group_by(timepoint,class,stype) %>%
  summarise(summarise_mean(percent), .groups = "drop")
write.csv(meanbreadth_ex20, "meanbreadth_ex20.csv", row.names=FALSE)

stats_breadth_ex20 <- prop_ex20 %>%
  filter (study_year ==1) %>%
  group_by(timepoint,stype)%>%
  select (percent, class) %>%
  wilcox_test( percent~ class, comparisons = list(c("prior A(H3N2)","prior vaccination")))
stats_breadth_ex20

###Table 2 breadth seroconvert -----
prop_conv <- data_extra %>%
  filter (study_year ==1, !is.na(conv),timepoint!="1", cohort!="naive") %>%
  group_by(stype,study_year, pid,class,timepoint) %>%
  summarise(percent = sum(conv)/sum(titre!="NA"))
prop_conv

mean_prop_conv <- prop_conv %>%
  group_by(timepoint,class,stype) %>%
  summarise(summarise_mean(percent), .groups = "drop")
write.csv(mean_prop_conv_r, "mean_prop_conv_r.csv", row.names=FALSE)

stats_breadthconv <-prop_conv %>%
  filter (study_year ==1) %>%
  group_by(timepoint,stype)%>%
  select (percent, class) %>%
  wilcox_test( percent~ class, comparisons = list(c("prior A(H3N2)","prior vaccination")))
stats_breadthconv

## Figure S1 S2 plot titres across H3N2 viruses--------------------------
indiv_plots <- data_extra %>%
  filter(study_year<3,stype=="H3N2")%>%
  mutate(abbr=fct_relevel(abbr,c("Sw13","NC14","Ncas16","Si16e","Ka17","Sw17","Sw17e","Br18","Sy18","Ncas18","SA19e",
                                 "HK19e" ))) %>%
  group_by(pid) %>%
  group_split() %>%
  map(function(onepid) {
    info <- onepid %>% select(pid, cohort,study_year, age) %>% distinct()
    info_str <- paste(paste0(colnames(info), ": ", info), collapse = " | ")
    plot <- onepid %>%
      ggplot(aes(abbr, titre, color = timepoint_lbl, shape = timepoint_lbl, lty = timepoint_lbl)) +
      theme_bw()+
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position =  "bottom"
      ) +
      scale_color_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                  "2~d7", "2~d28"),
                         labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                         values = c("turquoise","darkorange","orchid","#2e9d32","#927500","magenta"))+
      scale_shape_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                  "2~d7", "2~d28"),
                         labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                         values = c(1,2,3,4,6,5))+
      scale_linetype_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                     "2~d7", "2~d28"),
                            labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                            values = c("solid","solid","solid","dashed","dashed", "dashed"))+
      scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:20)) +
      coord_cartesian(ylim = c(5, 12000)) +
      labs(caption = info_str) +
      geom_line(aes(group = paste0(pid, timepoint_lbl))) +
      geom_point(size=2.5)
    attr(plot, "pidyear") <- paste(info$pid,info$study_year, sep = "_")
    plot
  })

walk(indiv_plots, function(pl) {
  pidyear <- attr(pl, "pidyear")
  ggsave(paste0("S:/Group/AF_LC_RT_share/WHO IMP CEIRS infant study/Write up papers conferences/comm med revision", pidyear, ".pdf"), pl, width = 12, height = 9, units = "cm")
})

## Figure S1 S2 plot titres across H1N1 viruses--------------------------
indiv_plotsH1 <- data_extra %>%
  filter(study_year<3,stype=="H1N1")%>%
  mutate(abbr=fct_relevel(abbr,c("Ca09","Mi15","Mi15e","Br18","Br18e","Vi19e","Sy21"))) %>%
  group_by(pid) %>%
  group_split() %>%
  map(function(onepid) {
    info <- onepid %>% select(pid, cohort,study_year, age) %>% distinct()
    info_str <- paste(paste0(colnames(info), ": ", info), collapse = " | ")
    plot <- onepid %>%
      ggplot(aes(abbr, titre, color = timepoint_lbl, shape = timepoint_lbl, lty = timepoint_lbl)) +
      theme_bw()+
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position =  "bottom"
      ) +
      scale_color_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                  "2~d7", "2~d28"),
                         labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                         values = c("turquoise","darkorange","orchid","#2e9d32","#927500","magenta"))+
      scale_shape_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                  "2~d7", "2~d28"),
                         labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                         values = c(1,2,3,4,6,5))+
      scale_linetype_manual(breaks=c("1-d0","1~d7","1~d28","2-d0",
                                     "2~d7", "2~d28"),
                            labels=c("Y1 d0","Y1 d7","Y1 d28","Y2 d0","Y2 d7","Y2 d28"),
                            values = c("solid","solid","solid","dashed","dashed", "dashed"))+
      scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:20)) +
      coord_cartesian(ylim = c(5, 12000)) +
      labs(caption = info_str) +
      geom_line(aes(group = paste0(pid, timepoint_lbl))) +
      geom_point(size=2.5)
    attr(plot, "pidyear") <- paste(info$pid,info$study_year, sep = "_")
    plot
  })

walk(indiv_plotsH1, function(pl) {
  pidyear <- attr(pl, "pidyear")
  ggsave(paste0("S:/Group/AF_LC_RT_share/WHO IMP CEIRS infant study/Write up papers conferences/comm med revisionH1", pidyear, ".pdf"), pl,width = 7, height = 9,  units = "cm")
})


#READ ELLA/NAI antibody data------------

setwd("S:/Group/AF_LC_RT_share/WHO IMP CEIRS infant study/Write up papers conferences/comm med revision")
# data <- read.csv("Ella_long_diff.csv",header = T, stringsAsFactors = F)
data_ella <- read_excel("WHOIMP_dataset.xlsx",sheet = "Ella_NI_titres")

data_ella_rn <- data_ella %>%
  mutate(
    pid = PID, timepoint = Visit2Y, visit = Visit, virus = Virus,
    abbr=Abbr,vax_strain = vax_strain,cohort = factor(Cohort),virus_year = VirusYear,
    , class = factor(Classification, levels =c("no/unclear prior","prior influenza A", "prior A(H3N2)", "prior vaccination","prior infection & vaccination")),
    study_year = StudyYear, sample_year = Year, cell = Egg_Cell, stype = factor(Subtype, levels=c("A/H3N2","A/H1N1"), labels =c("A(H3N2)","A(H1N1)")),titre = IC50_titre)

data_ella_rn$l2ni <- log(data_ella_rn$titre)/log(2)

#analyse prior effect by time
data_ella_ex <- data_ella_rn %>%
  mutate(
    timepoint_lbl = factor(
      timepoint, 1:8, c("1-d0", "1~d7", "1~d28","2-d0","2~d7","2~d28","3-d0","3~d28")
    )
  ) 
#exclude 9939-0007 who had both prior infection and prior vaccination
data_ella_ex <- data_ella_ex %>%
  filter (pid != "9-07")


##restructure ELLA to wide look at fold rise ------
ratios_ella <-subset(data_ella_ex,select = c("pid","class","cohort","study_year","visit","cell","titre","stype","abbr", "vax_strain"))

ratios_ella <- ratios_ella %>% 
  group_by(abbr,study_year, visit,pid) %>% 
  pivot_wider(names_from = visit, values_from = titre)
ratios_ella <- ratios_ella %>% 
  unnest (cols = c(`1`, `2`,`3`))
ratios_ella$FRd7<- ratios_ella$`2`/ratios_ella$`1`
ratios_ella$FRd28<- ratios_ella$`3`/ratios_ella$`1`
ratios_ella$conv7 [ratios_ella$FRd7 >=3] <- "Yes"
ratios_ella$conv7 [ratios_ella$FRd7 <3] <- "No"
ratios_ella$conv28 [ratios_ella$FRd28 >=3] <- "Yes"
ratio_ellas$conv28 [ratios_ella$FRd28 <3] <- "No"

#functions----
summarise_logmean <- function(vec, round_to = 0) {
  vec <- na.omit(vec)
  total <- length(vec)
  log_vec <- log(vec)
  logmean <- mean(log_vec)
  logse <- sd(log_vec) / sqrt(total)
  logmargin <- 1.96 * logse
  loglow <- logmean - logmargin
  loghigh <- logmean + logmargin
  mean <- exp(logmean)
  low <- exp(loglow)
  high <- exp(loghigh)
  tibble(total, mean, low, high)
}

summarise_prop <- function(vec) {
  vec <- na.omit(vec)
  success <- sum(vec)
  total <- length(vec)
  ci <- PropCIs::exactci(success, total, 0.95)
  prop <- success / total
  low <- ci$conf.int[[1]]
  high <- ci$conf.int[[2]]
  f <- function(x) round(x * 100)
  tibble(
    prop, low, high,
    comb = glue::glue("{f(prop)}% ({f(low)}%, {f(high)}%)")
  )
}

##NAI gmts gmrs conv vax strains-----------------

gmts_NI <- data_ella_ex %>% 
  filter(vax_strain == 1 & study_year ==1)%>%
  group_by(stype, class, timepoint_lbl) %>% 
  summarise(.groups = "drop", summarise_logmean(titre))
gmts_NI

gmrs_NI_7 <- ratios_ella %>% 
  filter(vax_strain ==1 &  study_year ==1)%>%
  group_by(stype, class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRd7))
gmrs_NI_7

gmrs_NI_28 <- ratios_ella %>% 
  filter(vax_strain ==1 &  study_year ==1)%>%
  group_by(stype, class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRd28))
gmrs_NI_28

###Figure 2B Vaccine Strain NAI titres -----
Figure2_NItitres_plot <- data_ella_ex %>%
  filter(timepoint <4 & vax_strain==1 ) %>%
  ggplot(aes(timepoint_lbl, titre, colour = class,fill= class, shape = pid ,linetype = pid )) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("NI antibody titre", breaks = 5 * 2^(0:15)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(1,29000)) +
  facet_grid(~ stype + class) +
  geom_line(aes(group = pid), alpha = 2,lwd = 0.3) +
  geom_jitter(width=0.1,height=0.02, size=2, alpha = 0.5) +
  scale_color_manual(values = c("#FF9933","#ff00ff","#ff00ff", "#69ba4c"))+
  scale_shape_manual(values = c(0,1,2,5,6,8,9,3,10,0,1,2,5,6,8,9,3,16,17,18)) +
  scale_linetype_manual(values = c("dotted","dashed","solid","solid","solid","dashed",
                                   "dotted","dashed","solid","dotted","dashed","dotted",
                                   "dashed","solid","dotted","dashed","dotted","solid"))
Figure2_NItitres_plot
ggsave("Figure2_NI_titres_exc907.pdf", Figure2_NItitres_plot, units="cm", width = 28, height = 10.6)


Figure2_2_NItitres_plot <-data_ella_ex %>%
  filter(timepoint <4 & vax_strain==1 ) %>%
  ggplot(aes(timepoint_lbl, titre, colour = class,fill= class)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("NI antibody titre", breaks = 5 * 2^(0:15)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(1,29000)) +
  facet_grid(~ stype + class) +
  geom_boxplot(aes(group = paste0(timepoint_lbl)),color = "grey", 
               width = 2.5, alpha = 0.2, fill = NA, outlier.alpha = 0.5, lwd = 0.3,outlier.shape = NA) 

Figure2_2_NItitres_plot
ggsave("Figure2_2_NI_titres_exc907.pdf", Figure2_2_NItitres_plot, units="cm", width = 28, height = 10.6)

statsNI_Y1 <- subset(data_ella_ex,vax_strain == 1 & class %in% c("prior A(H3N2)","prior vaccination") & study_year ==1)%>% 
  group_by(timepoint_lbl, stype)%>%
  select (class, l2ni) %>%
  wilcox_test(l2ni ~ class)
statsNI_Y1

#READ and re-structure stem antibody data -----
data_stem <- read_excel("WHOIMP_dataset.xlsx",sheet = "H3_HA_stem_binding")

data_stem_rn <- data_stem %>%
  mutate(
    pid = pid, timepoint = Visit2Y, visit = visit, cohort = factor(Cohort),
    class = factor(Classification, levels =c("no/unclear prior","prior influenza A", "prior A(H3N2)", "prior vaccination","prior infection & vaccination")),
    study_year = StudyYear, sample_year = Year,titre = Titre, vax_day = vaccine_day)


data_stem_rn$l2HA <- log(data_stem_rn$titre)/log(2)

#analyse prior effect by time
data_stem_ex <- data_stem_rn %>%
  mutate(
    timepoint_lbl = factor(
      timepoint, 1:8, c("1-d0", "1~d7", "1~d28","2-d0","2~d7","2~d28","3-d0","3~d28")
    )
  )

#exclude 9-07 who had both prior infection and prior vaccination
data_extra_ex <- data_stem_ex %>%
  filter (pid != "9-07")

## stem antibody ratios and conversion -----
ratios_stem <-subset(data_stem_ex ,!is.na (titre),select = c("pid","timepoint","titre",
                                                             "class","cohort"))
ratios_stem <- ratios_stem %>% 
  group_by(pid) %>% 
  pivot_wider(names_from = timepoint, values_from = titre)

ratios_stem <- ratios_stem %>% 
  unnest (cols = c(`1`, `2`, `3`,`4`,`5`,`6`,`7`,`8`))
ratios_stem$FRV2<- ratios_stem$`2`/ratios_stem$`1`
ratios_stem$FRV3 <-ratios_stem$`3`/ratios_stem$`1`
ratios_stem$FRV5<- ratios_stem$`5`/ratios_stem$`4`
ratios_stem$FRV6 <-ratios_stem$`6`/ratios_stem$`4`
ratios_stem$conv [ratios_stem$FRV2 >= 4] <- "Yes"
ratios_stem$conv [ratios_stem$FRV2 < 4] <- "No"
ratios_stem$conv3 [ratios_stem$FRV3 >= 4] <- "Yes"
ratios_stem$conv3 [ratios_stem$FRV3 < 4] <- "No"

##summarize stem antibody data----
summarise_logmean <- function(arr) {
  logarr <- log(arr)
  logmean <- mean(logarr)
  logse <- sd(logarr) / sqrt(length(arr))
  logerr_margin <- qnorm(0.975) * logse
  tibble(
    mean = exp(logmean),
    low = exp(logmean - logerr_margin),
    high = exp(logmean + logerr_margin)
  )
}

gmts_HAst <- data_stem_ex %>% 
  filter(! is.na(titre))%>%
  group_by(class, timepoint_lbl) %>% 
  summarise(.groups = "drop", summarise_logmean(titre))
gmts_HAst

gmrs_HAst_7 <- ratios_stem %>% 
  filter(! is.na(FRV2))%>%
  group_by(class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRV2))
gmrs_HAst_7

gmrs_HAst_28 <- ratios_stem %>% 
  filter(! is.na(FRV3))%>%
  group_by(class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRV3))
gmrs_HAst_28

gmrs_HAst_7Y2 <- ratios_stem %>% 
  filter(! is.na(FRV5))%>%
  group_by(class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRV5))
gmrs_HAst_7Y2

gmrs_HAst_28Y2 <- ratios_stem %>% 
  filter(! is.na(FRV6))%>%
  group_by(class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRV6))
gmrs_HAst_28Y2stats_stem_titre <- data_stem_ex %>%
  filter(timepoint == 2 & !is.na(titre))%>%
  select (titre, class) %>%
  wilcox_test(titre ~ class, comparisons = list (c("prior A(H3N2)","prior vaccination")))
stats_stem_titre

stats_stem_titre <- data_stem_ex %>%
  filter(timepoint == 3 & !is.na(titre))%>%
  select (titre, class) %>%
  wilcox_test(titre ~ class, comparisons = list (c("prior A(H3N2)","prior vaccination")))
stats_stem_titre

##Figure 3B HA stem TITRES -----
Figure_HAstem_titres_plot <- data_stem_ex %>%
  filter(! is.na(titre),timepoint<7, pid!= "9-07") %>%
  ggplot(aes(timepoint_lbl, titre, colour = class, linetype = pid, shape = pid, label = pid)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("HA stem antibody titre", breaks = 5 * 2^(0:15)) +
  coord_cartesian(ylim=c(1,1200)) +
  facet_grid(~ class) +
  geom_boxplot(aes(group = paste0(timepoint_lbl)), 
               width = 0.8, alpha = 0.1, fill = NA, outlier.alpha = 0.5, lwd = 0.6,outlier.shape = NA, color ="grey")+
  geom_line(aes(group = pid), alpha = 2,lwd = 0.3) +
  geom_point() +
  geom_text_repel (colour = "black",max.overlaps = 50, size = 2)+
  scale_linetype_manual(values = c("dotted","dashed","solid","dotted","solid","dashed",
                                   "solid","dashed","solid","solid","dashed","dotted",
                                   "dashed","solid","dotted","dashed","solid",
                                   "dotted","dashed","solid","dotted","dashed","solid"))+
  scale_shape_manual(values = c(0,0,16,1,2,5,6,9,8,3,16,10,1,17,5,6,8,9,3,16,17,18)) +
  scale_color_manual(values = c("#FF9933","#ff00ff","#ff00ff", "#69ba4c"))
Figure_HAstem_titres_plot

ggsave("Figure3B_H3HAstem_titres5.pdf", Figure_HAstem_titres_plot, units="cm", width = 18, height = 10)


Figure_HAstem_box <- data_stem_ex %>%
  filter(! is.na(titre),timepoint<7, pid!= "9-07") %>%
  ggplot(aes(timepoint_lbl, titre))+
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("HA stem antibody titre", breaks = 5 * 2^(0:15)) +
  coord_cartesian(ylim=c(1,1200)) +
  facet_grid(~ class) +
  geom_line(aes(group = pid), alpha = 2,lwd = 0.3) +
  geom_point() +
  geom_boxplot(aes(group = paste0(timepoint_lbl)), 
               width = 0.8, alpha = 0.1, fill = NA, outlier.alpha = 0.5, lwd = 0.6,outlier.shape = NA, color ="grey")

Figure_HAstem_box 
ggsave("Figure3B_H3HAstem_box.pdf", Figure_HAstem_box, units="cm", width = 18, height = 10)


# READ data HI antibodies against influenza B and restructure----
data_B <- read_excel("WHOIMP_dataset.xlsx",sheet = "B_HI")

data_B$l2hi <- log(data_B$Titer)/log(2)
data_Brenamed <- data_B %>%
  mutate(
    pid = PID, timepoint = Visit2Y, visit = Visit, sex = factor(`Sex (M=0,F=1)`, levels = c("0","1"), labels =c("female","male")),
    age = `Age at Vaccination`, virus = Virus, cohort = factor(Cohort), class = factor(Classification, levels =c("no/unclear prior","prior influenza A", "prior A(H3N2)", "prior vaccination", "prior infection & vaccination")),
    virus_year = VirusYear,vax_strain = vax_strain,study_year = StudyYear, sample_year = SampleYear, virus_year = VirusYear,cell = Egg_Cell, lineage = factor(Subtype, levels=c("B Vic","B Yam")),
    titre = Titer,l2hi = l2hi,gmr = ratio, vaccine_day = vaccine_day
  )

data_Brenamed$l2gmr <- log(data_Brenamed$gmr)/log(2)

data_Bextra <- data_Brenamed %>%
  mutate(
    timepoint_lbl = factor(
      timepoint, 1:8, c("1-d0", "1~d7", "1~d28","2-d0","2~d7","2~d28","3-d0","3~d28")
    )
  ) 

##restructure to wide to check fold rise ------
ratiosB <-subset(data_Bextra,select = c("pid","cohort","class","study_year","visit","cell","titre","lineage","virus","vax_strain"))
ratiosB <- ratiosB %>% 
  group_by(virus,study_year, visit,pid) %>% 
  pivot_wider(names_from = visit, values_from = titre)
ratiosB <- ratiosB %>% 
  unnest (cols = c(`1`, `2`,`3`))
ratiosB$FRd7<- ratiosB$`2`/ratiosB$`1`
ratiosB$FRd28<- ratiosB$`3`/ratiosB$`1`
ratiosB$conv7 [ratiosB$FRd7 >=4] <- "Yes"
ratiosB$conv7 [ratiosB$FRd7 <4] <- "No"
ratiosB$conv28 [ratiosB$FRd28 >=4] <- "Yes"
ratiosB$conv28 [ratiosB$FRd28 <4] <- "No"

##B HI data summaries-----
statsB_vax_strain_gmt_Y1 <- subset(data_Bextra, vax_strain == 1 & timepoint < 4 &  class %in% c("prior A(H3N2)","prior vaccination"))%>% 
  group_by(timepoint_lbl,lineage) %>%
  select (class, l2hi) %>%
  wilcox_test(l2hi ~ class)
statsB_vax_strain_gmt_Y1

gmts_HI_B <- data_Bextra %>% 
  filter(vax_strain ==1 & study_year ==1)%>%
  group_by(lineage, class, timepoint_lbl) %>% 
  summarise(.groups = "drop", summarise_logmean(titre))
gmts_HI_B

gmrs_B_28 <- ratiosB %>% 
  filter(vax_strain ==1 &  study_year ==1)%>%
  group_by(lineage, class) %>% 
  summarise(.groups = "drop", summarise_logmean(FRd28))
gmrs_B_28

##Figure S3 titer plots B vax strains-----
B_titres_plot <- data_Bextra %>%
  filter(timepoint <4 & vax_strain==1 & class != "prior infection & vaccination") %>%
  ggplot(aes(timepoint_lbl, titre, colour = class, fill= class,shape = pid ,linetype = pid, label = pid)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("HI antibody titre", breaks = 5 * 2^(0:13)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(3,5120)) +
  facet_grid(~ lineage + class) +
  geom_line(aes(group = pid), alpha = 2,lwd = 0.3) +
  geom_jitter(width=0.1,height=0.00, size=2, alpha = 0.5) +
  # geom_text_repel (aes(label=pid,segment.colour = NA,max.overlaps = 50))+
  scale_color_manual(values = c("#FF9933","#ff00ff","#ff00ff", "#69ba4c"))+
  scale_shape_manual(values = c(0,1,2,5,6,8,9,3,10,0,1,2,5,6,8,9,3,16,17,18)) +
  scale_linetype_manual(values = c("dotted","dashed","solid","solid","solid","dashed",
                                   "dotted","dashed","solid","dotted","dashed","dotted",
                                   "dashed","solid","dotted","dashed","dotted","solid"))
B_titres_plot


B_titres_boxplot <- data_Bextra %>%
  filter(timepoint <4 & vax_strain==1 & class != "prior infection & vaccination") %>%
  ggplot(aes(timepoint_lbl, titre, colour = class, fill= class)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position="none")+
  scale_y_log10("B Vic reactive HI antibody titre", breaks = 5 * 2^(0:13)) +
  scale_x_discrete ("Day", breaks = c("1-d0", "1~d7","1~d28"), labels = c("0", "~7","~14")) +
  coord_cartesian(ylim=c(3,5120)) +
  facet_grid(~ lineage + class) +
  geom_boxplot(aes(group = paste0(timepoint_lbl)), color = "grey",
               width = 2.5, alpha = 0.2, fill = NA, outlier.alpha = 0.5, lwd = 0.5,outlier.shape = NA)+
  stat_summary(aes(label=round(..y..,0)),fun.y="mean", geom="text", 
               position=position_nudge(x = 0, y = 0), fontface="bold", color = "black",size=3.5) 
B_titres_boxplot

#READ Supplementary data for comparison of titres across subtypes, viruses and assays-----
data_comp <- read_excel("WHOIMP_dataset.xlsx",sheet = "cross-assay comparison")

data_comp_rn <- data_comp %>%
  mutate(timepoint = Visit2Y, visit = Visit,  age = AgeAtVax, cohort = factor(Cohort),
         class = factor(Classification, levels =c("no/unclear prior","prior influenza A",
                                                  "prior A(H3N2)", "prior vaccination","prior infection & vaccination")),
         study_year = StudyYear,
         timepoint_lbl = factor(
           timepoint, 1:8, c("1-d0", "1~d7", "1~d28","2-d0","2~d7","2~d28","3-d0","3~d28")
         )
  ) 

##HI_NI ratios----
data_comp_rn$HI_NI_Vi19 <- data_comp_rn$HI_Vi19/data_comp_rn$NI_Vi19
data_comp_rn$HI_NI_SA19 <- data_comp_rn$HI_SA19/data_comp_rn$NI_SA19

stats_HI_NI_Vi19 <- subset(data_comp_rn,StudyYear==1 & class %in% c("prior A(H3N2)","prior vaccination"))%>% 
  group_by(visit) %>%
  select (class, HI_NI_Vi19) %>%
  wilcox_test(HI_NI_Vi19 ~ class)
stats_HI_NI_Vi19

box_HI_NI_ratio_Vi19 <- ggplot(subset(data_comp_rn,StudyYear<3 &  class %in% c("prior influenza A","prior A(H3N2)","prior vaccination")),
                               aes(timepoint_lbl, HI_NI_Vi19,color = class, shape = class,label = "")) +
  facet_grid(~ class) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = NA, alpha = 0.2, show.legend = FALSE, linewidth=0.6, coef=0, color = "grey")+
  geom_jitter(aes(group = pid,shape = pid),size=1.8,width=0.2,height = 0, show.legend = FALSE) +
  geom_line (aes(group = pid,linetype = pid))+
  geom_text (hjust=0, vjust=0)+
  scale_color_manual(values = c("#ff00ff","#ff00ff", "#69ba4c"))+
  scale_shape_manual(values = c(0,1,2,5,6,8,9,3,10,0,1,2,5,6,8,9,3,16,17,18)) +
  scale_linetype_manual(values = c("dotted","dashed","solid","solid","solid","dashed",
                                   "dotted","dashed","solid","dotted","dashed","dotted",
                                   "dashed","solid","dotted","dashed","dotted","solid"))+
  scale_y_log10("HI NI antibody ratio",limits = c(0.01,100), breaks = 1 * 2^(-5:10)) +
  geom_abline(intercept= log10(1), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  # stat_summary(aes(label=round(..y..,1)),fun ="mean", geom="text", 
  #              position=position_nudge(x = 0, y = 0), fontface="bold", color = "black",size =3.5)+
  theme_classic() + 
  ggtitle("A/Victoria/2570/2019 (H1N1") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=9, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")

box_HI_NI_ratio_Vi19

stats_HI_NI_SA19 <- subset(data_comp_rn,StudyYear==1 & class %in% c("prior A(H3N2)","prior vaccination"))%>% 
  group_by(visit) %>%
  select (class, HI_NI_SA19) %>%
  wilcox_test(HI_NI_SA19 ~ class)
stats_HI_NI_SA19

box_HI_NI_ratio_SA19 <- ggplot(subset(data_comp_rn,StudyYear<3 &  class %in% c("prior influenza A","prior A(H3N2)","prior vaccination")),
                               aes(timepoint_lbl, HI_NI_SA19,color = class, shape = class,label = "")) +
  facet_grid(~ class) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = NA, alpha = 0.2, show.legend = FALSE, linewidth=0.6, coef=0, color = "grey")+
  geom_jitter(aes(group = pid,shape = pid),size=1.8,width=0.2,height = 0, show.legend = FALSE) +
  geom_line (aes(group = pid,linetype = pid))+
  geom_text (hjust=0, vjust=0)+
  scale_color_manual(values = c("#ff00ff","#ff00ff", "#69ba4c"))+
  scale_shape_manual(values = c(0,1,2,5,6,8,9,3,10,0,1,2,5,6,8,9,3,16,17,18)) +
  scale_linetype_manual(values = c("dotted","dashed","solid","solid","solid","dashed",
                                   "dotted","dashed","solid","dotted","dashed","dotted",
                                   "dashed","solid","dotted","dashed","dotted","solid"))+
  scale_y_log10("HI NI antibody ratio",limits = c(0.01,100), breaks = 1 * 2^(-5:10)) +
  geom_abline(intercept= log10(1), slope=0, linetype="dotted", col="gray30", lwd=0.5) +
  theme_classic() + 
  ggtitle("A/South Australia/34/2019 (H3N2)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=9, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")

box_HI_NI_ratio_SA19

###Figure S6 HI to NI titre ratio -----
SuppFig_HI_NI_ratio_20vax <- ggarrange( box_HI_NI_ratio_Vi19, box_HI_NI_ratio_SA19, ncol=2)
SuppFig_HI_NI_ratio_20vax 
ggsave("SuppFig_HI_NI_ratio_20vax.pdf", SuppFig_HI_NI_ratio_20vax, units = "cm", width = 30, height = 11)
ggsave("FigureS3_line_BHI_titres.pdf", B_titres_plot, units="cm", width = 28, height = 10.6)


##HIvNI scatter----
scatterMi15_HI_NI <- ggplot(data_comp_rn, aes(HI_Mi15, NI_Mi15)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("HI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12)) +
  scale_y_log10("NI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12))+
  ggtitle("A/Michigan/45/2015 (H1N1)")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside =c(0.95,0.1)) +
  stat_cor()
scatterMi15_HI_NI

scatterSA19_HI_NI <- ggplot(data_comp_rn, aes(HI_SA19, NI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("HI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12)) +
  scale_y_log10("NI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12))+
  ggtitle("A/South Australia/34/2019 (H3N2)")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterSA19_HI_NI

scatterSw17_HI_NI <- ggplot(data_comp_rn, aes(HI_Sw17, NI_Sw17)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("HI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12)) +
  scale_y_log10("NI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12))+
  ggtitle("A/Switzerland/8060/2017 (H3N2)")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterSw17_HI_NI

scatterVi19_HI_NI <- ggplot(data_comp_rn, aes(HI_Vi19, NI_Vi19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("HI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12)) +
  scale_y_log10("NI antibody titre", limits = c(2.5, 15000),breaks = 5 * 2^(0:12))+
  ggtitle("A/Victoria/2570/2019 (H1N1)")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterVi19_HI_NI

scatterMi15NI_Sw17NI <- ggplot(data_comp_rn, aes(NI_Mi15, NI_Sw17)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("N1 (Mi15) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("N2 (Sw17) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("N1 versus N2")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterMi15NI_Sw17NI

scatterVi19NI_SA19NI <- ggplot(data_comp_rn, aes(NI_Vi19, NI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("N1 (Vi19) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("N2 (SA19) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("N1 versus N2 (2019)")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterVi19NI_SA19NI

scatterMi15HI_Sw17HI <- ggplot(data_comp_rn, aes(HI_Mi15, HI_Sw17)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H1 (Mi15) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("H3 (Sw17) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("H1 versus H3")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterMi15HI_Sw17HI

scatterVi19HI_SA19HI <- ggplot(data_comp_rn, aes(HI_Vi19, HI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H1 (Vi19) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("H3 (SA19) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("H1 versus H3")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterVi19HI_SA19HI

scatterMi15NI_Vi19NI <- ggplot(data_comp_rn, aes(NI_Mi15, NI_Vi19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("N1 (Mi15) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("N1 (Vi19) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("N1 cross strain")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterMi15NI_Vi19NI

scatterSw17NI_SA19NI <- ggplot(data_comp_rn, aes(NI_Sw17, NI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("N2 (Sw17) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("N2 (SA19) NI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("N2 cross strain")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterSw17NI_SA19NI

scatterMi15HI_Vi19HI <- ggplot(data_comp_rn, aes(HI_Mi15, HI_Vi19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H1 (Mi15) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("H1 (Vi19) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("H1 cross strain")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterMi15HI_Vi19HI 

scatterSw17HI_SA19HI <- ggplot(data_comp_rn, aes(HI_Sw17, HI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H3 (Sw17) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12)) +
  scale_y_log10("H3 (SA19) HI antibody titre", limits = c(2.5, 15000), breaks = 5 * 2^(0:12))+
  ggtitle("H3 cross strain")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterSw17HI_SA19HI 

###Figure S5------
Correl_matrix_HI_NI_H3_H1 <- ggarrange(scatterSw17_HI_NI,scatterSA19_HI_NI,scatterMi15_HI_NI,scatterVi19_HI_NI,
                                       scatterSw17HI_SA19HI,scatterMi15HI_Vi19HI,scatterSw17NI_SA19NI,scatterMi15NI_Vi19NI,
                                       scatterMi15HI_Sw17HI,scatterVi19HI_SA19HI,scatterMi15NI_Sw17NI,scatterVi19NI_SA19NI,
                                       nrow=3, ncol=4)
Correl_matrix_HI_NI_H3_H1
ggsave("Figure S5 NI HI correlationsV2.pdf", Correl_matrix_HI_NI_H3_H1, unit = "cm", width = 40, height = 30)

##NI HI versus H6 binding----
scatterH6_HA_Sw17NI <- ggplot(subset(data_comp_rn,!is.na(H6_HA_binding),!is.na(NI_Sw17)), 
                              aes(H6_HA_binding, NI_Sw17)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Sw17 (N2)) NI titre", limits = c(1, 1500), breaks = 5 * 2^(-2:8))+
  ggtitle("H6HA versus N2 NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Sw17NI

scatterH6_HA_SA19NI <- ggplot(data_comp_rn, aes(H6_HA_binding, NI_SA19)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("SA19 (N2)) NI titre", limits = c(1, 1500), breaks = 5 * 2^(-2:8))+
  ggtitle("H6HA versus N2 NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_SA19NI 

scatterH6_HA_Mi15NI <- ggplot(data_comp_rn, aes(H6_HA_binding, NI_Mi15)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Mi15 (N1)) NI titre", limits = c(1, 1500), breaks = 5 * 2^(-2:8))+
  ggtitle("H6HA versus N1 NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Mi15NI

scatterH6_HA_Vi19NI <- ggplot(data_comp_rn, aes(H6_HA_binding, NI_Vi19)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Vi19 (N1)) NI titre", limits = c(1, 1500), breaks = 5 * 2^(-2:8))+
  ggtitle("H6HA versus N1 NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Vi19NI

scatterH6_HA_Sw17HI <- ggplot(data_comp_rn, aes(H6_HA_binding, HI_Sw17)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Sw17 (H3) HI titre", limits = c(1, 15000), breaks = 5 * 2^(0:11))+
  ggtitle("H6HA versus H3 HI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Sw17HI

scatterH6_HA_SA19HI <- ggplot(data_comp_rn, aes(H6_HA_binding, HI_SA19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("SA19 (H3) HI titre", limits = c(1, 15000), breaks = 5 * 2^(0:11))+
  ggtitle("H6HA versus H3 HI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_SA19HI

scatterH6_HA_Mi15HI <- ggplot(data_comp_rn, aes(H6_HA_binding, HI_Mi15)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Mi15 (H1) HI titre", limits = c(1, 15000), breaks = 5 * 2^(0:11))+
  ggtitle("H6HA versus H1 HI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Mi15HI

scatterH6_HA_Vi19HI <- ggplot(data_comp_rn, aes(H6_HA_binding, HI_Vi19)) +
  geom_jitter(width=0.1, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("H6-HA binding titre", limits = c(0.01, 100),breaks = 1 * 2^(-3:6)) +
  scale_y_log10("Vi19 (H1) HI titre",limits = c(1, 15000),breaks = 5 * 2^(0:11))+
  ggtitle("H6HA versus H1 HI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=8, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatterH6_HA_Vi19HI
###Figure S7----
Correl_matrix_H6Binding_Ab <- ggarrange(scatterH6_HA_Sw17NI,scatterH6_HA_SA19NI,scatterH6_HA_Mi15NI,scatterH6_HA_Vi19NI,
                                        scatterH6_HA_Sw17HI,scatterH6_HA_SA19HI,scatterH6_HA_Mi15HI,scatterH6_HA_Vi19HI,
                                        nrow=2, ncol=4)
Correl_matrix_H6Binding_Ab 
ggsave("Figure S7 H6 binding Ab correlations.pdf", Correl_matrix_H6Binding_Ab , unit = "cm", width = 40, height = 20)


data$L2Mi15 <- log10(data$HI_Mi15)/log10(2)
scatter_H1HI_HAst <- ggplot(subset(data,!is.na(L2stem)), aes(L2Mi15, L2stem)) +
  geom_jitter(width=0.05, height=0.0, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_continuous("HI antibody titre (Mi15)", breaks=c(2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3),labels=c(5,10,20,40,80,160,320,640,1280,2560,5120,10240)) +
  scale_y_continuous("HA stem antibody titre", breaks=c(0,1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3),labels=c(1,2.5,5,10,20,40,80,160,320,640,1280)) +
  ggtitle("H1N1")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position=c(0.95,0.1)) +
  stat_cor()
scatter_H1HI_HAst
ggsave("HIMi15_HAstem_corr.pdf",scatter_H1HI_HAst, unit = "cm", width = 11, height = 10)

data$L2H6 <- log10(data$H6_HA_binding)/log10(2)

scatter_H6_HAst <- ggplot(subset(data,!is.na(L2H6)&!is.na(L2stem)),aes(L2H6,L2stem)) +
  geom_jitter(width=0.05, height=0.0, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_continuous("HI antibody titre (Mi15)", breaks=c(2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3),labels=c(5,10,20,40,80,160,320,640,1280,2560,5120,10240)) +
  scale_y_continuous("HA stem antibody titre", breaks=c(0,1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3),labels=c(1,2.5,5,10,20,40,80,160,320,640,1280)) +
  ggtitle("")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position=c(0.95,0.1)) +
  stat_cor()
scatter_H6_HAst

##Figure S8 NI N1 N2 versus N2 binding----
scatterN2_NA_SA19NI <- ggplot(data_comp_rn, aes(N2_NA_SA19_binding, NI_SA19)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("SA19 N2 NA binding titre", limits = c(2.5,320), breaks = 5 * 2^(-1:6)) +
  scale_y_log10("SA19 N2 NI titre", limits = c(2,5120), breaks = 5 * 2^(0:10))+
  ggtitle("N2 NA binding versus NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position= "none") +
  stat_cor()
scatterN2_NA_SA19NI

scatterN2_NA_Vi19NI <- ggplot(data_comp_rn, aes(N2_NA_SA19_binding, NI_Vi19)) +
  geom_jitter(width=0.05, height=0.05, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("SA19 N2 NA binding titre", limits = c(2.5,320), breaks = 5 * 2^(-1:6)) +
  scale_y_log10("Vi19 N1 NI titre", limits = c(2,5120), breaks = 5 * 2^(0:10))+
  ggtitle("N2 NA binding versus N1 NI")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position= "none") +
  stat_cor()
scatterN2_NA_Vi19NI

N2Binding_Ab <- ggarrange(scatterN2_NA_SA19NI,scatterN2_NA_Vi19NI,nrow=1, ncol=2)
N2Binding_Ab 
ggsave("Figure S8 N2 binding Ab correlations.pdf", N2Binding_Ab, unit = "cm", width = 20, height = 10)

##Figure S9 H3 HI v HA stem ------
scatter_HI_HAst <- ggplot(subset(data_comp_rn,!is.na(HA_stem_titre)), aes(HI_Sw17,HA_stem_titre)) +
  geom_jitter(width=0.08, height=0.0, size=2, shape =1) +
  geom_smooth (method ="lm", se =TRUE, fullrange=TRUE, level=0.95, alpha =0.15)+
  scale_x_log10("HI antibody titre (Sw17)",breaks = 5 * 2^(0:11) ) +
  scale_y_log10("HA stem antibody titre",breaks = 5 * 2^(0:11) ) +
  ggtitle("H3N2")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, margin = margin(2,0,0,0)),
        axis.title.y = element_text(size = 11, margin = margin(0,8,0,0)),
        axis.text.y = element_text(size=10, margin = margin(0,0,0,0)),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=12),
        legend.position.inside=c(0.95,0.1)) +
  stat_cor()
scatter_HI_HAst
ggsave("Figure S9 HISw17_HAstem_corr.pdf",scatter_HI_HAst, unit = "cm", width = 12, height = 10)
