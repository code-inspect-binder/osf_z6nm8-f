#R script used in Bracic et al (2024)

# Author: Marko Bracic (@markobracic)
# Profile: https://orcid.org/0000-0001-6528-3572
# Email: mbracic192@gmail.com


# Title: Data analysis of mice response in the JBT

################################################################################
# Description of script and instructions
################################################################################

# Used to compare reactions to different cues and produce response curve figure in JBT
# Output: 
  # Figure 4
  # Supplementary Figure S1
  # Supplementary Tables S3/S4

# Names
  # Judgement bias test (JBT) ~ CJB (cognitive judgement bias)
  # JBT paradigms: Touchscreen (TS) and Tunnel
  # JBT Cues: 
    # "P" = "Positive", 
    # "NP" = "Near Positive", 
    # "M" = "Middle", 
    # "NN" = "Near Negative", 
    # "N" = "Negative"))



################################################################################
# Packages used
################################################################################

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RVAideMemoire)# for Wilcox paired test
library(flextable) # for creating final table
library(officer)


################################################################################
# Custom functions used
################################################################################

# Creates table with summary statistics from cjb data set
table_summary <- function (data) {
  data %>% 
    group_by(cue) %>% 
    summarise(Mean = mean(opt_score),
              SD = sd(opt_score),
              Median=median(opt_score),
              "Range (min, max)" = paste (round(min(opt_score),2), 
                                          round(max(opt_score), 2), sep=", ")) %>% 
    mutate(across(where(is.numeric), function (.x) round(.x, digits=2))) %>% 
    rename(Cue=cue)
}


# Merges two stats tables 
bindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix("", mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}



################################################################################
# Loading and preparing the data
################################################################################

jbt_all <- read_csv("cleaned_data/JBT_all.csv")

jbt_all$cue <- factor(jbt_all$cue, levels = c("P","NP","M","NN", "N"))


################################################################################
# Figure 4: Optimism scores for each cue in the JBT (touchscreen paradigm)
################################################################################

ggplot(data = subset(jbt_all, cjb_test == "ts" &
                       !is.na(jbt_all$opt_score)),
                     aes(x=cue, y=opt_score)) +
  geom_line(aes(group = cjb_test, linetype = cjb_test), 
            stat = 'summary', 
            fun = mean,
            position=position_dodge(width=0.1)) +
  stat_summary(aes(linetype=cjb_test),
               fun.data = ggpubr::mean_sd, 
               geom = "errorbar", 
               width=0.2, 
               position=position_dodge(width=0.1)) +
  stat_summary(aes(group=cjb_test, shape=cjb_test), 
               fun=mean,
               geom="point", 
               size=2, 
               position=position_dodge(width=0.1)) +
  # geom_line(aes(group=id), alpha = 0.05) + #for individual lines
  theme_classic() +
  ylab("Optimism score") +
  xlab("") +
  ylim(-1.1, 1.1) +
  scale_x_discrete(labels=c("P" = "Positive", "NP" = "Near Positive",
                            "M" = "Middle", "NN" = "Near Negative", "N" = "Negative")) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14))+
  # Manually add lines with significance indicators
  annotate("segment", x = "P", y = 1.09, xend = "NP", yend = 1.09, linewidth = 0.1) +
  annotate("segment", x = "NP", y = 1.05, xend = "M", yend = 1.05, linewidth = 0.1) +
  annotate("segment", x = "M", y = 1, xend = "NN", yend = 1, linewidth = 0.1) +
  annotate("text", x = c(1.5, 2.5, 3.5), y = c(1.1, 1.06, 1.01), 
           label = c("*", "*", "*"), hjust = 0.5, size = 7)
# ggsave("response curve_ts.png", width = 22, height = 15, units = "cm")


################################################################################
# Supplementary Figure S1: Optimism score for each cue in the tunnel JBT paradigm
################################################################################

# Tunnel
ggplot(data = subset(jbt_all, cjb_test == "tunnel" &
                       !is.na(jbt_all$opt_score)),
       aes(x=cue, y=opt_score)) +
  geom_line(aes(group = cjb_test, linetype = cjb_test), 
            stat = 'summary', 
            fun = mean,
            position=position_dodge(width=0.1)) +
  stat_summary(aes(linetype=cjb_test),
               fun.data = ggpubr::mean_sd, 
               geom = "errorbar", 
               width=0.2, 
               position=position_dodge(width=0.1)) +
  stat_summary(aes(group=cjb_test, shape=cjb_test), 
               fun=mean,
               geom="point", 
               size=2, 
               position=position_dodge(width=0.1)) +
  #geom_line(aes(group=id), alpha = 0.05) + #for individual lines
  theme_classic() +
  ylab("Optimism score") +
  xlab("") +
  ylim(-1.1, 1.1) +
  scale_x_discrete(labels=c("P" = "Positive", "NP" = "Near Positive",
                            "M" = "Middle", "NN" = "Near Negative", "N" = "Negative")) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14)
  )
# ggsave("response curve_tun.png", width = 22, height = 15, units = "cm")


################################################################################
# Supplementary Table S3: Stats for each cue in the touchscreen JBT paradigm
################################################################################

# Summary statistics table

# Preparing a summary table for TS (using custom function defined above)
raw_table_ts <-  table_summary (subset(jbt_all, 
                                       jbt_all$cjb_test == "ts"))


################################################################################
# Pairwise comparison of optimism scores for each cue and table

# Wilcoxon signed-rank test
cue_diff_wilcox_ts = wilcox.paired.multcomp( opt_score ~ cue | id, 
                                             subset(jbt_all, 
                                                    jbt_all$cjb_test == "ts"),
                                             p.method = "holm")

################################################################################
# preparing a nice table

cue_diff_wilcox_ts = as.data.frame (cue_diff_wilcox_ts$p.value)
cue_diff_wilcox_ts = cue_diff_wilcox_ts %>% 
  tibble::rownames_to_column( "Cue") %>% 
  pivot_longer(values_to = "p.value", names_to = "Cue_comapred", cols=(P:NN))

cue_diff_wilcox_ts$Cue_comapred = factor(cue_diff_wilcox_ts$Cue_comapred, 
                                         levels = c("P", "NP", "M", "NN", "N"))

cue_diff_wilcox_ts = cue_diff_wilcox_ts %>% 
  arrange(Cue_comapred) %>% 
  unite("cue",Cue_comapred, Cue, sep = "-", remove =TRUE) %>%
  na.omit()

cue_diff_wilcox_ts$p.value = ifelse(cue_diff_wilcox_ts$p.value < 0.001,"< 0.001", 
                                    cue_diff_wilcox_ts$p.value)


numeric_values <- !is.na(as.numeric(cue_diff_wilcox_ts$p.value, 2))#takes only numeric values

# round the numeric values (p that is not < 0.001)
cue_diff_wilcox_ts$p.value[numeric_values] <-
  round(as.numeric(cue_diff_wilcox_ts$p.value[numeric_values]), 3)

names(cue_diff_wilcox_ts) <- c("Cues compared","p-value")


################################################################################
# Merging two tables

bigtable_ts = bindPad(raw_table_ts,cue_diff_wilcox_ts)

#Final table aesthetics
bigtable_ts %>% 
  qflextable() %>%
  bold(part = "header") %>%
  fontsize(size = 11, part = "all") %>%
  flextable::font(fontname = "Calibri", part = "all") %>%
  vline(j = 5) %>%
  hline_bottom(border = fp_border_default(width = 0))



################################################################################
# Supplementary Table S4: Stats for each cue in the tunnel JBT paradigm
################################################################################

# Summary statistics table tunnel paradigm

raw_table_tunnel = table_summary (subset(jbt_all, 
                                         jbt_all$cjb_test == "tunnel"))


################################################################################
# Pairwise comparison of optimism scores for each cue and table

# Wilcoxon signed-rank test
cue_diff_wilcox_tunnel = wilcox.paired.multcomp( opt_score ~ cue | id,
                                                 subset(jbt_all, 
                                                        jbt_all$cjb_test == "tunnel"),
                                                 p.method = "holm")

################################################################################
# preparing a nice table

cue_diff_wilcox_tunnel = as.data.frame (cue_diff_wilcox_tunnel$p.value)
cue_diff_wilcox_tunnel = cue_diff_wilcox_tunnel %>% 
  tibble::rownames_to_column( "Cue") %>% 
  pivot_longer(values_to = "p.value", names_to = "Cue_comapred", cols=(P:NN))

cue_diff_wilcox_tunnel$Cue_comapred = factor(cue_diff_wilcox_tunnel$Cue_comapred, 
                                             levels = c("P", "NP", "M", "NN", "N"))

cue_diff_wilcox_tunnel = cue_diff_wilcox_tunnel %>% 
  arrange(Cue_comapred) %>% 
  unite("cue",Cue_comapred, Cue, sep = "-", remove =TRUE) %>%
  na.omit()

cue_diff_wilcox_tunnel$p.value = ifelse(cue_diff_wilcox_tunnel$p.value < 0.001,
                                        "< 0.001", 
                                        cue_diff_wilcox_tunnel$p.value)


numeric_values <- !is.na(as.numeric(cue_diff_wilcox_tunnel$p.value, 2)) #takes only numeric values

cue_diff_wilcox_tunnel$p.value[numeric_values] <-
  round(as.numeric(cue_diff_wilcox_tunnel$p.value[numeric_values]), 3)

names(cue_diff_wilcox_tunnel) <- c("Cues compared","p-value")


################################################################################
# Merging two tables

bigtable_tunnel = bindPad(raw_table_tunnel,cue_diff_wilcox_tunnel)

# Final table aesthetics
bigtable_tunnel %>% 
  qflextable() %>%
  bold(part = "header") %>%
  fontsize(size = 11, part = "all") %>%
  flextable::font(fontname = "Calibri", part = "all") %>%
  vline(j = 5) %>%
  hline_bottom(border = fp_border_default(width = 0))





# Exporting the table into a word file
#read_docx() %>% 
#  body_add_par("Table S4") %>% 
#  body_add_flextable(value = bigtable_tunnel) %>% 
#  print(target = "CJB_cues.docx")
