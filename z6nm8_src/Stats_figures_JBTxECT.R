#R script used in Bracic et al (2023)

# Author: Marko Bracic (@markobracic)
# Profile: https://orcid.org/0000-0001-6528-3572
# Email: mbracic192@gmail.com


# Title: Individual variation in ECT and relationship between JBT and ECT

################################################################################
# Description of script and Instructions
################################################################################

# Script used to 
# 1) analyse the repeatability of the choice in ECT
# 2) analyse the relationship between JBT and ECT

# Output:
  # In-text values in the manuscript Results section: individual differences
  # Figure 3:  Individual differences and repeatability in ECT 
  # In-text values in the manuscript Results section: JBT vs. ECT
  # Figure 5:  Relationship between JBT and ECT
  # Supplementary Figure S3: Relationship between JBT and ECT for all cues



################################################################################
# Packages used
################################################################################

library(tidyverse)
library(lme4) # mixed models
library(lmerTest) # gives p-values from lme4 models
library(rptR) # extraction of variance components for multiple grouping levels for estimating repeatabilities
library(ggpubr) #for ggarange function
library(ggeffects) # estimate CI from model and use in ggplot
library(RColorBrewer) # colors for figure on individual differences
library(ggdist) # half-eye plot
library(report) # stats table



################################################################################
# Custom functions used
################################################################################

# Function to round p-values
pvalr <-
  function(pvals,
           sig.limit = .001,
           digits = 3,
           digits_big = 3,
           html = FALSE) {
    roundr <- function(x, digits = 1) {
      res <- sprintf(paste0("%.", digits, "f"), x)
      zzz <- paste0("0.", paste(rep("0", digits), collapse = ""))
      res[res == paste0("-", zzz)] <- zzz
      res
    }
    
    sapply(pvals, function(x, sig.limit) {
      if (x < sig.limit) {
        if (html) {
          return(sprintf("&lt; %s", format(sig.limit)))
        } else {
          return(sprintf("< %s", format(sig.limit)))
        }
      }
      if (x > .1) {
        return(roundr(x, digits = digits_big))
      } else {
        return(roundr(x, digits = digits))
      }
    }, sig.limit = sig.limit)
  }

# code adapted based on: https://stackoverflow.com/a/23018806



# Function to change p-values in the stat table with the ones from the model
update_table <- function(mymod) {
  modsum <- summary(mymod)
  p <- modsum$coefficients[, "Pr(>|t|)"]
  p <- as.vector(p)
  p <- pvalr(p, sig.limit = .001, digits = 3)
  p
}



################################################################################
# Loading the data
################################################################################

# Loading the data from repository
md_id <- read_csv("cleaned_data/JBTxECT_id.csv") #used for figures
md_day <- read_csv("cleaned_data/JBTxECT_day.csv")


# Or run the script for data preparation (not provided)
#source("01_all_data_prep.R")

################################################################################
# Preparing the data for statistical analysis
################################################################################

# Variables into Factors and setting reference level

# Individual identification
md_day$id <- factor(md_day$id)

# Two different sources of rat bedding
md_day$rat_bedding <- factor(md_day$rat_bedding)
md_day$rat_bedding <- relevel(md_day$rat_bedding, ref = "our")
# internal bedding as reference level


################################################################################
# Centering the factors -> increases model interpretation
  # see Schielzeth 2010 for details
  # https://doi.org/10.1111/j.2041-210X.2010.00012.x

# Rat bedding
md_day$rat_beddingC <- as.numeric(md_day$rat_bedding == "external") - 0.5
# external=1 and our=0, making our reference level and - 0.5 centers it (external=0.5, our=-0.5)

# Testing session (day)
md_day$ect_sessionC <- md_day$ect_session - 2



################################################################################
# Model for calculating REPEATABITY of choice score from ECT
################################################################################
# Includes individual tested in JBT test but also ind not tested in JBT tests


# Underling model
mod_repeatablity <- lmer(
  choice_score ~ rat_beddingC + ect_sessionC + (1 | id),
  data = subset(md_day, md_day$rat_bedding != "mixed")
) # centered to test

# Extracting model estimates
model_estimates_repeatablity <- summary(mod_repeatablity)$coefficient %>%
  as.data.frame() %>%
  rename("p" = "Pr(>|t|)") %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3)) %>%
  mutate(across(where(is.numeric), round, digits = 2))


################################################################################

# calculating R

set.seed(1902); # makes it reproducible, change the number for new bootstrapping
R_ect <- rptGaussian(
  choice_score ~ rat_beddingC + ect_sessionC + (1 | id),
  grname = "id",
  data = subset(md_day, md_day$rat_bedding != "mixed"),
  nboot = 1000
)



################################################################################
# FINAL MODELS (3 models from each of the JBT paradigms for 3 ambiguous cues)
################################################################################

# TOUCHSCEEN JBT

# Main model with middle ambiguous cue (Figure 4)
modM_ts <- lmer(
  choice_score ~ M + rat_beddingC + ect_sessionC + (1 | id),
  data = subset(
    md_day,
    cjb_type == "ts" & md_day$rat_bedding != "mixed" & !is.na(md_day$M)
  )
)

# Additional model for near negative cue (Supplementary Figure S3)
modNN_ts <- lmer(
  choice_score ~ NN + rat_beddingC + ect_sessionC + (1 | id),
  data = subset(
    md_day,
    cjb_type == "ts" & md_day$rat_bedding != "mixed" & !is.na(md_day$M)
  )
)
# to remove the influential data point (id = 3528), add this subset:
  # data = subset(
  #   md_day,
  #   cjb_type == "ts" & md_day$rat_bedding != "mixed" & !is.na(md_day$M) & id != 3528
# to correct for multiple testing, use this code
# p.adjust(summary(modNN_ts)$coefficients["NN", "Pr(>|t|)"], method = "holm", n = 3)

# Additional model for near positive cue (Supplementary Figure S3)
modNP_ts <- lmer(
  choice_score ~ NP + rat_beddingC + ect_sessionC + (1 | id),
  data = subset(
    md_day,
    cjb_type == "ts" & md_day$rat_bedding != "mixed" & !is.na(md_day$M)
  )
)



################################################################################
# TUNNEL JBT

# Main model with middle ambiguous cue (Figure 4)
modM_tunnel <- lmer(
  choice_score ~ M + ect_sessionC + (1 | id),
  data = subset(md_day, cjb_type == "tunnel" &
    md_day$rat_bedding != "mixed")
)

# Additional model for near negative cue (Supplementary Figure S3)
modNN_tunnel <- lmer(
  choice_score ~ NN + ect_sessionC + (1 | id),
  data = subset(md_day, cjb_type == "tunnel" &
    md_day$rat_bedding != "mixed")
)

# Additional model for near positive cue (Supplementary Figure S3)
modNP_tunnel <- lmer(
  choice_score ~ NP + ect_sessionC + (1 | id),
  data = subset(md_day, cjb_type == "tunnel" &
    md_day$rat_bedding != "mixed")
)




################################################################################
# Figure 3  Individual differences in ECT 
################################################################################

# Figure 3A) REPEATABITY estimate

# extracting bootstrapping data from rpt model
halfeye <- data.frame(R_ect$R_boot$id) %>%
  pivot_longer(., everything(), names_to = "Cue", values_to = "R_boot")

# Half eye plot with CI

  plot_halfeye <- ggplot(halfeye, aes(x = Cue, y = R_boot)) +
    # geom_point(
    #   size = 0.5,
    #   alpha = .05,
    #   position = position_jitter(seed = 1, width = .05)
    # ) +
    stat_halfeye(
      aes(slab_alpha = 0),
      show_point = F,
      adjust = 1,
      width = .2,
      .width = 1, # controls how long the line is (removes it)
      # justification = -.3, # controls position on the y axis
      ymin = R_ect$CI_emp$`2.5%`,
      ymax = R_ect$CI_emp$`97.5%`,
      interval_colour = "red",
      slab_type = "pdf",
    ) + 
    geom_point(aes(x= 1,
                   y = as.numeric(R_ect$R)),
               #colour = "red",
               size = 2.4
    ) +
    ylab("Repeatablity") +
    xlab("") +
    theme_classic() +
    # scale_x_discrete(expand = expansion(mult = c(0.2, 0))) +
    theme(
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 16),
      text = element_text(size = 14),
      legend.position = "none"
    ) +
    scale_x_discrete(expand = c(0, 0.2)) +
    scale_y_continuous(limits = c(0,1), 
                       expand = c(0,0))

################################################################################

# Figure 3B) Individual choice scores
  plot_ect_variation <- ggplot(md_day, aes(x = as.factor(ect_session), y = choice_score)) +
    geom_line(aes(group = as.factor(id), colour = as.factor(id)),
              linewidth = 0.6
    ) +
    ylab("Choice score") +
    xlab("Test session") +
    #scale_colour_grey(guide = "none") +
    scale_fill_manual(values =  colorRampPalette(brewer.pal(9, "Set3"), space = "rgb")(length(unique(md_day$id)))) +                                                              
    theme_classic() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none"
    ) +
    ylim(-1, 1) +
    scale_x_discrete(expand = c(0, 0.02))
# ggsave("Consistancy.png", width = 15, height = 12, units = "cm")

  
################################################################################

# Final figure 3)

ggarrange(
  plot_halfeye,
  plot_ect_variation,
  ncol = 2,
  nrow = 1,
  labels = "AUTO",
  widths = c(1, 2)
)
# ggsave("Ind_var_and_Rpt.png", width = 15, height = 12, units = "cm")



################################################################################
# Figure 5:  Relationship between JBT and ECT (only middle cue)
################################################################################

# Touchscreen paradigm only middle cue  

# Extracting and formatting p-values from models to be used in the figures
M_ts_table <- summary(modM_ts)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))
  # rounding bigger p-values to 2 digits using custom function  
  
# Preparing figure p-value annotation
plot_text_M_ts <- paste(
  "Optimism score: p = ",
  M_ts_table$p[M_ts_table$effect == "M"],
  "\nRat bedding: p = ",
  M_ts_table$p[M_ts_table$effect == "rat_beddingC"],
  "\nSession: p ",
  M_ts_table$p[M_ts_table$effect == "ect_sessionC"],
  sep = ""
)
  
# Figure 
set.seed(3451); # makes jitter reproducible, 
M_only_ts <- ggpredict(modM_ts, terms = "M") %>%
  #from ggeffects package: predicts CI
  ggplot(aes(x, predicted)) +
  geom_ribbon(
    aes(x, predicted,
        ymin = conf.low,
        ymax = conf.high),
    fill = "gray",
    alpha = 0.15
  ) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, cjb_type == "ts" &
                    rat_bedding != "mixed" & !is.na(md_id$M)),
    aes(y = choice_score, x = M),
    size = 2,
    colour = "gray",
    alpha = 0.8,
    position=position_jitter(h=0, w=0.01)
  ) +
  annotate( 
    "text",
    label = plot_text_M_ts,
    x = -0.2,
    y = 1,
    size = 3
    ) +
  ylab("Choice score") +
  xlab("Optimism score") +
  ylim(-1.0, 1.0) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")
# ggsave("JBT_M_ts vs ECT.png", width = 15, height = 12, units = "cm")



################################################################################
# Supplementary Table 5,6,7: Output from Statistical models 
################################################################################

# Output from summary() function (used for estimates, t-statistic, and p-values)
  
# Touchscreen JBT
model_estimates_ts <- rbind(
  summary(modNP_ts)$coefficients,
  summary(modM_ts)$coefficients,
  summary(modNN_ts)$coefficients
  ) %>%
  as.data.frame() %>%
  rename("p" = "Pr(>|t|)") %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3)) %>%
  mutate(across(where(is.numeric), round, digits = 2))

# Tunnel JBT
model_estimates_tun <- rbind(
  summary(modNP_tunnel)$coefficients,
  summary(modM_tunnel)$coefficients,
  summary(modNN_tunnel)$coefficients
) %>%
  as.data.frame() %>%
  rename("p" = "Pr(>|t|)") %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3)) %>%
  mutate(across(where(is.numeric), round, digits = 2))

# Using Report package (used to estimate confidence intervals and R squered)

# Touchscreen JBT
table_ts_NP <- report_table(modNP_ts)
table_ts_NP$p[1:4] <- update_table(modNP_ts)

table_ts_M <- report_table(modM_ts)
table_ts_M$p[1:4] <- update_table(modM_ts)

table_ts_NN <- report_table(modNN_ts)
table_ts_NN$p[1:4] <- update_table(modNN_ts)

# Tunnel JBT
table_tun_NP <- report_table(modNP_tunnel)
table_tun_NP$p[1:3] <- update_table(modNP_tunnel)

table_tun_M <- report_table(modM_tunnel)
table_tun_M$p[1:3] <- update_table(modM_tunnel)

table_tun_NN <- report_table(modNN_tunnel)
table_tun_NN$p[1:3] <- update_table(modNN_tunnel)

# Repeatability
table_repeatablity <- report_table(mod_repeatablity)
table_repeatablity$p[1:3] <- update_table(mod_repeatablity)


# CI computed using a Wald method t-distribution with Satterthwaite DF
# See https://cran.revolutionanalytics.com/web/packages/parameters/parameters.pdf

# To exporting the tables, display() function from insight package was used
  # Table then modified manually
  # Additional data removed, DF added, and table restructured



################################################################################
# Supplementary Figure S2: Relationship between JBT and ECT for all cues
################################################################################

# Extracting and formatting p-values from models to be used in the figures

# TOUCHSCEEN JBT (NP and NN ambiguous cue)

plot_text_NP_ts <- summary(modNP_ts)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))
  # rounding bigger p-values to 2 digits using custom function  

plot_text_NN_ts <- summary(modNN_ts)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))



################################################################################
# Tunnel JBT (same as for TOUCHSCEEN)

#Extracting and formatting p-values (for each ambiguous cue)

M_tun_table <- summary(modM_tunnel)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))
  # rounding bigger p-values to 2 digits using custom function

plot_text_NP_tun <- summary(modNP_tunnel)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))

plot_text_NN_tun <- summary(modNN_tunnel)$coefficients %>%
  subset(select = "Pr(>|t|)") %>%
  data.frame() %>%
  rownames_to_column(., var = "effect") %>%
  subset(., .$effect != "(Intercept)") %>%
  rename("p" = "Pr...t..") %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  mutate(p = pvalr(.$p, sig.limit = .001, digits = 3, digits_big = 2))


################################################################################
# Preparing figure p-value annotation (for each ambiguous cue)

# TOUCHSCEEN JBT

plot_text_NP_ts <- paste(
  "Optimism score: p = ",
  plot_text_NP_ts$p[plot_text_NP_ts$effect == "NP"],
  "\nRat bedding: p = ",
  plot_text_NP_ts$p[plot_text_NP_ts$effect == "rat_beddingC"],
  "\nSession: p ",
  plot_text_NP_ts$p[plot_text_NP_ts$effect == "ect_sessionC"],
  sep = ""
)

plot_text_NN_ts <- paste(
  "Optimism score: p = ",
  plot_text_NN_ts$p[plot_text_NN_ts$effect == "NN"],
  "\nRat bedding: p = ",
  plot_text_NN_ts$p[plot_text_NN_ts$effect == "rat_beddingC"],
  "\nSession: p ",
  plot_text_NN_ts$p[plot_text_NN_ts$effect == "ect_sessionC"],
  sep = ""
)

# Tunnel JBT (same procedure as  for TOUCHSCEEN)
plot_text_M_tun <- paste("Optimism score: p = ",
                         M_tun_table$p[M_tun_table$effect == "M"],
                         "\nSession: p = ",
                         M_tun_table$p[M_tun_table$effect == "ect_sessionC"],
                         sep = "")


plot_text_NP_tun <- paste("Optimism score: p = ",
                          plot_text_NP_tun$p[plot_text_NP_tun$effect == "NP"],
                          "\nSession: p = ",
                          plot_text_NP_tun$p[plot_text_NP_tun$effect == "ect_sessionC"],
                          sep = "")

plot_text_NN_tun <- paste("Optimism score: p = ",
                          plot_text_NN_tun$p[plot_text_NN_tun$effect == "NN"],
                          "\nSession: p = ",
                          plot_text_NN_tun$p[plot_text_NN_tun$effect == "ect_sessionC"],
                          sep = "")

################################################################################

# Figures TOUCHSCEEN JBT all ambiguous cues

NP_ts <- ggpredict(modNP_ts, terms = "NP") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, cjb_type == "ts" &
                    rat_bedding != "mixed" & !is.na(md_id$M)),
    aes(y = choice_score, x = NP),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_NP_ts,
    x = 0.5,
    y = Inf,
    size = 2
  ) +
  ylab("") +
  xlab("") +
  ylim(-1.1, 1.1) +
  ggtitle("Near positive") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")


M_ts <- ggpredict(modM_ts, terms = "M") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, cjb_type == "ts" &
                    rat_bedding != "mixed" & !is.na(md_id$M)),
    aes(y = choice_score, x = M),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_M_ts,
    x = -0.2,
    y = Inf,
    size = 2
  ) +
  ylab("Choice score") +
  xlab("") +
  ylim(-1.1, 1.1) +
  ggtitle("Middle") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")


NN_ts <- ggpredict(modNN_ts, terms = "NN") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, cjb_type == "ts" &
                    rat_bedding != "mixed" & !is.na(md_id$M))
    ,
    aes(y = choice_score, x = NN),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_NN_ts,
    x = -0.5,
    y = Inf,
    size = 2
  ) +
  ylab("") +
  xlab("Optimism score") +
  ylim(-1.1, 1.1) +
  ggtitle("Near negative") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")


# Figure TUNNEL JBT all ambiguous cues

M_tun <- ggpredict(modM_tunnel, terms = "M") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, md_id$M != "NA" & cjb_type == "tunnel"),
    aes(y = choice_score, x = M),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_M_tun,
    x = 0.2,
    y = Inf,
    size = 2
  ) +
  ylab("Choice score") +
  xlab("") +
  ylim(-1.1, 1.1) +
  ggtitle("Middle") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")


NN_tun <- ggpredict(modNN_tunnel, terms = "NN") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, md_id$M != "NA" & cjb_type == "tunnel"),
    aes(y = choice_score, x = NN),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_NN_tun,
    x = -0.2,
    y = Inf,
    size = 2
  ) +
  ylab("") +
  xlab("Optimism score") +
  ylim(-1.1, 1.1) +
  ggtitle("Near negative") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")


NP_tun <- ggpredict(modNP_tunnel, terms = "NP") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = "gray",
              alpha = 0.2) +
  geom_line() + #intercept and slope from the model
  geom_point(
    data = subset(md_id, md_id$M != "NA" & cjb_type == "tunnel"),
    aes(y = choice_score, x = NP),
    size = 1,
    colour = "gray"
  ) +
  annotate(
    "text",
    label = plot_text_NP_tun,
    x = 0.6,
    y = Inf,
    size = 2
  ) +
  ylab("") +
  xlab("") +
  ylim(-1.1, 1.1) +
  ggtitle("Near positive") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = -0.05)) +
  coord_cartesian(clip = "off")

################################################################################
# Final figure with both Touchscreen and tunnel JBT used in the supplements

ggarrange(NP_ts,
          NP_tun,
          M_ts,
          M_tun,
          NN_ts,
          NN_tun,
          labels = c("Touchscreen", "Tunnel"),
          vjust = 0,
          hjust = 0,
          nrow = 3, ncol = 2
) +
  theme(plot.margin = margin(20, 5.5, 5.5, 5.5))
# ggsave("optimism vs risk_merged_man.png", width = 27, height = 13, units = "cm")
