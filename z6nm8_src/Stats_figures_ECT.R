#R script used in Bracic et al (2023)

# Author: Marko Bracic (@markobracic)
# Profile: https://orcid.org/0000-0001-6528-3572
# Email: mbracic192@gmail.com


# Title: Behaviour in the ECT

################################################################################
# Description of script and Instructions
################################################################################

# Script used to analyse the behaviour in the ECT

# Output:
  # In-text values in the manuscript Results section: Safe vs predator (dangerous) chamber
  # Figure 2: Comparison of safe versus predator chamber 
  # Supplementary Figure S3: Scores in ECT compared with the training criteria

# For comparing ECT choice score vs training criteria  dataset with both ect and 
  #  jbt data was used (md_id: on the individual level)



################################################################################
# Packages used
################################################################################

library(tidyverse)
library(ggpubr) #for ggarange function
library(ggbeeswarm) #spreading data points in a nice way



################################################################################
# Loading the data
################################################################################

# File with combined data (ECT choice score vs traning criteria) 
md_id <- read_csv("cleaned_data/JBTxECT_id.csv")

# File with the all behaviors from ECT on the id level
ect_chambers <- read_csv("cleaned_data/ECT_chambers.csv")



################################################################################
# Statistical tests: ECT Choice score vs Training criteria 
################################################################################

# Did mice reduce their choice score in the test compared to 80% criteria?
  #   using  JBTxECT_id.csv dataset

wilcox.test(md_id$choice_score[md_id$rat_bedding != "mixed"],
            mu = 0.6) # 80% big = NCT score 0.6

#Result
  # yes



################################################################################
# Supplementary Figure S3: Comparison between the safe versus predator chamber
################################################################################
ggplot(subset(md_id,rat_bedding != "mixed"),
       aes(x=factor(0), y=choice_score)) +
  geom_boxplot( width = 0.2) +
  geom_beeswarm(alpha = 0.1) +
  ylab("Choice score") +
  ylim(-1.1,1) +
  geom_hline(yintercept=0.6, linetype="dashed", colour = "red4") +
  annotate(
    "text",
    label = "Training criteria",
    colour = "red4",
    x = 0.55,
    y = 0.66,
    size = 4
  ) +
  annotate("text", x = 1, y = 1, vjust = "top", label = "*", size = 7) +
  theme_classic(base_size = 14) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#ggsave("ect_response.png", width = 15, height = 12, units = "cm")




################################################################################
# Statistical tests: Comparison between the safe versus predator chamber
################################################################################

# Analysis of four behavioural measures using wilcox test (ect_chambers dataset)

# Entries into the chamber per trial
wilcox.test(
  data = ect_chambers,
  entries ~ chamber,
  paired = TRUE,
  alternative = "two.sided"
)


################################################################################
# Time spent in the chamber per trial
wilcox.test(
  data = ect_chambers %>%
    filter(!is.na(duration)) %>% #if duration not scored (=NA), remove row
    subset(., id %in% id[duplicated(id)]), #if one chamber was NA also remove the other chamber
  duration ~ chamber,
  paired = TRUE,
  alternative = "two.sided"
)
#if raw duration compared than the difference is significant


################################################################################
# Digging in the chamber per trial
wilcox.test(
  data = ect_chambers %>%
    filter(!is.na(digging)) %>% #remove row if did not enter chamber
    subset(., id %in% id[duplicated(id)]), #if one chamber was NA also remove the other chamber
  digging ~ chamber,
  paired = TRUE,
  alternative = "two.sided"
)


################################################################################
# Retractions from the chamber per trial
wilcox.test(
  data = ect_chambers %>%
    filter(!is.na(retraction)) %>% #remove if did not enter tunnel
    subset(., id %in% id[duplicated(id)]), #if one chamber was NA also remove the other chamber
  retraction ~ chamber,
  paired = TRUE,
  alternative = "two.sided"
)



################################################################################
# Figure 2: Comparison between the safe versus predator chamber
################################################################################

# difference in retraction between predator and safe (dangerous chamber
  # remove ind that did not enter dangerous tunnel
retraction = ggplot(
  ect_chambers %>%
    filter(!is.na(retraction)) %>% #remove if did not enter tunnel
    subset(., id %in% id[duplicated(id)]), #remove for other chamber
  aes(x = factor(chamber, levels = c("safe", "dangerous")), 
      y = retraction)
) +
  geom_boxplot(aes(group = chamber),
               width = 0.2,
               outlier.shape = NA,
               lwd=0.6 ) +
  geom_beeswarm(alpha = 0.1) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    vjust = "top",
    label = "*",
    size = 7
  ) +
  xlab("Chamber") +
  ylab("Retractions per entry") +
  scale_x_discrete(labels = c("dangerous" = "predator",
                              "safe" = "safe")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14))
#ggsave("retraction.png", width = 15, height = 12, units = "cm")


#ggplot(hesi_dif, aes(x=chamber, y=count)) +
# stat_summary(geom = "bar", fun = mean)+
# geom_jitter()


################################################################################
# Digging
  # remove ind that did not enter predator chamber

digging = ggplot(
  ect_chambers %>%
    filter(!is.na(digging)) %>% #remove if did not enter chamber
    subset(., id %in% id[duplicated(id)]), #remove for other chamber
  aes(x = factor(chamber, levels = c("safe", "dangerous")),
      y = digging)) +
  geom_boxplot(
    aes(group = chamber),
    width = 0.2,
    outlier.shape = NA,
    lwd = 0.6
  ) +
  geom_beeswarm(alpha = 0.1) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    vjust = "top",
    label = "*",
    size = 7
  ) +
  xlab("Chamber") +
  ylab("Digging per entry") +
  scale_x_discrete(
    labels = c("dangerous" = "predator",
               "safe" = "safe"
    )) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14))


################################################################################
# Duration

duration = ggplot(
  ect_chambers %>%
    filter(!is.na(duration)) %>% #remove if did not enter chamber
    subset(., id %in% id[duplicated(id)]), #remove for other chamber
  aes(x = factor(chamber, levels = c("safe", "dangerous")),
      y = duration)) +
  geom_boxplot(
    aes(group = chamber),
    width = 0.2,
    outlier.shape = NA,
    lwd = 0.6
  ) +
  geom_beeswarm(alpha = 0.1) +
  xlab("Chamber") +
  ylab("Duration per entry (s)") +
  scale_x_discrete(labels = c("dangerous" = "predator",
                              "safe" = "safe")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14))


################################################################################
# Entries

entries = ggplot(ect_chambers, 
                  aes(x = factor(chamber, levels = c("safe", 
                                                     "dangerous")), 
                      y = entries)) +
   geom_boxplot(aes(group = chamber),
                width = 0.2,
                outlier.shape = NA,
                lwd=0.6 ) +
   geom_beeswarm(alpha = 0.1) +
   annotate(
     "text",
     x = 1.5,
     y = Inf,
     vjust = "top",
     label = "*",
     size = 7
   ) +
   xlab("Chamber") +
   ylab("Entries per trial") +
   scale_x_discrete(labels = c("dangerous" = "predator",
                               "safe" = "safe")) +
   theme_classic() +
   theme(axis.title = element_text(size = 14),
         axis.text  = element_text(size = 14))


################################################################################
# Merging the figures 

ggarrange(
  entries + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ),
  digging + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ),
  retraction,
  duration,
  hjust = -2,
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  label.x	= -0.05
)
#ggsave("behaviours_man.png", width = 19, height = 16, units = "cm")
