########################################
# Script Name: Drug_dose_assays_script
# Description: Analysis of comprehensive results from drug dose assays
# Author: Megan Sandoval-Powers
# Last Updated: 8-2026
# R Version: 4.5.1
########################################

# =================== R PACKAGES =============================================
library(dplyr) #mutate/group_by/_summarize
library(ggplot2) #plotting
library(tidyr) #pivot_longer
library(readr) #read_csv
library(cowplot) #theme_cowplot
library(ggpattern) #geom_col_pattern, the hatched bars in plot


# =================== SET UP DATA FOR PLOTTING ============================================= 

#WHAT THIS DOES
#Reads the 24-well dose-response summary, averages replicates at 24 h and 48 h for each drug x dose, and draws one faceted bar chart with hatching
#to distinguish the two timepoints.

  #### Read in data ####

    #Read in summary dataframe with OD600 values
    data <- read_csv("Input_files/OD600_sumdata_24well_drugdose.csv")

  #### Build the dataframe #### 
    
    #Split drug from rep to isolate Drug in single column
      data <- data %>%
        mutate(Drug = sub("_rep.*", "", Name))
    
    #Get mean OD600 values and SD at 24 hr and 48 hr per drug dose
      summary_data <- data %>%
        group_by(Drug, Dose) %>%
        summarise(
          mean_24hr = mean(Relative_OD_24hr),
          sd_24hr   = sd(Relative_OD_24hr),
          mean_48hr = mean(Relative_OD_48hr),
          sd_48hr   = sd(Relative_OD_48hr),
          .groups = "drop")
      
    #Shift to long format for plotting
      summary_long <- summary_data %>%
        pivot_longer(
          cols = -c(Drug, Dose),
          names_to = c(".value", "Time"),
          names_pattern = "(mean|sd)_(.*)")
      
    #Change the timepoint names
      summary_long$Time <- recode(summary_long$Time,
                                  "24hr" = "24 hr",
                                  "48hr" = "48 hr")
      
    #Change the order of facets that will be plotted
      summary_long$Drug <- factor(summary_long$Drug,
                                levels = c("CAS", "CLO", "CASCLO"))
    
    #Change dose to be a factor
      summary_long$Dose <- factor(summary_long$Dose,
                                  levels = unique(summary_long$Dose))
      
  #### Plot mean relative OD600 ####

    #Plot with facets by drug. 
      p1 <- ggplot(summary_long,aes(x = Dose,y = mean,fill = Drug,pattern = Time)) +
        geom_bar_pattern(stat = "identity",position = position_dodge(width = 0.9),
          width = 0.9, colour = "black",pattern_fill = "black",pattern_density = 0.2,pattern_spacing = 0.03) +
        geom_errorbar(aes(ymin = mean - sd,ymax = mean + sd),position = position_dodge(width = 0.9),width = 0.25,
          linewidth = 0.3) +
        scale_y_continuous(limits=c(0,100),breaks=c(seq(0,100,by=25)), expand=c(0,0))+
        facet_wrap(~ Drug, scales = "free_x") +
        scale_fill_manual(values = c("CAS" = "gold","CLO" = "hotpink","CASCLO" = "orangered")) +
        scale_pattern_manual(values = c("24 hr" = "none","48 hr" = "stripe")) +
        theme_cowplot() +
        theme(
          strip.background = element_rect(fill = "grey85", colour = "black"),
          strip.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle=45,hjust=1,vjust=1)) +
        guides(fill=guide_legend(override.aes = list(pattern="none")),
               pattern = guide_legend(override.aes = list(fill="white"))) +
        labs(x = expression("Dose (" * mu * "M)"),y = expression("Mean relative OD"[600] ~ "(% of control)"))
      p1

    #Export plot
      #ggsave("Figures/Drug_dose_response_assays_sumfig.png", p1, width=8,height=5,units="in",dpi=300)
      