########################################
# Script Name: Growth_phenotyping_script
# Description: Plotting and analysis of growth data for phenotyping evolved pops
# Author: Megan Sandoval-Powers
# Last Updated: 8-2026
# R Version: 4.5.1
########################################


# =================== R PACKAGES =============================================
library(dplyr) #for the select function
library(cowplot) #Make clean graphs
library(patchwork) #Combine plots quickly
library(reshape2) #To melt things
library(stringr) #for str_split_fixed, parsing sample IDs
library(tibble) #used to hand-build some results tables
library(emmeans) #Pairwise comparisons of lm results with post hoc corrections. 
library(growthcurver) #for running the growthcurver package to get DT and CCs
library(ggh4x) #for plotting facet wrap
library(ggplot2) #all plotting
library(svglite) #for exporting figures
library(lme4)   #for lmer mixed effects model
library(lmerTest) #adds p values to lmer output

set.seed(1)           #fixes the random jitter so points land identically each run

   
# =================== SET UP LABELS, PLATE MAP, AND GLOBAL PLOTTING ELEMENTS =============================================  
      
  
    #### Plate well identifiers ####
      
      #Each assay used a different subset of the 96-well plate. These vectors say
      #which wells actually held samples, in the order the sample-name CSV lists them.
      
      #Build well names by crossing row letters with column numbers.
      #CAUTION: outer() varies its FIRST argument fastest, so rows cycle within each
      #column -> B2, C2, D2, E2, F2, G2, B3, C3, ... That column-major order is what
      #your sample-name CSVs assume.
      
      wells_grid <- function(rows, cols) as.vector(outer(rows, cols, paste0))
      
      WELLS_STD <- c(                                    #CAS and CLO assays:
        wells_grid(c("B","C","D","E","F","G"), 2:7),     #full block, rows B-G x cols 2-7
        wells_grid(c("B","C","D"), 8:10))                #partial block, rows B-D x cols 8-10
      
      WELLS_10COL <- wells_grid(c("B","C","D","E","F","G"), 2:11)  #CASCLO assay: rows B-G, cols 2-11
      WELLS_FULL  <- wells_grid(LETTERS[1:8], 1:12)                #cross-tolerance: whole plate
      
      
    #### Treatment labels ####  
      
      #The 5th field of each sample ID is a terse treatment code. These maps turn
      #those codes into the labels that appear on the figures. R
    
      #names  = raw codes as they appear in the sample IDs
      #values = display labels used on axes, strips, and in the models
      
      MAP_CAS      <- c(ND = "Control", Low = "Low CAS", High = "High CAS") #Change ND to Control, Low to Low CAS, High to High CAS and so on
      MAP_CLO      <- c(ND = "Control", Low = "Low CLO", High = "High CLO")
      MAP_CASCLO   <- c(ND = "Control", HighCAS = "High CAS",HighCLO = "High CLO", CASCLO = "CASCLO")
      MAP_CROSSTOL <- c(ND = "Control", HighCAS = "High CAS",CASCLO = "CASCLO", AmphoB = "Ampho B")     
      
      
    #### Color palette ####  
      
      #One palette for every population across every figure, so a given population is the same color everywhere and the legends stay consistent.
      POP_COLS <- c("ANC"      = "grey60",         
                    "Low CAS"  = "gold",           
                    "High CAS" = "darkgoldenrod",  
                    "Low CLO"  = "hotpink",        
                    "High CLO" = "firebrick4",     
                    "CASCLO"   = "orangered")      
      
      
    #### Fixed axis settings #### 
      
      #Fixed (not free) axes so panels are comparable by eye across the figure.
      DT_YLIM <- c(0, 6); DT_YBRK <- c(0, 2, 4, 6)   #DT plots will run from 0 to 6 hours; tick marks will appear at 0,2,4,6
      K_YLIM  <- c(0, 2); K_YBRK  <- c(0, 1, 2)      #K will run from 0 to 2; tickmarks at 0,1,2
      KLAB    <- expression(italic(K))               #italic K for the y-axis title
      
      
    #### Assay-environment x-axis labels ####
      
      #Left-to-right order of the assay environments on each panel's x-axis.
      #Control always first so the drug panels read as departures from it.
      
      CAS_ENV      <- c("Control", "Low CAS", "High CAS") 
      CLO_ENV      <- c("Control", "Low CLO", "High CLO")
      CASCLO_ENV   <- c("Control", "High CAS", "High CLO", "CASCLO")
      CROSSTOL_ENV <- c("Control", "High CAS", "CASCLO", "Ampho B")
      
# =================== HELPER FUNCTIONS =============================================
      
    #### Process raw OD files function ####
      
      #process_OD600_data()          RUN ONCE, THEN LEAVE ALONE ***
      #Turns a raw plate-reader CSV export into the tab-delimited table the rest of the pipeline reads. Strips the instrument metadata header, keeps the
      #per-well Mean rows, transposes so samples become columns, and adds a time column.
      
      #input_file    raw CSV straight off the plate reader
      #output_file   where to write the formatted table
      #header_lines  rows of instrument metadata to discard before the data
      #interval      minutes between reads (converted to hours later, in prep_plate)
      #step          rows per sample block in the export, used to locate the sample-ID rows
      
      process_OD600_data <- function(input_file, #path/name of the raw plate reader csv
                                     output_file, #path name where processed data will be saved
                                     header_lines = 90,  #number of metadata rows to remove
                                     interval = 30, #time between plate reader measurments, in minutes
                                     step = 12) #number of rows between successive sample ID rows
        { 
        
        dat <- read.csv(file = input_file, header = TRUE) #read the raw csv file, treat first row as column names
        
        dat2 <- dat[-c(1:header_lines), ]          #drop the instrument metadata block
        rownames(dat2) <- NULL                     #reset row numbers after the deletion
        colnames(dat2)[1] <- "SampleIDs"          #rename first column to sample IDs
        
        howmany <- length(which(dat2$SampleIDs == "Mean"))   #one "Mean" row per sample
        nums <- seq(2, by = step, length = howmany)          #rows holding the sample IDs
        data_labels <- dat2[nums, 1]                         #extract the sample IDs from the first column at those rows
        
        dat3 <- dat2 %>% filter(SampleIDs == "Mean")   #keep only the averaged reads
        dat3$SampleIDs <- data_labels                  #attach the IDs to them
        
        dat4 <- as.data.frame(t(dat3))             #transpose: samples become columns
        colnames(dat4) <- as.character(unlist(dat4[1, ]))   #unlist(): the first row is a data-frame row, not a vector
        dat4 <- dat4[-1, ]                         #that row is now the header, so drop it
        rownames(dat4) <- NULL
        
        time <- seq(0, by = interval, length = nrow(dat4))   #create the time sequence. starts at 0 minutes. 
        dat5 <- cbind(time, dat4)                          #add the time vector as the first column of the dataset             
        
        dat5[] <- lapply(dat5, function(x) as.numeric(as.character(x))) #convert every column to numeric
        
        write.table(dat5, file = output_file, quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE) 
        
        dat5
      }
      
    #### Relabel function ####
      
      #relabel()
      #Recode a character vector through a named lookup map.
      
      #x character vector (or factor) to recode map named vector; names = old values, values = new values
      #Values not found in the map pass through UNCHANGED. That makes the function safe to call twice on the same column
      
      relabel <- function(x, map) {
        x <- as.character(x)              #factors would index the map by level number, not label
        y <- map[x]                       #look up every value in x inside the named map. Example map["ND"] -> "Control"
        y[is.na(y)] <- x[is.na(y)]        #restore the originals where there was no match
        unname(y)                         #drop the names map[] attached, leaving a plain vector
      }
      
      
    #### Split the sample ID names function ####
      
      #split_ids()
      #Split the 5-part underscore-delimited sample ID into separate columns.

      #df   data frame containing the ID column
      #col  name of that column ("Sample" for long data, "sample" for growthcurver)
      #ID format:  Population_Strength_Replicate_Timepoint_Treatment
      #Returns df with five new columns plus ID_Treatment, the grouping key used to average replicates for the growth curves.
      
      split_ids <- function(df, col = "Sample") {
        p <- str_split_fixed(as.character(df[[col]]), #extract the ID column and make sure it is character
                             "_", #split each ID wherever there is an underscore
                             5)   #split into exactly 5 columns
        df %>% mutate(             #Add the split ID information as new columns
          Population = p[, 1],            #ANC, CAS, CLO, CASCLO
          Strength   = p[, 2],            #selection strength, e.g. S01 / S03
          Replicate  = p[, 3],            #rep01 ... rep12
          Timepoint  = p[, 4],            #evolution timepoint
          Treatment  = p[, 5],            #assay environment, raw code (recoded later)
          ID_Treatment = paste(Population, Strength, Timepoint, Treatment, sep = "_")) #combine fields except replicate into grouping ID
      }
      
      
    #### Organize the data for Growthcurver function ####
     
      #prep_plate()
      #Read one plate file and return it in both shapes the pipeline needs.
  
      #data_file   tab-delimited plate-reader export (time + one column per well)
      #names_file  CSV whose Sample column lists sample IDs in plate order
      #wells       character vector of wells that held samples
      
      #Returns list(wide, long):
      #wide  time + one column per SAMPLE -> what growthcurver expects
      #long  one row per sample per timepoint -> what ggplot expects
     
      prep_plate <- function(data_file, names_file, wells) { #data file of OD reads, sample name file, and expected wells
        
        raw <- read.table(data_file, header = TRUE)   #read the plate export
        raw$time <- raw$time / 60                     #convert minutes to hours
        
        wells <- wells[wells %in% names(raw)]         #keep only wells actually present

        wide  <- raw[, c("time", wells)]              #keep the time plus the selected sample wells; drop unused wells
        
        samp <- read.csv(names_file, header = TRUE)$Sample   #read sample IDs from the names file in plate order
        stopifnot(length(samp) == ncol(wide) - 1)     #names must match well count exactly.

        colnames(wide)[-1] <- as.character(samp)      #well names -> sample IDs ([-1] skips time)
        
        long <- melt(wide,                            #one row per well per timepoint
                     id            = "time",          #time stays as the identifier column
                     variable.name = "Sample",        #former column names become this
                     value.name    = "OD") %>%        #cell values become this
          split_ids("Sample")                   #then break the ID into components
        
        list(wide = wide, long = long) #return both wide and long versions of the processed data
      }
      
      
    #### Run growthcurver package to get DT and K function ####
      
     
      #summarise_plate()
      #Fit logistic growth curves and return one row per sample with DT and K.

      #wide        the wide data frame from prep_plate()
      #pop_rename  optional map to make Population specific, e.g. c(CAS = "Low CAS")
      #treat_map   optional map to turn raw treatment codes into display labels
     
      summarise_plate <- function(wide, pop_rename = NULL, treat_map = NULL) {
        
        out <- SummarizeGrowthByPlate(wide) %>%       #fit one logistic curve per sample column
          select(sample, k, t_gen) %>%                #keep sample ID, carrying capacity, doubling time
          split_ids("sample")                         #split sample ids into pop, strength, rep, timepoint, treatment
        
        if (!is.null(pop_rename))                     #if a population rename map was supplied, apply it 
          out <- out %>% mutate(Population = relabel(Population, pop_rename))   #recode pop labels like CAS to Low CAS
        
        if (!is.null(treat_map))                      
          out <- out %>% mutate(Treatment = relabel(Treatment, treat_map)) #convert treatment labels to figure labels 
        
        out %>% mutate(Population_Treatment = paste(Population, Treatment, sep = "_")) #create combined population-treatment key; to filter outliers later
  
      }
      
    #### Recode the ancestors replicates into single random-effect level function ####
      
      #collapse_anc_reps()
      #Pool the ancestor's replicates into a single random-effect level.
  
      #The ancestor's "replicates" are technical, not independent evolved lines, so treating them as separate levels would misstate the random effect.
    
      collapse_anc_reps <- function(df, pop = "ANC") {
        df %>% mutate(Replicate = factor(              #recode rep to a factor for lmer
          ifelse(Population == pop, "A", Replicate)))  #ancestor -> "A"; everyone else untouched
      }
      
      
    #### Renumber populations replicates function ####
      
      #offset_reps()
      #Renumber one population's replicates so they don't collide with the other's.

      #When two populations are combined into one model, "rep04" appears in both. lmer would read those as the SAME replicate and pair observations that
      #have nothing to do with each other. Shifting one set by 12 keeps them separate. 
    
      #pop  the population to shift
      #by   how far to shift (12 = the number of evolved replicate lines)
   
      offset_reps <- function(df, #dataframe
                              pop, #population to renumber
                              by = 12) { #offset amount
        df %>% mutate(Replicate = ifelse(             #recalculate replicate only for the specific population
          Population == pop,                         #check which rows belong to the population being shifted
          paste0("rep", as.integer(sub("rep", "", Replicate)) + by),  #strip "rep", add offset, re-prefix
          Replicate))                                                 #other population unchanged
      }
      
      
    #### Attach conventional significance labels to results table function ####
     
      #add_sig()
      #Attach conventional significance stars to a results table.
      
      add_sig <- function(df, #dataframe
                          p_col = "p.value") {      #name of p value column
        df %>% mutate(sig_label = case_when(          #create a new column called sig_label based on the pvalue
          .data[[p_col]] < 0.001 ~ "***",     #highly sig p < 0.001; .data[[]] allows p_col to be a variable
          .data[[p_col]] < 0.01  ~ "**",      #p < 0.01; very sig
          .data[[p_col]] < 0.05  ~ "*",        #p < 0.05; conventionally sig
          TRUE ~ NA_character_))              #NA, not "", so filtering is unambiguous
      }
      
    #### Run statistical analysis function ####
      
      #run_stats()
     
      #df        a growthcurver results frame (outliers already dropped)
      #response  "t_gen" for doubling time, "k" for carrying capacity
      #rhs       right-hand side of the model; the default is your standard one
  
      #Model:  log(response) ~ Population * Treatment + (1 | Replicate)
      #log()             improves normality of both responses
      #Population *      main effect + its interaction with environment
      #Treatment
      #(1 | Replicate)   random intercept: replicates are paired across environments, so they aren't independent observations
      
      #Returns a list so one call replaces five separate assignments:
      #$model  the fitted lmer
      #$anova  significance of the fixed effects
      #$emm    estimated marginal means object
      #$means  back-transformed means (readable on the original scale)
      #$pairs  pairwise contrasts within each environment, with sig_label
    
      run_stats <- function(df, response, rhs = "Population * Treatment + (1 | Replicate)") {
        
        f <- as.formula(sprintf("log(%s) ~ %s", response, rhs))   #build the formula usign the response and rhs supplied as text
        m <- lmer(f, data = df)                                   #fit it
        emm <- emmeans(m, ~ Population | Treatment)               #marginal means, split BY environment
 
        list(
          model = m,                                              #store the fitted mixed effects model
          anova = anova(m),                                       #test significance of fixed effects
          emm   = emm,                                            #store estimated marginal means object
          means = summary(emmeans(m, ~ Population | Treatment, type = "response")),            #undoes the log for reporting
          pairs = as.data.frame(contrast(emm, method = "pairwise")) %>%   #compare pops pairwise within each treatment
            add_sig())                                              #add sig stars
      }
      
    #### Summarize data for growth curve plots function ####
      
      #gc_summary()
      #Average replicates at each timepoint for the growth-curve plots.
      #Grouping by ID_Treatment (which excludes Replicate) is what collapses the replicates. Rounding to 2 dp 

      gc_summary <- function(long) {
        long %>%
          group_by(ID_Treatment, time) %>%                        #one group per curve per timepoint
          summarise(Sum_means = round(mean(OD), 2),               #plotted line
                    Sum_sd    = round(sd(OD), 2),                 #kept for reference
                    n         = n(),                              #replicates contributing
                    Sum_se    = round(sd(OD) / sqrt(n()), 2),     #plotted error bars
                    .groups   = "drop")                           #ungroup, else later joins misbehave
      }
      
    #### Plot the growth curves function ####
      
      #gc_plot()
      #One growth-curve figure. Called seven times with different palettes.

      #sum_df          output of gc_summary(), joined to Population + Treatment
      #pop_levels      the two populations, in legend order
      #strip_fill      facet strip colours, ONE PER ENVIRONMENT (3 or 4 here)
      #strip_text_col  "black" on light strips, "white" on dark ones
     
      gc_plot <- function(sum_df, pop_levels, strip_fill, strip_text_col = "black") {
        
        d <- sum_df %>% mutate(Population = factor(Population, levels = pop_levels)) #factor order drives legend order
        
        ggplot(d, aes(time, Sum_means, col = Population, shape = Population, group = interaction(Population, Treatment))) +
          geom_line(linewidth = 0.5) +
          geom_point(size = 1.5) +
          geom_errorbar(aes(ymin = Sum_means - Sum_se,            #mean +/- 1 SE
                            ymax = Sum_means + Sum_se),
                        width = 0.5, linewidth = 0.1) +
          scale_color_manual(values = POP_COLS) +                 #shared palette
          scale_shape_manual(values = setNames(c(16, 17), pop_levels)) + #^ circle for the first population, triangle for the second
          facet_wrap2(~Treatment, nrow = 1, strip.position = "top",
                      strip = strip_themed(                       #ggh4x: per-facet strip styling
                        background_x = elem_list_rect(fill = strip_fill),
                        text_x = elem_list_text(face = "bold", color = strip_text_col, size = 18))) +
          labs(y = expression(OD[600]), x = "Time (hr)",          #OD with subscript 600
               color = "Population", shape = "Population") +
          scale_x_continuous(breaks = seq(0, 48, 12),             #0, 12, 24, 36, 48
                             limits = c(0, 48), expand = c(0, 0)) +
          scale_y_continuous(breaks = c(0.1, 0.5, 1, 1.5), expand = c(0, 0)) +
          coord_cartesian(ylim = c(0, 1.8), xlim = c(0, 50)) +    #zoom, not clip: coord_cartesian
          guides(color = guide_legend(override.aes = list(size = 4, shape = c(16, 17))),
                 shape = "none") +                                #one merged legend, not two
          theme_cowplot() +
          theme(axis.text = element_text(size = 16),
                axis.title = element_text(size = 18),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 12),
                strip.text = element_text(size = 14, face = "bold"),
                strip.placement = "outside",
                strip.background = element_rect(fill = ""),
                legend.position = c(0.15, 0.5),                   #legend inside the first panel
                legend.justification = c(0.01, 1),
                legend.box.background = element_rect(color = "black", linewidth = 0.8),
                legend.margin = margin(2, 2, 2, 2),
                panel.spacing = unit(2, "lines"))                 #gap between facets
      }
      
            
      
      
      
# =================== CONFIGURE SETTINGS TO RUN PER ASSAY =============================================      
  
    #Everything that differs between assays lives in ASSAYS. Everything that doesn't lives in fit_assay(). To add an assay, add a config entry.
  
      #Config fields:
      #data        plate-reader export path
      #names       sample-names CSV path
      #wells       which wells held samples (one of the WELLS_* constants)
      #pop_rename  make Population specific, or NULL if it's already unambiguous
      #treat_map   raw treatment codes -> display labels
      #envs        x-axis order of assay environments
      #pops        populations in this assay, in legend/dodge order
      #strip       facet strip fill colours, one per environment
      #strip_col   strip text colour
      #offset_pop  population whose replicates get renumbered (cross-tol only)
      #raw         raw plate-reader CSV (input to process_OD600_data)
      #data        formatted table it writes, and what fit_assay reads
      #drop_DT     Population_Treatment values to drop before the DT model
      #drop_K      same, before the K model
  
      #WHY THESE ROWS ARE DROPPED: cases where the ancestor either didn't grow in drug media or never reached stationary phase, so growthcurver can't
      #estimate DT or K meaningfully.
    
      #THE FOUR CAS FILTERS ARE NOT SYMMETRIC. 
      #Low CAS assay,  High CAS env -> dropped from DT and K (no growth)
      #Low CAS assay,  Low CAS env  -> dropped from K only.
      #High CAS assay, High CAS env -> dropped from DT and K (no growth)
      #High CAS assay, Low CAS env  -> dropped from DT and K. 
      
      
      ASSAYS <- list(
        
        Low_CAS = list(
          raw   = "Input_files/Low_CAS_pops_Assay1.csv",
          data  = "Output_files/Low_CAS_pops_Assay1_formatted.txt",
          names = "Input_files/Sample_names_for_Low_CAS_assay1.csv",
          wells = WELLS_STD,
          pop_rename = c(CAS = "Low CAS"),              #this plate's "CAS" means low-selection
          treat_map  = MAP_CAS,
          envs  = CAS_ENV,
          pops  = c("ANC", "Low CAS"),
          strip = c("grey50", "gold", "darkgoldenrod"), strip_col = "black",
          drop_DT = "ANC_High CAS",                     #ancestor didn't grow in high drug
          drop_K  = c("ANC_High CAS", "ANC_Low CAS")),  #K also unusable in low drug (no plateau)
        
        High_CAS = list(
          raw   = "Input_files/High_CAS_pops_Assay1.csv",
          data  = "Output_files/High_CAS_pops_Assay1_formatted.txt",
          names = "Input_files/Sample_names_for_High_CAS_assay1.csv",
          wells = WELLS_STD,
          pop_rename = c(CAS = "High CAS"),             #same code, different meaning on this plate
          treat_map  = MAP_CAS,
          envs  = CAS_ENV,
          pops  = c("ANC", "High CAS"),
          strip = c("grey50", "gold", "darkgoldenrod"), strip_col = "black",
          drop_DT = c("ANC_High CAS", "ANC_Low CAS"),   
          drop_K  = c("ANC_High CAS", "ANC_Low CAS")),
        
        Low_CLO = list(
          raw   = "Input_files/Low_CLO_pops_Assay1.csv",
          data  = "Output_files/Low_CLO_pops_Assay1_formatted.txt",
          names = "Input_files/Sample_names_for_Low_CLO_assay1.csv",
          wells = WELLS_STD,
          pop_rename = c(CLO = "Low CLO"),
          treat_map  = MAP_CLO,
          envs  = CLO_ENV,
          pops  = c("ANC", "Low CLO"),
          strip = c("grey50", "hotpink", "firebrick4"), strip_col = "white",
          drop_DT = character(0),                       #no outliers: ancestor grew everywhere
          drop_K  = character(0)),
        
        High_CLO = list(
          raw   = "Input_files/High_CLO_pops_Assay1.csv",
          data  = "Output_files/High_CLO_pops_Assay1_formatted.txt",
          names = "Input_files/Sample_names_for_High_CLO_assay1.csv",
          wells = WELLS_STD,
          pop_rename = c(CLO = "High CLO"),
          treat_map  = MAP_CLO,
          envs  = CLO_ENV,
          pops  = c("ANC", "High CLO"),
          strip = c("grey50", "hotpink", "firebrick4"), strip_col = "white",
          drop_DT = character(0),
          drop_K  = character(0)),
        
        CASCLO = list(
          raw   = "Input_files/CASCLO_pops_Assay1.csv",
          data  = "Output_files/CASCLO_pops_Assay1_formatted.txt",
          names = "Input_files/Sample_names_for_CASCLO_assay1.csv",
          wells = WELLS_10COL,                          #wider plate layout than the CAS/CLO assays
          pop_rename = NULL,                            #"CASCLO" is already unambiguous
          treat_map  = MAP_CASCLO,
          envs  = CASCLO_ENV,
          pops  = c("ANC", "CASCLO"),
          strip = c("grey50", "darkgoldenrod", "firebrick4", "orangered"), strip_col = "white",
          drop_DT = c("ANC_High CAS", "ANC_CASCLO"),    #ancestor failed in both of these
          drop_K  = c("ANC_High CAS", "ANC_CASCLO")),   #same rows for K here (unlike the CAS assays)
        
        #Cross-tolerance assay: no ancestor at all. Two EVOLVED populations compared
        #head to head, including against a drug neither was selected on (Ampho B).
        Crosstol = list(
          raw   = "Input_files/CAS_CASCLO_pops_Crosstol_Assay.csv",
          data  = "Output_files/CAS_CASCLO_pops_Crosstol_Assay_formatted.txt",
          names = "Input_files/Sample_names_for_CAS_CASCLO_Crosstol_assay.csv",
          wells = WELLS_FULL,                           #whole plate used
          pop_rename = c(CAS = "High CAS"),
          treat_map  = MAP_CROSSTOL,
          envs  = CROSSTOL_ENV,
          pops  = c("High CAS", "CASCLO"),
          strip = c("grey50", "darkgoldenrod", "orangered", "royalblue"), strip_col = "white",
          offset_pop = "High CAS",                      #renumber so its reps don't pair with CASCLO's
          drop_DT = character(0),
          drop_K  = character(0))
      )
      
      
      #raw CSV -> formatted table  RUN ONCE, THEN RE-COMMENT 
      #Only needed when the raw exports change. The formatted .txt files are read back in by fit_assay() below, so day to day this stays commented out.
    
      #Driven from ASSAYS rather than a separate pair of file_list / output_list vectors: those duplicated every path, so adding an assay meant editing
      #three lists and keeping their ORDER aligned. Here each assay's raw and formatted paths sit next to each other in one place.
    
      # for (a in ASSAYS) process_OD600_data(a$raw, a$data)

      #To redo a single assay:
      # with(ASSAYS$Low_CAS, process_OD600_data(raw, data))
      
      
      
      
# =================== RUN FUNCTIONS FOR ANALYSIS =============================================       
      
      #fit_assay()
      #Run one assay end to end. Called once per config entry.
    
      #Order matters here:
      #1. read and reshape
      #2. fit growth curves and relabel
      #3. renumber replicates (BEFORE they become a factor)
      #4. set factor levels for plotting order
      #5. apply the response-specific outlier filters
      #6. collapse ancestor replicates, then fit the models
     
      fit_assay <- function(a) {
        
        plate <- prep_plate(a$data, a$names, a$wells)              #step 1: read plate data, assign sample IDs, and create wide/long data
        
        gc <- summarise_plate(plate$wide, a$pop_rename, a$treat_map)   #step 2: fit growth curves and calculate dt/k then relabel pops/treatments
        if (!is.null(a$offset_pop))                                #step 3: check whether a pop needs its replicate numbers shifted 
          gc <- offset_reps(gc, a$offset_pop)
        
        gc <- gc %>% mutate(                                       #step 4: convert population and treatment to factors 
          Treatment  = factor(Treatment,  levels = a$envs),        #x-axis order
          Population = factor(Population, levels = a$pops))        #set population order for legends and grouped/dodged plots
        
        sum_df <- gc_summary(plate$long) %>%                       #average replicate od values at each timepoint for plotting
          left_join(distinct(gc, ID_Treatment, Population, Treatment), by = "ID_Treatment") #join labels back onto the summarized growth curves; keep one unique population/treatment label for each ID_Treatment
        
        dt_df <- gc %>% filter(!Population_Treatment %in% a$drop_DT) %>%   #step 5: remove samples/groups excluded from the doubling time analysis
          collapse_anc_reps()                                             
        k_df  <- gc %>% filter(!Population_Treatment %in% a$drop_K) %>%    #step 6: combine ancestor replicates into one random-effect lvel
          collapse_anc_reps()
        
        list(
          long    = plate$long,                                    #raw long data, if you need it
          gcurver = gc,                                            #UNFILTERED DT/K per sample
          summary = sum_df,                                        #replicate-averaged curves
          dt_df   = dt_df,                                         #what the DT model actually saw
          k_df    = k_df,                                          #what the K model actually saw
          gc_plot = gc_plot(sum_df, a$pops, a$strip, a$strip_col), #create the growth curve plot
          DT      = run_stats(dt_df, "t_gen"),                     #run model for doubling time 
          K       = run_stats(k_df,  "k"))                        #run model for carrying capacity
      }
      
      
      #Run all six. lapply keeps the list names, so RES is indexed by assay name.
      RES <- lapply(ASSAYS, fit_assay)     #runs the complete analysis pipeline for every assay 
      
      
# =================== ANALYSIS AND PLOTTING MAIN FIGURE (EVOLVED POPS VERSUS EVOLVED) =============================================       
      
  #These pool two assays to compare low- vs high-selection populations directly. The ancestor is dropped: it's a within-assay reference, and it
  #isn't the comparison being made here.
     
    #### Set up data ####
      CAS_combined <- bind_rows(
        RES$Low_CAS$gcurver  %>% mutate(Assay = "Low CAS Assay"),   #add low cas assay data and identify its source assay
        RES$High_CAS$gcurver %>% mutate(Assay = "High CAS Assay")) %>%  #add high cas assay data and identify its source assay
        filter(Population != "ANC") %>%                            #drop ancestor (this also removes every row the outlier filters would have removed)
        offset_reps("High CAS") %>%                                # shift high cas rep numbers so they dont overlap low cas reps 
        mutate(Population = factor(as.character(Population), levels = c("Low CAS", "High CAS")), #set population order
               Treatment  = factor(as.character(Treatment),  levels = CAS_ENV),        #set treatment/x axis order
               Assay      = factor(Assay, levels = c("Low CAS Assay", "High CAS Assay"))) #set assay order
      
      CLO_combined <- bind_rows(
        RES$Low_CLO$gcurver  %>% mutate(Assay = "Low CLO Assay"),
        RES$High_CLO$gcurver %>% mutate(Assay = "High CLO Assay")) %>%
        filter(Population != "ANC") %>%
        offset_reps("High CLO") %>%
        mutate(Population = factor(as.character(Population), levels = c("Low CLO", "High CLO")),
               Treatment  = factor(as.character(Treatment),  levels = CLO_ENV),
               Assay      = factor(Assay, levels = c("Low CLO Assay", "High CLO Assay")))
      
    #### Plot growth curves #### 
      CAS_combined_GC_df <- bind_rows(RES$Low_CAS$summary, RES$High_CAS$summary) %>%  #combine summarized growth curves from both assays
        filter(Population %in% c("Low CAS", "High CAS")) %>%       #evolved only, no ancestor curves
        mutate(Population = factor(as.character(Population), levels = c("Low CAS", "High CAS")), #set population plotting order
               Treatment  = factor(as.character(Treatment),  levels = CAS_ENV))           #set treatment/x-axis order
      
      CLO_combined_GC_df <- bind_rows(RES$Low_CLO$summary, RES$High_CLO$summary) %>%
        filter(Population %in% c("Low CLO", "High CLO")) %>%
        mutate(Population = factor(as.character(Population), levels = c("Low CLO", "High CLO")),
               Treatment  = factor(as.character(Treatment),  levels = CLO_ENV))
      
      CAS_combined_GCPLOT <- gc_plot(CAS_combined_GC_df, c("Low CAS", "High CAS"), c("grey50", "gold", "darkgoldenrod"), "black")
      CLO_combined_GCPLOT <- gc_plot(CLO_combined_GC_df, c("Low CLO", "High CLO"), c("grey50", "hotpink", "firebrick4"), "white")
      
      
    #### Batch effect pre-tests ####

      #Each combined dataset spans two assays, so check whether assay identity alone explains variation before pooling. Three came back non-significant;
      #the fourth did not, which is why CLO K gets its own model below.
      
      anova(lm(log(t_gen) ~ Assay, data = CAS_combined))   #ns      -> standard model
      anova(lm(log(k)     ~ Assay, data = CAS_combined))   #ns      -> standard model
      anova(lm(log(t_gen) ~ Assay, data = CLO_combined))   #ns      -> standard model
      anova(lm(log(k)     ~ Assay, data = CLO_combined))   #p ~1e-08 -> nested model, see below
      
    #### Run standard models ####
      
      CAS_combined_DT <- run_stats(CAS_combined, "t_gen") #run mixed model for CAS DT
      CAS_combined_K  <- run_stats(CAS_combined, "k")  #run mixed model for CAS K
      CLO_combined_DT <- run_stats(CLO_combined, "t_gen") #run mixed model for CLO DT
      
      
    #### Run non-standard model for CLO K ####
      
      #Assay can't go in as a fixed effect: each population appears in only one
      #assay, so Assay and Population are perfectly confounded and the model is
      #not identifiable. It goes in as a nested random effect instead - (1 | Assay/Replicate) means "replicate nested within assay".
  
      #Consequence: emmeans can't produce marginal means for this fit, so the
      #three within-environment comparisons can't come from contrast(). Instead
      #the model is refitted three times, each with a different environment as the
      #reference level, and the Population coefficient is read at each.
    
      clo_k_p <- function(ref) {
        m <- lmer(log(k) ~ Population * relevel(Treatment, ref = ref) + (1 | Assay/Replicate),
                  data = CLO_combined)                       #relevel makes `ref` the baseline, so the Population coefficient is the Low-vs-High contrast IN that environment
        coef(summary(m))["PopulationHigh CLO", "Pr(>|t|)"]   #pull that one p-value.
      }
      
      CLO_combined_K_pairs <- tibble(                               #build a results table containing the three comparisons
        contrast  = "Low CLO - High CLO",                          #same pop comparison for every assay environment
        Treatment = factor(CLO_ENV, levels = CLO_ENV),            #store the three treatment environments in the correct order
        p.value   = vapply(CLO_ENV, clo_k_p, numeric(1))) %>%      #refit the model for each environment and collect its p value
        add_sig()                                                    #same star thresholds as everything else
      
      CLO_combined_K_pairs   
      
      
      
      
    #### Make DT and K barplots ####
      
      #Grouped bars, 3 rows x 2 columns: DT down the left, K down the right.
      #x    = assay environment
      #bars = dodged by population, filled from the shared palette
      #bar height = mean, error bars = mean +/- SE, points = individual replicates
      #fixed y-axis (DT 0-6, K 0-2); K > 2 dropped at PLOT time only
  
      
      #text sizes, all in one place
      #geom_text() measures in MILLIMETRES while theme sizes are in POINTS, so
      #TXT_SIG is on a different scale from the rest and needs its own knob.
      TXT_BASE     <- 14   #theme_cowplot base size
      TXT_AXIS_TTL <- 14   #axis titles
      TXT_AXIS_TXT <- 12   #axis tick labels
      TXT_LEG_TTL  <- 13   #legend title
      TXT_LEG_TXT  <- 12   #legend entries
      TXT_SIG      <- 6    #significance stars (mm; roughly 17pt)
      TXT_TAG      <- 18   #panel letters A, B, ...
      
      
      
      #build_panel()
      #One main-figure panel.
    
      #df          data for this panel (ancestor already excluded)
      #yvar        "t_gen" or "k"
      #sig_df      a results table with Treatment + sig_label
      #treat_lvls  x-axis order
      #pop_levels  the two populations, in dodge order
      #ylab/ylim/ybrk  y-axis title, fixed limits, tick positions
      #cut_above   drop values above this before plotting (NULL = keep all)
      #show_xlab   TRUE only on the bottom row, so the label isn't repeated
    
      build_panel <- function(df, yvar, sig_df, treat_lvls, pop_levels,
                              ylab, ylim, ybrk, cut_above = NULL, show_xlab = FALSE) {
        
        d <- df %>%
          mutate(Treatment  = factor(as.character(Treatment),  levels = treat_lvls), #set treatment/x-axis order
                 Population = factor(as.character(Population), levels = pop_levels),   #set population order
                 yval = .data[[yvar]]) %>%                   #generic name so the rest is response-agnostic
          filter(!is.na(yval))
        
        if (!is.null(cut_above)) d <- d %>% filter(yval <= cut_above)   #drops the K outlier. removes values above cutoff from plotting only
        
        gmax <- d %>% group_by(Treatment) %>%                #find highest observed value within each treatment
          summarise(gmax = max(yval, na.rm = TRUE),.groups = "drop")
        
        brack <- sig_df %>%
          filter(!is.na(sig_label)) %>%                      #significant comparisons only
          mutate(Treatment = factor(as.character(Treatment), levels = treat_lvls)) %>% #match treatment factor order
          left_join(gmax, by = "Treatment") %>%               #add maximum plotted value for each treatment
          filter(!is.na(gmax)) %>%                           #guards environments with no data
          mutate(xnum = as.numeric(Treatment),               #factor level -> x position
                 x_start = xnum - 0.2, x_end = xnum + 0.2,   #bracket spans the two dodged bars
                 y_position = pmin(gmax * 1.05,              #5% above the tallest point, but
                                   ylim[2] * 0.94))          #never off the top of the panel
        
        dodge <- position_dodge(width = 0.8)                 #ONE dodge object, reused below so bars, points, and error bars stay aligned
        
        p <- ggplot(d, aes(Treatment, yval, fill = Population)) +
          stat_summary(fun = mean, geom = "bar", position = dodge,      #bar = group mean
                       width = 0.7, color = "black", linewidth = 0.4, alpha = 0.85) +
          geom_point(shape = 21, color = "black", stroke = 0.2, size = 2, alpha = 0.7,show.legend = FALSE,  #bars already carry the legend
                     position = position_jitterdodge(jitter.width = 0.12,dodge.width = 0.8)) +
          stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge,
                       width = 0.1, color = "black", linewidth = 0.8) +
          scale_fill_manual(values = POP_COLS, name = "Population", drop = TRUE) +
          scale_y_continuous(name = ylab, breaks = ybrk,
                             expand = expansion(mult = c(0, 0.04))) +   #bars sit flush on the axis
          labs(x = "Assay environment") +
          coord_cartesian(ylim = ylim) +                                #zoom, don't drop rows
          theme_cowplot(font_size = TXT_BASE) +
          theme(axis.title.x = if (show_xlab) element_text(size = TXT_AXIS_TTL) else element_blank(),
                axis.text.x  = element_text(size = TXT_AXIS_TXT),
                axis.title.y = element_text(size = TXT_AXIS_TTL),
                axis.text.y  = element_text(size = TXT_AXIS_TXT),
                legend.position = "right",
                legend.title = element_text(size = TXT_LEG_TTL, face = "bold"),
                legend.text  = element_text(size = TXT_LEG_TXT),
                plot.background = element_blank())                      #transparent, so patchwork panels don't overpaint
        
        if (nrow(brack) > 0) {                                          #only add sig brackets if at least one exists
          p <- p +
            geom_segment(data = brack, inherit.aes = FALSE,             #inherit.aes = FALSE: this layer
                         aes(x = x_start, xend = x_end,                 #has its own data, unrelated
                             y = y_position, yend = y_position),        #to the main aes()
                         linewidth = 0.7) +
            geom_text(data = brack, inherit.aes = FALSE,
                      aes(x = xnum,                                     #centered over the bracket
                          y = y_position + diff(ylim) * 0.03,           #nudged up, scaled to the axis
                          label = sig_label),
                      size = TXT_SIG, fontface = "bold")
        }
        p
      }
      
      
      
      #left column: doubling time 
      pCAS_DT  <- build_panel(CAS_combined, "t_gen", CAS_combined_DT$pairs, #Data, response, and significance results
                              CAS_ENV, c("Low CAS", "High CAS"),           #Treatment and population order
                              "DT (hours)", DT_YLIM, DT_YBRK)              #Axis label, limits, and breaks
      
      pCLO_DT  <- build_panel(CLO_combined, "t_gen", CLO_combined_DT$pairs,
                              CLO_ENV, c("Low CLO", "High CLO"),
                              "DT (hours)", DT_YLIM, DT_YBRK)
      
      pComb_DT <- build_panel(RES$Crosstol$dt_df, "t_gen", RES$Crosstol$DT$pairs,
                              CROSSTOL_ENV, c("High CAS", "CASCLO"),
                              "DT (hours)", DT_YLIM, DT_YBRK,
                              show_xlab = TRUE)                     #bottom row gets the x-axis title
      
      
      #right column: carrying capacity 
      pCAS_K  <- build_panel(CAS_combined, "k", CAS_combined_K$pairs,
                             CAS_ENV, c("Low CAS", "High CAS"),
                             KLAB, K_YLIM, K_YBRK, cut_above = 2)
      
      pCLO_K  <- build_panel(CLO_combined, "k", CLO_combined_K_pairs,   #note: the hand-built table
                             CLO_ENV, c("Low CLO", "High CLO"),
                             KLAB, K_YLIM, K_YBRK, cut_above = 2)
      
      pComb_K <- build_panel(RES$Crosstol$k_df, "k", RES$Crosstol$K$pairs,
                             CROSSTOL_ENV, c("High CAS", "CASCLO"),
                             KLAB, K_YLIM, K_YBRK, cut_above = 2, show_xlab = TRUE)
      
      
      
      #Assemble combined figure 
      #Rows, not columns: each row is one comparison, DT left and K right, so the panel letters run ACROSS (A B / C D / E F). Assembling as
      #(DT column) | (K column) instead would letter them A,B,C down the left and D,E,F down the right.
      
      #Legend suppressed on the DT panels: both panels in a row show the same two populations, so the K panel's legend serves the whole row. This is the
      #duplicate-legend deletion you were doing by hand afterwards.
      
      #The trailing & applies only plot.tag to every panel; it does not touch the legend settings above.
      
      growth_mainfig <- ((pCAS_DT  + theme(legend.position = "none")) | pCAS_K)  /
        ((pCLO_DT  + theme(legend.position = "none")) | pCLO_K)  /
        ((pComb_DT + theme(legend.position = "none")) | pComb_K) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = TXT_TAG, face = "bold"))
      
      growth_mainfig
      
      #ggsave("Figures/growth_mainfig.svg", growth_mainfig, width = 11, height = 10, bg = "white", units = "in", dpi = 300)
      
      
      
# =================== ANALYSIS AND PLOTTING SUPPLEMENTARY FIGURE (EVOLVED VERSUS ANCESTOR) ============================================     

    #### Functions to set up data as needed and plot ####
      #Evolved vs ancestor, one panel per assay.
    
      #build_supp_panel()
      #Arguments mirror build_panel(), except pop_levels must be
      #c(ancestor, evolved) — that order sets which side each bar lands on.
     
      build_supp_panel <- function(df, yvar, sig_df, treat_lvls, pop_levels,
                                   ylab, ylim, ybrk, cut_above = NULL, show_xlab = FALSE) { #Set labels, axes, cutoff, and x-label option
        
        off <- c(-0.2, 0.2); names(off) <- pop_levels    #set fixed offsets: left = ANC, right = evolved
        
        d <- df %>%
          mutate(Treatment  = factor(as.character(Treatment),  levels = treat_lvls), #set treatment order
                 Population = factor(as.character(Population), levels = pop_levels), #set pop order
                 yval = .data[[yvar]]) %>%                                    #put seelcted response DT or K into yval
          filter(!is.na(yval), !is.na(Treatment), !is.na(Population))        #remove rows missing response or grouping variable
        if (!is.null(cut_above)) d <- d %>% filter(yval <= cut_above)
        
        d <- d %>% mutate(                                      
          xnum = as.numeric(Treatment),                              #environment -> x position
          xpt  = xnum + off[as.character(Population)] +              #+ the population's fixed offset
            runif(n(), -0.055, 0.055))                          #+ small jitter 
        
        summ <- d %>%                                                #bar heights and error bars
          group_by(Treatment, Population) %>%
          summarise(m  = mean(yval),                                 #calculate mean response
                    se = sd(yval) / sqrt(n()),                       #calculcate standard error
                    .groups = "drop") %>%                           #remove grouping after summarizing
          mutate(xnum = as.numeric(Treatment),                        #convert treatment to numeric x position
                 xpos = xnum + off[as.character(Population)])         #same offsets, no jitter
        
        gmax <- d %>% group_by(Treatment) %>%                          #find tallest point in each treatment
          summarise(gmax = max(yval, na.rm = TRUE), .groups = "drop")  #store max value for significance bracket placement
        
        brack <- sig_df %>%
          filter(!is.na(sig_label)) %>%                       #keep only significant comparisons
          mutate(Treatment = factor(as.character(Treatment), levels = treat_lvls)) %>% #set treatment order
          left_join(gmax, by = "Treatment") %>%                      #add tallest value for each treatment
          filter(!is.na(gmax)) %>%                                   #skip environments with no data
          mutate(xnum = as.numeric(Treatment),
                 x_start = xnum - 0.2, x_end = xnum + 0.2,           #matches the two bar positions
                 y_position = pmin(gmax * 1.05, ylim[2] * 0.94))
        
      
        p <- ggplot() +
          geom_col(data = summ, aes(xpos, m, fill = Population),     #geom_col: heights already computed
                   width = 0.36, color = "black", linewidth = 0.4, alpha = 0.85) +
          geom_point(data = d, aes(xpt, yval, fill = Population),    #individual replicates
                     shape = 21, color = "black", stroke = 0.2, size = 2, alpha = 0.7,
                     show.legend = FALSE) +
          geom_errorbar(data = summ, aes(x = xpos, ymin = m - se, ymax = m + se),
                        width = 0.12, color = "black", linewidth = 0.8) +
          scale_fill_manual(values = POP_COLS, name = "Population", drop = FALSE) +
          scale_x_continuous(breaks = seq_along(treat_lvls),          #x is numeric now, so the tick
                             labels = treat_lvls,                     #labels must be supplied by hand
                             expand = expansion(add = 0.55)) +        #padding so edge bars aren't clipped
          scale_y_continuous(name = ylab, breaks = ybrk,
                             expand = expansion(mult = c(0, 0.08))) + #extra headroom for the brackets
          labs(x = "Assay environment") +
          coord_cartesian(ylim = ylim) +
          theme_cowplot() +
          theme(axis.title.x = if (show_xlab) element_text(size = 12) else element_blank(),
                axis.text.x  = element_text(size = 9),                #smaller than the main figure:
                axis.title.y = element_text(size = 11),               #more panels, less room each
                axis.text.y  = element_text(size = 9),
                legend.position = "right",
                legend.title = element_text(size = 10, face = "bold"),
                legend.text  = element_text(size = 9),
                plot.background = element_blank())
        
        if (nrow(brack) > 0) {
          p <- p +
            geom_segment(data = brack, inherit.aes = FALSE,
                         aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
                         linewidth = 0.7) +
            geom_text(data = brack, inherit.aes = FALSE,
                      aes(x = xnum, y = y_position + diff(ylim) * 0.03, label = sig_label),
                      size = 4.5, fontface = "bold")
        }
        p
      }
      
      
      #One entry per assay: the x-axis environments and the two populations, in dodge order (ancestor first, so it lands on the left of each pair).
      supp_spec <- list(
        Low_CAS  = list(env = CAS_ENV,    pops = c("ANC", "Low CAS")),
        High_CAS = list(env = CAS_ENV,    pops = c("ANC", "High CAS")),
        Low_CLO  = list(env = CLO_ENV,    pops = c("ANC", "Low CLO")),
        High_CLO = list(env = CLO_ENV,    pops = c("ANC", "High CLO")),
        CASCLO   = list(env = CASCLO_ENV, pops = c("ANC", "CASCLO")))
      
      
     
      #supp_row()
      #One ROW of the figure: this assay's DT panel on the left, its K panel on the right.
      
      #assay      name of an entry in supp_spec / RES
      #show_xlab  TRUE only for the bottom row of the figure
      
      #Note the two panels draw from DIFFERENT data frames: dt_df and k_df have
      #different outlier filters, which is exactly why the ancestor appears in
      #the Low CAS slot of the DT panel but not the K one.
    
      supp_row <- function(assay, show_xlab = FALSE) {
        s <- supp_spec[[assay]]
        
        dt <- build_supp_panel(RES[[assay]]$dt_df, "t_gen", RES[[assay]]$DT$pairs,
                               s$env, s$pops, "DT (hours)", DT_YLIM, DT_YBRK,
                               show_xlab = show_xlab) +
          theme(legend.position = "none")     #the row's K panel carries the legend for both, since they share populations
        
        kk <- build_supp_panel(RES[[assay]]$k_df, "k", RES[[assay]]$K$pairs,
                               s$env, s$pops, KLAB, K_YLIM, K_YBRK,
                               cut_above = 2,     #hides K > 2 in the PLOT; models keep those rows
                               show_xlab = show_xlab)
        
        dt | kk
      }
      
      
      #assemble 
      #Five rows, one per assay, in the order they appear in your figure.
      #Panel letters run ACROSS (A B / C D / ...) because each row is nested as (dt | kk) and patchwork tags in the order panels are added.
      DT_K_supp_fig <- supp_row("Low_CAS")  /
        supp_row("High_CAS") /
        supp_row("Low_CLO")  /
        supp_row("High_CLO") /
        supp_row("CASCLO", show_xlab = TRUE) +   # bottom row only
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 14, face = "bold"))
      
      DT_K_supp_fig
      
      #ggsave("Figures/DT_K_supp_fig.svg", DT_K_supp_fig, width = 10, height = 12, bg = "white", units = "in", dpi = 300)
      
      
      
   
# =================== ANALYSIS AND PLOTTING ADDITIONAL FIGURES ============================================       
   
    #### Growth curves #### 
      
      #Not part of the two figures above. Kept because the curves are the only
      #place you can SEE why a row was filtered out — a flat ancestor trace in high drug
    
      #Nothing below this line is needed to produce growth_mainfig or
      #DT_K_supp_fig. Skip it, or delete it, without consequence.
      #plot_spacer() with a small height weight is just a gap between panels.
     
      
      #evolved vs evolved (matches the main figure's comparisons) 
      DS_GC_supp <- CAS_combined_GCPLOT / plot_spacer() / CLO_combined_GCPLOT +
        plot_layout(heights = c(1, 0.1, 1)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 30, face = "bold"))
      #ggsave("DS_GC_supp.png", DS_GC_supp, width = 12, height = 8, dpi = 300)
      
      DC_GC_supp <- RES$CASCLO$gc_plot / plot_spacer() / RES$Crosstol$gc_plot +
        plot_layout(heights = c(1, 0.1, 1)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 20, face = "bold"))
      #ggsave("DC_GC_supp.png", DC_GC_supp, width = 12, height = 8, dpi = 300)
      
      #evolved vs ancestor (matches the supplementary figure's comparisons)
      CAS_ANC_GC_supp <- RES$Low_CAS$gc_plot / plot_spacer() / RES$High_CAS$gc_plot +
        plot_layout(heights = c(1, 0.1, 1)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 30, face = "bold"))
      #ggsave("CAS_ANC_GC_supp.png", CAS_ANC_GC_supp, width = 12, height = 8, dpi = 300)
      
      CLO_ANC_GC_supp <- RES$Low_CLO$gc_plot / plot_spacer() / RES$High_CLO$gc_plot +
        plot_layout(heights = c(1, 0.1, 1)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 30, face = "bold"))
      #ggsave("CLO_ANC_GC_supp.png", CLO_ANC_GC_supp, width = 12, height = 8, dpi = 300)
      
      #Individual curves, if you only need one:
      #   RES$Low_CAS$gc_plot, RES$High_CAS$gc_plot, RES$Low_CLO$gc_plot,
      #   RES$High_CLO$gc_plot, RES$CASCLO$gc_plot, RES$Crosstol$gc_plot      
     
