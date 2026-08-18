########################################
# Script Name: SNP_QC_Workflow_Script.R
# Description: SNP processing pipeline for coverage QC, filtering, PCA, SFS, and de novo variant identification
# Author: Megan Sandoval-Powers
# Last Updated: 8-2026
# R Version: 4.5.1
########################################

# ================= R PACKAGES ==============================================

#Load packages that are needed
library(tidyverse)     # ggplot2, dplyr, tidyr, stringr, readr — the backbone of the whole script
library(cowplot)       # theme_cowplot() global theme + plot_grid() for the SFS multipanels
library(patchwork)     # the / and + operators, plot_layout(), plot_annotation(tag_levels)
library(shiny)         # fluidPage/renderPlot/shinyApp for the interactive coverage browser
library(factoextra)    # fviz_screeplot() for the PCA scree plot
library(ggiraph)       # geom_point_interactive(), girafe(), opts_hover_inv() for hover PCA
library(glue)          # glue('{Population.full}') tooltip labels in the ggiraph plot
library(ggvenn)        # ggvenn() for the fixed-SNP overlap diagrams
library(svglite)       # backend for ggsave(".svg") calls
library(ggrastr)       # for geom_point_rast(), rasterizes points so exported images stay small when needed


# ================= R ENVIRONMENT ==============================================

#Set theme globally for plots
theme_set(theme_cowplot())


# ================= DATAFILES ==============================================

#read in SNP datatable from all libraries
snps=read.table("Input_files/filtered_snps.txt",header=T)  #126,215 raw SNPs including the ancestor and founders and mito SNPs



# ================= ASSESSING GENOME-WIDE COVERAGE (PRE FILTERING) ==============================================

  ##### Step 1. Preparing SNP table ############

  #1. Replace any "." in the dataset with NAs
    snps_edited <- as.data.frame(lapply(snps,function(x) gsub("\\.","NA",x))) #swap the period with NA
    test <- apply(snps_edited,2, function(x) which(x =="NA")) #Look at where all the "NAs" occurred in the dataframe. I think all in N_hap_B4
    snps_edited<- as.data.frame(lapply(snps_edited,function(x) gsub("NA","0",x))) #Now change the NAs to 0
    rm(test)
    
  #2. Exclude the mito SNPs (can do this at the start, I end up doing later after filtering)
    #snps_edited <- snps_edited %>%
      #filter(!str_detect(CHROM, "mitochondrion")) 
    
  #3. Rename column identifiers 
    colnames(snps_edited) <- colnames(snps_edited) %>%
      #Special case FIRST
      str_replace("CTR_S01", "Control_ND") %>%
      
      #General replacements
      str_replace_all("CTR", "Control") %>%
      str_replace_all("CAC", "CASCLO") %>%
      str_replace_all("CAP", "CAS") %>%
      str_replace_all("S01", "Low") %>%
      str_replace_all("S03", "High") %>%
      
      #Timepoint rename
      str_replace_all("D01", "T01") %>%
      str_replace_all("D07", "T07") %>%
      str_replace_all("D14", "T14")
    
    
  ##### Step 2. Assessing Genome-Wide Coverage Across All Experimental Populations (PRE-FILTERING) ############
    
    #1. Subset just the experimental populations (ignore ancestor and founders)
      snps_exp_pops <- snps_edited %>% select(-starts_with(c("alt_EEA_anc","N_EEA_anc","alt_EEA_hap","N_EEA_hap")))


    #2. Extract the coverage columns (this is the "N_" columns)
      cov <- snps_exp_pops[,c(seq(7,ncol(snps_exp_pops),by=2))] #Start at 7th column because we don't need the POS info, go by the number of columns in the dataframe by 2 to grab the N columns.
      cov[] <- lapply(cov,as.numeric) #Now let's change the columns to be numeric. Need this to run colMeans


    #3. Get the mean total coverage for each population (column)
      meancov <- data.frame(colMeans(cov)) #make sure to ignore NA instances

    #4. Calculate the standard deviation in coverage
      sdcov <- as.data.frame(apply(cov,2,sd))

    #5. Combine mean coverage and standard deviation in a new summary dataframe
      summary_cov <- data.frame(meancov,sdcov)
      colnames(summary_cov)[1] <- ("Mean_Coverage") #rename coverage column
      colnames(summary_cov)[2] <- ("Standard_deviation")

    #6. Move the population names from the row name to its own column in the dataframe
      summary_cov$Population <- rownames(summary_cov) 
      rownames(summary_cov) <- NULL #Set the row numbers back to default

    #7. Look at other summary stats
      range(summary_cov$Mean_Coverage) #What is the range of coverage across all pops? 27X to 279X
      median(summary_cov$Mean_Coverage) #What is the median coverage across all pops? 83X

    #8. Look at populations below a certain threshold of coverage.
      cov_below_30 <- summary_cov[which(summary_cov$Mean_Coverage <30),] #Which are below 30X? N_EEA_CAC_S01_rep12_D07 (27X)
      nrow(summary_cov[which(summary_cov$Mean_Coverage <30),]) #How many are below 30X? 1.
      
      cov_below_40 <- summary_cov[which(summary_cov$Mean_Coverage <40),] #Which are below 40X? 
      nrow(summary_cov[which(summary_cov$Mean_Coverage <40),]) #How many are below 40X? 25 total. 


    #9. Make a quick histogram to plot the distribution of mean coverages across the genome for each sample
      cov_hist_plot <- summary_cov |>
        ggplot(aes(x=Mean_Coverage, fill = after_stat(x))) +
        geom_histogram(binwidth=12, show.legend=FALSE,color="black") +
        scale_fill_viridis_c(option = "C", alpha=0.8) +
        scale_x_continuous(name = "Mean Coverage", breaks=seq(0,300,by=40)) +
        scale_y_continuous(name = "Frequency", breaks=seq(0,30,by=5)) +
        theme(
          axis.text.y = element_text(size =12),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14))
      cov_hist_plot
      
    #10. Remove any df not needed
      rm(meancov,sdcov,cov_below_30,cov_below_40)

  ##### Step 3. Plotting coverage across the genome for each population (PRE-FILTERING) ##############
      
  #1. Setting up the dataframe 
      
    #Extract the required columns
      start_cols <- snps_edited[, c(2,3), drop=FALSE] #grabbing the POS and CHROM columns
      end_cols <- snps_edited[,c(439), drop=FALSE] #grabbing the ancestor coverage column
      
    #Combine these columns into the cov dataframe set up earlier which is just the coverages for all experimental samples
      combined_cov_df <- cbind(start_cols,cov,end_cols)
      
      
    #Change the POS column to be numeric 
      combined_cov_df$POS <- as.numeric(combined_cov_df$POS) 
      
    #for now, we will exclude the mito snps since it would interfere with the plot
      combined_cov_df <- combined_cov_df %>%
        filter(!str_detect(CHROM, "mitochondrion")) #now at 126,085 SNPs
      
      
    #Make a new cumulative position column since the current one restarts at each chromosome
      chr_sizes <-  c(I = 230218, II = 813184, III = 316620, IV = 1531933, V = 576874, VI = 270161,
                      VII = 1090940, VIII = 562643, IX = 439888, X = 745751, XI = 666816,
                      XII = 1078177, XIII = 924431, XIV = 784333, XV = 1091291, XVI = 948066) #store the chromosome sizes for our ref genome version R64-1-1
      
    #Calculate cumulative start positions for each chromosome. For example, chrII starts right after chr1 ends, and so on. 
      chr_starts <- c(0, cumsum(chr_sizes[-length(chr_sizes)]))  #exclude last chr because it doesnt need a start value for the next one
      names(chr_starts) <- names(chr_sizes) #assign chromosome names so we can use CHROM as an index to look up start positions
      
    #Add cumulative position to your dataframe. For each snp, we find the start position of its chromosome and add the snps position to get its genome-wide position
      combined_cov_df$Gaxis <- with(combined_cov_df, POS + chr_starts[CHROM])
      
    #convert cumulative position from base pairs to megabases for easier plotting
      combined_cov_df$MB <- combined_cov_df$Gaxis / 1e6
      
    #Reorder columns to place the MB column right after POS and remove the gaxis column
      combined_cov_df <- combined_cov_df %>%
        relocate(MB, .after = POS) %>%
        select(-Gaxis)  
      
      
  #2. Pivot to long format and normalize coverage to each sample's own median
      
      #Flipping the data vertically for coverage plots. 
      combined_cov_df <- combined_cov_df |>
        mutate(across(starts_with("N_"), ~ as.numeric(as.character(.)))) #first transform the columns from a character to numeric
      
      combined_cov_df_long <- combined_cov_df |>
        select(CHROM,MB, starts_with("N_")) |> #keep the CHROM and MB columns, and the coverage columns for each sample (they all start with N_)
        pivot_longer(
          cols = starts_with("N_"),
          names_to = "Sample",
          values_to = "Coverage"
        ) |>   #pivot the data from horizontal to long format. each row becomes a single cov measurement for a given snp, sample, and chromosome
        mutate(CHROM = as.factor(CHROM)) #convert chromosome to a factor for plotting
      
      
      #2. Normalize coverage to the genome average per sample
      #this code is extracting the timepoint and sample group from the sample name and saving in new columns
      combined_cov_df_long <- combined_cov_df_long |> 
        mutate(
          Timepoint = str_extract(Sample, "T\\d+"),
          SampleGroup = str_remove(Sample, "_T\\d+$"))
      
      # joins the per-sample mean onto each row by matching Sample to Population, 
      # then divides the position-level coverage by that samples genome-wide mean
      combined_cov_df_long_normalized <- combined_cov_df_long |>
        left_join(summary_cov |> select(Population, Mean_Coverage),
                  by = c("Sample" = "Population")) |>
        mutate(Normalized_Coverage = Coverage / Mean_Coverage)
      
      rm(combined_cov_df_long) # Remove intermediate frame
      
    
      
    #3. Export one plot per sample to a single PDF
      
      #stores a list of all unique sample names for later use
      samples <- unique(combined_cov_df_long_normalized$Sample)
      
      ## Export plots as a pdf instead of shiny app
      # This takes 3-5 minutes to complete with the raster step and dpi set to 200
      #pdf("Figures/Coverage/normalized_coverage_plots.pdf", width = 10, height = 6)
      
      for (s in samples) {
        df <- combined_cov_df_long_normalized |>
          filter(Sample == s)
        
        df$CHROM <- as.factor(df$CHROM)
        
        p <- ggplot(df, aes(x = MB, y = Normalized_Coverage, color = CHROM)) +
          geom_point_rast(size = 1, alpha = 0.7, raster.dpi = 100) + # Instead of geom_point, make rasterized version for faster viewing in pdf
          scale_color_manual(values = rep(c("black", "gray60"), length.out = length(levels(df$CHROM)))) +
          scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
          theme_cowplot() +
          theme(legend.position = "none",
                axis.title = element_text(size = 16, face = "bold"),
                axis.text = element_text(size = 14)) +
          coord_cartesian(ylim = c(0, 10)) +
          labs(title = s,
               x = "Genomic Position (MB)",
               y = "Normalized Coverage")
        
        print(p)
      }
      
      dev.off()
      
   
# ================= FILTERING FOR QC AND POLY SNPS ==============================================

    ##### Step 1. Filtering for candidate polymorphic SNPs ############   

    #1. Filter SNP table to remove sites with low minimum coverage.
      mincov <- apply(cov,1,min) #Find the minimum coverage across each row aka SNP (from dataframe excluding ancestor and founders, only doing this for experimental populations)
      
      snps_edited_temp <- snps_edited %>% #Bind the minimum coverage to the snp dataset, also exclude the nmiss,alt,and ref columns
        select(2, 3, 6:ncol(.)) %>%
        mutate(mincov = mincov)
     
      test_cov_filt_10X <- snps_edited_temp %>% #check how many SNPs if 10X was threshold
        filter(mincov > 10) #123,612 with 10X filter, after filtering for polymorphism in the ancestor this leaves us with 70K SNPs. Filtering with 20X would put at less than 40K SNPs. 
    
      snps_edited_cov_filt_10X <- test_cov_filt_10X #using 10X for now. Just making it a new dataframe called what i want. 
      rm(test_cov_filt_10X)
      
      snps_edited_cov_filt_10X[,3:ncol(snps_edited_cov_filt_10X)] <- 
      lapply(snps_edited_cov_filt_10X[,3:ncol(snps_edited_cov_filt_10X)], as.numeric) #changing the af and N columns to be numeric since they are characters previously


    #2. Filter for polymorphic sites based on the ancestor (those likely linked to SGV and not de novo mutations or sequencing errors)
      
      #a. Convert allele counts to allele frequencies. (divide the alt count columns by the N coverage columns at each SNP)
        snps_edited_cov_filt_10X[,c(seq(3,ncol(snps_edited_cov_filt_10X)-1,by=2))] <- #this will make sure it just overwrites the previous count columns in the same position
          snps_edited_cov_filt_10X[,c(seq(3,ncol(snps_edited_cov_filt_10X)-1,by=2))]/snps_edited_cov_filt_10X[,c(seq(4,ncol(snps_edited_cov_filt_10X),by=2))]
      
      
        colnames(snps_edited_cov_filt_10X) <- gsub("alt", "af", colnames(snps_edited_cov_filt_10X)) #change the alt_ to be af_ since they are now frequencies
      
      #b. Filter based on the ancestor. Remove sites that are not polymorphic in the ancestor aka these sites are fixed. 
        grep("^af_EEA_anc", colnames(snps_edited_cov_filt_10X)) #Figure out what column number the allele freq of the ancestor is 
        
        snps_edited_cov_poly_filt1 <- snps_edited_cov_filt_10X %>%
          filter(snps_edited_cov_filt_10X[,435]<0.98 & snps_edited_cov_filt_10X[,435]>0.02) #filter out any sites where ancestor freq is less than 2% or greater than 98%. Now have 70,137 SNPs 
      
        
    #3. Filter based on the founders 
      
        #a. Filter based on the founders. Two steps. Step 1, filter out SNPs that appear to be polymorphic within a founder
          #This is because our founders are haploid, isogenic. So we do not expect there to be variation within a founder. 
            snps_edited_cov_poly_filt1 |> 
              filter(af_EEA_hap_A1_00==0 | af_EEA_hap_A1_00==1) |>
              filter(af_EEA_hap_A2_00==0 | af_EEA_hap_A2_00==1) |>
              filter(af_EEA_hap_B3_00==0 | af_EEA_hap_B3_00==1) |>
              filter(af_EEA_hap_B4_00==0 | af_EEA_hap_B4_00==1) -> snps_edited_cov_poly_filt2 #now down to 64,809 SNPs. 
        
        #b. Step 2, filter sites that are fixed across ALL founders.These are distinct founders so we expect there to be differences across founders. 
            snps_edited_cov_poly_filt2 <- snps_edited_cov_poly_filt2 %>% 
              rowwise() %>%
              filter(!all(af_EEA_hap_A1_00==af_EEA_hap_A2_00,af_EEA_hap_A1_00==af_EEA_hap_B3_00,af_EEA_hap_A1_00==af_EEA_hap_B4_00)) #now down to 64,781 SNPs
        
    
    #4. Count how many of the filtered polymorphic SNPs are from nuclear genome vs mitochondrial genome
        mito_snps_count <- sum(snps_edited_cov_poly_filt2$CHROM == "mitochondrion")  #17 snps are from mitochondrion
            
        #Count number of non-mitochondrial SNPs
        non_mito_snps_count <- sum(snps_edited_cov_poly_filt2$CHROM != "mitochondrion") #64,764 from nuclear genome
        
    #5. Exclude mito SNPs for downstream analysis
        snps_edited_cov_poly_filt2 <- snps_edited_cov_poly_filt2[snps_edited_cov_poly_filt2$CHROM != "mitochondrion",] #down to 64,764 nuclear SNPs
        
        
    #6. Export this filtered SNP table for downstream analysis
      #write_tsv(snps_edited_cov_poly_filt2,"Output_files/snps_filtered_poly.txt")

        
    ##### Step 2. Assess genome-wide coverage post-filtering #####
        
        #Subset just the experimental populations (ignore ancestor and founders)
        filter_snps_pops_cov <- snps_edited_cov_poly_filt2 %>% select(-starts_with(c("af_EEA_anc","N_EEA_anc","af_EEA_hap","N_EEA_hap","mincov")))
        
        
        #Extract the coverage columns (this is the "N_" columns)
        filter_cov <- filter_snps_pops_cov[,c(seq(4,ncol(filter_snps_pops_cov),by=2))] 
        filter_cov[] <- lapply(filter_cov,as.numeric) #Now let's change the columns to be numeric. Need this to run colMeans
        
        
        #Get the mean total coverage for each population (column)
        filter_meancov <- data.frame(colMeans(filter_cov)) #make sure to ignore NA instances
        
        #Calculate the standard deviation in coverage
        filter_sdcov <- as.data.frame(apply(filter_cov,2,sd))
        
        #Combine mean coverage and standard deviation in a new summary dataframe
        filter_summary_cov <- data.frame(filter_meancov,filter_sdcov)
        colnames(filter_summary_cov)[1] <- ("Mean_Coverage") #rename coverage column
        colnames(filter_summary_cov)[2] <- ("Standard_deviation")
        
        #Move the population names from the row name to its own column in the dataframe
        filter_summary_cov$Population <- rownames(filter_summary_cov) 
        rownames(filter_summary_cov) <- NULL #Set the row numbers back to default
        
        #Look at other summary stats
        range(filter_summary_cov$Mean_Coverage) #What is the range of coverage across all pops? 27X to 290X
        median(filter_summary_cov$Mean_Coverage) #What is the median coverage across all pops? 86X
        
        #edit the table so that it can be exported for manuscript with conventional naming scheme
          filter_summary_cov <- filter_summary_cov %>%
              separate(Population,
                       into = c("N","Exp", "Population", "Dose", "rep", "Timepoint"),
                       sep = "_") %>%
              mutate(
                rep = gsub("rep", "", rep),
                Timepoint = gsub("T", "", Timepoint),
                Treatment = paste(Population, Dose, sep="_")) #separate population names
          
          filter_summary_cov <- filter_summary_cov %>%
            select(Treatment = Treatment,
                   Replicate = rep,
                   Timepoint = Timepoint,
                   Avg_Genome_Cov = Mean_Coverage) #grab just the columns we want and rename
          
          filter_summary_cov <- filter_summary_cov %>%
            mutate(Timepoint = as.numeric(Timepoint)) #remove the 0 before timepoint
          
          
          filter_summary_cov <- filter_summary_cov %>% #Rename populations 
            mutate(
              Treatment = recode(Treatment,
                                 "CASCLO_Low" = "CASCLO",
                                 "CAS_Low" = "Low CAS",
                                 "CAS_High" = "High CAS",
                                 "CLO_Low" = "Low CLO",
                                 "CLO_High" = "High CLO",
                                 "Control_ND" = "Control"),
              Treatment = factor(Treatment,
                                 levels = c("Low CAS", "High CAS", "Low CLO", "High CLO", "CASCLO", "Control")))
          
          #Export table
          #write.csv(filter_summary_cov,file="Output_files/Coverage_Summary_Filtered.csv",row.names=FALSE)
          
          #remove df not needed
          rm(snps_edited_cov_filt_10X,snps_edited_cov_poly_filt1,filter_cov,filter_meancov,filter_sdcov)
        

          
# ================= PCA OF FILTERED POLY SNPS ==============================================         

          
    ##### Step 1. Set up datatable for PCA (all treatments) ########
    
      #1. Pull out just the allele frequency columns and exclude the founder columns 
      poly_snps_af_only <- snps_edited_cov_poly_filt2[,c(seq(3,435, by=2))]
      
      #2. Transform the data into a matrix
      pca_data <- t(as.matrix(poly_snps_af_only))
      
      #3. Run PCA and standardize. 
      pca <- prcomp(pca_data,scale=TRUE) #scales = TRUE will standardize so that all variables have stdev of 1
      
      names(pca) #can look at column names from output of pca to find what you are interested in
      
      pca.var <- pca$sdev^2 #calculate the variance explained by each principal component
      
      pca.var.per <- round(pca.var/sum(pca.var)*100,1) #compute the proportion of variance explained and change to percentage. 
      
      ##looks like the most variance is explained by PC1 and PC2. But, still quite a few above a value of 1. 
      
      summary(pca) #here is the output of pca
      
      
    #4. Make a screeplot to see the eigenvalues for each individual PC. This can help figure out how much variation each PC captures
    scree <- fviz_screeplot(pca,ncp=15,hjust=-0.1,
                            barfill="darkturquoise")+
      ylim(0,25)+
      theme_classic()+
      labs(title= "Variances - PCA", x="Principal Components",y="% of variances")
    scree 
    
    
    #5. Now pull out the principal components and store so that we can easily plot any of them
    PC1 <- pca$x[,1]
    PC2 <- pca$x[,2]
    PC3 <- pca$x[,3]
    PC4 <- pca$x[,4]
    PC5 <- pca$x[,5]
    
    
    pca.data <- data.frame(Sample=rownames(pca$x), #Make dataframe with all of the PCs, we can choose which to plot later
                           PC1=pca$x[,1],
                           PC2=pca$x[,2],
                           PC3=pca$x[,3],
                           PC4=pca$x[,4],
                           PC5=pca$x[,5])
    rownames(pca.data) <- NULL #set row names back to default
    
    
    #6. Set up the dataframe for plotting 
    pca.data$Pop.edit <- gsub("af_EEA*_","",pca.data$Sample) #first remove the af_EEA before the population names
    
    #Now create new columns that split the identifiers up
    pca.data[c('Treat','Dose','Replicate','Timepoint')] <- str_split_fixed(pca.data$Pop.edit, '[_]',4)
    
    #create a sample identifier with treatment and strength combined
    pca.data$Population <- paste0(pca.data$Treat,"_",pca.data$Dose)
    
    #Reset row numbers to default
    rownames(pca.data) <- NULL
    
    #Change the names of the ancestor to be how I want
    pca.data[217,8:12] <- c("ANC","","rep01","T00","Ancestor") 
    
    #Change the control column to be named to just "control" instead of "control_S01"
    pca.data$Population <- gsub("Control_ND","Control",pca.data$Population)
    
    #Change CASCLO to just be named CASCLO instead of CASCLO_S01
    pca.data$Population <- gsub("CASCLO_Low","CASCLO",pca.data$Population)
    
    
    
    #Rearrange the order for legend and relabel for easier interpretation
    pca.data$Population <- factor(pca.data$Population,
                                  levels = c("Ancestor",
                                             "CAS_Low","CAS_High",
                                             "CLO_Low","CLO_High",
                                             "CASCLO","Control"),
                                  labels = c("Ancestor", 
                                             "Low CAS", "High CAS", 
                                             "Low CLO", "High CLO", 
                                             "CASCLO", "Control"))
    #Rename the timepoints
    pca.data$Timepoint <- factor(pca.data$Timepoint,
                                 levels = c("T00", "T01", "T07", "T14"),
                                 labels = c("0", "1", "7", "14"))
    
    #7. PCA plot (non-interactive version)
    
    #color palette for treatments
    pop.palette=c("black","gold","darkgoldenrod","hotpink","firebrick4","orangered","grey50")
    
    
    #For PC1 vs PC2
    pca.p1p2 <- pca.data |>
      ggplot(aes(x=PC1, y=PC2,color=Population,shape=Timepoint)) +
      geom_text(size=0,label="none")+
      geom_point(size=3.5,alpha=0.75) +   
      scale_shape_manual(values=c(18,16,17,15))+
      scale_color_manual(values=pop.palette)+
      xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
      ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
      theme(legend.position="right")
    pca.p1p2
    
  
    
    #8. PCA plot (interactive version)
    
    #Population full name for interactive plot
    # pca.data$Population.full <- paste0(pca.data$Population,"_",pca.data$Replicate,"_",pca.data$Timepoint)
    # 
    # #Make interactive plot now
    # pca.interac <- pca.data |> 
    #   ggplot(aes(x=PC1, y=PC2,color=Population))+ 
    #   geom_point_interactive(aes(tooltip = glue('{Population.full}'),shape=Timepoint), size=4.5, alpha=0.75)+
    #   geom_text(size=0,label="none")+
    #   scale_shape_manual(values=c(18,16,17,15))+
    #   scale_color_manual(values=pop.palette)+
    #   xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
    #   ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
    #   theme_bw()+
    #   theme(legend.position="right")+
    #   ggtitle("")
    # pca.interac
    # 
    # 
    # girafe(ggobj = pca.interac,
    #        options = list(
    #          opts_hover_inv(css = "opacity:0.1;"),
    #          opts_tooltip(css = "padding:3px;background-color:#333333;color:white;font-size:0.8rem",offx = 7, offy = 7)
    #        ),
    #        height_svg = 4.5,
    #        width_svg = 9
    # )
    
    
    
    ##### Step 2. PC1 by timepoint trajectory figure ########
    
    
      #1. Filter out the ancestor timepoint
      # pca.data.filtered <- pca.data %>%
      #   filter(Timepoint != "0")
      # 
      # #2. Make sure Timepoint is an ordered factor (important for correct x-axis order)
      # pca.data.filtered$Timepoint <- factor(
      #   pca.data.filtered$Timepoint,
      #   levels = c("1", "7", "14"),  # change levels if needed
      #   ordered = TRUE)
      # 
      # 
      # #3. Average the replicates per timepoint and population. Get the mean and SD
      # pca.summary_PC1 <- pca.data.filtered %>%
      #   group_by(Population, Timepoint) %>%
      #   summarise(
      #     mean_PC1 = mean(PC1, na.rm = TRUE),
      #     se_PC1 = sd(PC1, na.rm = TRUE) / sqrt(n()))
      # 
      # #4. Build the plot
      # 
      # #Order the populations according to how I want to show them on the figure
      # pca.summary_PC1$Population <- factor(
      #   pca.summary_PC1$Population,
      #   levels = c("Low CAS", "High CAS", "Low CLO", "High CLO", "CASCLO", "Control"))
      # 
      # #Color palette
      # pop.palette2 <- c(
      #   "gold",           # Low CAS
      #   "darkgoldenrod",  # High CAS
      #   "hotpink",        # Low CLO
      #   "firebrick4",     # High CLO
      #   "orangered",      # CASCLO
      #   "grey50")          # Control
      # 
      # 
      # #Plot PC1 trajectory per replicate line
      # pc1_time_plot <- 
      #   ggplot(pca.summary_PC1,aes(x = Timepoint,y = mean_PC1,color = Population,group = Population)) +
      #   geom_point(size = 3) +                             # mean points
      #   geom_line(size = 0.8) +                              # line through means
      #   geom_errorbar(aes(ymin = mean_PC1 - se_PC1,ymax = mean_PC1 + se_PC1),width = 0.1, linewidth = 0.3) +           # error bars (SE)
      #   scale_color_manual(values = pop.palette2) +
      #   scale_x_discrete(labels = c("D01" = "1", "D07" = "7", "D14" = "14")) +   # <- relabel here
      #   labs(x = "Timepoint",y = paste("PC1 - ", pca.var.per[1], "% variance explained", sep = "")) +
      #   theme(
      #     legend.position = "right")
      # pc1_time_plot
    
    
    ##### Step 3. PCA for individual treatments ########    
      
      ### CAS ###
      
        #1. Get sample names for CAS treatments, control, and ancestor
        cas_samples <- pca.data %>%
          filter(Population %in% c("Low CAS", "High CAS", "Control", "Ancestor")) %>%
          pull(Sample)
        
        #2. Subset AF matrix to CAS samples only
        pca_data_CAS <- t(as.matrix(poly_snps_af_only[, cas_samples]))  
        
        #3 Run PCA
        pca_CAS <- prcomp(pca_data_CAS, scale=TRUE)
        
        #4 Calculate variance explained
        pca_CAS.var <- pca_CAS$sdev^2
        pca_CAS.var.per <- round(pca_CAS.var/sum(pca_CAS.var)*100, 1)
        
        #5 Build plotting dataframe
        pca_CAS.data <- data.frame(
          Sample = rownames(pca_CAS$x),
          PC1 = pca_CAS$x[,1],
          PC2 = pca_CAS$x[,2]) %>%
          left_join(pca.data %>% select(Sample, Population, Timepoint, Replicate),
                    by="Sample")
        
        #6 Keep T0 only for ancestor, filter other populations to T1, T7, T14
        pca_CAS.data.plot <- pca_CAS.data %>%
          filter(Timepoint != "0" | Population == "Ancestor") %>%
          mutate(Timepoint_plot = ifelse(Population == "Ancestor", "0", as.character(Timepoint)),
                 Timepoint_plot = factor(Timepoint_plot, levels=c("0","1","7","14")))
        
        #7 Set factor levels for legend order
        pca_CAS.data.plot$Population <- factor(pca_CAS.data.plot$Population,
                                               levels=c("Ancestor","Control","Low CAS","High CAS"))
        
        #8. Colors and shapes
        pop.colors.CAS <- c("Ancestor"="black", "Control"="grey50",
                            "Low CAS"="gold", "High CAS"="darkgoldenrod")
        pop.shapes.CAS <- c("Ancestor"=16, "Control"=15, "Low CAS"=17, "High CAS"=17)
        
        #9. Plot
        pca_CAS_SCATTER <- pca_CAS.data.plot %>%
          ggplot(aes(x=PC1, y=PC2, color=Population, shape=Population)) +
          geom_point(size=5, alpha=1) +
          scale_color_manual(values=pop.colors.CAS) +
          scale_shape_manual(values=pop.shapes.CAS) +
          facet_wrap(~Timepoint_plot, nrow=1,
                     labeller=labeller(Timepoint_plot=c("0"="T0","1"="T1",
                                                        "7"="T7","14"="T14"))) +
          xlab(paste("PC1 - ", pca_CAS.var.per[1], "%", sep="")) +
          ylab(paste("PC2 - ", pca_CAS.var.per[2], "%", sep="")) +
          theme_cowplot() +
          theme(legend.position="right",
                strip.text=element_text(size=16, face="bold"),
                legend.text = element_text(size=14),
                axis.text = element_text(size=14),
                axis.title = element_text(size=16),
                legend.title = element_text(size=16),
                panel.spacing=unit(1,"lines"),
                panel.border=element_rect(color="black",fill=NA,linewidth=0.8))
        pca_CAS_SCATTER
        
        
      ### CLO ###
    
        #1. Get sample names for CLO treatments, control, and ancestor
        clo_samples <- pca.data %>%
          filter(Population %in% c("Low CLO", "High CLO", "Control", "Ancestor")) %>%
          pull(Sample)
        
        #2. Subset AF matrix to CLO samples only
        pca_data_CLO <- t(as.matrix(poly_snps_af_only[, clo_samples]))
        
        #3. Run PCA
        pca_CLO <- prcomp(pca_data_CLO, scale=TRUE)
        
        #4. Calculate variance explained
        pca_CLO.var <- pca_CLO$sdev^2
        pca_CLO.var.per <- round(pca_CLO.var/sum(pca_CLO.var)*100, 1)
        
        #5. Build plotting dataframe
        pca_CLO.data <- data.frame(
          Sample = rownames(pca_CLO$x),
          PC1 = pca_CLO$x[,1],
          PC2 = pca_CLO$x[,2]) %>%
          left_join(pca.data %>% select(Sample, Population, Timepoint, Replicate),
                    by="Sample")
        
        #6. Keep T0 only for ancestor
        pca_CLO.data.plot <- pca_CLO.data %>%
          filter(Timepoint != "0" | Population == "Ancestor") %>%
          mutate(Timepoint_plot = ifelse(Population == "Ancestor", "0", as.character(Timepoint)),
                 Timepoint_plot = factor(Timepoint_plot, levels=c("0","1","7","14")))
        
        #7. Set factor levels for legend order
        pca_CLO.data.plot$Population <- factor(pca_CLO.data.plot$Population,
                                               levels=c("Ancestor","Control","Low CLO","High CLO"))
        
        #8. Colors and shapes
        pop.colors.CLO <- c("Ancestor"="black", "Control"="grey50",
                            "Low CLO"="hotpink", "High CLO"="firebrick4")
        pop.shapes.CLO <- c("Ancestor"=16, "Control"=15, "Low CLO"=17, "High CLO"=17)
        
        #9. Plot
        pca_CLO_SCATTER <- pca_CLO.data.plot %>%
          ggplot(aes(x=PC1, y=PC2, color=Population, shape=Population)) +
          geom_point(size=5, alpha=1) +
          scale_color_manual(values=pop.colors.CLO) +
          scale_shape_manual(values=pop.shapes.CLO) +
          facet_wrap(~Timepoint_plot, nrow=1,
                     labeller=labeller(Timepoint_plot=c("0"="T0","1"="T1",
                                                        "7"="T7","14"="T14"))) +
          xlab(paste("PC1 - ", pca_CLO.var.per[1], "%", sep="")) +
          ylab(paste("PC2 - ", pca_CLO.var.per[2], "%", sep="")) +
          theme_cowplot() +
          theme(legend.position="right",
                strip.text=element_text(size=16, face="bold"),
                legend.text = element_text(size=14),
                axis.text = element_text(size=14),
                axis.title = element_text(size=16),
                legend.title = element_text(size=16),
                panel.spacing=unit(1, "lines"),
                panel.border=element_rect(color="black", fill=NA, linewidth=0.8))
        pca_CLO_SCATTER
        
        
      ### CASCLO ###
    
        #1. Get sample names for CASCLO treatments,High CLO, High CAS, control, and ancestor
        casclo_samples <- pca.data %>%
          filter(Population %in% c("CASCLO", "High CAS", "High CLO", "Control", "Ancestor")) %>%
          pull(Sample)
        
        #2. Subset AF matrix to CLO samples only
        pca_data_CASCLO <- t(as.matrix(poly_snps_af_only[, casclo_samples]))
        
        #3. Run PCA
        pca_CASCLO <- prcomp(pca_data_CASCLO, scale=TRUE)
        
        #4. Calculate variance explained
        pca_CASCLO.var <- pca_CASCLO$sdev^2
        pca_CASCLO.var.per <- round(pca_CASCLO.var/sum(pca_CASCLO.var)*100, 1)
        
        #5. Build plotting dataframe
        pca_CASCLO.data <- data.frame(
          Sample = rownames(pca_CASCLO$x),
          PC1 = pca_CASCLO$x[,1],
          PC2 = pca_CASCLO$x[,2]) %>%
          left_join(pca.data %>% select(Sample, Population, Timepoint, Replicate),
                    by="Sample")
        
        #6. Keep T0 only for ancestor
        pca_CASCLO.data.plot <- pca_CASCLO.data %>%
          filter(Timepoint != "0" | Population == "Ancestor") %>%
          mutate(Timepoint_plot = ifelse(Population == "Ancestor", "0", as.character(Timepoint)),
                 Timepoint_plot = factor(Timepoint_plot, levels=c("0","1","7","14")))
        
        #7. Set factor levels for legend order
        pca_CASCLO.data.plot$Population <- factor(pca_CASCLO.data.plot$Population,
                                                  levels=c("Ancestor","Control","CASCLO","High CAS", "High CLO"))
        
        #8. Colors and shapes
        pop.colors.CASCLO <- c("Ancestor"="black", "Control"="grey50", "CASCLO" = "orangered",
                               "High CAS"="darkgoldenrod", "High CLO"="firebrick4")
        pop.shapes.CASCLO <- c("Ancestor"=16, "Control"=15,"CASCLO"=17, "High CAS"=17, "High CLO"=17)
        
        #9. Plot
        pca_CASCLO_SCATTER <- pca_CASCLO.data.plot %>%
          ggplot(aes(x=PC1, y=PC2, color=Population, shape=Population)) +
          geom_point(size=5, alpha=1) +
          scale_color_manual(values=pop.colors.CASCLO) +
          scale_shape_manual(values=pop.shapes.CASCLO) +
          facet_wrap(~Timepoint_plot, nrow=1,
                     labeller=labeller(Timepoint_plot=c("0"="T0","1"="T1",
                                                        "7"="T7","14"="T14"))) +
          xlab(paste("PC1 - ", pca_CASCLO.var.per[1], "%", sep="")) +
          ylab(paste("PC2 - ", pca_CASCLO.var.per[2], "%", sep="")) +
          theme_cowplot() +
          theme(legend.position="right",
                strip.text=element_text(size=16, face="bold"),
                legend.text = element_text(size=14),
                axis.text = element_text(size=14),
                axis.title = element_text(size=16),
                legend.title = element_text(size=16),
                panel.spacing=unit(1, "lines"),
                panel.border=element_rect(color="black", fill=NA, linewidth=0.8))
        pca_CASCLO_SCATTER   
        
        
    ##### Step 4. Export multipanel PCA figure #####
    
      #1. Three rows each for the individual comparisons i want to make
        combined_PCA_plot_supp <- pca_CAS_SCATTER /
          pca_CLO_SCATTER /
          pca_CASCLO_SCATTER +
          plot_annotation(tag_levels = 'A') &
          theme(plot.tag = element_text(size=30,face="bold")) 
        
      
      #ggsave("Figures/PCA/combined_PCA_plot_supp.svg", combined_PCA_plot_supp, width = 14, height = 10, dpi = 300)
        
        rm(pca_data_CAS,pca_data_CASCLO,pca_data_CLO,pca_data,pca_CAS_SCATTER,pca_CLO_SCATTER,pca_CASCLO_SCATTER)
  
     
        
# ================= SITE FREQUENCY SPECTRA ==============================================

    ##### Step 1. Format filtered SNP dataframe and pull individual treatments ######
        
        poly_snps <- snps_edited_cov_poly_filt2  #just renaming the filtered poly snp df
        
      #1. Remove the "T" in front of the timepoint number so that this can be recognized as an integer for modeling 
        colnames(poly_snps) <-gsub("_T","_",colnames(poly_snps))
        
      #2. Rename the ancestor so it matches the format of the other samples and has the identifier "00" at the end to represent the time 0 
        colnames(poly_snps) <- gsub("af_EEA_anc_4SH_12","af_EEA_ANC_S00_rep00_00",colnames(poly_snps))
        colnames(poly_snps) <- gsub("N_EEA_anc_4SH_12","N_EEA_ANC_S00_rep00_00",colnames(poly_snps))
        
      #3. Ignore the mincov column right now
        poly_snps <- poly_snps[, -ncol(poly_snps)] #-ncol says to ignore the last column
        
      #4. Subset treatments individually to run for single treatment models
        
        #Each individual treatment with ancestor columns
        Low_CAS_only <- poly_snps[, c(1,2,435,436, grep("_CAS_Low_", colnames(poly_snps)))]
        High_CAS_only <- poly_snps[, c(1,2,435,436, grep("_CAS_High_", colnames(poly_snps)))]
        
        Low_CLO_only <- poly_snps[, c(1,2,435,436, grep("_CLO_Low_", colnames(poly_snps)))]
        High_CLO_only <- poly_snps[, c(1,2,435,436, grep("_CLO_High_", colnames(poly_snps)))]
        
        CASCLO_only <- poly_snps[, c(1,2,435,436, grep("_CASCLO_", colnames(poly_snps)))]
        
        Control_only <- poly_snps[, c(1,2,435,436, grep("_Control_", colnames(poly_snps)))]
        
        
      #5. Rename ancestor columns to match the other pops for modeling
        
        #Helper function for individual treatment subsets
        rename_ANC_cols <- function(df, new_tag) {
          colnames(df) <- gsub(
            "ANC_S00",
            new_tag,
            colnames(df))
          df
        }
        
        #Apply function.. Renaming ANC_S00 to the same treatment identifier.
        Low_CAS_only <- rename_ANC_cols(Low_CAS_only, "CAS_Low")
        High_CAS_only <- rename_ANC_cols(High_CAS_only, "CAS_High")
        
        Low_CLO_only <- rename_ANC_cols(Low_CLO_only, "CLO_Low")
        High_CLO_only <- rename_ANC_cols(High_CLO_only, "CLO_High")
        
        CASCLO_only <- rename_ANC_cols(CASCLO_only, "CASCLO_Low")
        
        Control_only <- rename_ANC_cols(Control_only, "Control_ND")
        
        
        
    ##### Step 2. Build SFS for each population at T14 #######
        
        #Set up the function so it does not average the reps (look just at final timepoint T14)
        plot_SFS_replicates_T14 <- function(df, treatment_name, binwidth = 0.02, color_high="white") {
          
          
          #1. Identify the final timepoint T14 allele frequency columns (and exclude ancestor columns)
          cols_T14 <- grep("^af_.*_rep[0-9]{2}_14$", colnames(df), value = TRUE)
          
          #remove ancestor 
          cols_T14 <- cols_T14[!grepl("rep00_00", cols_T14)]
          
          #2. Reshape from wide to long format. Each row = SNP x replicate
          SFS_df <- df %>%
            dplyr::select(all_of(cols_T14)) %>%
            
            pivot_longer(
              cols = everything(),
              names_to = "Replicate",
              values_to = "AF") %>%
            
            #3. Extract metadata from column names
            mutate(
              Timepoint = "T14",
              Replicate = gsub("^.*_(rep[0-9]{2})_14$", "\\1", Replicate)) #etract replicate ID
          
          #4. Summarize fixation per replicate. (So this will need to be done before the fixed SNPs filtering) 
          fixed_stats <- SFS_df %>%
            group_by(Replicate) %>%
            summarise(
              n_total = n(), #total SNP observations in replicate
              n_fixed = sum(AF == 0 | AF == 1, na.rm=TRUE), #number of fixed SNPs
              pct_fixed = 100 * n_fixed / n_total, #percent fixed SNPs
              .groups = "drop")
          
          #5. Ensure replicate order is consistent in plots
          SFS_df$Replicate <- factor(
            SFS_df$Replicate,
            levels = sort(unique(SFS_df$Replicate)))
          
          #6. Build histograms of AF distribution for each population at T14
          p <- ggplot(SFS_df, aes(x = AF)) +
            geom_histogram(binwidth = binwidth,fill = color_high,color = NA) +
            facet_wrap(~Replicate, nrow = 2, ncol = 6, axes = "all_x") +
            #geom_text(data=fixed_stats, aes(x=0.5,y=20000, label = paste0(round(pct_fixed,1), "% fixed")), vjust=1.3,size=5,inherit.aes=FALSE) +
            xlab("SNP Frequency") +
            ylab("Count") +
            scale_y_continuous(limits = c(0,30000),expand = c(0,0)) +
            scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),labels = c("0", "0.25", "0.5", "0.75", "1")) +
            theme_cowplot() +
            ggtitle(paste0(treatment_name, " (T14)")) +
            theme(
              axis.text = element_text(size=14),
              axis.title = element_text(size=16),
              strip.text = element_text(size=14),
              plot.title = element_text(size=16,face="bold"))
          
          return(p)
        }
        
        #Example usage for each treatment
        SFS_plots_Low_CAS <- plot_SFS_replicates_T14(Low_CAS_only,treatment_name = "Low CAS",color_high = "gold")
        
        SFS_plots_High_CAS <- plot_SFS_replicates_T14(High_CAS_only,treatment_name = "High CAS",color_high = "darkgoldenrod")
        
        SFS_plots_Low_CLO <- plot_SFS_replicates_T14(Low_CLO_only,treatment_name = "Low CLO",color_high = "hotpink")
        
        SFS_plots_High_CLO <- plot_SFS_replicates_T14(High_CLO_only,treatment_name = "High CLO",color_high = "firebrick4")
        
        SFS_plots_CASCLO <- plot_SFS_replicates_T14(CASCLO_only,treatment_name = "CASCLO",color_high = "orangered")
        
        SFS_plots_Control <- plot_SFS_replicates_T14(Control_only,treatment_name = "Control",color_high = "grey50")
        
        
        #Combine plots into multipanel figures
        SFS_multi_panel_CAS <- plot_grid(
          SFS_plots_Low_CAS,
          SFS_plots_High_CAS,
          ncol = 1,  
          align = 'v',
          label_size = 30,
          labels = c("A", "B"))
        
        SFS_multi_panel_CLO <- plot_grid(
          SFS_plots_Low_CLO,
          SFS_plots_High_CLO,
          ncol = 1,  
          align = 'v',
          labels = c("C", "D"),
          label_size = 30)
        
        SFS_multi_panel_CASCLO_CTR <- plot_grid(
          SFS_plots_CASCLO,
          SFS_plots_Control,
          ncol = 1,  
          align = 'v',
          label_size = 30,
          labels = c("E", "F"))
        
        
        #Save combined figures
        #ggsave("Figures/SFS/SFS_multi_panel_CAS.png",SFS_multi_panel_CAS,width = 12, height = 10, dpi = 300)
        #ggsave("Figures/SFS/SFS_multi_panel_CLO.png",SFS_multi_panel_CLO,width = 12, height = 10, dpi = 300)
        #ggsave("Figures/SFS/SFS_multi_panel_CASCLO_Control.png",SFS_multi_panel_CASCLO_CTR,width = 12, height = 10, dpi = 300)
        
       
         rm(SFS_plots_CASCLO,SFS_plots_Low_CAS,SFS_plots_High_CAS,SFS_plots_Low_CLO,SFS_plots_High_CLO,SFS_plots_Control)
        
        
        
# ================= IDENTIFYING "FIXED" CANDIDATE SNPS ==============================================     
        
    ##### Step 1. Filter and exclude fixed SNPs ######
        
        #1. Put treatments into a named list
        treatments <- list(
          Low_CAS = Low_CAS_only,
          High_CAS = High_CAS_only,
          Low_CLO = Low_CLO_only,
          High_CLO = High_CLO_only,
          CASCLO  = CASCLO_only,
          Control = Control_only)
        
        
        #2. Function to identify fixed SNPs (ancestor at some freq in a SNP but all evolved pops (all reps, times) == 0 or 1)
        
        get_fixed_snps <- function(df) {
          
          # 1. Identify all allele frequency columns start with _af
          af_cols <- grep("^af_", colnames(df), value = TRUE)
          
          # 2. Identify the ancestor allele frequency column
          ancestor_col <- grep("af_.*rep00_00$", af_cols, value = TRUE)[1]
          
          # 3. Define evolved population columns. Exclude the ancestor column
          evo_cols <- setdiff(af_cols, ancestor_col)
          
          # 4. Define the ancestral state for each SNP.
          ancestor_af <- df[[ancestor_col]] 
          
          # 5. Check fixation status across ALL evolved reps/timepoints
          all_zero <- apply(df[, evo_cols, drop = FALSE], 1, function(x) all(x == 0)) #true if every evolved pop af is exactly 0
          all_one  <- apply(df[, evo_cols, drop = FALSE], 1, function(x) all(x == 1)) #true if every evolved pop af is exactly 1
          
          # 6. Evolved mean AF and shift from ancestor (still useful biologically)
          evo_af <- rowMeans(df[, evo_cols, drop = FALSE], na.rm = TRUE) #mean allele frequency across all evolved reps/timepoints
          delta_af <- evo_af - ancestor_af #difference between evolved mean and ancestor
          
          # 7. Create SNP_ID identifier
          SNP_ID <- paste(df$CHROM, df$POS, sep = "_")
          
          # 8. Initialize fixation state vector (same length as number of SNPs)
          FixationState <- rep(NA_character_, nrow(df))
          FixationState[all_zero] <- "loss"
          FixationState[all_one]  <- "gain"
          
          #9 Fixed vs not fixed flag
          Fixed <- ifelse(!is.na(FixationState), "fixed", "not_fixed")
          
          #10. Output
          out <- data.frame(
            SNP_ID = SNP_ID,
            Fixed = Fixed,
            FixationState = FixationState,
            Ancestor_AF = ancestor_af,
            Evo_AF = evo_af,
            Delta_AF = delta_af
          )
          
          return(out)
        }
        
      
        #Apply function to each treatment and combine results into one dataframe
        fixed_snps <- lapply(names(treatments), function(trt) {
          res <- get_fixed_snps(treatments[[trt]])
          res$Treatment <- trt
          res
        }) %>% dplyr::bind_rows()
        
        
    ##### Step 2. Clean treatment names and count how many fixed #####
        
        treatment_names <- c(
          Low_CAS  = "Low CAS",
          High_CAS = "High CAS",
          Low_CLO  = "Low CLO",
          High_CLO = "High CLO",
          CASCLO   = "CASCLO",
          Control  = "Control")
        
        fixed_snps <- fixed_snps %>%
          dplyr::mutate(Treatment = treatment_names[Treatment]) #overwrite treatment names
        
        #Get the fixed candidates
        all_fixed_snps <- fixed_snps %>%
          dplyr::filter(Fixed == "fixed") %>%
          dplyr::pull(SNP_ID) %>%
          unique() 
        
        length(all_fixed_snps) #2159 fixed candidates total across all treatments
        
        #Count per treatment
        fixed_snps %>%
          dplyr::filter(Fixed == "fixed") %>%
          dplyr::group_by(Treatment, FixationState) %>%
          dplyr::summarise(n = dplyr::n_distinct(SNP_ID), .groups = "drop")

        
    ##### Step 3. Remove fixed snps from each df so that the rest can be used for modeling and export #####
        
        #Remove from all treatment dfs
          Low_CAS_only_filt  <- Low_CAS_only  %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          High_CAS_only_filt <- High_CAS_only %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          Low_CLO_only_filt  <- Low_CLO_only  %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          High_CLO_only_filt <- High_CLO_only %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          CASCLO_only_filt   <- CASCLO_only   %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          Control_only_filt  <- Control_only  %>% dplyr::filter(!paste(CHROM, POS, sep = "_") %in% all_fixed_snps)
          
          
          # write_tsv(Low_CAS_only_filt, "Output_files/Low_CAS_only_filt.txt")
          # write_tsv(High_CAS_only_filt, "Output_files/High_CAS_only_filt.txt")
          # write_tsv(Low_CLO_only_filt, "Output_files/Low_CLO_only_filt.txt")
          # write_tsv(High_CLO_only_filt, "Output_files/High_CLO_only_filt.txt")
          # write_tsv(CASCLO_only_filt, "Output_files/CASCLO_only_filt.txt")
          # write_tsv(Control_only_filt, "Output_files/Control_only_filt.txt")
        
        
        
    ##### Step 4. Further filter "fixed" candidates to make a list for separate analysis (since these cant be modeled) ######
        
        #1. remove any that also fixed in the control, and any that barely shifted from the ancestor starting af
          fixed_snps_filt <- fixed_snps %>%
            dplyr::filter(Fixed == "fixed") %>%
            #Remove any SNP that also fixed in Control
            dplyr::filter(!SNP_ID %in% (fixed_snps %>%
                                          dplyr::filter(Fixed == "fixed", Treatment == "Control") %>%
                                          dplyr::pull(SNP_ID) %>%
                                          unique())) %>%
            # Filter for meaningful shift from ancestor
            dplyr::filter(abs(Delta_AF) > 0.10)
    
        #2.Look at counts and overlap among treatments 
          
          #get total count
            length(unique(fixed_snps_filt$SNP_ID)) #755 total across all treatments
            
          #get count per treatment
            fixed_snps_filt %>%
              dplyr::group_by(Treatment) %>%
              dplyr::summarise(n = dplyr::n_distinct(SNP_ID), .groups = "drop") 
      
          
      #3. Export df with filtered fixed snps list
          #write_tsv(fixed_snps_filt, "Output_files/fixed_snps_filt.txt")
     
          
        

          
# ================= ASSESSING COVERAGE OF CHROM VIII (ERG11 region) FOR ANEUPLOIDY ============================================== 
    ##### Step 1. Set up dataframe with coverages and samples #####
          
          #Reshape the coverage data to long format and normalize read depth within each sample.
          #Normalizing to each sample's own genome-wide median puts all populations on a common
          #scale, so differences in library depth don't drive comparisons across samples.
          
          cov_long <- combined_cov_df %>%
            
            #Keep position info (chromosome, bp position, cumulative Mb) and the coverage columns.
            #Coverage columns all start with "N_" (see earlier formatting steps).
            select(CHROM, POS, MB, starts_with("N_")) %>%
            
            #Pivot from wide to long. Each row becomes one coverage measurement for a single
            #SNP position in a single sample.
            pivot_longer(cols = starts_with("N_"),
                         names_to = "Sample", values_to = "Coverage") %>%
            
            #Pull metadata out of the sample names into their own columns
            mutate(
              Timepoint = str_extract(Sample, "T\\d+"),      # T01, T07, or T14
              Replicate = str_extract(Sample, "rep\\d+"),    # rep01 through rep12
              Treatment = case_when(
                # NOTE: CASCLO must come first. Its sample names contain both "CAS" and "CLO",
                # so if it were placed lower it would be miscalled as a single-drug treatment.
                str_detect(Sample, "CASCLO")                           ~ "CASCLO",
                str_detect(Sample, "CAS") & str_detect(Sample, "Low")  ~ "Low CAS",
                str_detect(Sample, "CAS") & str_detect(Sample, "High") ~ "High CAS",
                str_detect(Sample, "CLO") & str_detect(Sample, "Low")  ~ "Low CLO",
                str_detect(Sample, "CLO") & str_detect(Sample, "High") ~ "High CLO",
                str_detect(Sample, "Control")                          ~ "Control")) %>%
            
            #Drop anything that didn't match a treatment pattern (e.g. the ancestor column),
            #since case_when returns NA for unmatched rows.
            filter(!is.na(Treatment)) %>%
            
            #Calculate the median coverage across all SNPs within each sample. This is the
            #per-sample normalization factor.
            group_by(Sample) %>%
            mutate(sample_median = median(Coverage, na.rm = TRUE)) %>%
            ungroup() %>%
            
            #Relative coverage: a value of 1 means the position is covered at the sample's
            mutate(rel_cov = Coverage / sample_median)
          
          
    ##### Step 2. Plot normalized coverage for chrom VIII #####    
          
          
          
          #Panel A: chr VIII per replicate 
          chromVIII_panelA_plot <- cov_long %>%
            filter(CHROM == "VIII") %>%                                          #ERG11 is on chr VIII
            group_by(Treatment, Replicate, Timepoint) %>%
            summarise(rel_cov = median(rel_cov, na.rm = TRUE), .groups = "drop") %>%  #one value per replicate per timepoint
            mutate(Treatment = factor(Treatment,
                                      levels = c("Low CAS","High CAS","Low CLO","High CLO","CASCLO","Control"))) %>%  #set facet order
            ggplot(aes(Timepoint, rel_cov, group = Replicate, color = Treatment)) +
            geom_hline(yintercept = c(1, 1.5), linetype = "dashed", color = "grey60") +  #1 = ancestral, 1.5 = one extra copy (diploid)
            geom_line(alpha = 0.7) + geom_point(size = 1.5) +                    #trajectory for each of the 12 replicates
            facet_wrap(~Treatment, nrow = 2) +                                   #one panel per treatment
            scale_color_manual(values = c("gold","darkgoldenrod","hotpink",
                                          "firebrick4","orangered","grey50")) +  #treatment palette used throughout
            scale_y_continuous(limits = c(0.9, 1.6), breaks = c(1.0, 1.25, 1.5)) +  #zoom to the informative range
            labs(x = NULL, y = "Normalized coverage (chr VIII)") +
            theme(legend.position = "none",  
                  axis.title = element_text(size = 15),
                  axis.text = element_text(size = 13),
                  strip.text = element_text(size = 13, face = "bold"))                 
          
          #Panel B: chr VIII along its length, High CLO vs Control
          chromVIII_panelB_plot <- cov_long %>%
            filter(CHROM == "VIII",
                   Treatment %in% c("Control", "Low CLO", "High CLO", "CASCLO"), #CLO comparisons only
                   POS > 20000, POS < 540000) %>%                                #trim subtelomeric mapping artifacts
            mutate(window = floor(POS / 20000) * 20000,                          #bin positions into 20 kb windows
                   Treatment = factor(Treatment,
                                      levels = c("Control", "Low CLO", "High CLO", "CASCLO"))) %>%  #set facet order
            group_by(Treatment, Timepoint, window) %>%
            summarise(rel_cov = median(rel_cov, na.rm = TRUE), .groups = "drop") %>%  #median across replicates per window
            ggplot(aes(window / 1000, rel_cov, color = Treatment,                #convert bp to kb for the x-axis
                       alpha = Timepoint, linewidth = Timepoint)) +              #timepoint shown by transparency and weight
            annotate("rect", xmin = 120.091, xmax = 121.683, ymin = -Inf, ymax = Inf,alpha = 0.25, fill = "blue") +
            annotate("text", x = 120.1, y = 1.63, label = "ERG11 region",size = 4, hjust = -0.15, color = "blue") +
            geom_line() +
            facet_wrap(~Treatment, nrow = 1) +
            scale_color_manual(values = c("Control"  = "grey50",
                                          "Low CLO"  = "hotpink",
                                          "High CLO" = "firebrick4",
                                          "CASCLO"   = "orangered")) +
            scale_alpha_manual(values = c("T01" = 0.35, "T07" = 0.65, "T14" = 1)) +      #later timepoints more opaque
            scale_discrete_manual("linewidth",
                                  values = c("T01" = 0.7, "T07" = 1.0, "T14" = 1.4)) +   #and drawn thicker
            scale_y_continuous(limits = c(0.8, 1.7)) +
            labs(x = "Position on chromosome VIII (kb)", y = "Normalized coverage") +
            guides(color = "none",                                               #treatment shown by facet, not legend
                   alpha = guide_legend(title = "Timepoint",
                                        override.aes = list(linewidth = 1.2, color = "black")),  #readable legend keys
                   linewidth = "none") +                                         #avoid duplicate timepoint legend
            theme(legend.position = "right",
                  axis.title = element_text(size = 15),
                  axis.text = element_text(size = 13),
                  strip.text = element_text(size = 13, face = "bold"),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 13))
          
  
          
          #Stack the two panels and label A and B
          chromVIII_panelA_plot / chromVIII_panelB_plot + plot_annotation(tag_levels = "A") &
            theme(plot.tag = element_text(size = 30, face = "bold"))
          
          #ggsave("Figures/Coverage/chrVIII_coverage.svg", width = 12, height = 8, dpi = 300)
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
        
          
          
# ================= IDENTIFYING DE NOVO CANDIDATES ==============================================


          
    ##### Step 1. Filtering for candidate de novo mutations and plotting ###############
  
  ## This is code to apply across ALL populations (NOT treatment specific)  
    #1. Make a copy of the coverage filtered SNPs
#         snps_de_novos <- snps_edited_cov_filt_10X
#      
#     #2. Define columns for different groups of samples
#         ancestor_col <- "af_EEA_anc_4SH_12" #store the ancestor af column
#             
#         founder_cols <- grep("^af_EEA_hap", colnames(snps_de_novos), value = TRUE) #store the founder af columns
#         
#         exp_pop_cols <- setdiff(grep("^af_", colnames(snps_de_novos), value = TRUE), c(ancestor_col, founder_cols)) #store the experimental pops af columns and exclude the ancestor and founder
#         
#     #3. Calculate some useful metrics on experimental populations
#         snps_de_novos <- snps_de_novos %>%
#           rowwise() %>%
#           mutate(
#             SNP_max = max(c_across(all_of(c(exp_pop_cols, ancestor_col))), na.rm = TRUE), #calculate max AF across exp pops and the ancestor
#             SNP_var = var(c_across(all_of(exp_pop_cols)), na.rm = TRUE) #calculate variance of AF across exp pops only. High variance = polymorphic. Low variance = sequencing noise or drift or consistently fixed
#           ) %>%
#           ungroup() %>%
#           
#           
#     #4. Keep SNPs fixed in the ancestor, we want to look for new mutations so exclude sites polymorphic in the ancestor
#     
#         filter(!!sym(ancestor_col) == 0 | !!sym(ancestor_col) == 1) %>% #!!sym converts the string in to a column reference. Keeps only SNPs where ancestor AF is exactly 0 or 1
#         
#     #5. Remove SNPs that are polymorphic in all populations
#         #Keeps SNPs that are fixed in some pops and polymorphic in others
#         
#         filter(!if_all(all_of(exp_pop_cols), ~ .x > 0 & .x < 1)) %>% #removes snps where every evolved pop has a freq between 0 and 1. At least one pop is fixed. 
#         
#     #6. Calculate the maximum AF change from the ancestor
#           
#         mutate(maxdiff = SNP_max - !!sym(ancestor_col)) %>% #difference between the highest allele freq observed and the ancestors AF. If maxdiff is large, it suggests that the AF changed dramatically in at least one population (stronger evidence of selection rather than just noise)
#         
#     #7. Keep SNPs with strong allele frequency shifts
#           
#         filter(maxdiff > 0.3) %>% #filters out snps that changed <30% in frequency compared to the ancestor. This is an empirical threshold. These are ones that change big in one or more pops since others could be drift or noise
#           
#     #8. Keep only SNPs with substantial variance
#         
#         filter(SNP_var > 0.1) #Removes SNPs where variance across exp pops is very low (needs to be greater than 5%). Low variance means the snp is the same across all pops so likely just an artifact.  
#         
#         
#     ## This chunk gives candidate de novo mutations that are likely the result of selection not standing variation or noise based on the criteria:
#         #SNPs that are fixed in the ancestor
#         #Remove any that are polymorphic across all experimental populations 
#         #with a large AF shift (> 0.3) relative to the ancestor
#         #and meaningful variance across pops (>0.05). 
#         
#         #The following chunk of code was used to see intermediates and troubleshoot thresholds to use. 
#           intermediate <- snps_edited_cov_filt_10X %>%
#             filter((!!sym(ancestor_col) == 0) | (!!sym(ancestor_col) == 1)) %>%
#             filter(!if_all(all_of(exp_pop_cols), ~ .x > 0 & .x < 1)) %>%
#           
#             rowwise() %>%
#             mutate(
#               SNP_max = max(c_across(all_of(c(exp_pop_cols, ancestor_col))), na.rm = TRUE),
#               SNP_var = var(c_across(all_of(exp_pop_cols)), na.rm = TRUE)
#             ) %>%
#             ungroup() %>%
#           
#             mutate(maxdiff = SNP_max - !!sym(ancestor_col)) #before applying the maxdiff and var filters 30,217 candidates. 
#         
#           summary(intermediate$maxdiff)
#           hist(intermediate$maxdiff, breaks = 50) #most fall below 0.2
#         
#           summary(intermediate$SNP_var)
#           hist(intermediate$SNP_var, breaks = 50) #most fall below 0.2
# 
#           test_thres <- intermediate %>%
#             filter(maxdiff > 0.3) %>% #just this filter will give 7K candidates but then goes to 0 after applying the next filter
#             filter(SNP_var > 0.1)
#           
#         
#         
#         
#     #9. Export this filtered SNP table of candidate de novos for downstream analysis
#         #write.table(snps_de_novos,"snps_filtered_denovos.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")
#         
#     
#           
#   ##This is code to look for treatment-specific de novos. Since the attempt above may have been too restrictive and these drugs could lead to differences in likelihood of de novos.
#     ##Note: this code will look separately at each drug so we can look independently to identify dose-specific or combo specific adaptive mutations. 
#      
#                
#     #1. Define a helper function to process each treatment group
#           ancestor_col <- "af_EEA_anc_4SH_12"
#           
#           process_de_novos <- function(data, treatment_prefix, ancestor_col, maxdiff_thresh = 0.3, var_thresh = 0.1) {
#             
#             af_cols <- grep(paste0("^af_EEA_", treatment_prefix, "(_|$)"), colnames(data), value = TRUE) #Get allele frequency columns for this treatment group (assumes 'af_' prefix)
#             
#             cov_cols <- grep(paste0("^N_EEA_", treatment_prefix, "(_|$)"), colnames(data), value = TRUE) #Get corresponding coverage columns for the treatment group 
#             
#             treatment_cols <- colnames(data)[colnames(data) %in% c(af_cols, cov_cols)] #Combine AF and coverage columns, preserve original order as they appear in the data
#             
#             
#             
#             data %>%
#               rowwise() %>%
#               mutate(
#                 SNP_max = max(c_across(all_of(c(af_cols, ancestor_col))), na.rm = TRUE),
#                 SNP_var = var(c_across(all_of(af_cols)), na.rm = TRUE)
#               ) %>%  #Calculate max AF across dose replicates + ancestor, and variance across dose replicates only
#               ungroup() %>%
#               
#               filter((!!sym(ancestor_col) == 0) | (!!sym(ancestor_col) == 1)) %>% #keep SNPs fixed in ancestor
#               
#               filter(!if_all(all_of(af_cols), ~ .x > 0 & .x < 1)) %>% #Remove snps polymorphic in all populations
#               
#               mutate(maxdiff = SNP_max - !!sym(ancestor_col)) %>% #Calculate max difference from ancestor AF
#               
#               filter(maxdiff > maxdiff_thresh) %>% #filters out snps that changed <30% in frequency compared to the ancestor.
#               filter(SNP_var > var_thresh) %>% #Removes SNPs where variance across exp pops is very low (needs to be greater than 10%). 
#               
#               select(any_of(c("CHROM", "POS")), #now select just the columns I want to keep, otherwise it will keep all treatment groups for each outputted df
#                      all_of(treatment_cols),
#                      SNP_max,SNP_var,maxdiff)
#             }
#           
#     #2. Now run this function for each treatment group:
#           cas_s01_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CAP_S01", ancestor_col) #4 candidates de novos based on above filters
#           cas_s03_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CAP_S03", ancestor_col) #22 candidate de novos
#           
#           clo_s01_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CLO_S01", ancestor_col) #0 candidate de novos
#           clo_s03_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CLO_S03", ancestor_col) #7 candidate de novos
#           
#           cac_s01_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CAC_S01", ancestor_col) #13 candidate de novos
#           
#           ctr_s01_denovos <- process_de_novos(snps_edited_cov_filt_10X, "CTR_S01", ancestor_col) #0 candidate de novos
#           
#         
#           
#   
#     #3. Add a column identifier to be used for upset plot
#           
#           denovos_dfs <- c("cas_s01_denovos", "cas_s03_denovos", "clo_s01_denovos", 
#                 "clo_s03_denovos", "cac_s01_denovos", "ctr_s01_denovos") #store all the dfs
# 
#           for (df_name in denovos_dfs) {
#             df <- get(df_name)
#             df <- df %>% mutate(snp_id = paste0(CHROM, "_", POS))
#             assign(df_name, df)
#           } #this will make a new column in each df called snp_id which is the CHROM_POS. Output will be in original denovos df. 
#           
#           
#     #4. Extract snp_id vectors for each treatment group
#           
#           snp_list <- list(
#             CAP_S01 = cas_s01_denovos$snp_id,
#             CAP_S03 = cas_s03_denovos$snp_id,
#             CLO_S01 = clo_s01_denovos$snp_id,
#             CLO_S03 = clo_s03_denovos$snp_id,
#             CAC_S01 = cac_s01_denovos$snp_id,
#             CTR_S01 = ctr_s01_denovos$snp_id)
#    
#       
#   ##Building the upset plot for de novo candidate SNPs       
#           
#     #1. Get all unique real SNPs
#           denovo_snps <- unique(unlist(snp_list, use.names = FALSE))
#           denovo_snps <- denovo_snps[!is.na(denovo_snps)]
#           
#     #2. Build a logical matrix: rows = SNPs, columns = samples
#           denovo_snps_matrix <- sapply(snp_list, function(snps){
#             snps <- snps[!is.na(snps)]
#             denovo_snps %in% snps
#           })
#           
#           #Convert to tibble, keeping SNP names as a column
#             denovo_snps_df <- as_tibble(denovo_snps_matrix)
#             denovo_snps_df <- denovo_snps_df %>%
#              mutate(SNP = denovo_snps) %>%     # add SNP names as a column
#               select(SNP, everything()) 
#           
#           
#     #3. Identify sample columns and set your order
#           sample_cols <- c("CAP_S01", "CAP_S03", "CLO_S01", "CLO_S03", "CAC_S01", "CTR_S01")
#           
#           
#     #4. Create a comma-separated combination column
#           denovo_snps_df <- denovo_snps_df %>%
#             mutate(combination = pmap_chr(
#               select(., all_of(sample_cols)),
#               function(...) {
#                 vals <- list(...)
#                 present_samples <- sample_cols[unlist(vals)]
#                 if(length(present_samples) == 0) "" else paste(present_samples, collapse = ", ")
#               }
#             )
#             )
#           
#           
#     #5. Suppose you want to order combinations by decreasing frequency:
#           comb_order <- denovo_snps_df %>%
#             dplyr::count(combination) %>%          # Count how many times each combination occurs
#             arrange(desc(n)) %>%            # Sort combinations by decreasing frequency
#             pull(combination)               # Extract the combination names as a vector
#           
#           #Then convert 'combination' to a factor with that order
#           denovo_snps_df <- denovo_snps_df %>%
#             mutate(combination = factor(
#               combination,                 # The column to convert
#               levels = comb_order          # Use the order we just defined
#             )
#             )
#           
#           
#     #6. Building a basic upset plot 
#           
#         #At this point, we can make a basic upset plot minus the horizontal bars with this code:
#           
#           #Plot with ggupset
#           denovo_snps_df %>%
#             ggplot(aes(x = combination)) +          # Map combinations to the x-axis
#             geom_bar(stat = "count") +             # Count occurrences of each combination
#             axis_combmatrix(
#               sep = ", ",                           # Separator used in the combination strings
#               levels = sample_cols                   # The order of samples on the axes
#             ) +
#             theme(
#               axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
#             )
#           
#           ### However, if we want to add the bars and make a fully customizable plot, we can make each element separately 
#           # and then patch them together with patchwork. 
#           # We just need to make some summary data frame that count how many instances we have of each SNP in each treatment
#           # 1. We will make a count of the combinations of treatments for each snp (counts_combinations)
#           # 2. and, vice versa, a count of the number of SNPs associated with each treatment 
#           # Basically have to do this because bar charts need raw totals for groupings to plot WITHOUT having to rely on stat = "count" like above
#           # Avoiding stat = "count" by manually counting things gives us more control to make the plot pretty. 
#           # Then we can make the three elements we want: the bar chart, the upset dumbbells, and the horizontal chart
#          
#            
#     #7. Building a more sophisticated upset plot
#           
#         #Creating initial summary data frames
#           
#           #Create a new dataframe where we count the combinations for each treatment
#           counts_combinations <- denovo_snps_df |>  
#             #Create a combination column: which treatments have this SNP
#             mutate(combination = pmap_chr( # pmap_chr will iterate row-wise across these 6 columns in snp_df:
#               list(CAP_S01, CAP_S03, CLO_S01, CLO_S03, CAC_S01, CTR_S01),
#               \(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6) { # For each row, assign the values in the 6 columns to lgl1–lgl6
#                 # Make a vector of group labels, one for each sample column
#                 # Use the logical values from this row (TRUE/FALSE) 
#                 # to select which groups are active (TRUE = keep, FALSE = drop)
#                 c('CAS_Low', 'CAS_High', 'CLO_Low', 'CLO_High', 'CASCLO', 'Control')[c(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6)] |> 
#                   paste(collapse = ',')
#               }
#             )
#             ) |> 
#             #Count how many SNPs fall into each combination
#             dplyr::count(combination) |> 
#             #Ensure that each individual treatment is included even if it has no SNPs
#             complete(combination = c('CAS_Low', 'CAS_High', 'CLO_Low', 'CLO_High', 'CASCLO', 'Control'), fill = list(n = 0))
#           
#           
#         #To make sure that our bar charts are sorted later on in the order of largest to smallest counts, 
#           #we convert the combination column into a factor.
#           counts_combinations <- counts_combinations |> 
#             mutate(combination = fct_reorder(combination, n, .desc = TRUE))
#           
#         #Similarly, we need to count how many SNPs are TRUE in each treatment 
#           #This requires first reordering the data so that the treatments are in a single column.
#           #Then we can count the number of SNPs in each treatment 
#           #Notice that this uses the fact that a TRUE value is treated as 1 and a FALSE value as 0 when we sum them up.
#           #Send it to denovo_snps_df2 so we can preserve denovo_snps_df if we want to use it for other plotting
#           
#           denovo_snps_df |> 
#             # Remove the 'combination' column if it exists
#             select(-combination) |>
#             pivot_longer(
#               cols = -SNP,
#               names_to = 'Treatment',
#               values_to = 'Present'
#             ) |> 
#             # Rename treatments after we pivot longer
#             mutate(Treatment = recode(
#               Treatment,
#               "CAP_S01" = "CAS_Low",
#               "CAP_S03" = "CAS_High",
#               "CLO_S01" = "CLO_Low",
#               "CLO_S03" = "CLO_High",
#               "CAC_S01" = "CASCLO",
#               "CTR_S01" = "Control"
#             )) |> 
#             summarize(
#               counts = sum(Present),
#               .by = Treatment,
#             ) -> denovo_snps_df2
#           
#         #Set color fill palette
#           my_pal = c("gold","darkgoldenrod","hotpink", "firebrick4", "orangered","black")
#           my_pal_reversed = c("orangered", "firebrick4", "darkgoldenrod", "gold")
#           
#           
#         #Make main bar chart
#           
#           #NOTE: We remove the grid expansion with coord_cartesian because that will later get in our way when we assemble the subplots.
#             bar_chart <- counts_combinations  |> 
#               filter(!combination %in% c("CLO_Low", "Control")) |>  # remove unwanted treatments
#               ggplot(aes(x = combination, y = n)) +
#               geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
#               geom_col(width = 0.6, fill = 'grey25') +
#               coord_cartesian(ylim = c(0,17 * 1.1), expand = FALSE) +
#               labs(x = element_blank(), y = element_blank())
#             bar_chart
#           
#           
#         #Next, we deal with the points. To do so, we need to split up the combinations into individual treatments
#           points_data <- counts_combinations |> 
#             mutate(Treatment = map(combination, 
#                                    ~str_split_1(as.character(.), ','))) |> 
#             unnest(Treatment) |> 
#             # Transform to factors with explicit order
#             mutate( Treatment = factor(Treatment, 
#                                        levels = c("Control", "CASCLO", "CLO_High", "CLO_Low", "CAS_High", "CAS_Low")  # specify desired order for plotting
#             ),
#             combination = factor(combination,levels = levels(counts_combinations$combination)
#             )
#             ) |> 
#             filter(!is.na(Treatment))  # Remove any missing values
#           points_data
#           
#           
#         #Make the dumbbell chart 
#           
#           point_chart <- points_data |> 
#             ggplot(aes(y = Treatment)) +  # use full Treatment factor for y-axis
#             #plot only rows with valid x-values
#             geom_line(data = subset(points_data, !Treatment %in% c("CLO_Low", "Control")),
#                       aes(x = combination, group = combination)) +
#             geom_point(data = subset(points_data, !Treatment %in% c("CLO_Low", "Control")),
#                        aes(x = combination, col = Treatment),size = 7) +
#             scale_color_manual(values = my_pal_reversed) +
#             scale_y_discrete(drop = FALSE) +  # keep all y-axis levels
#             theme(axis.text.y = element_text(hjust = 0.5),legend.position = "none") +
#             labs(x = element_blank(), y = element_blank()) +
#             scale_x_discrete(drop = TRUE)  # x-axis only includes plotted combinations
#           
#           point_chart
#           
# 
#         #create the bar charts for each of the three subjects.  
#           #Reorder the bars by the counts and use the negative value of the counts so that the bars go to the left.
#           treatment_total_bars <- denovo_snps_df2 |> 
#             mutate(Treatment = factor(Treatment, levels = c("Control", "CASCLO", "CLO_High", "CLO_Low", "CAS_High", "CAS_Low"))) |>
#             ggplot(aes(x = counts, y = Treatment)) +
#             geom_col(width = 0.6, fill = my_pal) +
#             scale_x_reverse(expand = c(0,0))+
#             coord_cartesian(expand = FALSE) +
#             labs(x = element_blank(), y = element_blank())
#           treatment_total_bars
#           
#           
#           
#         #Now we can put everything together.
#           #We’ll use the patchwork package for that. All we have to do is to specify the layout and then add the plots together.
# layout <- '
# ##AAAA
# BBCCCC'
#           
#           
#         #Patch each above the plots together
#           upset_plot <- (bar_chart + scale_x_discrete(labels = NULL)) + 
#             (treatment_total_bars + scale_y_discrete(labels = NULL) +
#                theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())) + 
#             (point_chart + 
#                scale_x_discrete(labels = NULL,drop = TRUE,expand = expansion(add = 0.3)) + # Add a bit of space 
#                scale_y_discrete(expand = expansion(add = 0.3),drop = FALSE, labels = c("Control" = "Control",
#                                                                                        "CASCLO" = "CASCLO",
#                                                                                        "CLO_High" = "High CLO",
#                                                                                        "CLO_Low" = "Low CLO",
#                                                                                        "CAS_High" = "High CAS",
#                                                                                        "CAS_Low" = "Low CAS"))) + 
#             plot_layout(design = layout) 
#           #plot_annotation(title = 'Title here', caption = 'This data is really cool. Wow science.',
#           #theme = theme(title = element_text(size = rel(2)),
#           #plot.title = element_text(face = 'bold', lineheight = 1.1),
#           #plot.caption = element_text(size = rel(0.5))))
#           upset_plot
#           
#         #Export graph
#           # jpeg(filename = "Denovo_snps_upset_plot.jpg", units= "in", width = 9, height = 6, res = 300)
#           # upset_plot
#           # dev.off()
#           
#        
#         ## Based on the plot above we see the following:
#           # 1 denovo variant occurs in both CAS_Low, CAS_High and CASCLO suggesting that they are not exclusive to one condition. 
#           # 4 denovo variant is shared between CAS_High and CASCLO
#           # 4 denovo variant is shared between CLO_High and CASCLO
#           # CAS_High had the most unique denovo variants. 
#           
#           #keep in mind that this analysis is not broken down by replicate. Meaning they needed to appear in at least one replicate or timepoint but likely aren't present in all replicates. 
#           #So the one shared between cas_low and cas_high hints these mutations might confer a general advantage under drug exposure, regardless of strength,
#           #but, the higher the strength the more likely for denovo mutations to occur. 
#           #some denovos unique to only high drug or low drug treatments may represent adaptations specific to the drug concentration or stress level. 
#         
          
         
          
      

      
      
       
          
          
          
      
          
    