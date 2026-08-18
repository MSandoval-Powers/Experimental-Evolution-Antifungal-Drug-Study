########################################
# Script Name: GLM_Script.R
# Description: Statistical analysis of frequency change from polymorphic SNP candidates
# Author: Megan Sandoval-Powers
# Last Updated: 8-2026
# R Version: 4.5.1
########################################


# =============== R PACKAGES ======================================= 
library(dplyr) #for the select function
library(MuMIn) #for multi model inference and to get the AICc statistic to compare models
library(data.table) #to read and manipulate data tables
library(tidyverse) #Data wrangling
library(reshape2) #To melt things
library(emmeans) #Pairwise comparisons of lm results with post hoc corrections. 
library(pbapply) #to look at progress of functions
library(purrr) #reduce function


# =============== DATAFILES ======================================= 


#read in datatables containing filtered polymorphic SNP candidates for each treatment
Low_CAS_only_filt=read.table("Input_files/Low_CAS_only_filt.txt",header=T)
High_CAS_only_filt=read.table("Input_files/High_CAS_only_filt.txt",header=T)
Low_CLO_only_filt=read.table("Input_files/Low_CLO_only_filt.txt",header=T)
High_CLO_only_filt=read.table("Input_files/High_CLO_only_filt.txt",header=T)
CASCLO_only_filt=read.table("Input_files/CASCLO_only_filt.txt",header=T)
Control_only_filt=read.table("Input_files/Control_only_filt.txt",header=T)



# =============== SETTING UP DATA AND INDICES FOR MODELING ======================================= 

  ### Step 1. Combine dataframes for cross-treatment models ###
  
    #1. Combine dataframes as desired and set up column names
          
        #For CAS dose cross-treatment model
          #Combine Low CAS, High CAS, and Control
          Low_High_CAS_comb <- Low_CAS_only_filt %>%
            inner_join(High_CAS_only_filt, by = c("CHROM","POS")) %>%
            inner_join(Control_only_filt, by = c("CHROM","POS"))
          
          #Rename the control column to match drug name (so model recognizes dose)
          colnames(Low_High_CAS_comb) <- gsub("_Control_ND_","_CAS_ND_",colnames(Low_High_CAS_comb))
          
          
        #For CLO dose cross-treatment model
          #Combine Low CLO, High CLO, and Control
          Low_High_CLO_comb <- Low_CLO_only_filt %>%
            inner_join(High_CLO_only_filt, by = c("CHROM","POS")) %>%
            inner_join(Control_only_filt, by = c("CHROM","POS"))
          
          #Rename the control column to match drug name (so model recognizes dose)
          colnames(Low_High_CLO_comb) <- gsub("_Control_ND_","_CLO_ND_",colnames(Low_High_CLO_comb))
          
          
        #For drug single vs combined drug model 
          #Combine High CAS, High CLO, CASCLO, and Control
          CAS_CLO_CASCLO_comb <- High_CAS_only_filt %>%
            inner_join(High_CLO_only_filt, by = c("CHROM","POS")) %>%
            inner_join(CASCLO_only_filt, by = c("CHROM","POS")) %>%
            inner_join(Control_only_filt, by = c("CHROM","POS"))
          
    
  ### Step 2. Setting up indices for modeling ###        
    
    #Now that the dataframes are set up, lets make an index that the model will work on. 
    #This is so that we don't have to "melt" the dataframe vertically which would be messy and huge but the model will still know which columns and rows to act on
          
    #1. Build the function to go through columns and store their position and identifiers
        index_funx1 <- function(df) {
          index_af <- c(seq(3,ncol(df),by=2)) #store just the AF values
          index_cov <- c(seq(4,ncol(df),by=2)) #store just the cov values
          
          id <- colnames(df)[index_af] #store just the column names which represent the population ID
          
          treatment <- sapply(id, function(x) unlist(strsplit(x,"_"))[3]) #store just the drug treatment portion of the population ID
          dose <- sapply(id, function(x) unlist(strsplit(x,"_"))[4]) #store just the dose treatment portion of the population ID
          rep <- sapply(id, function(x) unlist(strsplit(x,"_"))[5]) #store just the replicate portion of the population ID
          time <- sapply(id, function(x) unlist(strsplit(x, "_"))[6]) #store just the timepoint 
          
          output <- data.frame(
            idaf = index_af,
            idcov = index_cov,
            treatment = as.factor(treatment),
            dose = as.factor(dose),
            rep = as.factor(rep),
            time = as.numeric(time))
          
          return(output)
        }
        
          
    #2. Run the function on the dataframes to build indices
          
      ### For individual treatment models ###
          Low_CAS_only_filt_index <- index_funx1(Low_CAS_only_filt)
          High_CAS_only_filt_index <- index_funx1(High_CAS_only_filt)
          Low_CLO_only_filt_index <- index_funx1(Low_CLO_only_filt)
          High_CLO_only_filt_index <- index_funx1(High_CLO_only_filt)
          CASCLO_only_filt_index <- index_funx1(CASCLO_only_filt)
          Control_only_filt_index <- index_funx1(Control_only_filt)
          
      ### For cross-treatment models ###
          Low_High_CAS_comb_index <- index_funx1(Low_High_CAS_comb)
          
          Low_High_CLO_comb_index <- index_funx1(Low_High_CLO_comb)
          
          CAS_CLO_CASCLO_comb_index <- index_funx1(CAS_CLO_CASCLO_comb)
          
          
      #Remove rep prefix in the rep column so its just identified as a number
        #Build the function
          remove_rep_prefix <- function(df) {
            df$rep <- factor(
              gsub("^rep(\\d{2})$", "\\1", as.character(df$rep)),
              levels = unique(gsub("^rep(\\d{2})$", "\\1", as.character(df$rep)))
            )
            df
          }
          
        #Run the function
          Low_CAS_only_filt_index <- remove_rep_prefix(Low_CAS_only_filt_index)
          High_CAS_only_filt_index <- remove_rep_prefix(High_CAS_only_filt_index)
          
          Low_CLO_only_filt_index <- remove_rep_prefix(Low_CLO_only_filt_index)
          High_CLO_only_filt_index <- remove_rep_prefix(High_CLO_only_filt_index)
          
          CASCLO_only_filt_index <- remove_rep_prefix(CASCLO_only_filt_index)
          Control_only_filt_index <- remove_rep_prefix(Control_only_filt_index)
          
          Low_High_CAS_comb_index <- remove_rep_prefix(Low_High_CAS_comb_index)
          
          Low_High_CLO_comb_index <- remove_rep_prefix(Low_High_CLO_comb_index)
          
          CAS_CLO_CASCLO_comb_index <- remove_rep_prefix(CAS_CLO_CASCLO_comb_index)
       
           
          
          
# =============== TREATMENT-SPECIFIC MODELS ======================================= 

          
  ######### Step 1. Run each model on individual treatments #######
    
      #1. Function for modeling
        funx.model1 <- function(snp_df, index_df) { #snp_df = ONE SNP(one row), passed in from pbapply later; index_df = FULL index table
  
          #Build modeling dataframe. Rows come from index_df, AF and cov values come from a single SNP each time
          freq <- cbind(as.numeric(snp_df[index_df$idaf]),
                        as.numeric(snp_df[index_df$idcov]),
                        index_df) #freq is basically pulling the af and cov data from the snp table and binding it to the index table, one snp at a time
          colnames(freq)[1:2] <- c("af","cov")
  
          #fit SNP-specific model
          mod1 <- glm(af~time,family="quasibinomial",weights=cov,data=freq)
  
          #Extract per SNP statistics
          beta <- coef(mod1)["time"]
          pval <- summary(mod1)$coefficients["time", "Pr(>|t|)"]
  
          #Name columns based on drug
          treatment_name <- unique(index_df$treatment)
          dose_name <- unique(index_df$dose)
          col_label <- paste0(treatment_name, "_", dose_name)
          names(beta) <- paste0("beta_", col_label)
          names(pval) <- paste0("pval_", col_label)
  
          #Return one row of results per SNP
          model_results <- as.data.frame(cbind(t(pval),t(beta))) #transpose the pvalue objects then bind in a single dataframe
  
        }
            
            
            
            # funx.model1 <- function(snp_df, index_df) {
            #   
            #   freq <- cbind(as.numeric(snp_df[index_df$idaf]),
            #                 as.numeric(snp_df[index_df$idcov]),
            #                 index_df)
            #   colnames(freq)[1:2] <- c("af","cov")
            #   
            #   # container for warnings
            #   warn_msg <- NULL
            #   
            #   mod1 <- withCallingHandlers(
            #     glm(af ~ time, family = "quasibinomial", weights = cov, data = freq),
            #     warning = function(w) {
            #       warn_msg <<- c(warn_msg, conditionMessage(w))
            #       invokeRestart("muffleWarning")
            #     }
            #   )
            #   
            #   beta <- coef(mod1)["time"]
            #   pval <- summary(mod1)$coefficients["time", "Pr(>|t|)"]
            #   
            #   treatment_name <- unique(index_df$treatment)
            #   dose_name <- unique(index_df$dose)
            #   col_label <- paste0(treatment_name, "_", dose_name)
            #   
            #   names(beta) <- paste0("beta_", col_label)
            #   names(pval) <- paste0("pval_", col_label)
            #   
            #   # collapse warnings into one string
            #   warn_msg <- if (is.null(warn_msg)) NA else paste(unique(warn_msg), collapse = " | ")
            #   
            #   model_results <- data.frame(
            #     t(pval),
            #     t(beta),
            #     warning = warn_msg
            #   )
            #   
            #   return(model_results)
            # }
      
      
      #2. Chunk level orchestrator (reduce peak memory usage)
        funx.model1_chunk <- function(df_treat, index_df, chunk_size = 1000) { ##df_treat = FULL SNP table, index_df = FULL index table (same as above)
          
          n_rows <- nrow(df_treat)
          chunk_starts <- seq(1, n_rows, by = chunk_size)
          results <- list()
          
          for (start in chunk_starts) {
            end <- min(start + chunk_size - 1, n_rows)
            cat("Processing rows", start, "to", end, "\n")
            
            chunk <- df_treat[start:end, , drop = FALSE]
            
            #Apply per SNP model row-by-row
            chunk_result <- pbapply::pbapply(chunk, 1,function(x) funx.model1(x, index_df)) #build its own freq df, fit its own glm per snp, run its own emtrends
          
          results[[length(results) + 1]] <- chunk_result
          gc()
        }
        
        #Combine all SNP results
        model_results <- do.call(c, results)
        model_results <- as.data.frame(do.call(rbind, model_results))
        
        #Add CHROM/POS columns
        model_results <- cbind(df_treat[,1:2], model_results)
        return(model_results)
      }
       
      
      #3. Run the models per treatment
        model.res.Low_CAS <- funx.model1_chunk(Low_CAS_only_filt, Low_CAS_only_filt_index)
        model.res.High_CAS <- funx.model1_chunk(High_CAS_only_filt, High_CAS_only_filt_index)
        
        model.res.Low_CLO <- funx.model1_chunk(Low_CLO_only_filt, Low_CLO_only_filt_index)
        model.res.High_CLO <- funx.model1_chunk(High_CLO_only_filt, High_CLO_only_filt_index)
        
        model.res.CASCLO <- funx.model1_chunk(CASCLO_only_filt, CASCLO_only_filt_index)
        
        model.res.Control <- funx.model1_chunk(Control_only_filt, Control_only_filt_index)
        
  
  ######### Step 2. P-value adjust for multiple testing #######      
       
    #1. Build function to fdr adjust the pvalues from each model results dataframe created above
        
        #this will fdr adjust the pvalues, and replace that column
        adjust_pvals_overwrite <- function(df) {
          # Find all columns starting with "pval_"
          pval_cols <- grep("^pval_", colnames(df), value = TRUE)
          
          for (col in pval_cols) {
            # Make the new column name
            new_col <- sub("^pval_", "pval_adj_", col)
            
            # Adjust p-values and assign to new column
            df[[new_col]] <- p.adjust(df[[col]], method = "BH")
            
            # Remove the old column
            df[[col]] <- NULL
          }
          
          return(df)
        }
        
    #2. Run it on each model dataframe
        model.res.Low_CAS <- adjust_pvals_overwrite(model.res.Low_CAS)
        model.res.High_CAS <- adjust_pvals_overwrite(model.res.High_CAS)
        
        model.res.Low_CLO <- adjust_pvals_overwrite(model.res.Low_CLO)
        model.res.High_CLO <- adjust_pvals_overwrite(model.res.High_CLO)
        
        model.res.CASCLO <- adjust_pvals_overwrite(model.res.CASCLO)
        model.res.Control <- adjust_pvals_overwrite(model.res.Control)
      
  
        

# =============== CROSS-TREATMENT MODELS ======================================= 
   
        
  ######## Step 1. Practice model to see how it is working #####

      #Practice 1: Madeup SNP to see how the loop will function (no loop yet)       
      
      #1. Simulate some af and cov data for a single SNP, bind this to the index table
      #     fake_snp <- Low_High_CAS_comb_index
      #     fake_snp$af <- runif(111,min=0,max=1) #needs to match the number of rows in the index we made
      #     fake_snp$cov <- runif(111, min=30,max=200)
      # 
      #     #2. Set up a model when treating time as a continuous variable. Weighted by coverage. Does the rate of change in allele frequency over time differ between treatments?
      #     fake_snp_model <- glm(af~time*dose, family="quasibinomial",weights=cov,data=fake_snp)
      # 
      #     fake_snp_model_trends <- emtrends(fake_snp_model, specs = "dose", var = "time")
      # 
      #     #What if we want pairwise comparisons of what the trend is for that snp between each treatment.
      #     fake_snp_model_trends_pairwise <- summary(pairs(fake_snp_model_trends, adjust = "fdr"))
      # 
      # 
      #     #Pull estimates and standard errors
      #       pvals_adj <- fake_snp_model_trends_pairwise$p.value
      #       beta <- fake_snp_model_trends_pairwise$estimate
      # 
      #     #Name columns consistently
      #       names(pvals_adj) <- paste0("pval_adj_CAS_", gsub(" ", "",fake_snp_model_trends_pairwise$contrast)) #pull the names of the contrasts for the pvalues
      #       names(beta) <- paste0("beta_CAS_", gsub(" ", "",fake_snp_model_trends_pairwise$contrast)) #pull the names of the contrasts for the beta
      # 
      #       dd <- as.data.frame(cbind(t(pvals_adj), t(beta)))
      # 
      # 
      # 
      # 
      # 
      # #Practice 2: Small subset of real SNPs from the dataset (with loop)
      # 
      #     #Pull out first 20000 SNPs
      #     df_20000snps <- Low_High_CAS_comb[c(1:20000),]
      # 
      # 
      #     #Set up the test function
      #     funx.model.test <- function(df_20000snps,Low_High_CAS_comb_index){
      # 
      #       freq <- cbind(as.numeric(df_20000snps[Low_High_CAS_comb_index$idaf]),
      #                     as.numeric(df_20000snps[Low_High_CAS_comb_index$idcov]),
      #                     Low_High_CAS_comb_index) #freq is basically pulling the af and cov data from the snp table and binding it to the index table, one snp at a time
      #       colnames(freq)[1:2] <- c("af","cov")
      # 
      #       mod1 <- glm(af~time*dose,family="quasibinomial",weights=cov,data=freq)
      # 
      #       snp_model_trends <- emtrends(mod1, specs = "dose", var = "time")
      #       snp_model_trends2 <- test(emtrends(mod1, specs = "dose", var = "time"))
      # 
      # 
      #       pvals <- snp_model_trends2$p.value #extract the pvalues from the emtrends output
      #       beta <- snp_model_trends2$time.trend #pull the estimate from the emtrends output
      #       names(pvals) <- paste0("pval_CAS_", gsub(" ", "",snp_model_trends2$dose))
      #       names(beta) <- paste0("beta_CAS_", gsub(" ", "",snp_model_trends2$dose))
      # 
      # 
      #       #What if we want pairwise comparisons of what the trend is for that snp between each treatment.
      #       snp_model_trends_pairwise <- summary(pairs(snp_model_trends, adjust = "fdr"))
      # 
      # 
      #       pvals_pairwise <- snp_model_trends_pairwise$p.value #extract the p values from emtrends pairwise comparison
      #       beta_pairwise <- snp_model_trends_pairwise$estimate
      #       names(pvals_pairwise) <- paste0("pval_adj_CAS_", gsub(" ", "",snp_model_trends_pairwise$contrast)) #pull the names of the contrasts for the pvalues
      #       names(beta_pairwise) <- paste0("beta_CAS_", gsub(" ", "",snp_model_trends_pairwise$contrast)) #pull the names of the contrasts for the beta
      # 
      # 
      #       model_results <- as.data.frame(cbind(t(pvals),t(beta),t(pvals_pairwise),t(beta_pairwise))) #transpose the pvalue objects then bind in a single dataframe
      # 
      #     }
      # 
      #     #set up the test model
      #     test.run <- pbapply::pbapply(df_20000snps,1,function(x) funx.model.test(x,Low_High_CAS_comb_index))
      # 
      # 
      #     #Now make a dataframe
      #     test.run.datum <-  do.call(rbind.data.frame,test.run) #combine the results into a dataframe
      #     test.run.datum <- cbind(df_20000snps[,1:2],test.run.datum) #add the CHROM and POS columns
      # 
      #     rm(df_20000snps, freq,bad_row,i,test.run,test.run.datum)

      
  ######## Step 2. Run model for CAS doses ######
     
      #1. Set up the function
      funx.model2 <- function(snp_df,index_df){ #snp_df = ONE SNP(one row), passed in from pbapply later; index_df = FULL index table
        
        #Build modeling dataframe. Rows come from index_df, AF and cov values come from a single SNP each time
        freq <- cbind(as.numeric(snp_df[index_df$idaf]), 
                      as.numeric(snp_df[index_df$idcov]),
                      index_df) #freq is basically pulling the af and cov data from the snp table and binding it to the index table, one snp at a time
        colnames(freq)[1:2] <- c("af","cov")
        
        #fit SNP-specific model
        mod1 <- glm(af~time*dose,family="quasibinomial",weights=cov,data=freq)
        
        #Estimate time slopes within each dose (still one SNP at a time)
        mod1_trends <- emtrends(mod1, specs = "dose", var = "time")
       
        #Pairwise comparisons of slopes (one SNP at a time)
        mod1_trends_pair <- summary(pairs(mod1_trends, adjust = "fdr"))
        
        #Extract per SNP statistics
          pvals_adj <- mod1_trends_pair$p.value
          beta <- mod1_trends_pair$estimate

        #Name columns consistently
          names(pvals_adj) <- paste0("pval_adj_CAS_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the pvalues
          names(beta) <- paste0("beta_CAS_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the beta

        #Return one row of results per SNP
          model_results <- as.data.frame(cbind(t(pvals_adj),t(beta)))
      }
      
      
      #2. Chunk level orchestrator (reduce peak memory usage)
      funx.model2_chunk <- function(df_treat, index_df, chunk_size = 1000) { #df_treat = FULL SNP table, index_df = FULL index table (same as above)
        
        n_rows <- nrow(df_treat)
        chunk_starts <- seq(1, n_rows, by = chunk_size)
        results <- list()
        
        for (start in chunk_starts) {
          end <- min(start + chunk_size - 1, n_rows)
          cat("Processing rows", start, "to", end, "\n")
          
          chunk <- df_treat[start:end, , drop = FALSE]
          
          #Apply per SNP model row-by-row
          chunk_result <- pbapply::pbapply(chunk, 1,function(x) funx.model2(x, index_df)) #build its own freq df, fit its own glm per snp, run its own emtrends
          
          results[[length(results) + 1]] <- chunk_result
          gc()
        }
        
        #Combine all SNP results
        model_results <- do.call(c, results)
        model_results <- as.data.frame(do.call(rbind, model_results))
        
        #Add CHROM/POS columns
        model_results <- cbind(df_treat[,1:2], model_results)
        return(model_results)
      }
      
      #3. Run the model
      model.res.Low_High_CAS_comb <- funx.model2_chunk(Low_High_CAS_comb, Low_High_CAS_comb_index)
     
    
      
  ######## Step 3. Run model for CLO doses  ####
      #1. Set up the function
      funx.model3 <- function(snp_df,index_df){ #snp_df = ONE SNP(one row), passed in from pbapply later; index_df = FULL index table
        
        #Build modeling dataframe. Rows come from index_df, AF and cov values come from a single SNP each time
        freq <- cbind(as.numeric(snp_df[index_df$idaf]), 
                      as.numeric(snp_df[index_df$idcov]),
                      index_df) #freq is basically pulling the af and cov data from the snp table and binding it to the index table, one snp at a time
        colnames(freq)[1:2] <- c("af","cov")
        
        #fit SNP-specific model
        mod1 <- glm(af~time*dose,family="quasibinomial",weights=cov,data=freq)
        
        #Estimate time slopes within each dose (still one SNP at a time)
        mod1_trends <- emtrends(mod1, specs = "dose", var = "time")
        
        #Pairwise comparisons of slopes (one SNP at a time)
        mod1_trends_pair <- summary(pairs(mod1_trends, adjust = "fdr"))
        
        #Extract per SNP statistics
        pvals_adj <- mod1_trends_pair$p.value
        beta <- mod1_trends_pair$estimate
        
        #Name columns consistently
        names(pvals_adj) <- paste0("pval_adj_CLO_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the pvalues
        names(beta) <- paste0("beta_CLO_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the beta
        
        #Return one row of results per SNP
        model_results <- as.data.frame(cbind(t(pvals_adj),t(beta)))
      }
      
      
      #2. Chunk level orchestrator (reduce peak memory usage)
      funx.model3_chunk <- function(df_treat, index_df, chunk_size = 1000) { #df_treat = FULL SNP table, index_df = FULL index table (same as above)
        
        n_rows <- nrow(df_treat)
        chunk_starts <- seq(1, n_rows, by = chunk_size)
        results <- list()
        
        for (start in chunk_starts) {
          end <- min(start + chunk_size - 1, n_rows)
          cat("Processing rows", start, "to", end, "\n")
          
          chunk <- df_treat[start:end, , drop = FALSE]
          
          #Apply per SNP model row-by-row
          chunk_result <- pbapply::pbapply(chunk, 1,function(x) funx.model3(x, index_df)) #build its own freq df, fit its own glm per snp, run its own emtrends
          
          results[[length(results) + 1]] <- chunk_result
          gc()
        }
        
        #Combine all SNP results
        model_results <- do.call(c, results)
        model_results <- as.data.frame(do.call(rbind, model_results))
        
        #Add CHROM/POS columns
        model_results <- cbind(df_treat[,1:2], model_results)
        return(model_results)
      }
      
      
      #3. Run the model
      model.res.Low_High_CLO_comb <- funx.model3_chunk(Low_High_CLO_comb, Low_High_CLO_comb_index)
     
      
      
  ######## Step 4. Run model for single vs combined drug ####
      #1. Set up the function
      funx.model4 <- function(snp_df,index_df){ #snp_df = ONE SNP(one row), passed in from pbapply later; index_df = FULL index table
        
        #Build modeling dataframe. Rows come from index_df, AF and cov values come from a single SNP each time
        freq <- cbind(as.numeric(snp_df[index_df$idaf]), 
                      as.numeric(snp_df[index_df$idcov]),
                      index_df) #freq is basically pulling the af and cov data from the snp table and binding it to the index table, one snp at a time
        colnames(freq)[1:2] <- c("af","cov")
        
        #fit SNP-specific model
        mod1 <- glm(af~time*treatment,family="quasibinomial",weights=cov,data=freq)
        
        #Estimate time slopes within each dose (still one SNP at a time)
        mod1_trends <- emtrends(mod1, specs = "treatment", var = "time")
        
        #Pairwise comparisons of slopes (one SNP at a time)
        mod1_trends_pair <- summary(pairs(mod1_trends, adjust = "fdr"))
        
        #Extract per SNP statistics
        pvals_adj <- mod1_trends_pair$p.value
        beta <- mod1_trends_pair$estimate
        
        #Name columns consistently
        names(pvals_adj) <- paste0("pval_adj_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the pvalues
        names(beta) <- paste0("beta_", gsub(" ", "",mod1_trends_pair$contrast)) #pull the names of the contrasts for the beta
        
        #Return one row of results per SNP
        model_results <- as.data.frame(cbind(t(pvals_adj),t(beta)))
      }
      
      
      #2. Chunk level orchestrator (reduce peak memory usage)
      funx.model4_chunk <- function(df_treat, index_df, chunk_size = 1000) { #df_treat = FULL SNP table, index_df = FULL index table (same as above)
        
        n_rows <- nrow(df_treat)
        chunk_starts <- seq(1, n_rows, by = chunk_size)
        results <- list()
        
        for (start in chunk_starts) {
          end <- min(start + chunk_size - 1, n_rows)
          cat("Processing rows", start, "to", end, "\n")
          
          chunk <- df_treat[start:end, , drop = FALSE]
          
          #Apply per SNP model row-by-row
          chunk_result <- pbapply::pbapply(chunk, 1,function(x) funx.model4(x, index_df)) #build its own freq df, fit its own glm per snp, run its own emtrends
          
          results[[length(results) + 1]] <- chunk_result
          gc()
        }
        
        #Combine all SNP results
        model_results <- do.call(c, results)
        model_results <- as.data.frame(do.call(rbind, model_results))
        
        #Add CHROM/POS columns
        model_results <- cbind(df_treat[,1:2], model_results)
        return(model_results)
      }
      
      #3. Run the model
      model.res.CAS_CLO_CASCLO_comb <- funx.model4_chunk(CAS_CLO_CASCLO_comb, CAS_CLO_CASCLO_comb_index)

  
# =============== COMBINE RESULTS INTO MASTER TABLES ======================================= #
   
  ######## Step 1. Combine into master tables and rename columns #######   
      
    #### Master table 1: DRUG STRENGTH (INDIVIDUAL + CROSS-TREATMENT MODELS) ####
      #1. Combine the output from individual models with cross treatment models for drug dose
        model.res.Dose <- model.res.Low_CAS %>%
          inner_join(model.res.High_CAS, by = c("CHROM","POS")) %>%
          inner_join(model.res.Low_CLO, by = c("CHROM","POS")) %>%
          inner_join(model.res.High_CLO, by = c("CHROM","POS")) %>%
          inner_join(model.res.Control, by = c("CHROM","POS")) %>%
          inner_join(model.res.Low_High_CAS_comb, by = c("CHROM","POS")) %>%
          inner_join(model.res.Low_High_CLO_comb, by = c("CHROM","POS"))
        
        
      #2. Add a column with cumulative genome position (for plotting downstream) 
        
        chr_sizes <-  c(I = 230218, II = 813184, III = 316620, IV = 1531933, V = 576874, VI = 270161,
                        VII = 1090940, VIII = 562643, IX = 439888, X = 745751, XI = 666816,
                        XII = 1078177, XIII = 924431, XIV = 784333, XV = 1091291, XVI = 948066) #store the chromosome sizes for our ref genome version R64-1-1
        
        #Calculate cumulative start positions for each chromosome. For example, chrII starts right after chr1 ends, and so on. 
        chr_starts <- c(0, cumsum(chr_sizes[-length(chr_sizes)]))  #exclude last chr because it doesnt need a start value for the next one
        names(chr_starts) <- names(chr_sizes) #assign chromosome names so we can use CHROM as an index to look up start positions
        
        
        #Add MB column in the third position
        model.res.Dose <- model.res.Dose %>%
          mutate(Gaxis = POS + chr_starts[CHROM],
                 MB = Gaxis / 1e6) %>% #convert to Mb
          relocate(MB, .after = POS) %>% #Move Mb column to after POS
          select(-Gaxis) #Remove temporary Gaxis column
        
      #3. Add a SNP_ID column for downstream analysis
        model.res.Dose <- model.res.Dose %>%
          mutate(SNP_ID = paste0(CHROM,"_",POS)) %>%
          relocate(SNP_ID, .after = MB)
        
        
        
    #### Master table 2: SINGLE VS COMBINED DRUG (INDIVIDUAL + CROSS-TREATMENT MODELS) ####
      #1. Combine the output from individual models with cross treatment model 
        model.res.Complex <- model.res.High_CAS %>%
          inner_join(model.res.High_CLO, by = c("CHROM","POS")) %>%
          inner_join(model.res.CASCLO, by = c("CHROM","POS")) %>%
          inner_join(model.res.Control, by = c("CHROM","POS")) %>%
          inner_join(model.res.CAS_CLO_CASCLO_comb, by = c("CHROM","POS"))
      
        
      #2. Add a column with cumulative genome position (for plotting downstream) 
        
        #Add MB column in the third position
        model.res.Complex <- model.res.Complex %>%
          mutate(Gaxis = POS + chr_starts[CHROM],
                 MB = Gaxis / 1e6) %>% #convert to Mb
          relocate(MB, .after = POS) %>% #Move Mb column to after POS
          select(-Gaxis) #Remove temporary Gaxis column
        
      #4. Add a SNP_ID column for downstream analysis
        model.res.Complex <- model.res.Complex %>%
          mutate(SNP_ID = paste0(CHROM,"_",POS)) %>%
          relocate(SNP_ID, .after = MB)
        
        
    #### Export tables for downstream analysis ####
        
      #Export master table for low vs high drug strengths
        #write_tsv(model.res.Dose, "Output_files/model.results.Dose.txt")
        
      #Export master table for single vs combined drug results
        #write_tsv(model.res.Complex, "Output_files/model.results.Complex.txt")
        
        