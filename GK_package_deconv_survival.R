## functions for survival analysis of deconcolution results

#' library / function dependencies:
#' tidyverse
#' extractFromBetween function
  
#'mergeMultiRefResults
#'
#' a function that takes in a vector of files, that are the results of a multi reference deconvolution
#' 
#' @param files results of a multi- reference cellular composition deconvolution
#' @param cohort_pattern character vector of lenght = 2, to be used to the extractFromBetween function, to get the name of the cohort (or mix) from the file name. If lenght is only one, then it will be used as cohort name.
#' @param ref_pattern character vector of lenght = 2, to be used to the extractFromBetween function, to get the name of the reference matrix from the file name
#' 
#' @return list of 2. first element is the merged results of the multi reference deconvolution in one single long dataframe, and the second dataframe contains the mean , sd and rel.sd. of the predictions for each sample
#' @export

# summarise multi-reference deconvolution results
mergeMultiRefResults <- function(files, cohort_pattern, ref_patterns ) {
  
  #get files by lapply
  deconv_results <- lapply(files, function(file){
    
    df <- read.table(file, header = T, sep = "\t")
    colnames(df)[1] <- "sample"
    
    df$reference <- extractFromBetween(ref_pattern[1], ref_pattern[2], file)
    
    if (length(cohort_pattern == 2)) {
      df$cohort <- extractFromBetween(cohort_pattern[1], cohort_pattern[2], file)
    } else {
      df$cohort <- cohort_pattern
    }
    
    
    return(df)
  })
  
  #bind the different results into one dataframe
  deconv_results <- Reduce(bind_rows, deconv_results) 
  
  #pivot file into long format
  deconv_results <- deconv_results %>%
    pivot_longer(cols = c(-sample, -reference, -cohort), names_to = "celltype", values_to = "pred") %>%
    drop_na()

  #calculate mean sd and relative sd of the multi-ref predictions
  deconv_results_sum <- deconv_results %>% 
    group_by(cohort, sample, celltype) %>%
    summarise(mean_prediction = mean(pred, na.rm =T),
            sd_prediction = sd(pred, na.rm = T),
            n_prediction = n()) %>%
    ungroup() %>%
    mutate(rel_sd = sd_prediction/mean_prediction)
  
  out <- list(
    merged_results = deconv_results, 
    sum_results = deconv_results_sum
  )

  return(out)
}

##  functions to read in clinical data from TCGA cohorts

#' BindSimilarFiles
#' 
#' function to simultaneously read a number of similar files, and collect them in a list or in a dataframe
#' 
#' @param files a vector of file paths
#' @sur_pattern a character vector of length 2, to extract the name of the file from the path. usnig the extractFromBetween function
#' @as.char logical, to transform the data frames into characters or not. sometimes, the data types can be conflictins, and transformation could be necessary
#' @do.reduce logical, to reduce the read files ionto a fingle dataframe, of to keep is as a list
#' 
#' @return a data frame or a list
#' 
#' @export


BindSimilarFiles <- function(files, sur_pattern, as.char = F, do.reduce = T) {
  
  #core
  res <- lapply(files, function(file) {
    
    dat <-  read.table(file, sep = "\t", header = T )
    dat$file_ID <- extractFromBetween(sur_pattern[1], sur_pattern[2], file)
    
    # if the rows cannot be bound together because conflicting data types, all can be converted to character.
    if (as.char) {
     dat <-  mutate(dat, across(everything(), as.character))
    }
    
    return(dat)
  })
  
  if (do.reduce){
    res <- Reduce(bind_rows, res)
  }
  
  return(res)  
}

## simplified functions that use the BindSimilarFiles function, for some commonly used TCGA metadata

#' bindTCGAindexed
#' 
#' read in indexed clinical files of the TCGA cohort
#' all parameters are pre set
#' 
#' @projects character vector of TCGA cohorts to be included. all cohorts by default
#' 
#' @return clinical survival data from the TCGA cohorts

getTCGAindexed <- function( directory = "TCGA_processed", match_pattern = "_indexed_clinical.tsv", projects = NA,
                             as.char = F, do.reduce = T ) {
  
  files <- list.files(path = directory, pattern = match_pattern, full.names = T)
  sur_pattern <- c("TCGA-", match_pattern)
  
  if(!is.na(projects)){
    files <- unlist(lapply(projects, function(project) {
              grep(project, files, value = T)
              }))
  }
  
 out <- BindSimilarFiles(files = files, sur_pattern = sur_pattern,  as.char = as.char, do.reduce = do.reduce)
  
 rownames(out) <- NULL
 
 return(out)
}


#' bindTCGAdrug
#' 
#' read in the ddrug clinical files of the TCGA cohort
#' all parameters are pre set
#' 
#' @projects character vector of TCGA cohorts to be included. all cohorts by default
#' 
#' @return clinical survival data from the TCGA cohorts

getTCGdrug <- function( directory = "TCGA_processed", match_pattern = "_drug_clinical.tsv", projects = NA,
                            as.char = F, do.reduce = T ) {
  
  files <- list.files(path = directory, pattern = match_pattern, full.names = T)
  sur_pattern <- c("TCGA-", match_pattern)
  
  if(!is.na(projects)){
    files <- unlist(lapply(projects, function(project) {
      grep(project, files, value = T)
    }))
  }
  
  BindSimilarFiles(files = files, sur_pattern = patt,  as.char = as.char, do.reduce = do.reduce)
  
}


#### analyzing survival data, with deconvolution results

#' tailorSurvivalData
#' 
#' TCGA cohort has the survival data in a format that is not ideal for downstream analysis. 
#' This function takes the TCGA indexed clinical data (output of the getTCGAindexed function) and formats it to be compatible with the survival package
#' all parameters are pre set
#' @return a dataframe that contains the TCGA patient id, the censor event and the time


tailorTCGAsurvival <- function(data, 
                               cols= c("bcr_patient_barcode", "vital_status", "days_to_last_follow_up", "days_to_death") , 
                               na.rm = T,
                               status_options = list("censor" = "Alive", "event" = "Dead")
  ) {
    
    if (length(cols) != 4) {
      stop("number of selected columns have to be 4")
    }
    
    #subset for columns of interest
    data <- data[ , cols]
    
    #standardize column names
    colnames(data) <- c("id", "status", "censor_time", "event_time")
    
    #filter for accepted events
    data <- data[data$status %in% unlist(status_options), ]
    
    #define combined time
    data <- mutate(
      data, 
      status = ifelse(status %in% status_options$event, 1, 0)) %>%
      mutate(time = ifelse(status == 1, event_time, censor_time))
    
    
    data <- data[, c("id", "status", "time")]
    
    if (na.rm) {
      data <- data[!is.na(data$time), ]
    }
    return(data)
}

#' makeTertiles
#' 
#' function to create tertiles based on cell populations for
#' 
#' @param data a data frame, can be grouped if needed
#' @param value column of the dataframe, on which the tertiles will be based
#' 
#' @return a data frame, with an additional tertile column
#' 

makeTertiles <- function(data, value) {
  res <- mutate(data, tertile = ntile({{value}} , 3)) %>%
    mutate(tertile = if_else(tertile == 1, 'Low', if_else(tertile == 2, 'Medium', 'High')))
  
  return(res)
}

#'getCoxResults
#'
#'function to easily access the results of a COX hazard ratio analysis. 
#'
#'@param x the result of a cox analyis
#'
#'@return a dataframe the cintains the most important parameters of tha cox analysis

getCoxResults <- function(x){ 
  x <- summary(x)
  p.value<-signif(x$coefficients[5], digits=3)
  wald.z.score<-signif(x$coefficients[4], digits=3)
  se.beta <- signif(x$coef[3], digits = 3)
  beta<-signif(x$coef[1], digits=3);#coeficient beta
  HR <-signif(x$coef[2], digits=3);#exp(beta)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
  
  
  res <- data.frame(beta = beta, exp.beta_HR = HR, se.beta= se.beta, HR_confint_lower = HR.confint.lower,
                    HR_confint_upper = HR.confint.upper, wald.test = wald.z.score, p.value = p.value)
  
  return(res)
}


#'uniCoxHR
#'
#'function to perform multiple univariant cox hazard ratio analysis on a set of covariants of a survival data.

uniCoxHR <- function(x, covariants, time_col = "time", status_col = "status") {
  
  #get the results of the cox regression for each covariant in a list
  res <- lapply(covariants, function(covar){
    
    model <- coxph(formula = as.formula(paste('Surv(', time_col, ',', status_col, ')~', covar)) , 
                   data = x)
    
    model_res <- getCoxResults(model)
    model_res$covariant <- covar
    
    return(model_res)
  })
  
  #combine the results froma list
  res <- Reduce(rbind, res)
  
  return(res)
}

#'KMplotCellEntity 
#'
#'Kaplan meier plot of high and low tertiles of a cell type within a cancer entity
#'
#'@param x data frame with survival data and multiple covarinats
#'@param entity named character - the name is column name and the value is the entity of interest
#'@param cellcolumn  the column name that contains the terile information for the cell type of interest
#'@param time_col,status_col names of the columns tontaining survival time and status
#'
#'@return ggplot KM curve
#'@export


KMplotCellEntity <- function(x, entity, cellcolumn, time_col, status_col, remove.medium = T) {
  
  plotdata <-  x[ x[ , names(entity)] == entity   , ]
  
  if(remove.medium){
    plotdata <- plotdata[ plotdata[ , cellcolumn ] != "Medium" , ]
  }
  
  
  plotformula <- as.formula(paste('Surv(', time_col, ',', status_col, ')~', cellcolumn))
  
  p<- survfit2( plotformula , data = plotdata) %>%
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability",
      title = paste(entity, "-", cellcolumn),
      subtitle = paste("Number of patients:", nrow(plotdata))
    )+
    add_pvalue("annotation", size = 5)
    
  
  return(p)  
}


