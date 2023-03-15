
#### functions and tables for the DeconvCustomRefR package - pre deconvolution part


# TISCH_studies table: data to be included in the package, contains information on the downloadable studies on the TISCH2 portal

      #TISCH_studies <- read.csv("TISCH/TISCH_dataset_table.csv")
      #saveRDS(TISCH_studies, file = "TISCH/TISCH_studies.rds")

TISCH_studies <- readRDS("TISCH/TISCH_studies.rds")


#' Download files from the TISCH2 portal (http://tisch.comp-genomics.org/home) 
#' 
#' TISCH_download()  accepts a Dataset.ID from the TISCH site, and downloads processed single cell expression data, cell meatdata, and the average expression profile of clusters
#' 
#' @param study single character. Dataset.ID of the single-cell experiment to be downloaded -- UPDATE FOR VECTORS
#' @param files either "all" or a character vector containing "expr", "meta" or "cluster". Defines which files to download for the study
#' @param directory path to the directory where the files will be downloaded. working directory by default
#' @param cluster_unzip logical. whether to unzip the zip file containig the average expression profiles of clusters
#' @param wait Sys.sleep betwen file downloads. 5 seconds by default
#' 
#' @return the function has no output, the files are ownly downloaded
#' 
#' @export

TISCH_download <- function(studies, files = "all", directory = NA , cluster_unzip = F, wait = 5) {
  
  # assigning download dirrectory
  if (is.na(directory)) {
    directory = getwd()
  }
  
  for (study in studies) {
    
    # checking if study is available in TISCH
    if (!study %in% TISCH_studies$Dataset.ID) {
      message(paste0("Study ",study, " is not found on TISCH"))
    } else {
      
      ## downloading files
      # downloading expression matrix
      if (files == "all" | "expr" %in% files) {
        message("Downloading expression matrix...")
        
        url <- paste0(
          "https://biostorage.s3.ap-northeast-2.amazonaws.com/TISCH_2022/",
          study, "/", study, 
          "_expression.h5" )
        
        filedest <- paste0(directory, "/", study, "_expression.h5")
        if (!file.exists(filedest)) {        
        download.file(url = url, destfile = filedest)
        message("Waiting for 5 seconds...")
        Sys.sleep(wait)
        }
      }
      # downloading cellular metadata
      
      if (files == "all" | "meta" %in% files) {
        message("Downloading cell metadata...")
        
        url <- paste0(
          "http://tisch.comp-genomics.org/static/data/",
          study, "/", study, 
          "_CellMetainfo_table.tsv" )
        
        filedest <- paste0(directory, "/", study, "_CellMetainfo_table.tsv")
        if (!file.exists(filedest)) { 
        download.file(url = url, destfile = filedest)
        message("Waiting for 5 seconds...")
        Sys.sleep(wait) 
        }
      }
      
      # downloading average cluster expression profiles
      if (files == "all" | "cluster" %in% files) {
        message("Downloading average expression profile of clusters...")
        
        url <- paste0(
          "http://tisch.comp-genomics.org/static/data/",
          study, "/", study, 
          "_Expression.zip" )
        
        filedest <- paste0(directory, "/", study, "_Expression.zip")
        
        if (!file.exists(filedest)) { 
        download.file(url = url, destfile = filedest)
        message("Waiting for 5 seconds...")
        Sys.sleep(wait)
        }
      }
      
      if (files %in% c("all", "cluster") & cluster_unzip) {
        unzip(zipfile = paste0(directory, "/", study, "_Expression.zip"), exdir = paste0(directory, "/Expression_unzip"))
      }
      
      message(paste0("Download of ", study, " DONE"))
      
    }
  }
}

#' Generating reference or pseudo mix matrices out of single cell RNASeq data
#' 
#'
#' @param matrix expresion matrix, sparse (.mtx) or normal (.tsv)
#' @param clusters cluster annotation of single cells - vector with an equal length to matrix width
#' @param features name of genes  - vector with a length equal to matrix rows
#' @param method "median" or "mean"  - metric to be used for calculating the expression signature. for pseudo mixes mean is highly recommended
#' 
#' @return The function returns a data frame that contains the mean or median values for each feature per cell cluster, of sample
#' 
#' @export

makeRefMatrix = function(matrix, clusters, features, method = "median") {
  
  # check method
  if (!method %in% c("median", "mean")) {
    stop("Method not known. Possible options: median, mean")
  }
  
  # assign an empty matrix 
  out_matrix = matrix(nrow = length(features), ncol = length(clusters))
  
  # find the list of individual clusters
  cluster_names = unique(clusters)
  
  group_means = matrix(nrow = length(features), ncol = length(cluster_names))
  
  #calculate median or mean of 
  for (i in 1:length(cluster_names)) {
    
    group_cell_index = which(cluster_names[i] == clusters)
    
    if ( method == "mean") {
      
      group_means[,i] = rowMeans(raw_matrix[, group_cell_index], na.rm = T)
    } else {
      
      group_means[,i] = apply(raw_matrix[, group_cell_index], 1, median, na.rm = TRUE)
    }
    
    message(paste0(method ,"expressions of group ", i, " of ", length(group_names)  ," were calculated."))
  }
  
  colnames(group_means) = cluster_names
  group_means = as.data.frame(group_means)
  group_means = cbind(data.frame(geneID = features), group_means)
  
  
  return(group_means)
  
}




#' selecting mos variable features of a dataframe
#'
#' @param data dataframe, rows are features, and columns aare observations
#' @param feature_no number of variable features to be selected
#' @param feature_col integer or character indicating the column to be used as features
#'
#'@return a character vector of the most variable features
#'
#'@export

#selection of most variable genes/features from an expression matrix

mostVariableFeatures <- function(data, feature_no = 2000, feature_col = NA) {
  #get feature names
  if (is.na(feature_col)) {
    feat <- rownames(data)
  } else if (feature_col <= ncol(data)) {
    feat <- data[ , feature_col]
  } else {
    stop("feature_col must be lower or equal to the ncol of data")
  }
  
  # calculate coefficient of variance for each feature
  numdata = data[ , unlist(lapply(data, is.numeric), use.names = F)]
  
  means = apply(numdata, 1, mean)
  vars = apply(numdata, 1, var)
  
  cv = vars/means
  
  # bind cv and feature name together
  feat_var = data.frame(
    feature = feat,
    cv = cv
  )
  
  #select most variabe features
  feat_var = arrange(feat_var, desc(cv))
  
  feat_var = feat_var[!is.na(feat_var$cv), ]
  feat_var = feat_var[1:feature_no , ]
  
  return(feat_var$feature)
}

#' extract the sring between the two patterns
#' 
#' this function uses regexp to extract a substring from in between two patterns.
#' 
#' @param p1 pattern before the string to be eextracted
#' @param p2 pattern after the string to be ectracted
#' @string character or character vector 
#' 
#' @return a vector of the extracted pattern
#' 
#' @export

extractFromBetween <- function(p1, p2, string){
  out<- gsub(paste0(".*", p1, "(.+)" , p2 ,".*"), "\\1", string)
  
  return(out)
}
 

### to be continued: update makeRefMatrix to be compatible with downsteram elements

