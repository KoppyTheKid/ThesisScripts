# custom plot fuctions based on ggplot

# load library dependencies:

library(tidyverse)
library(RColorBrewer)
library(ggpubr)

### VARIABLES AND ADATA
default_palette <- c(brewer.pal(n= 12, name = "Paired") , brewer.pal(n= 8, name = "Set2"))

### FUNCTIONS

#' RepeatColorPalette 
#' function to be used by other plots to reuse a color palette if the original palette is not long enouh
#' 
#' @return a RcolorBrewer palette
#' 
#' @export

RepeatColorPalette <- function(filldistinct, colpalette = default_palette)  {
  
  if (filldistinct > length(colpalette)) {
    message("Number of unique elemnts in fillvar is larger than number of available colors. Colors will be repeated.")
    rep = ceiling(filldistinct / length(colpalette))
    colpalette <- rep(colpalette, times = rep )
  }
  return(colpalette)
}

#'CorHeatmap - a heatmap to visualize correlation matrixes
#'
#' takes a correlation matrix and cerates a standardised ggplot heatmap
#' 
#' @param data dataframe, in ggplot compatible format
#' @param xvar aes x
#' @param yvar aes y 
#' @param cor.value correlation values
#' 
#' @return a ggplot tile plot
#' 
#' @export

CorHeatmap <- function(data , xvar, yvar, cor.value, lab = T) {
  
  p <- ggplot(data = data, aes(x = {{xvar}} , y= {{yvar}}, fill = {{cor.value}}, label = round({{cor.value}}, digits = 2)))+
    geom_tile()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_gradientn(colours = c("blue4", "white", "red4"), limits = c(-1,1))
  
  if (lab == T){
    p <- p +  geom_text( size = 3)
  }
  
  return(p)
}


#'StackedBarChart
#'
#'stacked barchart for visualizing proportions whitihin a sample
#' @param data ggplot compatible dataframe
#' @param yvar quantitative variable for the y axis
#' @param xvar categorical variable for x axis
#' @param fillvar categorical variable for coloring
#' @param flip_coord logical, to flip the axises or not
#' 
#' @return a ggplot stacked barchart
#' 
#' @export

StackedBarChart <- function(data, yvar, xvar, fillvar, flip_coord = T) {
  
  colVector<- select(data, {{xvar}})
  
  #repeat color palette if necessary
  colpalette <- RepeatColorPalette(filldistinct = n_distinct(colVector))
  
  p <- ggplot(
    data = data,
    aes(x={{xvar}}, y={{yvar}}, fill={{fillvar}}))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = colpalette)+
    theme(text = element_text(size = 6))
  
  if (flip_coord) {
    p <- p + 
      coord_flip()+
      theme(legend.position = "right")
  } else {
    p <- p +
      theme(legend.position = "left") 
  }
  
  return(p)
}

#'SimpleBarChart
#'
#' @param data ggplot compatible dataframe
#' @param yvar quantitative variable for the y axis
#' @param xvar categorical variable for x axis
#' @param unique_colors logical
#' 
#' @return a ggplot barchart
#' 
#' @export

SimpleBarChart <- function(data, xvar, yvar, unique_colors = F) {
  
  colVector<- select(data, {{xvar}})
  
  #define colors
  if (unique_colors){
    colpalette <- RepeatColorPalette(filldistinct = n_distinct(colVector))
  } else {
    colpalette <- rep("dodgerblue", times = n_distinct(colVector)) ### fix the bug in this line
  }
  
  p <- ggplot(
    data = data,
    aes(x = {{xvar}}, y = {{yvar}}, fill = {{xvar}}))+
    geom_bar(stat = "identity", width = 0.6)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none" )+
    scale_fill_manual(values = colpalette)
  
  return(p)
} 

#'BoxplotWithSignificance
#'
#' @param data ggplot compatible dataframe
#' @param yvar quantitative variable for the y axis
#' @param xvar categorical variable for x axis
#'
#' @return a ggplot boxplor
#' 
#' @export

BoxplotWithSignificance <- function(data, xvar, yvar) {
  
  p <- ggplot(
    data = data,
    aes(x = {{xvar}}, y = {{yvar}}))+
    geom_boxplot(outlier.size = 0.3, notch = F)+
    #geom_jitter(size = 0.4, color = "black")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #significance layer has to be added as well.
  
  return(p)
}

#'ViolinPlot
#'
#' @param data ggplot compatible dataframe
#' @param yvar quantitative variable for the y axis
#' @param xvar categorical variable for x axis
#'
#' @return a ggplot violinplot
#' 
#' @export

ViolinplotWithSignificance <- function(data, xvar, yvar) {
  
  p <- ggplot(
    data = data,
    aes(x = {{xvar}}, y = {{yvar}}))+
    geom_violin()+
    #geom_jitter(size = 0.4, color = "black")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #significance layer has to be added as well.
  
  return(p)
}

#' corLong - correlation matrix instantly transformed into a three column dataframe without repeats
#' 
#' for more ducumentation, see cor()
#' @param x a numeric vector, matrix or data frame, passed on to cor()
#' @param y NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method 	a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' 
#' @return a data frame containing correlation values in a long pivoted format
#' 
#' @export

corLong <-function(x, y =NULL, use = "everything", method = c("pearson", "kendall", "spearman"), rm.triangle = F) {
  res<-cor(x = x, y = y, use = use, method = method)
  
  if(rm.triangle){
    res[lower.tri(res, diag = T)] <- NA
  }
  
  res <- as.data.frame(res)%>%
    rownames_to_column("cat1") %>%
    pivot_longer(cols = -c("cat1"), values_to = "cor.p", names_to = "cat2") %>%
    drop_na() 
  
  return(res)
}

#'plotHeatSurv
#'
#'plot a heatmap showing the survival risks, the results of a Hazard Ratio regession
#'
#' @param data ggplot compatible dataframe
#' @param yvar categorical variable for the y axis
#' @param xvar categorical variable for x axis
#' @param fillvar numericat variable, pvalu or z-scores to indicate the risk score
#' @param mode pval or zscore  -  this setting just effects the color range. if none, then the defaazlut color gradieent is used

plotHeatSurv <- function (data, xvar, yvar, fillvar, mode = "zscore", show.label = T) {
  
  if (!mode %in% c("pval", "zscore", "none")) {
    stop("mode can only be pval, zscore or none")
  }
  
  p <- ggplot(
    data = data,
    aes(x = {{xvar}}, y = {{yvar}}, fill = {{fillvar}}, label = signif({{fillvar}}, digits = 2))
  ) +
    geom_tile(col = "black")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (show.label) {
    p <- p + geom_text(col = "white", size = 2)
  }
  
  if (mode == "pval") {
    p <- p +
      scale_fill_gradient(low = "red", high = "white")
  } 
  
  if (mode == "zscore") {
    p <- p+
      scale_fill_gradient2(low = "green", high = "red", mid = "black")
  }
  
  return(p)
}
