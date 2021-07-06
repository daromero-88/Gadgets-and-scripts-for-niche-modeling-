#####BINOMIAL TEST FOR THRESHOLDED MODELS######

--------------------------------------------------------
  #'  Author: Daniel Romero-Alvarez
  #'                
  #' The following function will allow the user to perform the 
  #' binomial test to assess if a model is different than random
  #' and therefore assess the statistic significance of a threhsolded 
  #' model. 
--------------------------------------------------------

####DETAILS#####
#' The following function has been adapted from the supplementary material of the publication:  
#' Novel Methods in Disease Biogeography: A Case Study with Heterosporosis
#' Escobar et al. 2017 
#' doi: 10.3389/fvets.2017.00105
#' More information regardint the in-built R function to develop the binomial distribution: 
#' binomial probability function explanation: https://www.tutorialspoint.com/r/r_binomial_distribution.htm
#' 
#' This test was presented for the first time in: https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1466-822X.2002.00275.x
#' And a great explanation can be found in the following talk: https://www.youtube.com/watch?v=BIfXB6McNTY, 
#' the latter in the context of the free Ecological Niche Modeling course 2020. 
#' 


####PACKAGES#####
library(raster)
library(sp)
library(rgeos)
library (rgdal)


#####FUNCTION#####

#' ARGUMENTS: 
#' @param bin_model a raster file with 0 for unsuitable regions and 1 for suitable regions
#' It is usually obtained after applying a binarization rule over a continuous modeling output. 
#' 
#' @param ind_pts dataframe with three columns, species name, longitude and latitude. This is the
#' typical arrange of an occurrence dataframe in ecological niche modeling experiments. The most important
#' detail is that longitude and latitude should be the columns 2 and 3 of the dataframe. 
#' 
#' @describeIn bin_test is a function that allows that applies the binomial test over a binarized (thresholded)
#' model output to assess for its statistical significance using independent occurrences. See section details
#' for further information. 

# CODE bin_test -------
# DEPENDENCIES: raster, sp, rgeos, rgdal 

#FUNCTION: 
bin_test = function (bin_model, ind_pts){
  #raw data
  ind_pts_vals = data.frame (extract (bin_model, ind_pts[,2:3]))
  mod_pxls = data.frame(rasterToPoints(bin_model))
  #counts
  pt_fail = length (ind_pts_vals[ind_pts_vals==0])
  pt_sucess = length (ind_pts_vals[ind_pts_vals==1])
  unsuit_pxls = dim (mod_pxls[which(mod_pxls[,3] == 0),])[1]
  suit_pxls = dim (mod_pxls[which(mod_pxls[,3] == 1),])[1]
  #calculation of binomial probabiliy 
  bn_prob = 1-pbinom(pt_sucess, #occurrences predicted as present 
                     pt_fail+pt_sucess, #total number of occurrences (trials)
                     (suit_pxls/(suit_pxls+unsuit_pxls))) #proportion of area predicted suitable 
  return (bn_prob)
}


#####EXAMPLE#####
