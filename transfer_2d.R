#####ENM TRANSFERENCES 2D######

--------------------------------------------------------
  #'  Author: Daniel Romero-Alvarez
  #'                
  #' The following function will allow the user to transfer the information 
  #' enclosed by a convex hull defined by occurrence points in a particular 
  #' environmental space to another environmental space and project it to 
  #' geography. The function will work with 2 environmental variables—hence
  #' the name—but allows a quick exploration of the environmental potential of 
  #' the interpolation of environments defined by occurrences in e-space. 
  
--------------------------------------------------------

####PACKAGES#####
library(raster)
library(sp)
library(rgeos)
library (rgdal)


#####OPTIMIZATION######: 
#' Add if we want or not individual plots with a pflag, check code in melioidoosis
#' add corresponding labels to big plots number 1
#' add legend? 
#' in general plot, change xlab and ylab 
# 
  
#####FUNCTION#####

#' ARGUMENTS: 
#' @param pts1 dataframe with longitude and latitude as columns,
#' they represent observations that are going to be transfered to other environments
#' 
#' @param bck1 rasterstack of environmental variables, 
#' values for the points defined in the first parameter. 
#' 
#' @param bck2 rasterstack of environmental variables, 
#' there variables represent environemnts for the transference region 
#'
#' @return 
#' \code{df_trans} dataframe with information regarding the overlap between the background of transference
#' convex hull and the points (observation) based convex hull 
#' 
#' @describeIn transfers_2d is a function that allows the user the explore the environmental space 
#' of any particular set of points in a studied area, create a convex hull in the environmental space 
#' and then obtain the overlap of points from other environmental conditions to that particular convex hull
#' this allows transferences in environmental space that can be seen then transfered back to geography
#' 

# CODE transfers_2d -------
# DEPENDENCIES: raster, sp

#FUNCTION: 
transfers_2d = function (pts1, bck1, bck2){
  
  #convex hull function 
  #from: https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
  Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
    hpts <- chull(x = xcoord, y = ycoord)
    hpts <- c(hpts, hpts[1])
    lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  }  
  
  #data from observations
  df = data.frame (pts1[,1:2], extract (bck1, pts1[,1:2])) #dataframe with points and values of environments
  
  #transfer data 
  df0 = rasterToPoints(bck1) #coordinates and values of ALL environments in the study area 
  
  #Using partial data: (if needed)
  #if less data is needed, we can obtain values from a shape file as random points: 
  # pt_random = spsample (plg2, n = 10000, 'random')
  # df2 = data.frame (pt_random, extract(bck2, pt_random@coords))
  # df2 = df2[!is.na(df2[,3] & df2[,4]),] #eliminate NAs for sure...otherwise convex hull function not working 
  # dim(df2)
  
  #Using all data 
  df2 = rasterToPoints(bck2) #coordinates and values of ALL environments in the TRASNFER area 
  
  #plot1: transfer environments  
  dev.new()
  x_max = max (df2[,3], df[,3])
  x_min = min (df2[,3], df[,3])
  y_max = max (df2[,4], df[,4])
  y_min = min (df2[,4], df[,4])
  plot (df2[,3], df2[,4], col = 'grey', pch = 3, cex = 0.5, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) # PC environments from the background transfer
  points (df[,3], df[,4], col = 'black', pch = 16, cex = 0.5) #outbreak localities
  Plot_ConvexHull(xcoord = df2[,3], ycoord = df2[,4], lcolor = "grey") #convex hull in background transfer
  Plot_ConvexHull(xcoord = df[,3], ycoord = df[,4], lcolor = "black") #convex hull in points transfer
  
  #Overlap dataframe for transference
  #from: https://www.rdocumentation.org/packages/grDevices/versions/3.6.1/topics/chull
  ch_index = chull (x = df[,3], y = df[,4]) #obtaining the indexes of the convex hulls for the overlap 
  
  #coordinates of convex hull in the corresponding variables
  ver_t = df[,3][c(ch_index)] 
  ver_h = df[,4][c(ch_index)]
  
  #overlap between both convex hulls
  #https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/point.in.polygon
  overlap = point.in.polygon(df2[,3], df2[,4], ver_t, ver_h, mode.checked = F) #find those values that are inside the convex hull defined by the observations
  
  #transference dataframe 
  df_trans = cbind (df2, overlap) #dataframe with the overlap vector, 0 = not present in transference background, 1 = present in transference background 
  
  #Plot 2: To geography  
  dev.new()
  plot (bck2[[1]])
  points (df_trans[which(df_trans[,dim(df_trans)[2]]== 1),1:2], col = 'black', pch = 16, cex = 0.5) #add only the points inside the polygon (==1)
  
  #Plot 3: all information in one plot 
  dev.new ()
  par(mfrow = c(2,2))
  plot (bck1[[1]], cex = 0.5, main = 'Geographical study area (ENV1)')
  points (df[,1], df[,2], col = 'black', pch = 16, cex = 0.5)
  
  plot (df0[,3], df0[,4], col = 'grey', pch = 3, cex = 0.5, xlim = c(x_min, x_max), ylim = c(y_min, y_max), 
        main = 'Environmental study area', xlab = 'ENV1', ylab = 'ENV2')
  points (df[,3], df[,4], col = 'black', pch = 16, cex = 0.5) 
  Plot_ConvexHull(xcoord = df0[,3], ycoord = df0[,4], lcolor = "grey")
  Plot_ConvexHull(xcoord = df[,3], ycoord = df[,4], lcolor = "black")
  
  plot (df2[,3], df2[,4], col = 'grey', pch = 3, cex = 0.5, xlim = c(x_min, x_max), ylim = c(y_min, y_max), 
        main = 'Environmental trasnfer area', xlab = 'ENV1', ylab = 'ENV2')
  points (df[,3], df[,4], col = 'black', pch = 16, cex = 0.5) 
  Plot_ConvexHull(xcoord = df2[,3], ycoord = df2[,4], lcolor = "grey")
  Plot_ConvexHull(xcoord = df[,3], ycoord = df[,4], lcolor = "black")
  
  plot (bck2[[1]], cex = 0.5, main = 'Geographical transfer area (ENV1)')
  points (df_trans[which(df_trans[,dim(df_trans)[2]]== 1),1:2], col = 'black', pch = 16, cex = 0.5) #add only the points inside the polygon (==1)
  
  return (df_trans)
}

