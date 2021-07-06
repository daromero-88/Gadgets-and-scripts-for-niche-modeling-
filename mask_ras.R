#####RASTER - RASTER CROPPING######

--------------------------------------------------------
  #'  Author: Daniel Romero-Alvarez
  #'                
  #' The following function will allow the user to crop and mask 
  #' raster files using other rasters. useful for cropping particular
  #' produts with the output of an species distribution model 
  --------------------------------------------------------
  
  ####DETAILS#####
#' The function upgrade the 
#' 
#' 


####PACKAGES#####
library(raster)
library(sp)
library(rgeos)
library (rgdal)


#####FUNCTION#####

#' ARGUMENTS: 
#' @param raster1 raster object to be cropped 
#' 
#' @param raster2 raster object that is going to be used to crop and mask the desired object. 
#' 
#' @param projection_def projection for portion of the function where a mask is applied 
#' 
#' @describeIn 
#' 

# CODE mask_ras -------
# DEPENDENCIES: raster, sp, rgeos, rgdal 

#FUNCTION: 
mask_ras= function (raster1, raster2, projection_def){
  r1 = rasterToPoints(raster2)
  rp1 = data.frame (r1[which (r1[,3]==1),])
  sp_r1 = SpatialPointsDataFrame(coords = rp1[,1:2], data = rp1,
                                 proj4string = projection_def)
  r2 = mask (raster1, sp_r1)
  return (r2)
}

#add this option: 
# if (write.output = T){
#   writeRaster(r2, filename = paste (getwd (), '/r2', sep  = ''), 
#               format = 'ascii', bylayer = T, suffix = 'r2', 
#               overwrite = T, NAflag = -9999)
# }



#####EXAMPLE#####


####OPTIMIZATION####
#' Add the section of WRITE OUPUT to automatically write the ouput in the folder of the 
#' working directory. 

