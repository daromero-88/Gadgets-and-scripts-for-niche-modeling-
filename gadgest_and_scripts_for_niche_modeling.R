#####Gadgets and scripts for niche modeling######

#By Daniel Romero-Alvarez 


#INSTALLING PACKAGES-----------------

install.packages ('devtools')
install.packages('raster')
install.packages('sp')
install.packages('rgeos')
install.packages ('rgdal')
install.packages ('maptools')
install.packages('mapdata')
install.packages ('dismo')
install.packages('spThin')
install.packages ('ENMeval')

if(!require(kuenm)){
  devtools::install_github("marlonecobos/kuenm")
}

#LIBRARIES----------------------------

library (devtools)
library(raster)
library(sp)
library(rgeos)
library (rgdal)
library (maptools)
library(mapdata)
library (dismo) #SDM algorithms plus other functions
library(spThin) #thinning of occurrences
library (ENMeval) #maxnet algorithm 
library (kuenm) #maxent algorithm dissected 
library (virtualspecies) #create virtual species
library(corrplot) #create correlation matrices
library (gatepoints) #for advance functions
library (rgl)
library (Rcpp)
library (hypervolume) #hypervolume based models 

#GREAT RESOURCES----------------------
#MYBLOG

#WORKING DIRECTORY--------------------
setwd ('/Users/daniel/Documents/GitHub/Gadgets-and-scripts-for-niche-modeling-')


#MANIPULATING SHAPES-------------------

#WORLD-MAP
data("wrld_simpl", package = "maptools")
WGS84 = crs(wrld_simpl) # geographic projection

plot (wrld_simpl)

#ZOOM
zoom(wrld_simpl, ext=drawExtent(), new=TRUE, useRaster=FALSE) #zooming in a particular part of the world

head (wrld_simpl@data) #reviewing first 6 rows
colnames(wrld_simpl@data) #reviewing column names

#subsetting by country: 
idd = subset (wrld_simpl, NAME == 'India')

plot (idd)


#DOWNLOADING NATURAL EARTH SHAPE FILES

#' Manual download: 
#' https://www.naturalearthdata.com/downloads/

#' Downloading with the script

dir.create('shapes') #create folder

shp1 = 'https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip' #website for direct download, obtained from 'redirecton' a the download manager

download.file(shp1, destfile = './shapes/countries.zip') 

unzip(zipfile = "./shapes/countries.zip", exdir = './shapes')

#reading downloaded file
world2 = readOGR ('./shapes/ne_10m_admin_0_countries.shp')

head (world2@data) #reviewing first 6 rows
colnames(world2@data) #reviewing column names

View (world2@data)
unique (world2@data$REGION_UN) #unique values in the REGION_UN column 
unique (world2@data$SUBREGION)

plot (world2[which(world2@data$SUBREGION=='Southern Asia'),])

#subsetting by subregion 

sa1 = subset (world2, SUBREGION == 'Southern Asia')
plot (sa1)


#SOURCES FOR ENVIRONMENTAL VARIABLES------------------
#bioclimatic variables
#https://www.worldclim.org/data/bioclim.html

#worldclim: 
#https://www.worldclim.org/data/index.html

#dowloading the data for the first time: 
wc1 = getData('worldclim', var='bio', res=10)
par (mfrow = c(2,1))
plot (wc1[[1]], main = 'bio1') #please notice the DOUBLE SQUARE BRACKET
plot (wc1[[12]], main = 'bio12')

wc2 = wc1[[c(1,4,12)]]

#dev.new()
plot (wc1)

#For more variables, please review: 
#https://www.romerostories.com/post/resources-for-ecological-niche-modeling-environmental-variables-1


#CROPING AND MASKING-----------------------------------

sa2 = crop (wc2, sa1) #cropping extracts the extent of the area studied
plot (sa2)
sa2 = mask (sa2, sa1) #masking delineates the shape file in the raster 
plot (sa2)

#writing raster data: 

writeRaster(sa2, filename = names (sa2), format = 'ascii', 
            bylayer = T, suffix = names (sa2))

#CREATING A VIRTUAL SPECIES-----------------------------

#' Crucial literature for this section: 
#' Recent review on the topic: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04385 
#' We are using Leroy et al 2016 virtualspecies package, 
#' The original manuscript: https://borea.mnhn.fr/en/virtualspecies-r-package-generate-virtual-species-distributions
#' A step by step super detailed tutorial: http://borisleroy.com/files/virtualspecies-tutorial.html


#If not aware of how to implement this portion:
#READ THE CODE UNTIL THE BEGINING OF CALIBRATION AREAS
#READ THE CODE UNTIL THE BEGINING OF CALIBRATION AREAS
#READ THE CODE UNTIL THE BEGINING OF CALIBRATION AREAS
#READ THE CODE UNTIL THE BEGINING OF CALIBRATION AREAS


#1. Define environmental suitability function

#Calculate descriptive statistics from the environmental variables if needed: 
mn = cellStats(sa2[[1]], mean) #obtaining the mean of BIO1 in the Southern Asia region
st_dev = cellStats(sa2[[1]], sd) #obtaining the standard deviation of BIO1 in the Southern Asia region

mn2 = cellStats(sa2[[3]], mean) #obtaining the mean of BIO12 in the Southern Asia region
st_dev2 = cellStats(sa2[[3]], sd) #obtaining the standard deviation of BIO12 in the Southern Asia region


# Here I am using only bio1 and bio12 to create the virtual species
my.parameters = formatFunctions(bio1 = c(fun = 'dnorm', mean = 300, sd = 50), 
                                 bio12 = c(fun = 'dnorm', mean = 600, sd = 200),
                                 rescale = T) #we transformed to TRUE because I am using two variables

#2. projection of the simulated suitability to the landscape

species1 = generateSpFromFun(raster.stack = sa2[[c(1,3)]], parameters = my.parameters, 
                              rescale.each.response = T, #we transformed to TRUE because I am using two variables
                              rescale = T, #we transformed to TRUE because I am using two variables
                              plot = TRUE) 

plotResponse(species1)

plot (species1$suitab.raster)

#Step 3: Conversion to presence-absence

cellStats(species1$suitab.raster, stat = "mean") #default prevalence, the mean of the suitability
## [1] 0.1095543

#default parameters: 
species1.prev1 <- convertToPA(species1, 
                              PA.method = "probability",
                              prob.method = "linear") #uses the default prevalence, the mean of the suit rasters! 

plot (species1.prev1$pa.raster)

#tuned parameters, user defined threshold: 
species2.prev2 <- convertToPA(species1, 
                              PA.method = "probability",
                              prob.method = "linear",
                              species.prevalence = 0.05) #user define prevalence, has to be lower than predefined prevalence! 

species2.prev2$probability.of.occurrence

par (mfrow = c(1,1))
plot (species2.prev2$pa.raster)
plot (species2.prev2$suitab.raster)
plot (species2.prev2$probability.of.occurrence)


#writing raster: 
writeRaster(species2.prev2$probability.of.occurrence, filename = 'v_species_niche1', format = 'ascii', 
            bylayer = T, suffix = 'v_species_niche1', overwrite = T)


plotSuitabilityToProba(species2.prev2)


#Step 4: sample occurrences

occs2 <- sampleOccurrences(species2.prev2, 200, type = "presence only", #can be changed to presence/absence
                           plot = TRUE, extract.probability = T)

occs2$sample.points
#occs2 here is a complex object and should be examined with $

#visualization: 
plot(species2.prev2$probability.of.occurrence)
points (occs2$sample.points[,1:2][which (occs2$sample.points[,4]==1),], 
        col = 'red', cex = 0.3, pch = 3)

sp_points = occs2$sample.points[,1:2][which (occs2$sample.points[,4]==1),] #vector of occurrences 


#' columns 3 and 4 contains 0s and 1s representing the points that were sampled from the 
#' virtual species raster


#writing occurrence data: 
write.csv (sp_points, './tp1_points.csv', row.names = F)


#' Create the virtual species using the three selected variables! 
#' EXPLORE THE PARAMETERS OF THE FUNCTIONS MORE IN DETAIL 
#' 


#SAMPLE POINTS FROM TRUE DISTRIBUTION-----------------------------------

#reading existing niche of vs: 
v_species1 = raster('v_species_niche1.asc')

#reading points
tp1 = read.csv ('tp1_points.csv', header= T)

#reference column: 
tp1$ind =  paste(tp1[,1], tp1[,2], sep = "_") #new column for reference 


#' from the true observations, we collect only the 35% that are available from  
#' a true atual sampling scheme, assuming non-biased excercises: 
#' 

rsp = tp1[sample(nrow(tp1),(0.35*nrow(tp1))),]
dim(rsp)

#collecting independent points:
ind_p =  tp1[!tp1[,3] %in% rsp[,3],] #points not sampled in the true occurrence data 
dim(ind_p)

#' we select again the 35% of the points not used in the distribution to obtain an independent set of points 
#' for future evaluation 

ind_p2 = ind_p[sample(nrow(ind_p),(0.25*nrow(ind_p))),] 
ind_p2$sp1 = rep ('sp1', length(ind_p2[,1]))
ind_p2 = ind_p2[,(c(4,1,2,3))]
dim(ind_p2)

#eliminating the reference column 
tp1$ind = NULL
rsp$ind=NULL
ind_p2$ind = NULL


#writing definitive points: 
write.csv(rsp, 'model_pts.csv', row.names = F)
write.csv (ind_p2, 'ext_eval_pts.csv', row.names = F)

#visualization 
plot (v_species1)
points (tp1, col = 'red', pch = 3, cex = 0.6) #points from the truth distribution 
points (rsp, col = 'blue', pch = 16, cex = 0.6) #points for modeling 
points (ind_p2, col = 'black', pch = 16, cex = 0.6) #points for evaluation 

#CLEANING AND THINNING------------------------
model_pts = read.csv ('model_pts.csv', header = T)

#adding spp column: 
model_pts$sp1 = rep('sp1', nrow (model_pts))


#eliminating duplicates
model_pts$dup = paste(model_pts$sp1, model_pts$x, model_pts$y, sep = '_')
dim(model_pts)
model_pts = model_pts[!duplicated(model_pts$dup),1:3]
dim(model_pts)


#thinning the data using 30 km datasets

thin(model_pts, lat.col = "y", long.col = "x", spec.col = "sp1",
     thin.par = 30, reps = 3, write.files = TRUE, max.files = 3, 
     out.dir = "thining", out.base = "model_pts",
     write.log.file = FALSE, verbose = TRUE)


#PARTITION OF POINTS FOR INTERNAL EVALUATION------------------------

#WORKING WITH MODEL PTS CURATED: 
model_pts2 = read.csv ('./thining/model_pts_thin1.csv', header = T)
dim(model_pts2) #64 points to work 

model_pts2$check = paste(model_pts2[,2], model_pts2[,3], sep = "_")
train = model_pts2[sample(nrow(model_pts2), round((length(model_pts2[,1])/2))),]
test = model_pts2[!model_pts2[,4] %in% train[,4],]

dim(train)
dim(test)

model_pts2$check = NULL
train$check = NULL
test$check = NULL

write.csv(model_pts2, "model_pts_def.csv", row.names = FALSE)
write.csv(train, "in_train.csv", row.names = FALSE)
write.csv(test, "in_test.csv", row.names = FALSE)

#visualization
plot (sa1)
points (train[,2:3], col = 'red', pch = 3, cex = 0.6) #training points 
points (test[,2:3], col = 'blue', pch = 16, cex = 0.6) #evaluation points  


#CALIBRATION AREAS-------------------------------------

#DATA -> Points
model_pts3 = read.csv ('model_pts_def.csv', header = T)
dim(model_pts3)

#creating spatial points dataframe 
occ_sp = SpatialPointsDataFrame(coords = model_pts3[,2:3], 
                                data = model_pts3,proj4string = WGS84)

dev.new()
plot (occ_sp)

#RANDOM EXTENT-------------------------------------------
plot(sa1)
plot (occ_sp, add = T, pch = 3, cex = 0.6, col = 'red')
e = drawExtent()

cp = as(extent(e), 'SpatialPolygons')
proj4string(cp) = WGS84 #project extent to the same projection as wrld_simpl map

mask1 = crop (sa1, cp)

#visualization
plot(mask1)
plot (occ_sp, add = T, pch = 3, cex = 0.6, col = 'red')


#CONVEX HULL AROUND THE POINTS OF INTEREST-----------------------

ch_index = chull (x = model_pts3[,2], y = model_pts3[,3]) #obtaining the indexes of the convex hulls for the overlap 

#coordinates of convex hull in the corresponding variables
ver_t = model_pts3[,2][c(ch_index)] #vetices on longitude
ver_h = model_pts3[,3][c(ch_index)] #vetices on latitude

ch = convHull(cbind (ver_t, ver_h)) #create the convex hull

#visualization
plot(sa1)
plot (occ_sp, add = T, pch = 3, cex = 0.6, col = 'red')
plot (ch, add = T, border = 'red')

mask2 = crop (sa1, ch@polygons)
plot (mask2)


#BUFFER WITH ANY PARTICULAR REASON--------------------------------------
dd = 500000 #Define the distance buffer in meters because it has a WGS84 projection! 

#the buffer itself: 
buff1 = buffer(occ_sp, width = dd, dissolve = T) 

#visualization: 
plot (sa1)
plot (occ_sp, add = T, col = 'blue')
plot (buff1, add = T, border = 'red')


#IMPARTIAL BUFFER-------------------------------
cnt = data.frame(mean(occ_sp$x),mean(occ_sp$y)) #obtaining the centroid of points 
dist = mean(pointDistance(occ_sp, cnt, longlat = T)) #obtaining the mean across the distances of all occurrences to the centroid

buff2 = buffer(occ_sp, width = dist, dissolve = T) #create the buffer

df = data.frame (matrix(nrow = 1, ncol = 1, 1)) #fake dataframe to convert the spatial polygon to a spatial polygon dataframe...only this way can be written... 
buff2 = SpatialPolygonsDataFrame(buff2, df, match.ID = F) #transforming in correct object

plot (buff2, add= T, border = 'red')
points (cnt, pch = 16, cex = 2)


writeOGR(buff2, layer= 'buff2',  
         dsn = './BUFF2', driver = 'ESRI Shapefile')


#FRONTIERS--------------------------------------

#Dispersal dynamics of the studied species
#Youtube video --> https://www.youtube.com/watch?v=uzCj5JHGrQs


#PROCESSING ENVIRONMENTAL VARIABLES--------------------
#material for this section: 

#environmental variables to be assessed 
wc1

#shape file representing study area 
buff2

#READING RASTERS IN A BACTH (STACK) AND CLIPPING USING A MASK----------------------

#Establish your own working directory, where you have all your data: 
#' all the data means: 
#' - folder with all your rasters
#' - folder with the shapefile, REMEMBER: shapefiles use different files and they have to be all together 

setwd ('/Users/daniel/Documents/CURSO_INDIA/Divakar1/')

#you can test that you are in the correct folder using the get working directory function: 
getwd()

#read all the rasters you want to RESAMPLE as a STACK, so you can work as a batch: 
#NOTICE: './z' = represent the subfolder where all the raster files are located 
modisd = stack (list.files(path = './z',full.names = T))

#Manually select the shape file you will use to clip (crop/mask) your files
cliping_lay1 = readOGR (file.choose()) #select the .shp file, it has to have the rest of the corresponding shape files around (i.e., .dbf, .prj)

#cropping/masking: 
crop_files = crop (modisd, cliping_lay1) #cropping = match the extend of the shape file 
mask_files = mask (crop_files, cliping_lay1) #masking = match the contour of the shape file 

#Write your rasters as an .ascii file:  
writeRaster(mask_files, filename = 'masked1', format = 'ascii', 
            bylayer = T, suffix = names (mask_files), overwrite = T,
            NAflag = -9999)


#RESAMPLING FOR OBTAINING THE SAME PIXEL SIZE AND EXTENT--------------------------

#Establish your own working directory, where you have all your data: 
#' all the data means: 
#' - folder with all your rasters
#' - folder with the raster that is going to be use as a base layer

setwd ('/Users/daniel/Documents/CURSO_INDIA/Divakar1/')

#you can test that you are in the correct folder using the get working directory function: 
getwd()

#read all the rasters you want to RESAMPLE as a STACK, so you can work as a bacth: 
modisd = stack (list.files(path = './z',full.names = T))

#Manually select the raster file that you want to use as the transformation raster
#this one will have the pixel size and the extent that you are interested to obtain: 
base_layer = raster (file.choose())


#Use the RESAMPLE approach to obtain the appropriate: 
results1 = resample (modisd, base_layer, method = 'bilinear')

#Visualize your results: 
dev.new()
plot (results1[[1]])
dev.new()
plot (base_layer)

#Write your rasters as an .ascii file:  
writeRaster(results1, filename = 'modis1', format = 'ascii', 
            bylayer = T, suffix = names (results1), overwrite = T,
            NAflag = -9999)


#CORRELATION MATRIX------------------------------------

#cropping to calibration area: 
buf_ras = crop(wc1, buff2)
buf_ras = mask (buf_ras, buff2)

plot (buf_ras[[1]])


#correlation matrix: 

rc = layerStats(buf_ras, 'pearson', na.rm = T)

class(rc) #list with two objects, the correlation values and the mean of each variable 

dat = data.frame (rc[2]) #transform the second element in a dataframe
names = row.names (dat) #obtain the names 
length(names)

#transform vector into a matrix, manually specify the number of columns ncol= 19 in this case 
mat = matrix (unlist (rc, recursive = F, use.names = T), ncol = 19, byrow = T)
colnames(mat) = names #adding the names to the columns
rownames(mat) = c(names, 'mean') #addig names to the rows
View (mat) #checking the correlation matrix 

#selecting variables 
mat[mat[,1] <0.8,]

#visualization 
#library(corrplot)
mat_def = mat[-nrow(mat),] #matrix without the mean 

#dev.new()
corrplot (mat_def, method = 'number', type = 'upper', 
          tl.col = 'black', is.corr = F)

#variales for modeling 
un_wc1 = buf_ras[[c(1,2,7,12, 14, 15)]]

writeRaster(un_wc1, filename = paste(names (un_wc1), 'unc', sep ='_'), format = 'ascii', 
            bylayer = T, suffix = paste(names (un_wc1), 'unc', sep ='_'))

#variables for projection 
sel_vars = wc1[[c(1,2,7,12, 14, 15)]]

un_project = crop (sel_vars, sa1) #using the original mask 
un_project = mask (un_project, sa1)

writeRaster(un_project, filename = paste(names (un_wc1), 'proj', sep ='_'), format = 'ascii', 
            bylayer = T, suffix = paste(names (un_wc1), 'proj', sep ='_'))

#PRINCIPAL COMPONENT ANALYSIS--------------------------

out_folder <- "PCA_vs" #name of the outside folder
n_pcs <- 5 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

kuenm_rpca(buf_ras, var.scale = TRUE, write.result = T, project = F,
           proj.vars = proj_folder, n.pcs = n_pcs, in.format = 'ascii', out.format = 'ascii', 
           out.dir = out_folder)

#read PCs 
pca_vs= stack(list.files('./PCA_vs/Initial', full.names = T, pattern = 'asc'))

#table 
pc_res1 = read.table ('./PCA_vs/Initial/pca_results.txt', header = T, sep = '\t')

#dev.new()
plot (pca_vs)


#VARIABLE COMBINATIONS--------------------------------

var.dir = './var_unc/'
out.dir = './var_comb/'
min.number = 3 #change accorgingly, >2, <total # vars
in.format = 'ascii'
out.format = 'ascii'
kuenm_varcomb(var.dir, out.dir, min.number, in.format = "ascii",
              out.format = "ascii")

#FUNCTION IS NOT WORKING IN MY VERSION, SHOULD WORK MOST RECENT VERSIONS
#PROBLEM WITH THE PROGRESSION BAR 
#USE fix (kuenm_varcomb) IN THE CONSOLE TO FIX IT... 

#MODELING EXCERCISES-----------------------------------------

#MAXENT MODEL KUENM - KUENM-----------------------------------
#objects required for this section
#occurrences divided in three sets: 
model_pts3 = read.csv ('model_pts_def.csv', header = T)
train = read.csv ('in_train.csv', header = T)
test = read.csv ('in_test.csv', header= T)

#selected environmental variables: 
un_wc1

#projection variables: 
un_project

#folder with all these items plus the maxent algorithm
#more info: https://github.com/marlonecobos/kuenm

#define folder 
setwd('./KUENM')
getwd()

#MODEL CALIBRATION - KUENM---------------------------------------

occ_joint <- "model_pts_def.csv"
occ_tra <- "in_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(0.1, 0.5, 1, 1.5, 2)
f_clas <- c("l", 'lq', 'lqp')
maxent_path <- '/Users/daniel/Documents/GitHub/Gadgets-and-scripts-for-niche-modeling-/KUENM'
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, 
          maxent.path = maxent_path, wait = wait, run = run)

#help (kuenm_cal)

#MODEL EVALUATION AND SELECTION - KUENM---------------------------------------
occ_test <- "in_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 100
kept <- FALSE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform MOP calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for previous function

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)
#help (kuenm_ceval)

#FINAL MODEL AND TRANSFERENCE- KUENM---------------------------------------

batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- TRUE
G_var_dir <- "G_variables"
out_format <- "logistic"
project <- TRUE 
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)

#help (kuenm_mod)

#ANALYSIS OF EXTRAPOLATION MESS/MOP- KUENM-------------------------

sets_var <- "Set_2" # a vector of various sets can be used
out_mop <- "MOP_results"
percent <- 5
paral <- FALSE

kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var, out.mop = out_mop,
           percent = percent, parallel = paral)


#VISUALIZATION OF MODELS---------------------------------
#median of bootstrap
ne_mod = raster ('./Final_Models/M_0.1_F_lq_Set_2_NE/sp1_se_asia_median.asc')


#uncertainty depicted as range: 
ne_max = raster ('./Final_Models/M_0.1_F_lq_Set_2_NE/sp1_se_asia_max.asc')
ne_min = raster ('./Final_Models/M_0.1_F_lq_Set_2_NE/sp1_se_asia_min.asc')
uncer_mod = ne_max - ne_min 

#color ramp for uncerntainty: 
color = colorRampPalette(c('#c90020', 'azure2', '#0571b1'))

#visualization
dev.new()
par (mfrow = c(1,1))
plot(species2.prev2$probability.of.occurrence, main = 'Pathogen HK')
plot (ne_mod, main = 'KUENM model')
plot (uncer_mod, col= color(100), main = 'uncertainty (range)')

par (mfrow = c(1,1))
plot (ne_mod, main = 'KUENM model')
points (model_pts3[,2:3], col = 'red', pch = 4, cex = 0.3)
points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

zoom(ne_mod, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

writeRaster(uncer_mod, filename = 'ne_uncertain', format = 'ascii', 
            bylayer = T, suffix = 'ne_uncertain')


#MODEL THRESHOLDS---------------------------------------

#extracting the values from the model using all the calibration points
suit_vals = extract (ne_mod, model_pts3[,2:3]) 

#using the value of 5% as threshold 
th1 = quantile(suit_vals, 0.05) 

#collecting only those values above that threhsold, automatically transforms in 1 presence 0 absence 
ne_th1 = ne_mod >= th1

plot (ne_th1)
points (model_pts3[,2:3], col = 'red', pch = 4, cex = 0.3) #calibration points
points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3) #independent points 

writeRaster(ne_th1, filename = 'ne_th5', format = 'ascii', 
            bylayer = T, suffix = 'ne_th5')



#MAXENT MODEL - DEFAULT - KUENM---------------------------------

setwd ('/Users/daniel/Documents/GitHub/Gadgets-and-scripts-for-niche-modeling-/KUENM_DEFAULT')

occ_joint <- "model_pts_def.csv"
occ_tra <- "in_train.csv"
M_var_dir <- "M_variables"
occ_test <- "in_test.csv"
out_eval <- "Calibration_results"
batch_fin <- "Final_models2"
mod_dir <- "Final_Models2"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- TRUE
G_var_dir <- "G_variables"
out_format <- "logistic"
project <- TRUE 
ext_type <- "no_ext"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


max_default = raster ('./Final_Models/M_1_F_lqpth_Set_1_NE/sp1_se_asia_median.asc')
plot (max_default)

#binarize the model 
suit_vals2 = extract (max_default, model_pts3[,2:3]) 

#using the value of 5% as threshold 
th3 = quantile(suit_vals2, 0.05) 

#collecting only those values above that threhsold, automatically transforms in 1 presence 0 absence 
default_th = max_default >= th3

plot (default_th)
points (model_pts3[,2:3], col = 'red', pch = 4, cex = 0.3) #calibration points
points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3) #independent points 

writeRaster(default_th, filename = 'default_th5', format = 'ascii', 
            bylayer = T, suffix = 'default_th5')


#MAXNET MODEL ENMEVAL---------------------------------
#more info: 
#https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html
#https://github.com/jamiemkass/ENMeval

#required material for this section
#occurrences
model_pts3

#environemntal variables 
un_wc1 #for calibration 
un_project #for transference, using this variables for obtaining similar output

#modeling using MAXNET algorithm, not MAXENT, taking the advantage of leave-one-out approach

clus1_model = ENMevaluate(model_pts3[,2:3], un_project, bg.coords = NULL, occ.grp = NULL,
                          bg.grp = NULL, RMvalues = c (0.1, 0.5, 1, 2),
                          fc = c("L", "LQ", "LQP"), 
                          categoricals = NULL, n.bg = 10000, method = 'randomkfold',
                          overlap = FALSE, aggregation.factor = c(2, 2),
                          kfolds = 2, bin.output = TRUE, clamp = FALSE,
                          rasterPreds = TRUE, parallel = FALSE, numCores = NULL)

#whie using randomkfold, kfolds parameters should be larger than 1..

#FOR JACKKNIFE PARTITION (small occurrences) NOTICE THE PARAMETERS MODIFIED
# clus1_model = ENMevaluate(model_pts3[,2:3], un_project, bg.coords = NULL, occ.grp = NULL,
#                           bg.grp = NULL, RMvalues = c (0.1, 0.5, 1, 2),
#                           fc = c("L", "LQ", "LQP"), 
#                           categoricals = NULL, n.bg = 10000, method = 'jackknife',
#                           overlap = FALSE, aggregation.factor = c(2, 2),
#                           kfolds = NA, bin.output = FALSE, clamp = FALSE,
#                           rasterPreds = TRUE, parallel = FALSE, numCores = NULL)



#writing results table
clus1_res = clus1_model@results
View (clus1_res) #checking the results in R
write.csv (clus1_res, file = 'clus1_res.csv', row.names = F)

#evaluation plot: 

#omission rates based on Minimum Training Presence 
eval.plot(clus1_model@results, 'avg.test.orMTP',
          legend = TRUE, legend.position = 'topright')


#based on Akaike Information Criteria corrected by sample size (AICc)
eval.plot(clus1_model@results, 'AICc',
          legend = TRUE, legend.position = 'topright')

#based on Delta AICc
eval.plot(clus1_model@results, 'delta.AICc',
          legend = TRUE, legend.position = 'topright')


#Selecting best model based on minimun training presence and AICc

min_MTP = min(clus1_res[,9]) #0.03125 --> minimum omission rate 

minMTP_mods = clus1_res[which(clus1_res[,9] == min_MTP),]
dim(minMTP_mods) #8 models 

minAICc_mods = minMTP_mods[which(minMTP_mods[,14] == min (minMTP_mods[,14])),]
dim(minAICc_mods)

#selected models: 
enmeval_def = minAICc[,1]

#' if more than one model is selected at this step it will be convenient to obtain the 
#' median of all the best models to maximize the information available


enmeval_sel = clus1_model@predictions[[which (clus1_model@results$settings == enmeval_def)]]
plot (enmeval_sel)

#binarize the model: 

vals_clus1 = extract(enmeval_sel,model_pts3[,2:3])
mod_th = quantile(vals_clus1, 0.05) 
#mod_th = sort (vals_clus1)[1] if minimum training presence required
enmeval_bin = enmeval_sel>=mod_th
plot(enmeval_bin)

#writing continuos and binary raster: 

writeRaster(enmeval_sel, filename = 'enmeval_cont', format = 'ascii', 
            bylayer = T, suffix = 'enmeval_cont', overwrite = T,
            NAflag = -9999)

writeRaster(enmeval_bin, filename = 'enmeval_bin', format = 'ascii', 
            bylayer = T, suffix = 'enmeval_bin', overwrite = T,
            NAflag = -9999)




#HYPERVOLUME MODEL---------------------------------

#based on the Supporting Text S3 file downloaded from: 
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12865

#material required for this section 
#values of the environemntal variales defining the study area 
env_vals = data.frame(cbind (model_pts3[,2:3], extract(un_wc1, model_pts3[,2:3])))

#geographical projection raster 
un_project

#calibration: 
#SVM: parameters = svm.nu (default = 0.01) and svm.gamma (default = 0.5)
svm_mod = hypervolume_svm(data = env_vals[3:8]) #environmental values as dataframe
plot(svm_mod)

#projection geography: 
svm_geo = hypervolume_project (hv = svm_mod, rasters = un_project, type = "inclusion")
plot(svm_geo)

plot (un_project[[1]])

#writing results: 
write.csv (svm_mod@RandomPoints, 'hypervolume_points.csv', row.names = F)

#writing the raster files (svm models):
writeRaster(svm_geo, filename = 'hypervolum_svm', format = 'ascii',
            bylayer = T, suffix = 'hypervolum_svm', overwrite = T,
            NAflag = -9999)


#MODEL EVALUATION WITH INDEPENDENT DATA--------------------------------------------

######I AM HERE --> READING RASTERS######

#' decide where to put the folder, inside KUENM??
#' decide if I should add folders for each modeling effort... maybe this is better
#' previous option is crucial because a lot of stuff that is not related with KUENM
#' is inside the KUENM folder... so better to rearrange
#' 


setwd ('/Users/daniel/Documents/GitHub/Gadgets-and-scripts-for-niche-modeling-/')

#Continuous models: 
ne_mod = raster ('./eval_diff_models/sp1_se_asia_median_default.asc')
enmeval_sel = raster ('./eval_diff_models/enmeval_cont.asc')
max_default = raster ('./eval_diff_models/sp1_se_asia_median_original.asc')


#using binary models: 
#KUENM selected parameters
ne_th1 = raster ('./eval_diff_models/ne_th5.asc')

#KUENM default
default_th = raster ('./eval_diff_models/default_th5.asc')

#ENMeval selected
enmeval_bin = raster ('./eval_diff_models/enmeval_bin.asc')

#Hypervolume MTP 
svm_geo = raster ('./eval_diff_models/hypervolum_svm.asc')


plot (enmeval_sel)
plot (max_default)
plot (ne_mod)

plot (ne_th1)
plot (default_th)
plot (enmeval_bin)
plot (svm_geo)

#MANUAL EVALUATIONS - PARTIAL ROC----------------------------------------------------


#material for this section: 
#independent occurrences for evaluation 
ind_p2 

#continuous output models: 
#KUENM selected parameters
ne_mod

#KUENM default
max_default

#ENMeval selected
enmeval_sel

#Parameter E, error in occurrences, Peterson et al. 2008, in this case --> E = 5

#PARTIAL ROC 
#KUENM selected parameters
proc_res1 = kuenm_proc(occ.test = ind_p2[,2:3], model = ne_mod, 
                       threshold = 5, rand.percent = 100, 
                       iterations = 100, parallel = F)

proc_res1$pROC_summary 

#KUENM default
proc_res2 = kuenm_proc(occ.test = ind_p2[,2:3], model = max_default, 
                       threshold = 5, rand.percent = 100, 
                       iterations = 100, parallel = F)

proc_res2$pROC_summary 


#ENMeval selected
proc_res3 = kuenm_proc(occ.test = ind_p2[,2:3], model = enmeval_sel, 
                       threshold = 5, rand.percent = 100, 
                       iterations = 100, parallel = F)

proc_res3$pROC_summary 


#PENDING VISUALIZATION


#MANUAL EVALUATIONS - OMISSION RATE---------------------------------------------------


#Independent set of points: 
ind_p2 

#using binary models: 
#KUENM selected parameters
ne_th1

#KUENM default
default_th

#ENMeval selected
enmeval_bin

#Hypervolume MTP 
svm_geo


#' we use the independent points to extract the 1/0 values from the binarized model 
#' and then we obtain the average of that information 

#KUENM selected parameters
sens1 = mean(extract(ne_th1, ind_p2[,2:3])) #sensitivity, number of points correctly predicted as present
omr1 = 1-sens1 #omission rate, 1-sensitivity


#KUENM function can be used to calculate the omission rates directly on the continuous output: 
#' 5% threshold using the continuous model: 
#' any threshold should specified 
or_res1 = kuenm_omrat (model = ne_mod, occ.tra = model_pts3[,2:3],
                       occ.test = ind_p2[,2:3], threshold = 5)

or_res1

#' notice that here, the omission rate =  0.1875
#' which is one unit larger than the one produced with the manual binarization
#' where omr1 = 0.15625, this has to do with the use of either > or >= during the 
#' binarization and the internal calculations of the kuenm_omrat function. 
#' If we reduce the threshold in the kuenm function to threshold = 4, we obtain the same result
#' 

#KUENM default
sens2 = mean(extract(default_th, ind_p2[,2:3])) #sensitivity, number of points correctly predicted as present
omr2 = 1-sens2 #omission rate, 1-sensitivity


#ENMeval selected
sens3 = mean(extract(enmeval_bin, ind_p2[,2:3])) #sensitivity, number of points correctly predicted as present
omr3 = 1-sens3 #omission rate, 1-sensitivity


#Hypervolume MTP 
sens4 = mean(extract(svm_geo, ind_p2[,2:3])) #sensitivity, number of points correctly predicted as present
omr4 = 1-sens4 #omission rate, 1-sensitivity


#MANUAL EVALUATIONS - BINOMIAL TEST----------------------------------------------------------------------

#material required for this section 
#please read the function 'bin_test' available here: 

#binarized models 

#using binary models: 
#KUENM selected parameters
ne_th1

#KUENM default
default_th

#ENMeval selected
enmeval_bin

#Hypervolume MTP 
svm_geo

#Applying the 'bin_test' function: 

#KUENM selected parameters
bin_test(ne_th1, ind_p2)

#KUENM default
bin_test (default_th, ind_p2)

#ENMeval selected
bin_test (enmeval_bin, ind_p2)

#Hypervolume MTP 
bin_test(svm_geo, ind_p2)

#all the models are statistically significant... 


#EVALUATION OF FINAL MODEL WITH INDEPENDENT DATA - KUENM------------------------------

#help(kuenm_feval)

occ_ind <- "ext_eval_pts.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for previous functions

fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)


#EVALUATION FRONTIERS --> NULL MODEL------------------------------------------------

#Paper depicting this method: 
#https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.13573

#example with default svm_hypervolume model: 

#independent points: 
ind_p2

#material required for this section 
#values of the environemntal variales defining the study area 
env_vals = data.frame(cbind (model_pts3[,2:3], extract(un_wc1, model_pts3[,2:3])))

#USING ONLY BIO1, BIO7 and BIO12: 

#geographical projection raster 
un_project2 = un_project[[c(1,3,4)]]

#calibration: 
#SVM: parameters = svm.nu (default = 0.01) and svm.gamma (default = 0.5)
svm_mod_obs = hypervolume_svm(data = env_vals[,c(3,5,6)]) #working with only two variables 
plot(svm_mod_obs)

#projection geography: 
svm_geo_obs = hypervolume_project (hv = svm_mod_obs, rasters = un_project2, type = "inclusion")
plot(svm_geo_obs)

#obtaining observed omission rates: 
sens5 = mean(extract(svm_geo_obs, ind_p2[,2:3])) #sensitivity, number of points correctly predicted as present
omr5 = 1-sens5 #omission rate, 1-sensitivity

#generating null models, n = 100:

#' for this particular algorithm is crucial to obtain only points with values and without NAs
#' thus, first step relates with obtaining only pixels with information across the rasters: 

px_values = which (!is.na (values(un_wc1[[1]]))) 
  
or_null = NULL #dataframe to be filled 

for (jj in 1:50){
  #generating random points only on pixels with values 
  rand_indx = sample (px_values, nrow(model_pts3), replace = TRUE) #obtain the indexes of sampled values
  #rand_vals = un_wc1[[1]][rand_indx] #obtain the values of those indexes if needed
  rand_coords = xyFromCell(un_wc1[[1]], rand_indx) #obtain the coordinates of the indexes
  #obtaining environmental dataframe
  envs = data.frame(extract(un_wc1,rand_coords)) #obtain working dataframe
  svms = hypervolume_svm(data = envs[,c(1,3,4)]) #model developed only with three env. variables
  svm_proj = hypervolume_project (hv = svms, rasters = un_project2, type = "inclusion") #projecting to geography on those three env. variables
  or_n = 1-(mean (extract(svm_proj, ind_p2[,2:3]))) #calculating omission rates using independent points
  pre_df = cbind (rep = jj, or_n = or_n) #joining a dataframe with the replicate and statistic
  or_null = rbind (or_null, pre_df) #writing the results
}

#visualization of one random model: 
plot (svm_proj)
points (rand_coords, pch = 3, col = 'red')

#dataframe with null statistics
or_null

dev.new()
hist (or_null[,2], main = 'Randomization test \n\ Algorithm: hypervolume', 
      xlab = 'Omission rate', col = 'grey')
abline (v = omr5, col = 'red')
qn975 = quantile (or_null[,2], probs = 0.975)
qn025 = quantile (or_null[,2], probs = 0.025)
abline (v = qn975, lty = 3, col = 'black') 
abline (v = qn025, lty = 3, col = 'black')

#model is statistically significant, observation is uncontained in the null distribution 

write.csv (or_null, 'or_null.csv', row.names = F)


#EVALUATION FRONTIERS --> HYPERTEST------------------------------------------------

#Paper discussing this method: 
#https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13479
#Jimenez and Soberón, et al 2020. 


#TRANSMISSION RISK---------------------------------------------------------------

#definitive model:

#ENMeval model has the lowest OR from all the models
enmeval_sel #continuous model 
enmeva_bin #binary model 

dev.new()
par(mfrow = c(2,1))
plot (enmeval_sel, main = 'Continuous')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)
plot (enmeva_bin, main = 'Binarized (5%)')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)

#checking the environemtnal space of suitable pixels (see advance toosl)
hh = data.frame (rasterToPoints(enmeval_bin))
hh2 = e_space_back (data.frame(hh[hh[,3]==1,1:2]), sa2, sa2, pflag = T)


#OTHER MODES: 
#KUENM SELECTED MODEL
ne_mod #continuous model 
ne_th1 #binary model 

par(mfrow = c(2,1))
plot (ne_mod, main = 'Continuous')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)
plot (ne_th1, main = 'Binarized')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)


#KUENM DEFAULT
max_default #continuous model 
default_th #binary model 

par(mfrow = c(2,1))
plot (max_default, main = 'Continuous')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)
plot (default_th, main ='Binarized')
points (model_pts3[,2:3], col = 'red', pch = 3, cex = 0.3)

#OVERLAP WITH HUMAN POPULATION------------------------------------------------
#SOURCES: 
#visit: MY BLOG... 

# Population of the world
# https://sedac.ciesin.columbia.edu/data/collection/gpw-v4

# 2020, 2.5 mins --> 5 km 

world_pop1 = raster ('./world_pop/gpw_v4_population_count_rev11_2020_2pt5_min.tif')
plot (world_pop1)
zoom(world_pop1, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#crop and aligning population raster to the same spatial resolution --> 10 min 

pop_sa = crop (world_pop1, sa1)
pop_sa = mask (pop_sa, sa1)

plot (pop_sa)

#' Resampling: 
#' Basically upscaling data from 2.5 to 10 min to match the pixel size,
#' the spatial resolution, of our model. You can approach this problem using either
#' the aggregate or resample function. I prefer the latter, in this case you need
#' to specify the method, in this case BILINEAR replicate the results obtained with 
#' aggregate. More on this discussion can be found here: 
#' https://gis.stackexchange.com/questions/255150/using-resample-vs-aggregate-extend-in-r-to-have-rasters-of-matching-resolutio
#' Downscaling is usually not recommended. 

popsa_res = resample (pop_sa, sa2[[1]], method = 'bilinear')
plot (popsa_res)
plot(sa1, add = T)

#' Cropping the population raster to match the raster of our model: 
#' using the function mask_ras available at: 

pop_mod = mask_ras (popsa_res, enmeval_bin, WGS84)
plot(pop_mod)
plot (sa1, add = T)

#' By obtaining summary statistics from this raster, 
#' we can count the amount of people at risk or infer the percentage 
#' with respect of the entire population: 

#total population obtained as a sum of all the values of the pixels available: 
pop_count = cellStats(pop_mod, sum)
pop_total = cellStats(popsa_res, sum)

#percentage in overlap with the model: 
pop_perc = (pop_count/pop_total)*100

#Thus, 34.37567% of the population is overlapping the suitable areas of the model... 


#OVERLAP WITH OTHER SPECIES------------------------------------------------

#Land mammals of the world 
#downloaded from: https://www.iucnredlist.org/resources/spatial-data-download
mmi = readOGR (file.choose()) #it is a 1GB file so not added in github 

#Subseeting database with an endemic Indian species: 
#' endemic mammals of India: http://lntreasures.com/indiam.html
# 'Semnopithecus hypoleucos'
# 'Black-footed gray langur'

pl = mmi[which(mmi$binomial == 'Semnopithecus hypoleucos'),]

writeOGR(pl, layer= 'primate1', #writing shapefiles 
         dsn = 'primate', driver = 'ESRI Shapefile')

#Visualization:
par (mfrow = c(1,1))
plot (enmeval_bin)
plot (sa1, add = T)
plot (pl, add = T, border = 'blue')


#calculating percentage occupied by the species:
#cropping the model to the respective area:
ovp = crop (enmeval_bin, pl)
ovp = mask (ovp, pl)

plot (ovp)

#obtaining values in the cropped raster 

ovp_vals= rasterToPoints(ovp) 

ovp_pos = ovp_vals[ovp_vals[,3]==1,] #pixels suitable 
ovp_neg = ovp_vals[ovp_vals[,3]==0,] #pixels unsuitable 

percent_inrisk = (length(ovp_pos[,1])/length(ovp_vals[,1]))*100

#' Thus, 68.30645% of the theoretical range of the Black-footed gray langur 
#' overlaps with the distribution of the studied species. 

#Exploring environmental space (see advance tools): 

#Checking the environmental space
f = e_space (data.frame(ovp_vals[ovp_vals[,3]==1,1:2]), sa2, pflag = T)
f2 = e_space (data.frame(ovp_vals[ovp_vals[,3]==0,1:2]), sa2, pflag = T)

#Contrasting the environmental space of the spss with the entire study environmental area
g = e_space_back (data.frame(ovp_vals[ovp_vals[,3]==1,1:2]), sa2, sa2, pflag = T)
g2 = e_space_back (data.frame(ovp_vals[ovp_vals[,3]==0,1:2]), sa2, sa2, pflag = T)


#checking environmental space used by the species in comparison with the entire model output

rrr = mask_ras(sa2, enmeval_bin, projection_def = WGS84) #cropping environemnts with model output
plot (rrr)

#occupation of the spp in the total model output
h1 = e_space_back (data.frame(ovp_vals[ovp_vals[,3]==1,1:2]), rrr, rrr, pflag = T)


#SPECIFIC LANDSCAPE FEATURES----------------------------------------------

#subset by a particular feature such as ALTITUDE: 
#downloaded from: 

alt_wrld = raster (file.choose()) #to large for adding it in github

#crop to the study area: 

alt_sa = crop (alt_wrld, sa1)
alt_sa = mask (alt_sa,sa1)

dev.new()
plot(alt_sa)

#resampling 
alt_res = resample (alt_sa, sa2[[1]], method = 'bilinear')
plot (alt_res)

#subsetting by the model
alt_ss = mask_ras (alt_res, enmeval_bin, WGS84)
plot (alt_ss)

#' Let's say the species can live only BELOW 500
#' we should threshold our raster file at that particular number 

alt_th = 500 #threshold

alt_def = alt_ss<alt_th

par (mfrow = c(2,1))
plot (alt_ss)
plot (sa1, add = T)

plot (alt_def)
plot (sa1, add = T)

#' Then we can calculate again percentage of population exposed 
#' percent of area occupied or other measures...
#' this step is also called POST-MODELING 




#post-modeling steps:
#altitude stuff... 
#other species? add the variable itself vs postmodeling? 
#leave all this stuff for the study area!!!! CROPPING STUFF!!! 


#overlap with populations: INCIDENCE 
#overlap with urban regions (nightlight times)
#urban rural areas 
#overlap with livestock 
#overlap with other distributions (ILCP, anthropozoonosis)




#ADVANCE TOOLS-----------------------------------------
#' subsetting principal components for model assessment 
#' visualization of environmental space: 
#' function to crop rasters using other rasters 
#' occurrence contribution: jackknife approach 
#' assessing response curves
#' environmental convex hulls and environmental transferences  
#' suitability/uncertainty plots 
#' 
#' 


#SPLITTING ENVIRONMENTAL VARIABLES for PCs-------------------------------- 

setwd ('/Users/daniel/Documents/MIELOIDOSIS/envs_sets_BI')

#variables
temp =  (list.files('./temp', full.names= T))
hum =  (list.files('./hum', full.names= T))
soil =  (list.files('./soil', full.names= T))

#list of names of rasters 
b = list()
for (i in 1:length(temp)) {
  for (j in 1:length(hum)) {
    for (k in 1:length(soil)) {
      a = list(paste0(c(temp[i], hum[j], soil[k])))
      b[length(b)+1] =a
    }
  }
}

#PC2 can not be alone without PC1...so: 
#https://www.dataquest.io/blog/control-structures-in-r-using-loops-and-if-else-statements/

#is paramount to add the [[length(cc[[tp]]+1)]], otherwise 
#it will replace the elements added... 

for (tp in 1:length(b)){
  if(b[[tp]][1] == './temp/temp_pc_2.asc'){
    b[[tp]][[length(b[[tp]])+1]] = './temp/temp_pc_1.asc'
  }
  if(b[[tp]][2] == './hum/hum_pc_2.asc'){
    b[[tp]][[length(b[[tp]])+1]] = './hum/hum_pc_1.asc'
  }
  if(b[[tp]][3] == './soil/soil_pc_2.asc'){
    b[[tp]][[length(b[[tp]])+1]] = './soil/soil_pc_1.asc'
  }
}

#same for PC3, can not be alone: 

for (tp in 1:length(b)){
  if(b[[tp]][1] == './temp/temp_pc_3.asc'){
    b[[tp]][[length(b[[tp]])+1]] = './temp/temp_pc_1.asc'
    b[[tp]][[length(b[[tp]])+1]] = './temp/temp_pc_2.asc'
  }
  if(b[[tp]][2] == './hum/hum_pc_3.asc'){
    b[[tp]][[length(b[[tp]])+1]]= './hum/hum_pc_1.asc'
    b[[tp]][[length(b[[tp]])+1]] = './hum/hum_pc_2.asc'
  }
  if(b[[tp]][3] == './soil/soil_pc_3.asc'){
    b[[tp]][[length(b[[tp]])+1]]= './soil/soil_pc_1.asc'
    b[[tp]][[length(b[[tp]])+1]]= './soil/soil_pc_2.asc'
  }
}


#folders
dir.create ('./folders3')
setwd('/Users/daniel/Documents/MIELOIDOSIS/envs_sets_BI/folders3')
fol = c(1:length(b))
for (ras in 1:length(b)){
  dir.create(paste('Set', fol[ras], sep = '_'))
}

#rasters
setwd('/Users/daniel/Documents/MIELOIDOSIS/envs_sets_BI')

r = list.files('./folders3', full.names = T) #list of folders to get names

for (i in 1:length(b)){ #the key was to use i in 1: length(b) otherwise is not looping... 
  tt = stack (b[[i]])
  writeRaster(tt, paste (r[i], names (tt), sep = '/'), 
              format = 'ascii', bylayer = T, suffix = names(tt), 
              overwrite = T, NAflag = -9999)
}


#Visualization of sets: 
mat3 = matrix (unlist (b, use.names = TRUE), 
               ncol = 18, byrow = F) #number of columns should be known and modified 
View (mat3)

mat4 = t(mat3) #transponse the matrix so better visualization
rownames(mat4) = paste('Set', seq(1:18), sep = '_') #add row names
colnames (mat4) = c('Temp', 'Hum', 'Soil') #add column names 

write.csv (mat4, './folders/Sets_env.csv')

#CONSENSUS MODELS FROM MORE THAN ONE SET AND INTERQUANTILE RANGES----------------------------------


#selecting only those belonging to the NE category
avg_model_NE = grep (pattern = '._NE/', avg_model_boots, value = T)

#an internal list with only with the true replicates (avoiding averages, sd, max, min, etc)
models_NE1 = list()

for (i in 0:9){
  rr = grep (pattern = paste ('./Australia_', i, '_world.asc', sep = ''), #collect the name with this pattern, notice the i, changing for each replicate
             avg_model_NE, value = T)
  models_NE1[[length(models_NE1)+1]] = rr #TWO SQUARE BRACKETS FOR ADDING TWO ELEMENTS!!! 
}

models_NE1_def = unlist (models_NE1) #transform the list of lists in one list of elements 

#creating a stack of these models: 
NE1_ras = stack (models_NE1_def) #read this lists as stacks, that's is why the full.names is crucial


NE1_ras[[1]] #definitely work 

plot (NE1_ras[[1]]) 

#calculate required summary statistics: 

NE_mean = overlay (NE1_ras, fun = mean, na.rm = T) #I am using overlay here for all the functions for descriptive stats
NE_sd = overlay (NE1_ras, fun = sd, na.rm = T)
NE_median = overlay(NE1_ras, fun = median, na.rm = T)
NE_max = overlay(NE1_ras, fun = max, na.rm = T)
NE_min = overlay(NE1_ras, fun = min, na.rm = T)
NE_range = NE_max-NE_min

NE_quantiles = calc(NE1_ras, fun= quantile, na.rm = T) #dividing the rasters in five files of quantiles
NE_75 = NE_quantiles[[4]] #quantile 75
NE_25 = NE_quantiles[[2]] #quantile 25

NE_iqr = NE_75 - NE_25 #obtaining interquantile range 

dir.create ('./Summaries_NE_1')

NE_sum_vector = stack (NE_mean, NE_sd, NE_median, 
                       NE_max, NE_min, NE_range, NE_75, NE_25, NE_iqr)

NE_names = c('NE_mean', 'NE_sd', 'NE_median', 
             'NE_max', 'NE_min', 'NE_range', 'NE_75', 'NE_25', 'NE_iqr') #create a vector of names 

writeRaster (NE_sum_vector, filename = paste ('./Summaries_NE_1', '/', NE_names, sep = ''), 
             #filename should include the path of the final file, 
             #before the .asc which is included in FORMAT 
             format = 'ascii', bylayer = T, suffix = NE_names,
             overwrite = T, NAflag = -9999)

#EXPLORING THE ENVIRONMENTAL SPACE-----------------------------------------
#read the following functions published in: 
#https://github.com/daromero-88/Environmental-sampling

#' 1. E-space-functions.R
#' 2. Hutchinson-functions.R
#' 3. Post-track-functions.R

#Points in the environmnetal space: 
a = e_space (model_pts3[,2:3], un_wc1, pflag = T, wrld_map = sa1)

#Points in the environmental background: 
b = e_space_back(model_pts3[,2:3], un_wc1, un_wc1, pflag = T)
bb = e_space_back(model_pts3[,2:3], un_wc1, un_project, pflag = T)

#selecting points from the environmental space
pt_sel = hutchinson(EtoG = T, b , c(3,6), sa1, 2)

#selecting points from the geographical space
hutchinson(EtoG = F, b , c(3,6), sa1, 2)

#selected points in a raster 
post_track(pt_sel, un_wc1[[3]], sa1)


#quick implementation for environmental space exploration: 
c = data.frame(rasterToPoints (un_wc1[[1]])) #obtaining points from BIO1
d = c[c[,3]<1,] #selecting only those with temperatures less than 1

#visualizing those points in the environmental space 
e = e_space_back (data.frame(d[,1:2]),un_wc1,un_wc1, pflag = T)

#selecting group of points from the environmental space 
select1 = hutchinson(EtoG = F, e , c(3,6), sa1, 2) #A lot of NAs because e, has data from the entire raster 

#check surface under the uncertainty model 
post_track(select1, uncer_mod, sa1)

#PLOTTING MAXENT RESPONSE CURVES-------------------------------------------



#OTHERS-------------------------------------------

#ECOREGIONS

#STUFF TO DO---------------------------------------
#' complete transference_2d function 
#' add the other functions 
#' calculate AIC from different models
#' came up with an idea to present all this ideas beautifully in GITHUB... :(
#' ADD MATERIAL FOR THIS SECTION + REFERENCES FOR THE CORRESPONDING SECTION
#' 
