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


#SETTING THE WORKING DIRECTORY----------------------------

#' define folder 
#' Please select a folder in your computer called KUENM and 
#' use it as the main folder with the data provided

setwd('/PATH TO THE FOLDER/KUENM') 

#MODEL CALIBRATION - KUENM---------------------------------------
occ_joint <- "model_pts_def.csv"
occ_tra <- "in_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(0.1, 0.5, 1, 1.5, 2)
f_clas <- c("l", 'lq', 'lqp')
maxent_path <- 'PATH TO THE FOLDER/KUENM'
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
