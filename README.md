# Gadgets-and-scripts-for-niche-modeling-
From 5th July to 9th July, 2021, we developed a workshop on ecological niche modeling for infectious diseases with A. Townsend Peterson, Marlon Cobos and Abdelghafar Alkishe, together with the American Society for Microbiology. In this repository you can find a complete implementation of the ideas discussed including: 
- Manipulation of occurrences.
- Manipulation of environmental variables.
- Model development based on Maxent (default, KUENM, ENMeval).
- Evaluation metrics for model selection. 
- Post-modeling (e.g., transmission risk). 

The main R script compiling all the information is called: gadgest_and_scripts_for_niche_modeling.R. You can open this script and follow all the instructions. 

You also need to open and read the information of other available scripts found in this repository. These scripts will allow you to develop some of the actions available inside the main script: 
- bin_test.R: function to implement the binomial test to evaluate thresholded models. 
- mask_ras.R: function to crop and mask raster variables using other raster variables (e.g., cropping population rasters with model output). 
- transfer_2d.R: function to transfer convex hulls defined in two dimensions from one area to another. 
- Functions for other manipulations of the environmental space can be found here: https://github.com/daromero-88/Environmental-sampling

Virtual species were not discussed in this course but allowed me to create dummy data to manipulate and implement the tools presented here. I used Leroy et al 2016 'virtualspecies' R package, crucial information to play with virtual species can be found via: 
- A recent review on the topic: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04385 
- The original manuscript presenting the package: https://borea.mnhn.fr/en/virtualspecies-r-package-generate-virtual-species-distributions
- A step by step super detailed tutorial for using the package: http://borisleroy.com/files/virtualspecies-tutorial.html

The virtual species was labeled 'Pathogen HK', the driver of the zombie outbreak from the game series Dead Island (https://deadisland.fandom.com/wiki/Pathogen_HK). For the present examples, we assumed: 
- An invasion of Pathogen HK to Southern Asia including regions of India.
- That abiotic conditions determine the presence of the pathogen in the environment.
