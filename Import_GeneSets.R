# Import Gene Sets that have been made elsewhere
# E. Lamont
# 7/29/25


###########################################################
###################### LOAD GENE SETS #####################

# To put them in a list of lists
rda_files <- list.files("Data/GeneSet_Data", pattern = "\\.rda$", full.names = TRUE)
allGeneSetList <- list()

for(file in rda_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  env <- new.env()
  load(file, envir = env)  # loads allGeneSets into env
  allGeneSetList[[file_name]] <- env$allGeneSets  # store it in our list
  rm(env)
}

# Now update the names for each gene set in the list of lists:
allGeneSetList <- lapply(allGeneSetList, function(gset) {
  if (!is.null(gset)) {
    names(gset) <- gsub("<.*", "", names(gset))
  }
  return(gset)
})


###########################################################
################# TO MAKE NEW GENE SETS ###################

# First! Make the .csv file by hand. Needs a "Gene" and a "GeneSet" column

# XXX_GeneSets <- read.csv("GeneSet_Data/XXX_GeneSets.csv")

# Create a list where each GeneSet is a named element
# Needs to be called allGeneSets so it is easier to load with all the others
# allGeneSets <- split(XXX_GeneSets$Gene, XXX_GeneSets$GeneSet)

# SAVE AS RDA FOR LATER 
# save(allGeneSets, file = "GeneSet_Data/XXXGeneSets.rda")
