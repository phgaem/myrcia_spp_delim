#Niche overlap analysis - Myrcia neoobscura species complex (Myrtaceae)
#Author: Giovani Carlos Andrella, PPG Biologia Vegetal, Unicamp (giovani.andrella@gmail.com)

#nstall.packages("phyloclim")
library(phyloclim)

#models
modelos <- "./ENMs"

# Transform .tiff files to .asc (requirement for analysis)

#Directory to save .asc files

Dirsave <- "./asc"
dir.create(Dirsave)

#Selecting ENMs

b.files <- list.files(pattern=".tif") 
names(b.files) <- gsub("_", " ", gsub(".tif$", "", basename(b.files)))
i = 1

for(i in 1:length(b.files)) {
  print(i)
  layer <- raster(b.files[i])
  writeRaster(layer, paste(Dirsave, paste0(names(b.files[i]), '.asc'), sep = '/'),
              format = 'ascii',
              overwrite = F)
}

#Selecting ENMs in .asc format
asc.files <- list.files(pattern=".asc") 


#Analysis with niche.overlap function

niche.overlap(asc.files) -> lyc.ro

#Rounding the values
round(lyc.ro, 3) -> lyc.ro

#Saving the results

write.table(lyc.ro, file="niche_overlap.csv", sep=",")