
library(raster)
library(rgdal)
library(rgeos)

land = raster('./NLCD_2016_Land_Cover_L48_20190424/NLCD2016.tif')
proj1 = proj4string(land)

coi.shp = spTransform(readOGR(dsn='./Insecta_COI_USA_deme_diversity_SpeciesSpecific.shp',layer='Insecta_COI_USA_deme_diversity_SpeciesSpecific',stringsAsFactors=F),CRS(proj1)) #shapefile with COI accessions
coi.shp2 = coi.shp[which(coi.shp@data$Scale==1 & (coi.shp@data$Order=='Ephemeroptera' | coi.shp@data$Order=='Plecoptera' | coi.shp@data$Order=='Trichoptera')),]

# Define raster reclassification scheme: Developed = 21-24 https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend
built.reclass = cbind(seq(0,95,1),rep(0,96))
built.reclass[match(21:24,built.reclass[,1]),2]= 1

coi.shp2@data$p_built = NA
for (i in 1:length(coi.shp2)){
	print(noquote(i))
	buffer1 = buffer(coi.shp2[i,],1000)
	land2 = crop(land,buffer1)
	land3 = mask(land2,buffer1)
	land4 = as.matrix(reclassify(land3,built.reclass))
	coi.shp2@data$p_built[i] = length(which(land4==1)) / length(which(!is.na(land4)))
}
writeOGR(coi.shp2,dsn='./EPT_COI_USAcont3_built.shp',layer='EPT_COI_USAcont3_built',driver='ESRI Shapefile',overwrite=T)
