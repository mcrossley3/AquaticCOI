

library(raster)
library(rgdal)
library(rgeos)

coi.shp = spTransform(readOGR(dsn='./Insecta_COI_USA_deme_diversity_SpeciesSpecific.shp',layer='Insecta_COI_USA_deme_diversity_SpeciesSpecific',stringsAsFactors=F),CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0')) #shapefile with COI accessions (transform CRS to match PRISM rasters)
coi.shp2 = coi.shp[which(coi.shp@data$Scale==1 & (coi.shp@data$Order=='Ephemeroptera' | coi.shp@data$Order=='Plecoptera' | coi.shp@data$Order=='Trichoptera')),]
D:\Landscape_Environmental_data\PRISM_data
usa.tmin = raster('/path/to/PRISM_tmin_30yr_normal_4kmM2_all_bil/PRISM_tmin_30yr_normal_4kmM2_annual_bil.bil')
usa.tmax = raster('/path/to/PRISM_tmax_30yr_normal_4kmM2_all_bil/PRISM_tmax_30yr_normal_4kmM2_annual_bil.bil')
usa.precip = raster('/path/to/PRISM_data/PRISM_ppt_30yr_normal_4kmM2_all_bil/PRISM_ppt_30yr_normal_4kmM2_annual_bil.bil')
usa.vpdmax = raster('/path/to/PRISM_vpdmax_30yr_normal_4kmM2_all_bil/PRISM_vpdmax_30yr_normal_4kmM2_annual_bil.bil')
x11()
par(mfrow=c(2,2))
plot(usa.tmin)
plot(usa.tmax)
plot(usa.precip)
plot(usa.vpdmax)

# Extract climate metrics for each accession
coi.shp2@data$tmin = coi.shp2@data$tmax = coi.shp2@data$precip = coi.shp2@data$vpdmax = NA
for (i in 1:nrow(coi.shp2)){
	print(noquote(i))
	p1 = buffer(coi.shp2[i,],1000)
	coi.shp2@data$tmin[i] = extract(usa.tmin,p1,fun=mean)
	coi.shp2@data$tmax[i] = extract(usa.tmax,p1,fun=mean)
	coi.shp2@data$precip[i] = extract(usa.precip,p1,fun=mean)
	coi.shp2@data$vpdmax[i] = extract(usa.vpdmax,p1,fun=mean)
}
writeOGR(coi.shp2,dsn='./EPT_COI_USAcont3_climate.shp',layer='EPT_COI_USAcont3_climate',driver='ESRI Shapefile',overwrite=T)

