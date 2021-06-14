
library(rgdal)

# Species-specific demes
built.shp = readOGR(dsn='./EPT_COI_USAcont3_built.shp',layer='EPT_COI_USAcont3_built',stringsAsFactors=F)
ag.shp = readOGR(dsn='./EPT_COI_USAcont3_ag.shp',layer='EPT_COI_USAcont3_ag',stringsAsFactors=F)
clim.shp = readOGR(dsn='./EPT_COI_USAcont3_climate.shp',layer='EPT_COI_USAcont3_climate',stringsAsFactors=F)
icide.shp = readOGR(dsn='./EPT_COI_USAcont3_insecticides.shp',layer='EPT_COI_USAcont3_insecticides',stringsAsFactors=F)

merged = ag.shp
merged@data$built = built.shp$p_built
merged@data$tmax = clim.shp$tmax
merged@data$tmin = clim.shp$tmin
merged@data$vpdmax = clim.shp$vpdmax
merged@data$precip = clim.shp$precip
merged@data$icide = icide.shp$kgh_nsc
colnames(merged@data)[match(c('cr_1940','cr_1950'),colnames(merged@data))] = c('ag_post','ag_pre')
merged@data$Lon = merged@coords[,1]
merged@data$Lat = merged@coords[,2]
dat = merged@data; dim(dat) #2,017 demes
thresh = sd(dat$Nuc_Div)*4 #0.07059178
dat2 = dat[which(dat$Nuc_Div<thresh),]; dim(dat2) #1,986 demes
merged2 = merged[which(dat$Nuc_Div<thresh),]

write.table(dat2,'./EPT_COI_USA_div.txt',sep='\t',quote=F,row.names=F)
writeOGR(merged2,dsn='./diversity/EPT_COI_USA_div.shp',layer='EPT_COI_USA_div',driver='ESRI Shapefile',overwrite=T)

