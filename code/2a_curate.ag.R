
library(raster)
library(rgdal)
library(rgeos) #for gArea
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' #Albers projection for gArea

source('C:/Users/mcros/Desktop/Postdoc UGA/AgChange/import.AgChange.R') #import shapefile containing full ag change dataset. Available from: https://www.openicpsr.org/openicpsr/project/115795/version/V3/view
crops = c('Brl','Bkw','Crn','Ctn','Flx','Hay','Oat','Pnt','Ptt','Pls','Ric','Rye','Sgm','Soy','Swt','Sgc','Tbc','Wht')
years = c(1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1959,1974,1982,1992,2002,2012,2017)

coi.shp = readOGR(dsn='./Insecta_COI_USA_deme_diversity_SpeciesSpecific.shp',layer='Insecta_COI_USA_deme_diversity_SpeciesSpecific',stringsAsFactors=F) #shapefile with COI accessions
coi.shp2 = coi.shp[which(coi.shp@data$Scale==1 & (coi.shp@data$Order=='Ephemeroptera' | coi.shp@data$Order=='Plecoptera' | coi.shp@data$Order=='Trichoptera')),]

# Calculate total proportion county area in crop lands for each year
total.crop = c()
for (y in 1:length(years)){
	year.add = apply(counties@data[,grep(years[y],colnames(counties@data))[1:18]],1,function(x){sum(as.numeric(x),na.rm=T)})
	total.crop = cbind(total.crop,year.add)
}
colnames(total.crop) = years
crop.totals = data.frame('FIPS'=counties@data$FIPS,total.crop,stringsAsFactors=F)

# Extract land cover metrics for each accession
coi.shp2@data$crop_pre1950 = coi.shp2@data$crop_post1940 = NA
for (i in 1:nrow(coi.shp2)){
	print(noquote(i))
	p1 = buffer(coi.shp2[i,],1000)
	p1.int = raster::intersect(p1,counties)
	fips1 = p1.int@data$FIPS
	if (length(p1.int)==1){
		#deme contained by one county, no need to use areal weighting to get average cropland
		c.totals = crop.totals[which(crop.totals$FIPS==fips1),] #Total crop lands (average proportion county area in crops)
		coi.shp2@data$crop_pre1950[i] = mean(as.numeric(c.totals[,grep(1840,colnames(c.totals)):grep(1940,colnames(c.totals))])) #average total crop land cover between 1840-1940 (pre-WWII)
		coi.shp2@data$crop_post1940[i] = mean(as.numeric(c.totals[,grep(1950,colnames(c.totals)):grep(2017,colnames(c.totals))])) #average total crop land cover between 1840-1940 (pre-WWII)
	} else {
		#multiple counties intersected - use areal weighting to get average cropland
		areas1 = gArea(p1.int,byid=T)
		weights1 = areas1 / sum(areas1) #proportions of deme area overlapping counties
		c.totals = crop.totals[match(fips1,crop.totals$FIPS),]
		for (j in 1:length(fips1)){ #weight cropland proportions by area of deme overlap
			c.totals[j,-1] = c.totals[j,-1] * weights1[j]
		}
		c.totals2 = colSums(c.totals[,-1]) #area-weighted average for each county*year
		coi.shp2@data$crop_pre1950[i] = mean(as.numeric(c.totals2[grep(1840,names(c.totals2)):grep(1940,names(c.totals2))])) #average total crop land cover between 1840-1940 (pre-WWII)
		coi.shp2@data$crop_post1940[i] = mean(as.numeric(c.totals2[grep(1950,names(c.totals2)):grep(2017,names(c.totals2))])) #average total crop land cover between 1840-1940 (pre-WWII)
	}
}
writeOGR(coi.shp2,dsn='./EPT_COI_USAcont3_ag.shp',layer='EPT_COI_USAcont3_ag',driver='ESRI Shapefile',overwrite=T)


####
# Create shapefile with cropland pre-WWII

source('C:/Users/mcros/Desktop/Postdoc UGA/AgChange/import.AgChange.R')

crops = c('Brl','Bkw','Crn','Ctn','Flx','Hay','Oat','Pnt','Ptt','Pls','Ric','Rye','Sgm','Soy','Swt','Sgc','Tbc','Wht')
years = c(1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1959,1974,1982,1992,2002,2012,2017)
cropyears = c()
for (i in 1:length(crops)){
	for (j in 1:11){
		cropyears = c(cropyears,paste0(crops[i],years[j]))
	}
}
cropyears2 = c()
for (i in 1:length(crops)){
	for (j in 12:19){
		cropyears2 = c(cropyears2,paste0(crops[i],years[j]))
	}
}

counties@data$cr1840 = counties@data$cr1850 = counties@data$cr1860 = counties@data$cr1870 = counties@data$cr1880 = counties@data$cr1890 = counties@data$cr1900 = counties@data$cr1910 = counties@data$cr1920 = counties@data$cr1930 = counties@data$cr1940 = -999
counties@data$cr1950 = counties@data$cr1959 = counties@data$cr1974 = counties@data$cr1982 = counties@data$cr1992 = counties@data$cr2002 = counties@data$cr2012 = counties@data$cr2017 = -999
counties@data$pre1950 = counties@data$pre1960 = counties@data$pre1970 = counties@data$post1940 = counties@data$post1950 = counties@data$post1960 = counties@data$post1970 = -999
for (i in 1:length(counties)){
	print(noquote(i))
	counties@data$cr1840[i] = sum(counties@data[i,match(cropyears[grep('1840',cropyears)],colnames(counties@data))])
	counties@data$cr1850[i] = sum(counties@data[i,match(cropyears[grep('1850',cropyears)],colnames(counties@data))])
	counties@data$cr1860[i] = sum(counties@data[i,match(cropyears[grep('1860',cropyears)],colnames(counties@data))])
	counties@data$cr1870[i] = sum(counties@data[i,match(cropyears[grep('1870',cropyears)],colnames(counties@data))])
	counties@data$cr1880[i] = sum(counties@data[i,match(cropyears[grep('1880',cropyears)],colnames(counties@data))])
	counties@data$cr1890[i] = sum(counties@data[i,match(cropyears[grep('1890',cropyears)],colnames(counties@data))])
	counties@data$cr1900[i] = sum(counties@data[i,match(cropyears[grep('1900',cropyears)],colnames(counties@data))])
	counties@data$cr1910[i] = sum(counties@data[i,match(cropyears[grep('1910',cropyears)],colnames(counties@data))])
	counties@data$cr1920[i] = sum(counties@data[i,match(cropyears[grep('1920',cropyears)],colnames(counties@data))])
	counties@data$cr1930[i] = sum(counties@data[i,match(cropyears[grep('1930',cropyears)],colnames(counties@data))])
	counties@data$cr1940[i] = sum(counties@data[i,match(cropyears[grep('1940',cropyears)],colnames(counties@data))])
	counties@data$cr1950[i] = sum(counties@data[i,match(cropyears2[grep('1950',cropyears2)],colnames(counties@data))])
	counties@data$cr1959[i] = sum(counties@data[i,match(cropyears2[grep('1959',cropyears2)],colnames(counties@data))])
	counties@data$cr1974[i] = sum(counties@data[i,match(cropyears2[grep('1974',cropyears2)],colnames(counties@data))])
	counties@data$cr1982[i] = sum(counties@data[i,match(cropyears2[grep('1982',cropyears2)],colnames(counties@data))])
	counties@data$cr1992[i] = sum(counties@data[i,match(cropyears2[grep('1992',cropyears2)],colnames(counties@data))])
	counties@data$cr2002[i] = sum(counties@data[i,match(cropyears2[grep('2002',cropyears2)],colnames(counties@data))])
	counties@data$cr2012[i] = sum(counties@data[i,match(cropyears2[grep('2012',cropyears2)],colnames(counties@data))])
	counties@data$cr2017[i] = sum(counties@data[i,match(cropyears2[grep('2017',cropyears2)],colnames(counties@data))])
	counties@data$pre1950[i] = mean(as.numeric(counties@data[i,match(c('cr1840','cr1850','cr1860','cr1870','cr1890','cr1900','cr1910','cr1920','cr1930','cr1940'),colnames(counties@data))]))
	counties@data$pre1960[i] = mean(as.numeric(counties@data[i,match(c('cr1840','cr1850','cr1860','cr1870','cr1890','cr1900','cr1910','cr1920','cr1930','cr1940','cr1950'),colnames(counties@data))]))
	counties@data$pre1970[i] = mean(as.numeric(counties@data[i,match(c('cr1840','cr1850','cr1860','cr1870','cr1890','cr1900','cr1910','cr1920','cr1930','cr1940','cr1950','cr1959'),colnames(counties@data))]))
	counties@data$pre1980[i] = mean(as.numeric(counties@data[i,match(c('cr1840','cr1850','cr1860','cr1870','cr1890','cr1900','cr1910','cr1920','cr1930','cr1940','cr1950','cr1959','cr1974'),colnames(counties@data))]))
	counties@data$post1940[i] = mean(as.numeric(counties@data[i,match(c('cr1950','cr1959','cr1974','cr1982','cr1992','cr2002','cr2012','cr2017'),colnames(counties@data))]))
	counties@data$post1950[i] = mean(as.numeric(counties@data[i,match(c('cr1959','cr1974','cr1982','cr1992','cr2002','cr2012','cr2017'),colnames(counties@data))]))
	counties@data$post1960[i] = mean(as.numeric(counties@data[i,match(c('cr1974','cr1982','cr1992','cr2002','cr2012','cr2017'),colnames(counties@data))]))
	counties@data$post1970[i] = mean(as.numeric(counties@data[i,match(c('cr1982','cr1992','cr2002','cr2012','cr2017'),colnames(counties@data))]))
}
counties@data$diff50 = counties@data$post1940 - counties@data$pre1950
counties@data$diff60 = counties@data$post1950 - counties@data$pre1960
counties@data$diff70 = counties@data$post1960 - counties@data$pre1970
counties@data$diff80 = counties@data$post1970 - counties@data$pre1980

counties2 = counties
counties2@data = counties2@data[,476:515]
writeOGR(counties2,dsn='./counties_ag.shp',layer='counties_ag',driver='ESRI Shapefile',overwrite=T,verbose=F)

