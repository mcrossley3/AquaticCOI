
# Curate cropland and insecticide use data from USDA census of agriculture 2017

source('C:/Users/mcros/Desktop/Postdoc UGA/AgChange/import.AgChange.R')

state.files = list.files('/path/to/ag/census/files/',full.names=T) #https://www.nass.usda.gov/Publications/AgCensus/2017/Full_Report/Census_by_State/index.php

# Define function that reformats census values
reformat.census = function(xv){
	xv = gsub(',','',xv) # remove "," & convert to numeric
	xv = gsub('-',0,xv) # replace "-" with 0
	xv = gsub(';D;',0,xv) # replace "(D)" with 0
	xv = gsub(';Z;',0,xv) # replace "(Z)" with 0 (Z implies < 1 acre of crop)
	return(xv)
}

# Define function that grabs crop values from counties
get.crops = function(crop.start,crop.stop,co.names=county.names,TXT2=txt2,st.name=state1){
	# crop.start = string that defines beginning of county-level data for this crop
	# crop.stop = string that defines end of county-level data for this crop
	# co.names = a vector of county names in this state
	# TXT2 = the txt2 object subsetted for this state
	crop.txt = TXT2[grep(crop.start,TXT2)[1]:grep(crop.stop,TXT2)[1]]
	crop.counties = gsub('\\s+','',gsub('[.]','',apply(array(crop.txt),1,function(x){strsplit(x,':')[[1]][1]})))
	out = rep(0,length(co.names))
	for (iz in 1:n.counties){
#		print(iz)
		c1 = gsub(' ','',gsub('[.]','',co.names[iz]))
		cpos1 = which(crop.counties==c1)
		if (c1==gsub(' ','',st.name)){
			cpos1 = cpos1[2] #use 2nd occurence of county name, because the 1st is actually numbers for the entire state
		} else {
		}
		if (length(cpos1)==0){
		} else if (length(cpos1)>1){
			print('          >1 match!.......................')
		} else {
			crop.txt2 = strsplit(strsplit(crop.txt[cpos1],':')[[1]][2],"\\s+")[[1]]
			out[iz] = crop.txt2[3]
		}
	}
	# Reformat
	out = reformat.census(out)
	return(as.numeric(out))
}

census.all = c()
for (s in 1:length(state.files)){
	state1 = gsub('.txt','',strsplit(state.files[s],'/')[[1]][3])
	print(noquote(state1))
	con = file(state.files[s],'r')
	txt = readLines(con)
	close(con)
	t11 = txt[grep('Table 11.',txt)]
	txt2 = txt[-c(1:(which(txt==t11[3])+2))] #Table trimmed of preceding state-level data
	txt2 = gsub('\\(',';',txt2)
	txt2 = gsub('\\)',';',txt2)

	# County names
	item = gsub(' ','',txt2[grep('Item',txt2)])
	item2 = item[1:((grep('ItemNonresponse',item)-1))] #lines with county names (headers of each table with county-level data)
	county.names = gsub('[.]','. ',unique(c(unlist(apply(array(item2),1,function(x){strsplit(x,':')[[1]]})))))
	if (state1=='Arkansas' | state1=='Idaho' | state1=='Iowa' | state1=='New York' | state1=='Oklahoma' | state1=='Utah'){ #Arkansas also has a county named "Arkansas"
		county.names = county.names[-which(county.names=='Item')] #names of counties
	} else {
		county.names = county.names[-which(county.names=='Item' | county.names==state1 | county.names==gsub(' ','',state1))] #names of counties
	}
	n.counties = length(county.names)

	# Harvested cropland
	txt.cropland = txt[grep('Table 8.  Farms, Land in Farms, Value of Land and Buildings, and Land Use:  2017 and 2012',txt)[1]:grep('Table 9.  Harvested Cropland by Size of Farm and Acres Harvested:  2017 and 2012',txt)[1]]
	txt.cropland = gsub('\\(',';',txt.cropland)
	txt.cropland = gsub('\\)',';',txt.cropland)
	cropland = gsub('acres, 2017:','',txt.cropland[grep('Total cropland ......................................farms, 2017:',txt.cropland)+2])
	cropland2 = unlist(apply(array(cropland),1,function(x){strsplit(x,"\\s+")[[1]]}))
	cropland3 = cropland2[which(nchar(cropland2)>0)][2:(n.counties+1)] #vector of no. cattle for each county
	cropland = as.numeric(reformat.census(cropland3))

	# Acres treated with insecticide
	txt.insecticide = txt[grep('Table 40.  Fertilizers and Chemicals Applied:  2017 and 2012',txt)[1]:grep('Table 41.  Land Use Practices:  2017 and 2012',txt)[1]]
	insecticide = gsub('acres, 2017:','',txt.insecticide[grep('Insects',txt.insecticide)+2])
	insecticide2 = unlist(apply(array(insecticide),1,function(x){strsplit(x,"\\s+")[[1]]}))
	insecticide3 = insecticide2[which(nchar(insecticide2)>0)][2:(n.counties+1)] #vector of acres treated to control insects for each county
	insecticide = as.numeric(reformat.census(insecticide3))

	add1 = data.frame('State'=rep(state1,length(county.names)),'County'=county.names,'Cropland'=cropland,'Insecticide'=insecticide,stringsAsFactors=F)
	census.all = data.frame(rbind(census.all,add1),stringsAsFactors=F)
}

# Match rows to 1840-2012 dataframe
stkey = read.table('./StateFIPS_STCTY.txt',sep='\t',as.is=T,check.names=F,header=T)
ctys = gsub(' City','',gsub(' Parish','',gsub(' County','',counties@data$COUNTY)))
sts = apply(array(counties@data$STCTY),1,function(x){y=strsplit(x,'_')[[1]][1];z=stkey$STATE[which(stkey$STCTY==y)];if(length(z)==0){return(NA)}else{return(z)}})
dat = data.frame('County'=ctys,'State'=sts,stringsAsFactors=F)

# Get row position of 2017 census records in counties (final data for 1840-2012)
con = file('county.match.key.txt','r')
match.key1 = readLines(con)
close(con)
match.key = data.frame(t(apply(array(match.key1),1,function(x){y=strsplit(x,'\t')[[1]];if(length(y<6)){c(y,rep(NA,6-length(y)))}else{y}})),stringsAsFactors=F)
colnames(match.key) = match.key[1,]
match.key = match.key[-1,]

new.pos = apply(dat,1,function(x){y=which(census.all$State==x[2] & census.all$County==x[1]);if(length(y)==0){return(NA)}else{return(y)}})
miss.pos = which(is.na(new.pos))
census.orphans = counties@data[miss.pos,]
for (i in 1:nrow(census.orphans)){
	orphan.state = census.orphans$STATE[i]
	orphan.county = gsub(' Parish','',gsub(' County','',census.orphans$COUNTY[i]))
#	print(noquote(c(i,orphan.county)))
	match1 = match.key$MatchID[which(match.key$Counties.State2==orphan.state & match.key$Counties.County==orphan.county)]
	if (length(match1)==0){
		if (orphan.state=='MO' & orphan.county=='Saint Louis City'){
			add.pos = NA
		} else if (orphan.state=='NV'){
			match1 = match.key$MatchID[which(match.key$Counties.State2==orphan.state & match.key$Counties.County==gsub(' City','',orphan.county))]			
			coco = match.key$Census2017.County[which(match.key$MatchID==match1)]
			cost = match.key$Census2017.State[which(match.key$MatchID==match1)]
			add.pos = which(census.all$State==cost & census.all$County==coco)
		} else if (orphan.state=='VA'){
			match1 = match.key$MatchID[which(match.key$Counties.State2==orphan.state & match.key$Counties.County==gsub(' City','',orphan.county))]
			if (is.na(match1)){
				print(noquote(paste(orphan.state,orphan.county)))
				add.pos = NA
			} else {
				coco = match.key$Census2017.County[which(match.key$MatchID==match1)]
				cost = match.key$Census2017.State[which(match.key$MatchID==match1)]
				add.pos = which(census.all$State==cost & census.all$County==coco)
			}
		} else {
			print(c(i,'exception'))
		}
	} else {
		if (orphan.state=='DC'){
			add.pos = NA
		} else {
			coco = match.key$Census2017.County[which(match.key$MatchID==match1)]
			cost = match.key$Census2017.State[which(match.key$MatchID==match1)]
			add.pos = which(census.all$State==cost & census.all$County==coco)
		}
	}
	new.pos[miss.pos[i]] = add.pos
}

add.census2017 = census.all[new.pos,] #re-arranged 2017 census data

check = cbind(counties@data$STATE,counties@data$COUNTY,add.census2017$County) #make sure county names align in difficult cases
check[which(check[,1]=='VA'),]
check[which(check[,1]=='NV'),]

# Convert animals to density & crops to prop. area
counties2 = counties
counties2@data = data.frame(cbind(counties2@data[,476:484],add.census2017[,-c(1:2)]),stringsAsFactors=F)
counties2@data$p.Cropland = (counties2@data$Cropland/247.105) / as.numeric(counties2@data$AREA_KM)
counties2@data$p.Insecticide = counties2@data$Insecticide / counties2@data$Cropland #11 cases where p.Cropland > 1
writeOGR(counties2,dsn='./US_county_2017_insecticides.shp',layer='US_county_2017_insecticides',driver='ESRI Shapefile',overwrite=T)


###############################################
# Get insecticide use data from NAQWA

library('rgdal')

#NOTE: this code takes a long time
pesticide.files = list.files('./USGS pesticide data/txt',full.names=T)
years = apply(array(pesticide.files),1,function(x){sub('.txt','',strsplit(x,'pesticides')[[1]][2])})
counties = readOGR(dsn='./US_counties_2012_fips.shp',layer='US_counties_2012_fips',stringsAsFactors=F)
colnames(counties@data) = 'FIPS'
for (i in 23:length(pesticide.files)){
	chem.data = read.table(pesticide.files[i],sep='\t',as.is=T,check.names=F,header=T)
	chem.data$STATE_FIPS_CODE = apply(array(chem.data$STATE_FIPS_CODE),1,function(x){if(nchar(x)==1){return(paste0('0',x))}else{return(x)}}) #ensure state FIPS has 2 digits
	chem.data$COUNTY_FIPS_CODE = apply(array(chem.data$COUNTY_FIPS_CODE),1,function(x){if(nchar(x)==1){return(paste0('00',x))}else if(nchar(x)==2){return(paste0('0',x))}else{return(x)}}) #ensure county FIPS has 3 digits
	chem.data$FIPS = paste(chem.data$STATE_FIPS_CODE,chem.data$COUNTY_FIPS_CODE,sep='') #generate FIPS to match counties shapefile
	chemicals = unique(chem.data$COMPOUND)
	counties.high = counties
	counties.low = counties
	lc = length(chemicals)
	high.add = low.add = counties@data #becomes a dataframe to merge with counties attribute table
	for (j in 1:length(chemicals)){
		print(paste0(j,' of ',lc,'   ',chemicals[j]))
		high.data = low.data = rep(-999,lc) #vector of high chemical values, to be added to high.add
		for (k in 1:length(counties)){
			cdat = chem.data[which(chem.data$COMPOUND==chemicals[j] & chem.data$FIPS==counties@data[k,]),]
			if (nrow(cdat)==0){ #deal with counties that lack data for chemical j
				high.data[k] = NA
				low.data[k] = NA
			} else {
				high.data[k] = cdat$EPEST_HIGH_KG
				low.data[k] = cdat$EPEST_LOW_KG
			}
		}
		high.add = cbind(high.add,high.data)
		low.add = cbind(low.add,low.data)
	}
	colnames(high.add) = colnames(low.add) = c('FIPS',chemicals)
	counties.high@data = high.add
	counties.low@data = low.add
	writeOGR(counties.high,dsn=paste0('./NAQWA/US_counties_pesticides_',years[i],'_high.shp'),layer=paste0('US_counties_pesticides_',years[i],'_high'),driver="ESRI Shapefile",overwrite=T)
	writeOGR(counties.low,dsn=paste0('./NAQWA/US_counties_pesticides_',years[i],'_low.shp'),layer=paste0('US_counties_pesticides_',years[i],'_low'),driver="ESRI Shapefile",overwrite=T)
}

# Create new shapefiles with totals for pesticide types, & totals for pesticide groups
pesticide.files = list.files('./NAQWA',full.names=T)
years = apply(array(pesticide.files),1,function(x){sub('.txt','',strsplit(x,'pesticides')[[1]][2])})
years2 = c('92','93','94','95','96','97','98','99','00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17') #WARNING: being lazy. Order must match years
counties = readOGR(dsn='./US_counties_2012_fips.shp',layer='US_counties_2012_fips',stringsAsFactors=F)
fips = counties@data[,1]
key = read.csv('./pesticide_key.csv',as.is=T,check.names=F,header=T)
key$Group[which(key$Group=='')] = 'other'
types = c('Insecticide','Herbicide','Fungicide','Nematicide')
types2 = c('I','H','F','N')
add.data = fips
add.header = 'FIPS'
# Add totals for pesticide types
for (p in 1:length(pesticide.files)){
	print(years[p])
	chem.data = read.table(pesticide.files[p],sep='\t',as.is=T,check.names=F,header=T)
	for (i in 1:length(types)){
		print(types[i])
		type.chems = key$Chemical[which(key$Type==types[i])]
		tcpos = match(chem.data$COMPOUND,type.chems)
		chem.data2 = chem.data[which(!is.na(tcpos)),]
		fips.chem.data2 = as.character(apply(chem.data2[,3:4],1,function(x){
			x1=as.character(x[1]); x2=as.character(x[2])
			if(nchar(x1)<2){y1=paste0(0,x1)}else{y1=x1};
			if(nchar(x2)<3){m1=3-nchar(x2); m2=paste0(rep(0,m1),collapse=''); y2=paste0(m2,x2)}else{y2=x2};
			return(paste0(y1,y2))})) #fix FIPS codes
		# Add chemical data
		add.chem = c()
		for (f in 1:length(fips)){
			fips.pos = which(fips.chem.data2==fips[f])
			if (length(fips.pos)>0){
				add.chem = c(add.chem, sum(chem.data2[which(fips.chem.data2==fips[f]),6],na.rm=T))
			} else {
				add.chem = c(add.chem,NA) #chem is not reported in this county in this year
			}
		}
		add.header = c(add.header,paste(types2[i],years2[p],sep='_'))
		add.data = cbind(add.data,add.chem)
		# Add chemical group data
		group.data = key[which(key$Type==types[i]),]
		groups = unique(group.data$Group)
		print('...adding group data...')
		for (j in 1:length(groups)){
			add.group=c()
			group.chems = key$Chemical[which(key$Type==types[i] & key$Group==groups[j])]
			gcpos = match(chem.data2$COMPOUND,group.chems)
			chem.data3 = chem.data2[which(!is.na(gcpos)),]
			fips.chem.data3 = fips.chem.data2[which(!is.na(gcpos))]
			for (f2 in 1:length(fips)){
				fips.pos2 = which(fips.chem.data3==fips[f2])
				if (length(fips.pos2)>0){
					add.group = c(add.group, sum(chem.data3[which(fips.chem.data3==fips[f2]),6],na.rm=T))
				} else {
					add.group = c(add.group,NA) #chem group is not reported in this county in this year
				}
			}
			add.header = c(add.header,paste0(types2[i],groups[j],'_',years2[p]))
			add.data = cbind(add.data,add.group)
		}
	}
}

# Write all data to one text file
counties2 = counties
counties2@data = data.frame(add.data,stringsAsFactors=F)
colnames(counties2@data) = add.header
# Ensure numbers are written as numbers
for (i in 2:ncol(counties2@data)){
	counties2@data[,i] = as.numeric(counties2@data[,i])
}
write.table(counties2@data,'./US_counties_pesticides_totals.txt',sep='\t',quote=F,row.names=F)


###############################################
# Join insecticide use data to demes

library(rgdal)
library(rgeos) #for spTransform()
library(raster) #for intersects()
library(sf)
proj1 = '+proj=longlat +datum=NAD83 +no_defs'

# Genetic diversity dataset
coi.shp = spTransform(readOGR(dsn='./Insecta_COI_USA_deme_diversity_SpeciesSpecific.shp',layer='Insecta_COI_USA_deme_diversity_SpeciesSpecific',stringsAsFactors=F),CRS(proj1)) #shapefile with COI accessions
coi.shp2 = coi.shp[which(coi.shp@data$Scale==1 & (coi.shp@data$Order=='Ephemeroptera' | coi.shp@data$Order=='Plecoptera' | coi.shp@data$Order=='Trichoptera')),]

# Link to HUC10 boundaries
huc10 = readOGR('./WBDHU10.shp',layer='WBDHU10',stringsAsFactors=F,verbose=F) #https://www.usgs.gov/core-science-systems/ngp/national-hydrography/watershed-boundary-dataset?qt-science_support_page_related_con=4#qt-science_support_page_related_con
over1 = over(coi.shp2,huc10,byid=T)
coi.shp2@data$huc10 = over1$huc10
coi.shp2$huc10.name = over1$name
coi.shp2$huc10.km2 = over1$areasqkm

# Insecticide use in kg/ha 1992-2017 https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/
USDA.insecticide = readOGR(dsn='./US_county_2017_insecticides.shp',layer='US_county_2017_insecticides',verbose=F,stringsAsFactors=F)
pesticides = read.table('./US_counties_pesticides_totals.txt',sep='\t',as.is=T,check.names=F,header=T)
pesticides$FIPS = apply(array(pesticides$FIPS),1,function(x){if(nchar(x)<5){return(paste0(0,x))}else{return(x)}})
merge1 = merge(USDA.insecticide,pesticides,by='FIPS')
# Convert kg to kg/ha cropland
cols = paste('I',c(92:99,'01','02','03','04','05','06','07','08','09',10:17),sep='_')
for (j in 1:length(cols)){
	merge1@data[,which(colnames(merge1@data)==cols[j])] = merge1@data[,which(colnames(merge1@data)==cols[j])] / (merge1@data$Croplnd*0.404686)
}
merge1@data$I.avg = apply(merge1@data[,match(cols,colnames(merge1@data))],1,mean,na.rm=T) #1992-2017 average insecticide use
merge1 = spTransform(merge1,CRS(proj1))

# Extract insecticide data for each deme
coi.shp2@data$kgha.insecticide = NA
for (i in 1:nrow(coi.shp2)){
	print(noquote(i))
	p1 = buffer(coi.shp2[i,],1000)
	p1.int = raster::intersect(p1,merge1)
	if (length(p1.int)==1){
		#deme contained by one county, no need to use areal weighting to get average
		coi.shp2@data$kgha.insecticide[i] = p1.int@data$I.avg
	} else {
		#multiple counties intersected - use areal weighting to get average
		areas1 = gArea(p1.int,byid=T)
		weights1 = areas1 / sum(areas1) #proportions of deme area overlapping counties
		coi.shp2@data$kgha.insecticide[i] = sum(p1.int@data$I.avg * weights1)
	}
}
writeOGR(coi.shp2,dsn='./EPT_COI_USAcont3_insecticides.shp',layer='EPT_COI_USAcont3_insecticides',driver='ESRI Shapefile',overwrite=T)
