
library('rentrez')

set_entrez_key('<your key here>')
Sys.getenv('ENTREZ_KEY')

orders = c('Protura', 'Collembola', 'Diplura', 'Microcoryphia', 'Thysanura', 'Ephemeroptera', 'Odonata', 'Orthoptera', 'Phasmatodea', 'Grylloblattodea', 'Mantophasmatodea', 'Dermaptera', 'Plecoptera', 'Embiidina', 'Zoraptera', 'Isoptera', 'Mantodea', 'Blattodea', 'Hemiptera', 'Heteroptera', 'Thysanoptera', 'Psocoptera', 'Phthiraptera', 'Coleoptera', 'Neuroptera', 'Hymenoptera', 'Trichoptera', 'Lepidoptera', 'Siphonaptera', 'Mecoptera', 'Strepsiptera', 'unclassified Diptera','Brachycera','Nematocera') #source: http://www.entnemdept.ufl.edu/choate/insect_orders.htm
# parse.taxon doesn't properly handle Diptera, so I've broken it down to suborder.

# Define recursive function that parses insect orders into lower taxonomic levels such that retmax (99,999 records) is not exceeded during nuccore search
parse.taxon = function(scientificname,tx.out=c(),tx.hits=c()){
# INPUT: scientificnames; a vector of scientific names directly below an upper taxonomic level (e.g. suborders within an order)
# VALUE: list of scientific names
	tx1 = scientificname
	NS = entrez_search(db='nuccore',term=paste0(tx1,'[ORGN] AND (COI[GENE] OR COX1[GENE])'),use_history=T,retmax=99999)
	nh = NS$count
	ni = length(NS$ids)
#	print(tx1) #troubleshooting
	if (nh > ni){ #search returned too many records; parse order down taxonomic levels until retmax is not exceeded
		TSearch = entrez_search(db='taxonomy',term=paste0(tx1,'[NXLV]'),use_history=T,retmax=99999)
		if (tx1=='unclassified Cecidomyiidae'){ #deal with special case: Skip recursion
			tx.out = c(tx.out,tx1)
			tx.hits = c(tx.hits,nh)
		} else {
			TSearch.sum = entrez_summary(db='taxonomy',web_history=TSearch$web_history)
			scientificnames2 = unlist(lapply(TSearch.sum,function(x){x$scientificname}))
			for (sn in 1:length(scientificnames2)){
				parse.taxon(scientificnames2[sn])
			}
		}
	} else { #hits = ids; return scientific name
		tx.out = c(tx.out,tx1)
		tx.hits = c(tx.hits,nh)
	}
	return(cbind(tx.out,tx.hits))
}

# Get taxonomic levels that can be converted to UIDs for capturing insect COI sequences and collection attributes from Genbank
taxa.table1 = c()
taxa.table2 = c()
for (o in 1:length(orders)){
	order1 = orders[o]
	print(order1)
	order.nuc.search = entrez_search(db='nuccore',term=paste0(order1,'[ORGN] AND (COI[GENE] OR COX1[GENE])'),use_history=T,retmax=99999)
	n.hits = order.nuc.search$count
	n.ids = length(order.nuc.search$ids)
	if (n.hits > n.ids){ #search returned too many records; parse order down taxonomic levels until retmax is not exceeded
		taxon.search = entrez_search(db='taxonomy',term=paste0(order1,'[NXLV]'),use_history=T,retmax=99999)
		taxon.search.sum = entrez_summary(db='taxonomy',web_history=taxon.search$web_history)
		scinames = unlist(lapply(taxon.search.sum,function(x){x$scientificname}))
		taxa.names = c()
		taxa.hits = c()
		for (s in 1:length(scinames)){
			taxa = parse.taxon(scinames[s])
			taxa.names = c(taxa.names,taxa[,1])
			taxa.hits = c(taxa.hits,taxa[,2])
		}
		taxa.table1 = c(taxa.table1,taxa.names)
		taxa.table2 = c(taxa.table2,taxa.hits)		
	} else {
		taxa.table1 = c(taxa.table1,order1)
		taxa.table2 = c(taxa.table2,n.hits)
	}
}
taxa.table = data.frame('Taxon'=taxa.table1,'N.accessions'=taxa.table2)
unique.tax = unique(as.character(taxa.table[,1]))
taxa.table2 = taxa.table[match(unique.tax,taxa.table[,1]),] #remove duplicate 'environmental samples'
write.table(taxa.table2,'Searchable_taxa_COI.txt',sep='\t',quote=F,row.names=F)


########################################################################################################
# Get metadata (split workload among many CPU)

taxa.table = read.table('Searchable_taxa_COI.txt',sep='\t',as.is=T,check.names=F,header=T)
#terminal.seq = seq(1,nrow(taxa.table),5)

library('rentrez')

set_entrez_key('<your key here>')
Sys.getenv('ENTREZ_KEY')

N=1
con.geo = file('Insecta_COI_geo_metadata.txt','w') #first file, use headers. All subsequent files (to eventually be merged with first), no headers.
writeLines(paste('UID','ORGN','Species','Lat','Lon','Date',sep='\t'),con.geo)
con.nogeo = file('Insecta_COI_nogeo.metadata.txt','w')
writeLines(paste('UID','ORGN','Species','Date',sep='\t'),con.nogeo)

N=2
con.geo = file('Insecta_COI_geo_metadata2.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata2.txt','w')

N=3
con.geo = file('Insecta_COI_geo_metadata3.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata3.txt','w')

N=4
con.geo = file('Insecta_COI_geo_metadata4.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata4.txt','w')

N=5
con.geo = file('Insecta_COI_geo_metadata5.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata5.txt','w')

N=6
con.geo = file('Insecta_COI_geo_metadata6.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata6.txt','w')

N=7
con.geo = file('Insecta_COI_geo_metadata7.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata7.txt','w')

N=8
con.geo = file('Insecta_COI_geo_metadata8.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata8.txt','w')

N=9
con.geo = file('Insecta_COI_geo_metadata9.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata9.txt','w')

N=10
con.geo = file('Insecta_COI_geo_metadata10.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata10.txt','w')

N=11
con.geo = file('Insecta_COI_geo_metadata11.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata11.txt','w')

N=12
con.geo = file('Insecta_COI_geo_metadata12.txt','w')
con.nogeo = file('Insecta_COI_nogeo.metadata12.txt','w')

# Spread process across multiple CPU
#for (i in 1:6){ #1
#for (i in 6:11){ #2
#for (i in 12:16){ #3
#for (i in 17:21){ #4
#for (i in 22:26){ #5
#for (i in 27:31){ #6
#for (i in 32:36){ #7
#for (i in 37:41){ #8
#for (i in 42:46){ #9
#for (i in 47:51){ #10
#for (i in 52:56){ #11
for (i in 57:61){ #12

con.geo = file('esamples_geo.txt','w')
con.nogeo = file('esamples_nogeo.txt','w')

	orgn = taxa.table$Taxon[i]
	print(orgn)
	if (taxa.table$N.accessions[i]==0){
		#ignore rows in taxa.table that have 0 records
	} else {
		orgn.search = entrez_search(db='nuccore',term=paste0(orgn,'[ORGN] AND (COI[GENE] OR COX1[GENE])'),use_history=F,retmax=99999)
		seq.ids = c(seq(1,length(orgn.search$ids),300),length(orgn.search$ids)+1) #add one to length because of code in for loop
		for (k in 1:(length(seq.ids)-1)){
			print(paste0('   ',k,' of ',length(seq.ids)-1))
			meta = entrez_summary(db='nuccore',id=orgn.search$ids[seq.ids[k]:(seq.ids[k+1]-1)])
			meta.subtype = extract_from_esummary(esummaries=meta,elements='subtype')
			meta.subname = extract_from_esummary(esummaries=meta,elements='subname')
			meta.organism = extract_from_esummary(esummaries=meta,elements='organism')
			llgrep = grep('lat_lon',meta.subtype)
			if (length(llgrep)>0){
				#at least some records have lat_lon info
				meta.subtype.geo = meta.subtype[llgrep] #use to get positions of lat_lon info
				meta.subname.geo = meta.subname[llgrep] #contains actual attribute data
				meta.organism.geo = meta.organism[llgrep] #contains species name
				for (l in 1:length(meta.subtype.geo)){
					split = strsplit(meta.subname.geo[l],'[|]')[[1]] #split fields of subname (for use in special cases)
					#Get lat/lon
					latlon.pos = which(strsplit(meta.subtype.geo[l],'[|]')[[1]]=='lat_lon')
					latlon = strsplit(meta.subname.geo[l],'[|]')[[1]][latlon.pos]
					latlon_split = strsplit(latlon,' ')[[1]]
					lat = as.numeric(latlon_split[1])
					lon = as.numeric(latlon_split[3])
					#Deal with special case where subtype and subname fields do not align
					if (is.na(lat) | is.na(lon)){
						#occurs if latlon is not actually lat and lon (subname fields do not match subtype fields)
						for (sp in 1:length(split)){
							split1 = strsplit(split[sp],' ')[[1]]
							if (length(split1) != 4){
								#ignore; lat_lon will have 4 strings upon splitting
							} else {
								if (split1[2]=='N' | split1[2]=='S'){
									latlon_split = split1
									lat = as.numeric(split1[1])
									lon = as.numeric(split1[3])
									if (latlon_split[2]=='S'){lat=lat*(-1)}else{temp=c()} #correct sign
									if (latlon_split[4]=='W'){lon=lon*(-1)}else{temp=c()}
								} else {
								}
							}
						}
					} else { #lat/lon looks good
						if (latlon_split[2]=='S'){lat=lat*(-1)}else{temp=c()} #correct sign
						if (latlon_split[4]=='W'){lon=lon*(-1)}else{temp=c()}
					}
					#Get date
					date.pos = which(strsplit(meta.subtype.geo[l],'[|]')[[1]]=='collection_date')
					if (length(date.pos)==0){
						date = NA #no date associated with accession
					} else {
						date = strsplit(meta.subname.geo[l],'[|]')[[1]][date.pos]
						date.split = strsplit(date,'-')[[1]]
						if (length(date.split) != 3){ #check that date is actually a date (should have 3 strings upon splitting)
							for (ds in 1:length(split)){ #if not, find it!
								split2 = strsplit(split[ds],'-')[[1]]
								if (length(split2) != 3){
									#ignore; date will have 3 strings upon splitting
								} else {
									if (is.numeric(split2[1]) & is.character(split2[2]) & is.numeric(split2[3])){
										date = split[ds] #found the date!
									} else {
									}
								}
							}
						} else {
						}
					}
					add.line1 = paste(names(meta.subtype.geo[l]),orgn,meta.organism.geo[l],lat,lon,date,sep='\t')
					writeLines(add.line1,con.geo)
				}	
				#Deal with special case
				if (length(llgrep)==length(meta)){
					#all records have lat_lon
				} else {
					#only some records lack lat_lon
					meta.subtype.nogeo = meta.subtype[-llgrep]
					meta.subname.nogeo = meta.subname[-llgrep]
					meta.organism.nogeo = meta.organism[-llgrep]
					for (j in 1:length(meta.subtype.nogeo)){
						split.subname.nogeo = strsplit(meta.subname.nogeo[j],'[|]')[[1]] #split fields of subname (for use in special cases)
						date.pos1 = which(strsplit(meta.subtype.nogeo[j],'[|]')[[1]]=='collection_date')
						if (length(date.pos1)==0){
							date1 = NA
						} else {
							date1 = strsplit(meta.subname.nogeo[j],'[|]')[[1]][date.pos1]
							date.split1 = strsplit(date1,'-')[[1]]
							if (length(date.split1) != 3){ #check that date1 is actually a date (should have 3 strings upon splitting)
								for (ds1 in 1:length(split.subname.nogeo)){ #if not, find it!
									split.subname.nogeo1 = strsplit(split.subname.nogeo[ds1],'-')[[1]]
									if (length(split.subname.nogeo1) != 3){
										#skip; this is not the date
									} else {
										if (is.numeric(split.subname.nogeo1[1]) & is.character(split.subname.nogeo1[2]) & is.numeric(split.subname.nogeo1[3])){
											date1 = split.subname.nogeo[ds1] #found the date!
										} else {
										}
									}
								}
							} else {
							}
						}
						add.line2 = paste(names(meta.subtype.nogeo[j]),orgn,meta.organism.nogeo[j],date1,sep='\t')
						writeLines(add.line2,con.nogeo)
					}
				}
			} else {
				#no lat_lon in any records
				for (m in 1:length(meta.subtype)){
					split.subname.nogeo2 = strsplit(meta.subname[m],'[|]')[[1]]
					date.pos2 = which(strsplit(meta.subtype[m],'[|]')[[1]]=='collection_date')
					if (length(date.pos2)==0){
						date2 = NA
					} else {
						date2 = strsplit(meta.subname[m],'[|]')[[1]][date.pos2]
						date.split2 = strsplit(date2,'-')[[1]]
						if (length(date.split2) != 3){ #check that date is actually a date (should have 3 strings upon splitting)
							for (ds2 in 1:length(split.subname.nogeo2)){ #if not, find it!
								split.subname.nogeo2.1 = strsplit(split.subname.nogeo2[ds2],'-')[[1]]
								if (length(split.subname.nogeo2.1) != 3){
									#skip; this is not the date
								} else {
									if (is.numeric(split.subname.nogeo2.1[1]) & is.character(split.subname.nogeo2.1[2]) & is.numeric(split.subname.nogeo2.1[3])){
										date2 = split.subname.nogeo2[ds2] #found the date!
									} else {
									}
								}
							}
						} else {
						}
					}
					add.line3 = paste(names(meta.subtype[m]),orgn,meta.organism[m],date2,sep='\t')
					writeLines(add.line3,con.nogeo)
				}
			}
		}
	}
}
close(con.geo)
close(con.nogeo)


# Merge all of the metadata files into one
con.geo = file('Insecta_COI_geo_metadata.txt','r')
geo = readLines(con.geo)
close(con.geo)
con.nogeo = file('Insecta_COI_nogeo.metadata.txt','r')
nogeo = readLines(con.nogeo)
close(con.nogeo)

geo.con = file('Insecta_COI_geo_metadata_ALL.txt','w')
writeLines(geo,geo.con)
nogeo.con = file('Insecta_COI_nogeo_metadata_ALL.txt','w')
writeLines(nogeo,nogeo.con)

for (f in 2:12){
	con1 = file(paste0('Insecta_COI_geo_metadata',f,'.txt'),'r')
	add.geo = readLines(con1)
	close(con1)
	con2 = file(paste0('Insecta_COI_nogeo.metadata',f,'.txt'),'r')
	add.nogeo = readLines(con2)
	close(con2)
	for (i in 1:length(add.geo)){
		writeLines(add.geo[i],geo.con)
	}
	for (j in 1:length(add.nogeo)){
		writeLines(add.nogeo[j],nogeo.con)
	}
}
close(geo.con)
close(nogeo.con)


# Check attributes/coordinates of samples appearing in ocean
# "environmental samples" includes mostly non-insects!
metadata = read.table('Insecta_COI_geo_metadata_ALL.txt',sep='\t',as.is=T,check.names=F,header=T)
check = metadata[which(metadata$Lat<41 & metadata$Lat>40 & metadata$Lon>(-68) & metadata$Lon<(-66)),]
str(check) #these are fish samples!
esamples = read.table('esamples_geo.txt',sep='\t',as.is=T,check.names=F)
ematch = match(esamples[,1],metadata$UID)
metadata2 = metadata[-ematch,]
write.table(metadata2,'Insecta_COI_geo_metadata_ALL2.txt',sep='\t',quote=F,col.names=T,row.names=F)


########
# Write shapefile

library('rgdal')
metadata = read.table('Insecta_COI_geo_metadata_ALL2.txt',sep='\t',as.is=T,check.names=F,header=T)
ll.proj = '+init=epsg:4326' #WGS84
shape = SpatialPointsDataFrame(coords=data.frame('lon'=metadata$Lon,'lat'=metadata$Lat),data=metadata[,c(1:3,6)],proj4string=CRS(ll.proj))
writeOGR(shape,dsn='./Insecta_COI.shp',layer='Insecta_COI',driver='ESRI Shapefile',overwrite=T)


###############################################################################
# Create fasta with sequences (order matches ~metadata_ALL)

library('rentrez')
library('rgdal')
library('rgeos')

set_entrez_key('<your key here>')
Sys.getenv('ENTREZ_KEY')

ll.proj = '+init=epsg:4326' #WGS84

coi = spTransform(readOGR(dsn='Insecta_COI.shp',layer='Insecta_COI',stringsAsFactors=F),CRS(ll.proj))
usa = spTransform(readOGR(dsn='/path/to/National_atlas/statesp010g.shp',layer='statesp010g'),CRS(ll.proj))
usa = usa[-which(usa@data[,1]=='U.S. Virgin Islands' | usa@data[,1]=='Hawaii' | usa@data[,1]=='Alaska' | usa@data[,1]=='Puerto Rico'),]

# Subset COI accessions to those which intersect the conterminous US
usa.intersect = gIntersects(coi,usa,byid=T) #rows=states,columns=accessions
intersect.count = apply(usa.intersect,2,function(x){length(which(x==TRUE))})
coi.usa = coi[which(intersect.count>0),]
writeOGR(coi.usa,dsn='Insecta_COI_USAcont.shp',layer='Insecta_COI_USAcont',driver='ESRI Shapefile',overwrite=T)

# Get fastas
coi.usa = readOGR(dsn='Insecta_COI_USAcont.shp',layer='Insecta_COI_USAcont',stringsAsFactors=F)

uids = coi.usa@data$UID
con1 = file('Insecta_COI_USA.fasta','w')
for (k in 1:length(uids)){
	print(paste0('   ',k,' of ',length(uids)))
	fastas = entrez_fetch(db='nuccore',id=uids[k],rettype='fasta')
	fastas1 = sub('>',paste0('>',uids[k],'; '),fastas)
	writeLines(fastas1,con1)
	Sys.sleep(1)
}
close(con1)


#####################################
# Summarize representation of taxa in US

library('rgdal')

coi.usa = readOGR(dsn='./Insecta_COI_USAcont.shp',layer='Insecta_COI_USAcont',stringsAsFactors=F)
coi.usa.data = coi.usa@data
orgn.key = read.table('ORGNkey.txt',sep='\t',as.is=T,check.names=F,header=T)
orders = unique(orgn.key[,1])
order.count = c()
for (i in 1:length(orders)){
	print(orders[i])
	o1 = orders[i]
	orgns = orgn.key[which(orgn.key[,1]==o1),2]
	o.count = 0
	for (j in 1:length(orgns)){
		orgn1 = orgns[j]
		orgn.count = length(which(coi.usa.data$ORGN==orgn1))
		o.count = o.count + orgn.count
	}
	order.count = c(order.count,o.count)
}
cbind(orders,order.count)

# Add a column to the coi.usa shapefile that contains insect order of each accession
# Add column with year
order.col = c()
for (i in 1:nrow(coi.usa.data)){
	orgn1 = coi.usa.data$ORGN[i]
	order.match = orgn.key[which(orgn.key[,2]==orgn1),1]
	order.col = c(order.col,order.match)
}
years = apply(array(coi.usa.data$Date),1,function(x){strsplit(x,'[-]')[[1]][3]})

coi.usa.data2 = data.frame(coi.usa.data,'Year'=as.numeric(years),'Order'=order.col,stringsAsFactors=F)

# Append to shapefile
coi.usa@data = coi.usa.data2
writeOGR(coi.usa,dsn='Insecta_COI_USAcont2.shp',layer='Insecta_COI_USAcont2',driver='ESRI Shapefile',overwrite=T)


# Quick visuals
ll.proj = '+init=epsg:4326' #WGS84
coi.usa = readOGR(dsn='Insecta_COI_USAcont2.shp',layer='Insecta_COI_USAcont2',stringsAsFactors=F)
usa = spTransform(readOGR(dsn='/path/to/National_atlas/statesp010g.shp',layer='statesp010g'),CRS(ll.proj))
usa = usa[-which(usa@data[,1]=='U.S. Virgin Islands' | usa@data[,1]=='Hawaii' | usa@data[,1]=='Alaska' | usa@data[,1]=='Puerto Rico'),]

orders = unique(coi.usa@data$Order)
for (i in 1:length(orders)){
	print(orders[i])
	order.pos = which(coi.usa@data$Order==orders[i])
	coi.subset = coi.usa[order.pos,]
	u.coords = unique(paste(coi.subset@coords[,1],coi.subset@coords[,2]))
	png(paste0('./COI_',orders[i],'_map.png'))
	plot(usa,main=paste0(orders[i],'\n',length(order.pos),' accessions, ',length(u.coords),' sites'))
	plot(coi.subset,add=T,pch=16)
	dev.off()
}

streamdwellers = which(coi.usa@data$Order=='Ephemeroptera' | coi.usa@data$Order=='Plecoptera' | coi.usa@data$Order=='Trichoptera' | coi.usa@data$Order=='Odonata')
coi.subset = coi.usa[streamdwellers,]
u.coords = unique(paste(coi.subset@coords[,1],coi.subset@coords[,2]))
png('./COI_streamdwellers_map.png')
plot(usa,main=paste0('Stream-dwellers\n',length(streamdwellers),' accessions, ',length(u.coords),' sites'))
plot(coi.subset,add=T,pch=16)
dev.off()


# Create taxonomy key for metadata
metadata = read.table('Insecta_COI_geo_metadata_ALL2.txt',sep='\t',as.is=T,check.names=F,header=T)
species.u = unique(metadata$Species)
species.split = apply(array(species.u),1,function(x){strsplit(x,' sp.')[[1]][1]})
species.u2 = unique(species.split)
species = c()
genera = c()
families = c()
for (i in 1:length(species.u2)){
	print(paste0(i,' of ',length(species.u2)))
	sp.split1 = strsplit(species.u2[i],' ')[[1]]
	if (length(grep('sp.',sp.split1))==1){ #get genus & family from ambiguous species records
		gen.fam = tax_name(query=sp.split1[1], get=c('species','genus','family'), db='ncbi')
		spec1 = gen.fam$species
		gen1 = gen.fam$genus
		fam1 = gen.fam$family
	} else { #normal case where accession already contains genus & species
		gen.fam2 = tax_name(query=species.u2[i], get=c('species','genus','family'), db='ncbi')
		spec1 = gen.fam2$species
		gen1 = gen.fam2$genus
		fam1 = gen.fam2$family
	}
	families = c(families,fam1)
	genera = c(genera,gen1)
	species = c(species,spec1)
	Sys.sleep(1) #pause for # sec
}
add.info = data.frame('Species'=species,'Genus'=genera,'Family'=families)
write.table(add.info,'add.info.temp.txt',sep='\t',row.names=F,col.names=T)

add.info[which(is.na(add.info$Family)),]
add.info[which(is.na(add.info$Genus)),]

# Match metadata entries with taxonomy info
add.info = read.table('add.info.temp.txt',sep='\t',as.is=T,check.names=F,header=T)
add.tax.u = data.frame('Species'=NA,'Genus'=NA,'Family'=NA)
problems = c()
for (i in 1:length(species.u2)){
	fam.match = which(add.info$Family==species.u2[i])
	if (length(fam.match)>0){ #entry is a family name
		add.tax.u[i,] = c(NA,NA,species.u2)
	} else {
		genus.match = which(add.info$Genus==species.u2[i])
		if (length(genus.match)>0){ #entry is a genus name
			add.info.g = add.info[genus.match,]
			fam1 = unique(add.info.g$Family)
			add.tax.u[i,] = c(NA,species.u2,fam1)
		} else {
			species.match = which(add.info$Species==species.u2[i])
			if (length(species.match)>0){ #entry is a full (genus &) species name
				add.tax.u[i,] = add.info[species.match,]
			} else { #entry has no match in taxonomy key
				problems = c(problems,species.u2[i])
			}
		}
	}
}

# Attempt to recover taxonomy for problem species (recovered several hundred)
species.u2 = problems
species = c()
genera = c()
families = c()
for (i in 35:length(species.u2)){
	print(paste0(i,' of ',length(species.u2)))
	gen.fam = tax_name(query=species.u2[i], get=c('species','genus','family'), db='ncbi')
	families = c(families,gen.fam$family)
	genera = c(genera,gen.fam$genus)
	species = c(species,gen.fam$species)
	Sys.sleep(1) #pause for # sec
}
add.info2 = data.frame('Species'=species,'Genus'=genera,'Family'=families,stringsAsFactors=F)[-35,]

species.u2[which(is.na(add.info2$Family))]


# Append taxonomy info to metadata
add.infos = rbind(add.info,add.info2) #combine tables
add.infos.merge = apply(add.infos,1,function(x){paste(x[1],x[2],x[3],sep='_')})
add.unique = unique(add.infos.merge)
keep1 = match(add.unique,add.infos.merge)
add.infos2 = add.infos[keep1,] #prune duplicates
write.table(add.info,'add.info.temp2.txt',sep='\t',row.names=F,col.names=T)

coi.usa = readOGR(dsn='Insecta_COI_USAcont2.shp',layer='Insecta_COI_USAcont2',stringsAsFactors=F)
metadata.usa = coi.usa@data

add.tax = data.frame('Species'=NA,'Genus'=NA,'Family'=NA,'Order'=NA)
problems = c()
for (i in 1:nrow(metadata.usa)){
	print(paste0(i,' of ',nrow(metadata.usa)))
	meta.species.split = strsplit(metadata.usa$Species[i],' sp. ')[[1]]
	meta.orgn = metadata.usa$ORGN[i]
	meta.order = orgn.key[which(orgn.key[,2]==meta.orgn),1]
	fam.match = which(add.infos2$Family==meta.species.split)
	if (length(fam.match)>0){ #entry is a family name
		add.tax[i,] = c(NA,NA,meta.species.split[1],meta.order)
	} else {
		genus.match = which(add.infos2$Genus==meta.species.split)
		if (length(genus.match)>0){ #entry is a genus name
			add.infos2.g = add.infos2[genus.match,]
			fam1 = unique(add.infos2.g$Family)
			add.tax[i,] = c(NA,meta.species.split[1],fam1,meta.order)
		} else {
			species.match = which(add.infos2$Species==meta.species.split)
			if (length(species.match)>0){ #entry is a full (genus &) species name
				add.tax[i,] = c(add.infos2[species.match,],meta.order)
			} else { #entry has no match in taxonomy key
				add.tax[i,] = rep(NA,4)
				problems = c(problems,i)
			}
		}
	}
}
coi.usa@data = data.frame(metadata.usa,add.tax[,-4],stringsAsFactors=F)
colnames(coi.usa@data)[3] = 'Name'
colnames(coi.usa@data)[7] = 'Species'
writeOGR(coi.usa,dsn='Insecta_COI_USAcont3.shp',layer='Insecta_COI_USAcont3',driver='ESRI Shapefile',overwrite=T)

length(which(is.na(coi.usa@data$Family))) / nrow(coi.usa@data)
length(which(is.na(coi.usa@data$Genus))) / nrow(coi.usa@data)
length(which(is.na(coi.usa@data$Species))) / nrow(coi.usa@data)


####################################################################
# Create MAFFT alignments for each insect order

library('rgdal')
library('rlang')

coi.shp = readOGR(dsn='./Insecta_COI_USAcont3.shp',layer='Insecta_COI_USAcont3',stringsAsFactors=F)
coi.data = coi.shp@data
orders = unique(coi.data$Order)
order.count = apply(array(orders),1,function(x){length(which(coi.data$Order==x))})
orders.counts = cbind(orders,order.count)

con = file('./Insecta_COI_USA_single.fasta','r') #import sequences
fasta = readLines(con)
close(con)
headers = fasta[seq(1,length(fasta),2)]
sequences = fasta[seq(2,length(fasta),2)]
uids = apply(array(headers),1,function(x){sub('>','',strsplit(x,';')[[1]][1])})

# Run MAFFT
for (o in 1:length(orders)){
	print(orders[o])
	system(paste0('mkdir ./mafft_tmp/',orders[o])) #create subdirectory for order o
	# Create order-specific fasta
	o.uids = coi.data$UID[which(coi.data$Order==orders[o])]
	o.pos = match(o.uids,uids)
	print('...creating fasta...')
	outfile = file(paste0('./mafft_tmp/',orders[o],'/COI_USA_',orders[o],'.fasta'),'w') #start new fasta file in order's directory
	for (l in 1:length(o.pos)){
		writeLines(headers[o.pos[l]],outfile)
		writeLines(sequences[o.pos[l]],outfile)
	}
	close(outfile)
	# Launch MAFFT
	print('...running MAFFT...')
	mafft.command = paste0('mafft --auto --thread 12 --threadtb 12 --threadit 12 ./mafft_tmp/',
		orders[o],'/COI_USA_',orders[o],'.fasta > ./mafft_tmp/',
		orders[o],'/COI_USA_',orders[o],'_aligned.fasta')
	system(mafft.command)
}


####################################################################
# In command terminal, linearize fastas
files="Blattodea
Coleoptera
Collembola
Dermaptera
Diptera
Ephemeroptera
Hemiptera
Heteroptera
Hymenoptera
Isoptera
Lepidoptera
Mantodea
Mecoptera
Microcoryphia
Neuroptera
Odonata
Orthoptera
Phasmatodea
Phthiraptera
Plecoptera
Psocoptera
Siphonaptera
Thysanoptera
Trichoptera"
i=1
for file in $files
do
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ./${file}/COI_USA_${file}_aligned.fasta > ./${file}/COI_USA_${file}_aligned_single.fasta
let "i+=1";
done


####################################################################
# Group sequences by spatial proximity; species-specific demes

library('rgdal')
library('geosphere')

coi.shp = readOGR(dsn='./Insecta_COI_USAcont3.shp',layer='Insecta_COI_USAcont3',stringsAsFactors=F)
coi.data = coi.shp@data
orders = unique(coi.shp@data$Order); orders = orders[-which(orders=='Mantodea')] #Mantodea has only 1 accession in USA
dists = c(1,10,100) #distance thresholds in km

log1 = data.frame('Order'=NA,'Species'=NA,'Scale'=-999,'N.demes'=-999,'N.seqs'=-999,'N.singletons'=-999,stringsAsFactors=F)
deme.div = data.frame('Order'=NA,'Species'=NA,'Scale'=-999,'Deme'=NA,'Nuc.Div'=-999,'U.Haps'=-999,'N.pairs'=-999,'N.seqs'=-999,'Median.date'=-999,'N.dates'=-999,'Median.dist'=-999)
deme.key = c()
lx = 1 #row index for log
dx = 1 #row index for deme.div
for (o in 1:length(orders)){
	print(orders[o])
	#import sequences for order o
	con = file(paste0('./mafft_tmp/',orders[o],'/COI_USA_',orders[o],'_aligned_single.fasta'),'r')
	o.fasta = readLines(con)
	close(con)
	o.fasta.uids = apply(array(o.fasta[seq(1,length(o.fasta),2)]),1,function(x){sub('>','',strsplit(x,';')[[1]][1])})
	o.fasta.seqs = o.fasta[seq(2,length(o.fasta),2)]
	seq.len = nchar(o.fasta.seqs)[1]
	#subset by species
	o.species = unique(coi.shp@data$Species[which(coi.shp@data$Order==orders[o])])
	n.species = length(o.species)
	for (s in 1:length(o.species)){
		print(paste0('...species ',s,' of ',n.species))
		if (is.na(o.species[s])){ #skip calculations for sequences that lack a species designation
			o.uids = coi.shp@data$UID[which(coi.shp@data$Order==orders[o] & is.na(coi.shp@data$Species))]
			log1[lx,1] = orders[o] #insect order
			log1[lx,2] = o.species[s] #species
			log1[lx,3] = NA #spatial scale
			log1[lx,4] = NA #no. demes (with >1 sequence)
			log1[lx,5] = NA #no. sequences kept
			log1[lx,6] = length(o.uids) #no. singletons
			lx = lx + 1
		} else if (length(which(coi.shp@data$Order==orders[o] & coi.shp@data$Species==o.species[s]))==1){ #skip species that only have 1 sequence represented
			log1[lx,1] = orders[o] #insect order
			log1[lx,2] = o.species[s] #species
			log1[lx,3] = NA #spatial scale
			log1[lx,4] = NA #no. demes (with >1 sequence)
			log1[lx,5] = NA #no. sequences kept
			log1[lx,6] = 1 #no. singletons
			lx = lx + 1
		} else {
			o.uids = coi.shp@data$UID[which(coi.shp@data$Order==orders[o] & coi.shp@data$Species==o.species[s])]
			o.latlon = coi.shp@coords[which(coi.shp@data$Order==orders[o] & coi.shp@data$Species==o.species[s]),]
			#get spatial data & calculate pairwise distances
			dist.mat = matrix(-999,nrow=nrow(o.latlon),ncol=nrow(o.latlon))
			print('......calculating geographic distance matrix')
			for (i in 1:nrow(o.latlon)){ #fill distance matrix with each pairwise distance
				for (j in 1:nrow(o.latlon)){
					dist.mat[i,j] = distGeo(p1=c(o.latlon[i,1],o.latlon[i,2]),p2=c(o.latlon[j,1],o.latlon[j,2])) / 1000 #distance between sample i and j in km
				}
			}
			clust = hclust(as.dist(dist.mat)) #hierarchical clustering based on distance (km)
			#get dates
			o.years = coi.data$Year[which(coi.data$Order==orders[o])]
			# Calculate genetic diversity at each distance threshold
			print('......calculating per-deme nucleotide diversity')
			for (d in 1:length(dists)){
				print(paste0('......',dists[d],'km'))
				demes = cbind(o.latlon,cutree(clust,h=dists[d])) #matrix with lon, lat, and spatial cluster
				deme.key = data.frame(rbind(deme.key,data.frame('UID'=o.uids,'Order'=rep(orders[o],length(o.uids)),'Species'=rep(o.species[s],length(o.uids)),'Scale'=rep(dists[d],length(o.uids)),'Lon'=demes[,1],'Lat'=demes[,2],'Deme'=demes[,3]))) #keep track of lat/lon of samples within each deme
				deme.sizes = tabulate(demes[,3])
				d.keep = which(deme.sizes>1)
				d.toss = which(deme.sizes==1)
				# Add details to log
				log1[lx,1] = orders[o] #insect order
				log1[lx,2] = o.species[s] #species
				log1[lx,3] = dists[d] #spatial scale
				log1[lx,4] = length(d.keep) #no. demes (with >1 sequence)
				log1[lx,5] = sum(deme.sizes[d.keep]) #no. sequences kept
				log1[lx,6] = sum(deme.sizes[d.toss]) #no. singletons
				lx = lx + 1
				# Estimate genetic diversity of each deme
				for (i2 in d.keep){ #for each appropriately large deme
					o.d.years =  o.years[which(demes[,3]==i2)]
					o.d.uids = o.uids[which(demes[,3]==i2)]
					o.d.taxonomy = coi.data[match(o.d.uids,coi.data$UID),7:9]
					o.d.seqs = o.fasta.seqs[match(o.d.uids,o.fasta.uids)]
					seq.matrix = t(apply(array(o.d.seqs),1,function(x){strsplit(x,'')[[1]]}))
					pairs1 = combn(length(o.d.uids),2)
					n.pairs = ncol(pairs1) #number of pairwise comparisons
					lost.pairs = 0
					n.seqs = length(o.d.uids)
					u.seq = length(unique(o.d.seqs))
					pi.sum = 0 #sum of kij/mij for all sequence pairs
					for (j2 in 1:n.pairs){ #calculate per site nucleotide diversity for each sequence pair (eq. 2 from Miraldo et al. 2016)
						mij = 0 #number of shared base pairs between sequence i and j
						kij = 0 #number of nucleotide differences between sequence i and j
						for (l in 1:seq.len){
							n1 = seq.matrix[pairs1[1,j2],l]
							n2 = seq.matrix[pairs1[2,j2],l]
							if (n1=='-' | n2=='-'){
								#ignore base pair (not shared between sequence i and j)
							} else {
								mij = mij + 1 #add to shared base pair counter
								if (n1==n2){
									#base pair is identical between sequence i and j
								} else {
									kij = kij + 1 #add to nucleotide difference counter
								}
							}
						} #end per base pair loop
						if ((mij/seq.len)<0.5){
							lost.pairs = lost.pairs + 1 #pair has <50% shared base pairs
						} else {
						}
						pi.sum = pi.sum + (kij / mij)
					} #end per pair loop
					#get median distance among sequences
					mat.pos = which(demes[,3]==i2)
					o.d.dist.mat = dist.mat[mat.pos,mat.pos]
					o.d.dist.mat.lower = o.d.dist.mat[lower.tri(o.d.dist.mat)]
					#add data to output table
					GD = pi.sum / n.pairs #divide by number of pairwise comparisons
					deme.div[dx,1] = orders[o] #order
					deme.div[dx,2] = o.species[s] #species
					deme.div[dx,3] = dists[d] #scale
					deme.div[dx,4] = i2 #deme
					deme.div[dx,5] = GD #nucleotide diversity
					deme.div[dx,6] = u.seq / n.seqs #fraction of haplotypes that are unique
					deme.div[dx,7] = n.pairs - lost.pairs #number of pairs
					deme.div[dx,8] = n.seqs #number of sequences
					deme.div[dx,9] = median(o.d.years,na.rm=T) #median date of sequences
					deme.div[dx,10] = length(unique(o.d.years)) #number of dates represented by sequences
					deme.div[dx,11] = median(o.d.dist.mat.lower,na.rm=T) #median distances among sequences
					dx = dx + 1
				} #end per deme loop
			} #end per scale loop	
		} #end if species is not NA
	} #end species loop
} #end per order loop
write.table(log1,'demes_attributes_SpeciesLevel.txt',sep='\t',row.names=F,quote=F)
write.table(deme.div,'demes_diversity_SpeciesLevel.txt',sep='\t',row.names=F,quote=F)
write.table(deme.key,'demes_key_SpeciesLevel.txt',sep='\t',row.names=F,quote=F)


# Create shapefiles for each order with genetic diversity of demes

library('rgdal')
library('rgeos')

coi.shp = readOGR(dsn='./Insecta_COI_USAcont3.shp',layer='Insecta_COI_USAcont3',stringsAsFactors=F)
orders = unique(coi.shp@data$Order); orders = orders[-which(orders=='Mantodea')] #Mantodea has only 1 accession in USA
dists = c(1,10,100) #distance thresholds in km

deme.div = read.table('demes_diversity_SpeciesLevel.txt',sep='\t',as.is=T,check.names=F,header=T)
deme.key = read.table('demes_key_SpeciesLevel.txt',sep='\t',as.is=T,check.names=F,header=T)

# Get centroids for each deme
deme.centroids = c()
for (i in 1:nrow(deme.div)){
	deme1 = deme.key[which(deme.key$Order==deme.div$Order[i] & deme.key$Deme==deme.div$Deme[i] & deme.key$Scale==deme.div$Scale[i] & deme.key$Species==deme.div$Species[i]),]
	deme.coords = deme1[,5:6] #WARNING: sensitive to removal/addition of variables
	d.centroid = gCentroid(SpatialPoints(coords=deme.coords,proj4string=CRS('+init=epsg:4326')))
	deme.centroids = rbind(deme.centroids,d.centroid@coords)	
}
deme.div2 = data.frame(deme.div,'Lon'=deme.centroids[,1],'Lat'=deme.centroids[,2],stringsAsFactors=F)
div.shp = SpatialPointsDataFrame(coords=deme.div2[,12:13],data=deme.div2[,1:11],proj4string=CRS('+init=epsg:4326')) #WARNING: sensitive to removal/addition of variables
writeOGR(div.shp,dsn='./Insecta_COI_USA_deme_diversity_SpeciesSpecific.shp',layer='Insecta_COI_USA_deme_diversity_SpeciesSpecific',driver='ESRI Shapefile',overwrite=T)

# What years were represented in EPT sequences?
orders = c('Ephemeroptera','Plecoptera','Trichoptera')
for (o in 1:length(orders)){
	print(orders[o])
	#import sequences for order o
	con = file(paste0('./mafft_tmp/',orders[o],'/COI_USA_',orders[o],'_aligned_single.fasta'),'r')
	o.fasta = readLines(con)
	close(con)
}

