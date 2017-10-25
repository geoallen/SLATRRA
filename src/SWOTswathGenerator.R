# SWOTswathGenerator_forRapidHub.R
# George H. Allen, October 2017

############################
# Description:
# Generates the SWOT dual swath paths from Aviso ephemeris table used in the Allen et al.
# "Global estimates of river flow wave travel times and implications for low-latency satellite data".
# Data output is a polygon shapefile with attributes corresponding to track orbit, swath side, direction,
# satellite location, and time along the full orbit cycle. See see the associated dataset descriptions section 
# for details of output. 


############################
# load libraries:
if (!"geosphere" %in% rownames(installed.packages())) {install.packages("geosphere")}; require(geosphere)
if (!"shapefiles" %in% rownames(installed.packages())) {install.packages("shapefiles")}; require(shapefiles)


############################
# path to ephemeris table (from https://www.aviso.altimetry.fr/en/missions/future-missions/swot/orbit.html):
ephPath = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/SWOTtracks/orig_ephemis_data/ephem_science_sept2015_ell.txt"

# polygon shapefile outpath:
#OutSFilePath = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/SWOTtracks_official/calval_tracks/calval_tracks"
OutSFilePath = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/SWOTtracksl/SWOTtracks_polygon/ephem_science_sept2015_track"


############################
# functions:

# insert a row into a table:
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1), ] = existingDF[seq(r,nrow(existingDF)), ]
  existingDF[r,] = newrow
  return(existingDF)
}

# Add an ID field to indicate orbit number:
orbitIDer <- function(DF){
  NArows = which(abs(diff(DF$lon)) > 100)
  print(paste(length(NArows), "orbits found"))
  NAbreaks = diff(c(0, NArows, nrow(DF))); 
  orbitID = rep(1:(length(NArows)+1), NAbreaks, each=T)
  DF = cbind(DF, orbitID)
  return(DF)
}

# Add ascending/descending rag and an corresponding ID field: 
dirIDer <- function(DF){
  dif = diff(DF$lat)
  ascInd = which(dif>0)
  dscInd = which(dif<=0)
  Track_type = rep("ascending", nrow(DF))
  Track_type[dscInd] = "descending"
  
  jumpInd = sort(c(1, ascInd[which(diff(ascInd)>1)], dscInd[which(diff(dscInd)>1)], nrow(DF)))
  breaks = diff(jumpInd)
  breaks[1] = breaks[1]+1
  dirID = rep(1:(length(jumpInd)-1), breaks, each=T)
  print(paste(length(jumpInd), "ascending/descending cycles"))
  DF = cbind(DF, Track_type, dirID)
  return(DF)
}

# create a combined ID for orbit and ascending/descending:
combineIDs <- function(DF){
  ID = paste0(DF$orbitID, "_", DF$dirID, "_", DF$swathSide)
  DF = cbind(DF, ID)
  return(DF)
}

# duplicate the first element of each asc/desc ID to ensure polygons touch each other:
addExtraCoordinate <- function(DF){
  firstUIDind = rev(match(unique(DF$dirID), DF$dirID))
  for (i in 1:length(firstUIDind)){
    DF = insertRow(DF, DF[firstUIDind[i],], firstUIDind[i])
  }
  # change indices of these inserted rows so that they match the previous index 
  # and remove the first index, since it doesn't have a previous index:
  firstUIDind = match(unique(DF$dirID), DF$dirID)[-1]
  # if the dirID matches up with the edge of an orbit, remove it from list:
  cornerFlag = which(DF$orbit[firstUIDind] != DF$orbit[firstUIDind-1])
  if(length(cornerFlag)>0){
    DF = DF[-firstUIDind[cornerFlag],]
    firstUIDind = firstUIDind[-cornerFlag]
  }
  DF[firstUIDind, 3:ncol(DF)] = DF[firstUIDind-1, 3:ncol(DF)]
  DF = DF[-1, ]
  return(DF)
}

# fill in shapefile attribute table 
ddTableFiller <- function(ddTable, swath, i, ind, mInd){
  ddTable$mean_lon_dd[i] = swath$lon[mInd]
  ddTable$mean_lat_dd[i] = swath$lat[mInd]
  ddTable$orbit_ID[i] = swath$orbitID[mInd]
  ddTable$swathSide[i] = swath$swathSide[mInd]
  ddTable$tracktype[i] = as.vector(swath$Track_type[mInd])
  ddTable$mean_time_s[i] = round(mean(swath$time[ind]))
  ddTable$mean_elev_m[i] = round(mean(swath$elev[ind]))
  return(ddTable)
}

############################
# main code:

# read in ephemeris table:
eph = read.table(ephPath, sep = ' ')
names(eph) = c("time", "lon", "lat", "elev")

# write out as a csv for Arc:
FID = 1:nrow(eph)
write.csv(cbind(FID,eph), sub('txt', 'csv', ephPath), row.names=F)

# convert lon from range of (0:360) to(-180:180):
eph$lon[eph$lon>180] = eph$lon[eph$lon>180] - 360

# insert ID field for ascention/descention:
eph = dirIDer(eph)

# insert ID field for each orbit:
eph = orbitIDer(eph)

# duplicate the first element of each unique ID to ensure polygons touch
# at tops and bottoms of peaks:
#eph = addExtraCoordinate(eph)

# concat. the first and last index to the last and first index:
nrows = nrow(eph)
lon = c(eph$lon[nrows], eph$lon, eph$lon[1])
lat = c(eph$lat[nrows], eph$lat, eph$lat[1])
len = length(lon)
# remove first two and last two indices respectively to create pair of coordinates for bearing calc:
p1 = cbind(lon[-c(len-1, len)], lat[-c(len-1, len)])
p2 = cbind(lon[-c(1,2)], lat[-c(1,2)])

# calculate bearing between each adjacent nadir path:
b = bearing(p1, p2)
b[b<0] = b[b<0]+180 # take care of the first and last index

# calcuate inner and outer SWOT wide swath coordinates:
nadirXY = as.data.frame(cbind(eph$lon, eph$lat)); names(nadirXY) = c("lon", "lat")
leftOuterXY = as.data.frame(destPoint(p=nadirXY, b=b-90, d=60e3))
leftInnerXY = as.data.frame(destPoint(p=nadirXY, b=b-90, d=10e3))
rightInnerXY = as.data.frame(destPoint(p=nadirXY, b=b+90, d=10e3))
rightOuterXY = as.data.frame(destPoint(p=nadirXY, b=b+90, d=60e3))

# # plot:
# plot(nadirXY$lon, nadirXY$lat, type='l', lwd=0.1)
# lines(leftOuterXY$lon, leftOuterXY$lat, col="blue", lwd=0.1)
# lines(leftInnerXY$lon, leftInnerXY$lat, col="light blue", lwd=0.1)
# lines(rightInnerXY$lon, rightInnerXY$lat, col="red", lwd=0.1)
# lines(rightOuterXY$lon, rightOuterXY$lat, col="pink", lwd=0.1)

# insert a field for left or right swoth:
leftOuterXY$swathSide = 1
leftInnerXY$swathSide = 1
rightOuterXY$swathSide = 2
rightInnerXY$swathSide = 2

# add eph fields to dataframes:
names(eph)[names(eph)=="lon"] = "centerLon"
names(eph)[names(eph)=="lat"] = "centerLat"

leftOuterXY = cbind(leftOuterXY, eph)
leftInnerXY = cbind(leftInnerXY, eph)
rightInnerXY = cbind(rightInnerXY, eph)
rightOuterXY = cbind(rightOuterXY, eph)

# create a combined ID for orbit, ascending/descending, left/right:
leftOuterXY = combineIDs(leftOuterXY)
leftInnerXY = combineIDs(leftInnerXY)
rightOuterXY = combineIDs(rightOuterXY)
rightInnerXY = combineIDs(rightInnerXY)

# reverse sort inner tracks:
leftInnerXY = leftInnerXY[nrows:1, ]
rightInnerXY = rightInnerXY[nrows:1, ]

# concatenate inner and outer swath boundaries:
swath = rbind(leftOuterXY, leftInnerXY, rightOuterXY, rightInnerXY)

# sort by the ID:
swath = swath[order(swath$ID),]
swath = swath[order(swath$orbitID),]
swath = swath[order(swath$dirID),]


# modify swath side attribute:
swath$swathSide[swath$swathSide==1] = "left"
swath$swathSide[swath$swathSide==2] = "right"

# calculate attribute table:
ID = as.vector(unique(swath$ID))
mean_lon_dd = mean_lat_dd = orbit_ID = swathSide = tracktype = mean_time_s = mean_elev_m = rep(NA, length(name))
ddTable = data.frame(ID, mean_lon_dd, mean_lat_dd, orbit_ID, swathSide, tracktype, mean_time_s, mean_elev_m)

for (i in 1:nrow(ddTable)){
  ind = which(swath$ID == ddTable$ID[i])
  mInd = round(mean(ind))
  
  # remove vertices that jump orbits:
  farWest = which(swath$lon[ind] < -170)
  farEast = which(swath$lon[ind] > 170)
  NfarWest = length(farWest)
  NfarEast = length(farEast)
  if (NfarWest==0 | NfarEast==0){
    ddTable = ddTableFiller(ddTable, swath, i, ind, mInd)
  }else{
    if (NfarWest > NfarEast){
      swath = swath[-ind[farEast], ] 
    }else{
      swath = swath[-ind[farWest], ] 
    }
    ind = which(swath$ID == ddTable$ID[i])
    mInd = round(mean(ind))
    ddTable = ddTableFiller(ddTable, swath, i, ind, mInd)
  }
}

# build shapefile and write to file:
dd = data.frame(ID=as.vector(swath$ID), X=swath$lon, Y=swath$lat)
shapefile = convert.to.shapefile(dd, ddTable, "ID", 5)
write.shapefile(shapefile, OutSFilePath, arcgis=T)


system('say script done runny honey')



##########
## Below are issues that could be fixed if these orbits are used for other purposes than
# what is used in Allen et al. However, for the purposes of the Allen et al. study, these issues
# did not pose problems (#2 & #3) or they were resolved by using Arc (#1). 

# 1. Aviso ephemeris data is for a 21 day period. The exact repeat cycle 
# duration is 20.86455 days. Delete the following redundant tracks:
# IDs2rm = c("1_1_1", "1_1_2", 
#        "272_586_1", "272_586_2", 
#        "272_587_1", "272_587_2", 
#        "273_587_1", "273_587_2", 
#        "273_588_1", "273_588_2")

# 2. FIXME: tops and bottoms of orbit need to be touching 
# (i.e. using the same coordiniates - almost got this working using the subfunction above)

# 3. FIXME: edges (at 0 & 360 longitude are jagged. Extend orbits out past these boundries,
# then clip in Arc to make clean shapefile)

