# waveCelerityCalculator.R
# George H. Allen, October 2017

#### DESCRIPTION ####

# The following script contains the data analysis and production of figures and tables 
# in the Allen et al., "Global estimates of river flow wave travel times and 
# implications for low-latency satellite data"  
# 
# The most basic sections of the code include: a calculation of river slope
# Joining points of intererests (POIs: cities, dams, and gauges) to the river network 
# Modeling flow wave celerity and calculating travel time
# Validating the model with a gauge-based estimate of celerity
# Generating figures and tables, and alculating the various statistics presented in the manuscript
# For input data, see the associated dataset description section 

#### LIBRARIES ####
if (!"foreign" %in% rownames(installed.packages())) {install.packages("foreign")}; require(foreign)
if (!"MASS" %in% rownames(installed.packages())) {install.packages("MASS")}; require(MASS)
if (!"geosphere" %in% rownames(installed.packages())) {install.packages("geosphere")}; require(geosphere)
if (!"shapefiles" %in% rownames(installed.packages())) {install.packages("shapefiles")}; require(shapefiles)

#### HARD-CODED CONSTANTS ####
rivEndPtTabFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/andreadisEtAl2013/hydrosheds_wqdles'
rivMidPtTabFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/andreadisEtAl2013/hydrosheds_wqdLCbasFswot'
rivNetworkConnectFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/andreadisEtAl2013/hydrosheds_connectivity'
origPolylinesFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/andreadisEtAl2013/hydrosheds_wqd'
cityFpath = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/cities/NaturalEarth/orig/ne_10m_populated_places.dbf'
damFpath = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/dams/dams-rev01-global-shp/GRanD_dams_v1_1.dbf'
gaugeFpath = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/gauges/gages_x010g.shp_nt00884/gages_x010g.dbf'
slopeOutFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/slopeflowlines'
cityDamGaugeFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/cityDamGaugeflowlines'
celTabOutFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/celerityTables'
obsOutFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/obsflowlines'
gaugeRecordFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/inputDatasets/vector/gauges/USGS/2014' # downloaded in 2014 :(
rivOutFdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/celerityFlowlines'

# get list of tables:
rivEndPtTabFnames = list.files(rivEndPtTabFdir, 'dbf', recursive=T)
rivEndPtTabNames = substr(rivEndPtTabFnames, 1, 2)
rivEndPtTabFpaths = paste0(rivEndPtTabFdir, '/', rivEndPtTabFnames)

rivMidPtTabFnames = list.files(rivMidPtTabFdir, 'dbf', recursive=T)
rivMidPtTabNames = substr(rivMidPtTabFnames, 1, 2)
rivMidPtTabFpaths = paste0(rivMidPtTabFdir, '/', rivMidPtTabFnames)
midPtMatch = match(rivEndPtTabNames, rivMidPtTabNames)

rivNetworkConnectFnames = list.files(rivNetworkConnectFdir)
rivNetworkConnectFpaths = list.files(rivNetworkConnectFdir, 'csv', recursive=T, full.names=T)

origPolylinesFpaths = paste0(origPolylinesFdir, '/', rivEndPtTabFnames)

slopeOutFpaths = paste0(slopeOutFdir, '/', rivEndPtTabFnames)
cityDamGaugeFpaths = paste0(cityDamGaugeFdir, '/', rivEndPtTabFnames)
obsOutFpaths = paste0(obsOutFdir, '/', rivEndPtTabFnames)
rivOutFpaths = paste0(rivOutFdir, '/', rivEndPtTabFnames)
celTabOutFpaths = paste0(rivOutFdir, '/', rivEndPtTabNames, 'riv_rangeGT05.csv')
# list gauge files:
gFnames = list.files(gaugeRecordFdir)
gNames = as.numeric(sub('.csv', '', gFnames))
gFpaths = paste0(gaugeRecordFdir, '/', gFnames)

# specify hydroBasins which have river networks that drain north
# into SRTM data gap, north of 60N lat:
basin2rm = c(3000001840,3000004740,3000009130,2000028310,2000041840,2000043720,8000009560)

# cubic feet per second to cubic peters per second conversion:
cfs2mfsConv = 0.028316847

#### FUNCTIONS #### 

# for a given river segment, returns the upstream segment with the greatest
# drainage area. Takes in connectivity table and info about drainage area:
largerBasin = function(cTab, mTab, UID, ncolCtab){
  upUIDs = cTab[UID, 4:ncolCtab]
  upUIDs = upUIDs[upUIDs != 0]
  # in the case of two basins being of equal size, the first listed seg is used:
  upUID = upUIDs[which.max(mTab$AREA[upUIDs])]
  return(upUID)
}
tableMaker = function(colNames, nrows, fill=NA){
  tab = as.data.frame(array(fill, c(nrows, length(colNames))))
  names(tab) = colNames
  return(tab)
}
slopeCalc <- function(tab, cTab, zeroSlope=1e-4){
  # calculate the slope of each segment and add it to the table
  slope = (tab$ELEV_M[oddInd] - tab$ELEV_M[evenInd])/(1e3*tab$LENGTH_KM[evenInd])
  
  #########
  # for zero or negative slopes, increase slope smoothing window 
  # so that slopes become greater than zero:
  
  # copy tab and remove every other row:
  mTab = tab[oddInd, ]
  if (length(which(mTab$ARCID != cTab$UID))>0){
    message("UIDs between connectivity table and Andreadis dataset do not line up 1:1")
  }
  if (length(which(1:nrow(cTab) != cTab$UID)) > 0){
    message("UIDs do not match table indices 1:1")
  }

  # for rivers north of 60N, set their slope equal to 
  mTab$SLOPE = slope
  ncolCtab = ncol(cTab)
  
  # for zero segs that are both a head and mouth, set their slopes to NA or a small value:
  mTab$SLOPE[mTab$SLOPE <= 0 & cTab$DNSTR == 0 & cTab$N_UPSTRSEGS==0] = zeroSlope
  
  # find zero or negative slope segs:
  #zSlopeInd = which(mTab$SLOPE <= 0)
  zSlopeBool = mTab$SLOPE <= 0
  
  # for each negative or zero slope segments:
  #for (j in 1:length(zSlopeInd)){
  while(T %in% zSlopeBool){
     
    # select zero-slope segment:
    #upUID = dnUID = zSlopeInd[j]
    upUID = dnUID = which(zSlopeBool)[1]
    
    # create slope table to store input from the crawler functions and fill in the first row:
    sTabNames = c("len", "uprElev", "lwrElev", "UID")
    sTab = as.data.frame(array(NA, c(5, length(sTabNames)))); names(sTab) = sTabNames
    sTab[1,] = cbind(mTab$LENGTH_KM[upUID],tab$ELEV_M[upUID*2-1],tab$ELEV_M[upUID*2],upUID) 
    
    # start out by crawling both up and downstream to expand slope window.
    # if we reach upstream or downstream network limit, switch to downstream
    # or upstream only crawl direction:
    k = 1
    crwlMode = c("crawlUp","crawlDown")
    elevDiff = max(sTab$uprElev, na.rm=T) - min(sTab$lwrElev, na.rm=T)
    
    # while there is no differnce in elevation between the top and bottom of the 
    # slope calculation window, expand the slope calc window:
    while(elevDiff <= 0){
      
      # upstream crawler function:
      if (crwlMode[1] == "crawlUp"){
        upUID = largerBasin(cTab, mTab, upUID, ncolCtab)
        # if there are no segments upstream, 
        # trigger break to while statement:
        if (length(upUID) == 0){
          crwlMode[1] = 0
        }else{
          k = k + 1 
          sTab[k,] = cbind(mTab$LENGTH_KM[upUID],tab$ELEV_M[upUID*2-1],tab$ELEV_M[upUID*2],upUID)
        }
      }
      
      # downstream crawler function:
      if (crwlMode[2] == "crawlDown"){
        dnUID = cTab$DNSTR[dnUID]
        # if there are no segments downstream, 
        # trigger break to while statement:
        if (dnUID == 0){
          crwlMode[2] = 0
        }else{ 
          k = k + 1 
          sTab[k,] = cbind(mTab$LENGTH_KM[dnUID],tab$ELEV_M[dnUID*2-1],tab$ELEV_M[dnUID*2],dnUID)
        }
      }
      elevDiff = max(sTab$uprElev, na.rm=T) - min(sTab$lwrElev, na.rm=T)
      
      # rarely there there are multiple segments strung together without any positive
      # upstream elevation gain that go from the head to mouth of a river network. 
      # Set their slopes equal to a small value because their slopes cannot made positive 
      # through increasing the slope calculation window:
      if (identical(crwlMode, c("0","0"))){ elevDiff = 1 }
    }
    
    # assign mean slope to all segments within slope calculation window:
    meanSlope = elevDiff/(1e3*sum(sTab$len, na.rm=T))
    mTab$SLOPE[na.omit(sTab$UID)] = meanSlope
    
    zSlopeBool = mTab$SLOPE <= 0
  }
  
  return(mTab)
  
}
POI2RivJoin <- function(POItab, POIXY, pLineTab, attrTab, inAttr, POIname, scoreFunc, sRad, smallSrad){
  # for each point of interest (POI), find the river segment that the 
  # POI is most likely to be located on.
  
  # extract lat lon from tables:
  rXY = cbind(pLineTab$X, pLineTab$Y)
  tXY = cbind(attrTab$POINT_X, attrTab$POINT_Y)
  
  inAttrColInd = match(inAttr, names(POItab))
  # add new columns to attrTab:
  # add join score (degree of certainty whether the given POI lies on the 
  # selected river segment). Can add other attributes here too (e.g. JOINDIST):
  colNames = paste(POIname, c("JOINSCORE", inAttr), sep="_")
  outTab = data.frame(array(0, c(nrow(attrTab), length(colNames))))
  names(outTab) = colNames
  
  # if output columns already exist in attrTab, remove them from attrTab:
  colMatch = match(names(outTab), names(attrTab))
  if (T %in% !is.na(colMatch)){ attrTab = attrTab[, -colMatch] }
  
  # subset POI table based on spatial extent of polylines:
  rangeXY = cbind(range(rXY[,1]), range(rXY[,2]))
  POIsubInd = which(POIXY[,1] > rangeXY[1,1] & POIXY[,1] < rangeXY[2,1] & 
                      POIXY[,2] > rangeXY[1,2] & POIXY[,2] < rangeXY[2,2])
  
  # find best matching river segment for each point of interest:
  for (j in POIsubInd){ #nrow(POIXY)){
    # subset river polylines with an endpoints within a square distance from each POI:
    NESW = destPoint(POIXY[j,], c(0,90,180,270), sRad)
    closeRivs = which(tXY[,1] > NESW[4,1] & tXY[,1] < NESW[2,1] & tXY[,2] > NESW[3,2] & tXY[,2] < NESW[1,2])
    # if there are no nearby rivers, skip to next POI: 
    if (length(closeRivs)==0){next}
    # get distance of nearby polyline verticies:
    # match pLineTab ID with closeRiv Index:
    IDmatch = which(!is.na(match(pLineTab$Id, closeRivs)))
    closeXY = rXY[IDmatch, ]
    closeD = distGeo(POIXY[j,], closeXY)
    closestRivs = which(closeD < smallSrad)
    # if there are no nearby rivers, skip to next POI: 
    if (length(closestRivs)==0){next}
    closestRivsOrd = order(closeD[closestRivs])
    closestD = closeD[closestRivs[closestRivsOrd]]
    # get drainage area from attrTab by backing out nearest line dists UIDs:
    closestTabInd = pLineTab$Id[IDmatch[closestRivs[closestRivsOrd]]]
    closestA = attrTab$AREA[closestTabInd]
    # score river segments by their drainage area and their proximity to the POI center:
    score = scoreFunc(closestA, closestD)
    #print(cbind(closestTabInd, closestD, closestA, score/max(score)))
    maxScoreInd = which.max(score)
    topScoreInd = closestTabInd[maxScoreInd]
    # add info to the attrTab:
    #print(attrTab$ARCID[topScoreInd])
    outTab[topScoreInd, ] = as.numeric(as.matrix(cbind(score[maxScoreInd], POItab[j, inAttrColInd])))
  }
  POItab[j, inAttrColInd]
  as.numeric(as.matrix(POItab[j, inAttrColInd]))
  
  # add output table to attrTab:
  attrTab = cbind(attrTab, outTab)
  
  return(attrTab)
  # # Plot:
  # plot(tXY[closeRivs, 1], tXY[closeRivs, 2])
  # points(POIXY[j,1], POIXY[j,2], col=2, pch=15)
  # lines(closeXY[, 1], closeXY[, 2])
  # points(closeXY[, 1], closeXY[, 2], cex=0.4, pch=16)
  # lines(closeXY[closestRivs, 1], closeXY[closestRivs, 2], col=4)
  # points(tXY[closestTabInd, 1], tXY[closestTabInd, 2], pch=16)
  # 
  # plot(tXY[closestTabInd, 1], tXY[closestTabInd, 2])
  # line(closeXY[closestRivs, 1], closeXY[closestRivs, 2], col=4)
  # points(POIXY[j,1], POIXY[j,2], col=2, pch=15)
  # x = which(pLineTab$Id == closestTabInd[closestRivsOrd[which.max(score)]])
  # line(pLineTab$X[x], pLineTab$Y[x], col=2)
  # 
}
dnStrGaugeCrawler <- function(tab, cTab){
  # for each gauge, find the downstream gauge and add that index to a table. 
  # The first column of this table lists the UID of the segments with gauges. 
  # the next column is the corresponding downstream segment, etc. 
  gCol = grep('GAUGE_Site|GAUGE_STAI', names(tab))
  gBoo = tab[ ,gCol] != 0
  gID = tab[gBoo ,gCol]
  segID = tab$ARCID[gBoo]
  
  # check to make sure that gauges are only assigned to one segment each:
  if (length(unique(gID)) != length(gID)){
    message("length of unique gauge list is not equal to length of non unique gauge list,
            indicating that single gauges assigned to multiple segments!")
  }
  
  segIDtab = segID
  thisID = segID
  nextID = cTab$DNSTR[thisID]
  
  while (T %in% c(nextID>0)){
    thisID = nextID
    zInd = thisID == 0
    thisID[zInd] = NA
    segIDtab = cbind(segIDtab, thisID)
    nextID = cTab$DNSTR[thisID]
  }
  
  return(segIDtab)
  
  }
lagRangeCalc = function(minCel, maxCel, dnStrDist){
  lagRange = ceiling(cbind(dnStrDist/maxCel, dnStrDist/minCel)*0.01157407) # km to m, and sec to days conversion
  lagRange[lagRange > 200] = 200
  return(lagRange)
}
QreadAndProc = function(qTab, quantile){
  # find which column contains dsicahrge records:
  colNames = names(qTab)
  #ind = grep("00060_00003|00060_00001|X01_00060_00011|00060_32400", colNames)
  ind = grep("00060", colNames)
  Qcol = ind[grep("_cd", colNames[ind], inv=T)][1]
  # only consider discharges > 90th percentile:
  Q = qTab[, Qcol]
  Q = suppressWarnings(as.numeric(levels(Q))[Q])
  highQ = quantile(Q, quantile, na.rm=T)[[1]]
  Q[Q<highQ] = NA
  return(Q)
}
lagCor = function(q1, q2, lagRange, k){
  ########
  # over a reasonable range of celerities (<10 mps), find 
  # what time lag corresponds to the maximum correlation:
  lagWin = 1:lagRange[k,2]
  corVec = rep(NA, lagRange[k,2])
  for (l in lagWin){
    q1 = c(q1, NA)
    q2 = c(NA, q2)
    corVec[l] = cor(q1, q2, use="pairwise.complete.obs")
  }
  return(corVec)
}
lagCorrPlot <- function(tab, corVec, lag_day, gDnStrmDist, modelTT, i, j, k, Q1, Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID){
  
  par(mfrow=c(2,1))
  par(mar=c(4.1,5.1,2.1,7.1))
  # plot lag-correlation:
  plot(corVec, type='l',
       ylim = c(0, 1),
       xlab = "Lag (days)",
       ylab = "Correlation, R")
  points(lag_day, corVec[lag_day], cex=0.8)
  abline(v=lagRangeCalc(0.5, 2.5, gDnStrmDist[k]), lty=2, col='gray')
  abline(v=modelTT)
  mtext(paste("i:", i,"   j:", j, "    k:", k))
  
  # plot hydrograph time series and lag correlations
  Q1[Q1==0] = 1
  Q2[Q2==0] = 1
  ylim = range(c(Q1, Q2, na.rm=T)*cfs2mfsConv, na.rm=T)
  plot(dates, Q1*cfs2mfsConv,
       ylim=ylim,
       xlab="Year",
       ylab="Flow (cms)",
       #main=paste(i, j, k)
       #log='y',
       type='l', col=1, lwd=0.4, las=1)
  # par(new=T)
  # plot(c(rep(NA,lag_day),dates), c(rep(NA,lag_day), Q2)*cfs2mfsConv,
  #      ylim=ylim,
  #      yaxt = "n",
  #      ylab = NA, xlab = NA,
  #      #log='y',
  #      type='l', col='green', lwd=0.3,
  #      axes=F)
  par(new=T)
  plot(dates, Q2*cfs2mfsConv,
       ylim=ylim,
       yaxt = "n",
       ylab = NA, xlab = NA,
       #log='y',
       type='l', col='blue', lwd=0.3,
       axes=F)
  axis(4, las=1, col=4, col.ticks=4, col.axis=4)
  corns = par("usr"); par(xpd=T)
  text(x=corns[2]+40, y=mean(corns[3:4]),
       'Flow (cms)',
       srt=270,
       col=4)
  
  mtext(paste("dist km:", round(gDnStrmDist[k]),
              "   cel:", round(cel_mps, 2), "m/s",
              "   R:", round(R, 3),
              "   Q ovrlp:", Qoverlap))
  mtext(paste("Area dif:", round(gAreaDif_per), "%   modCel:",
              round(sum(tab$LENGTH_KM[ID])/(modelTT)*0.01157407,2)), line=-1) #*1.3# convert from m/s to days and km to m
}
celerityModel_mann_rect <- function(tab, B){
  # flow wave celerity model:
  # segment length (m): --> hydrosheds
  L = tab$LENGTH_KM*1e3 * 1.3 # convert to meters and sinuosity correction
  # Manning's roughness coefficient (s/m^0.33):
  n = 0.035
  # bankful width (m): --> Andreadis et al
  w = tab$WIDTH
  # bankful depth (m): --> Andreadis et al
  h = tab$DEPTH
  # river slope (m/m): --> hydrosheds
  S = tab$SLOPE / 1.3 # sinuosity correction
  # hydraulic radius
  R = w * h / (2*h + w)
  # manning's equation for flow velocity (m/s):
  v = n^(-1) * R^(2/3) * S^(1/2)
  # Wave celerity equation from Lighthill and Whitham (1955) (m/s):
  celerity = v * B 
  # wave propagation time (days):
  propT = L / (celerity*86400) # convert from s to days
  # add information to tab:
  tab$CELER_MPS = celerity # (m/s)
  tab$PROPTIME_DAY = propT # (days)
  
  return(tab)
}
celerityModel_mann_tri <- function(tab, B){
  # flow wave celerity model:
  # segment length (m): --> hydrosheds
  L = tab$LENGTH_KM*1e3 * 1.3 # convert to meters and sinuosity correction
  # Manning's roughness coefficient (s/m^0.33):
  n = 0.035
  # bankful width (m): --> Andreadis et al
  w = tab$WIDTH
  # bankful depth (m): --> Andreadis et al
  h = tab$DEPTH
  # river slope (m/m): --> hydrosheds
  S = tab$SLOPE / 1.3 # sinuosity correction
  # hydraulic radius
  R = ((w/(2*h))*h^2)/(2*h*(1+(w/(2*h))^2))^0.5 # triangular cross section
  # manning's equation for flow velocity (m/s):
  v = n^(-1) * R^(2/3) * S^(1/2)
  # Wave celerity equation from Lighthill and Whitham (1955):
  celerity = v * B 
  # wave propagation time (days):
  propT = L / (celerity*86400) # convert from s to days
  # add information to tab:
  tab$CELER_MPS = celerity # (m/s)
  tab$PROPTIME_DAY = propT # (days)
  
  return(tab)
}
nextDnStrmPOIlength <- function(tab, cTab, mouthInd, POIbool, outField){
  outFieldInd = which(names(tab)==outField)
  mouthInd = mouthInd[!POIbool[mouthInd]]
  thisUID = cTab$UID[mouthInd]
  tM = match(thisUID, tab$ARCID)
  tab[tM, outFieldInd] = 0
  lastInd = mouthInd
  thisM = match(cTab$DNSTR, cTab$UID[lastInd])
  thisInd = which(!is.na(thisM))
  while (length(thisInd) > 0){
    # remove from iteration segments which have a dam located on them:
    thisInd = thisInd[!POIbool[thisInd]]
    thisUID = cTab$UID[thisInd]
    dnstrUID = cTab$UID[lastInd][thisM[thisInd]]
    tM = match(thisUID, tab$ARCID)
    dM = match(dnstrUID, tab$ARCID)
    # zero-out all cumulative times for segments that do not have a dam downstream:
    tab[tM, outFieldInd] = 0
    # set last index to current index:
    lastInd = thisInd
    thisM = match(cTab$DNSTR, cTab$UID[lastInd])
    thisInd = which(!is.na(thisM))
  }
  return(tab)
}
cumRivTime <- function(tab, cTab){
  # start with river segment at bottom of network and 
  # work upstream to calculate cumulative length and time:
  
  # identify which segments are at the bottom of networks:
  mouthInd = which(cTab$DNSTR == 0)
  thisUID = cTab$UID[mouthInd]
  
  # match current ID in cTab with ID in tab:
  tM = match(thisUID, tab$ARCID)
  tab$CITY_UPSTR_TIME_DAY[tM] = tab$DAM_UPSTR_TIME_DAY[tM] = 0
  
  # add current segment length & time to upstream length & time:
  tab$UPSTR_DIST_KM[tM] = tab$LENGTH_KM[tM]
  tab$UPSTR_TIME_DAY[tM] = tab$CITY_UPSTR_TIME_DAY[tM] = tab$DAM_UPSTR_TIME_DAY[tM] = tab$PROPTIME_DAY[tM]
  
  # determine which segments have cities and dams. Adjust thresholds
  # to limit cities and dams to certain sizes or join score quality:
  cityBool = tab$CITY_POP_M > 1e5 & tab$CITY_JOINS > 1.148072e-5 # 10%ile score
  damBool = tab$DAM_CAP_MC > 0 # & tab$DAM_JOINSC > 0 
  
  # set matching index to mouth index:
  lastInd = mouthInd
  thisM = match(cTab$DNSTR, cTab$UID[lastInd])
  thisInd = which(!is.na(thisM))
  
  # iterate:
  #  j = 0
  ptm = proc.time()
  
  # match downstream IDs to current IDs:
  while (length(thisInd) > 0){
    
    thisUID = cTab$UID[thisInd]
    dnstrUID = cTab$UID[lastInd][thisM[thisInd]]
    
    tM = match(thisUID, tab$ARCID)
    dM = match(dnstrUID, tab$ARCID)
    
    # add downstream data to data of current segment length:
    tab$UPSTR_DIST_KM[tM] = tab$UPSTR_DIST_KM[dM] + tab$LENGTH_KM[tM]
    tab$UPSTR_TIME_DAY[tM] = tab$UPSTR_TIME_DAY[dM] + tab$PROPTIME_DAY[tM]
    
    # zero-out all cumulative time for segments that don't have a 
    # city or dam downstream:
    tab$CITY_UPSTR_TIME_DAY[tM] = tab$CITY_UPSTR_TIME_DAY[dM] + tab$PROPTIME_DAY[tM]
    tab$DAM_UPSTR_TIME_DAY[tM] = tab$DAM_UPSTR_TIME_DAY[dM] + tab$PROPTIME_DAY[tM]
    tab$CITY_UPSTR_TIME_DAY[cityBool] = 0
    tab$DAM_UPSTR_TIME_DAY[damBool] = 0
    
    # set last index to current index:
    lastInd = thisInd
    thisM = match(cTab$DNSTR, cTab$UID[lastInd])
    thisInd = which(!is.na(thisM))
    #j = j+1
  }
  
  #### TIME TO THE NEXT POI
  # zero-out the cumulative time of river segments that are
  # not upstream from a city or dam (just upstream of the basin outlet):
  # TIME TO NEXT DOWNSTREAM CITY:
  # zero-out segments which have a city located on them:
  tab = nextDnStrmPOIlength(tab, cTab, mouthInd, POIbool=cityBool, outField="CITY_UPSTR_TIME_DAY")
  # TIME TO NEXT DOWNSTREAM DAM:
  # zero-out segments which have a dam located on them:
  tab = nextDnStrmPOIlength(tab, cTab, mouthInd, POIbool=damBool, outField="DAM_UPSTR_TIME_DAY")
  
  #print(paste("N upstream iterations:", j))
  print(proc.time() - ptm)
  return(tab)
}
distPlot <- function(x, weight, uprQuant, breaks=NA, j, xLab, yLab, makeTab=T, latencies = NA, latencyTab=NA){
  
  # split segments into 1 km segments and multiply length 
  # by the number of SWOT overpasses at that segment:
  x = rep(x, round(weight), each=T)
  
  # define upper quantile
  quant = quantile(x, uprQuant)
  
  if (is.na(breaks[1])){
    if (quant-min(x) > 10){
      breaks = seq(0, ceiling(quant), 1)
    }else{
      breaks = seq(0, ceiling(quant), length.out=20)
    }
  }else{
    quant = max(breaks)
  }
  
  xlim = c(min(x), quant)#range(x)#
  
  par(mar=c(5.1,5.1,2.1,4.1))
  h = hist(x[x<xlim[2]], breaks=breaks, freq=T,
           xlim=xlim, 
           main='', #'Global distribution of flow celerity',
           xlab=xLab,
           ylab="",
           las=1,
           border=0,
           col=rgb(0.7, 0.7, 0.7, 1))
  mtext(yLab, 2, line = 4)
  
  cdf = ecdf(x)
  cdfXseq = seq(xlim[1], xlim[2], length.out=100)
  par(new=T); plot(cdfXseq, 
                   cdf(cdfXseq),
                   xlim=xlim, 
                   yaxt = "n", 
                   ylab = NA,
                   xlab = NA,
                   type='l',
                   lwd=2,
                   col=4)
  axis(4,
       las=1,
       col=4,
       col.ticks=4,
       col.axis=4)
  corns = par("usr"); par(xpd=T)
  text(x=corns[2]+(corns[2]-corns[1])/4, y=mean(corns[3:4]), 
       'Probability', 
       srt=270,
       col=4)
  text(corns[1]-(corns[2]-corns[1])/16, corns[4]+0.05, 
       letters[j], 
       cex=1.5,
       font=2)
  
  # add median point to plot:
  med = median(x)
  points(med, 0.5, col=4)
  text(med, 0.5, paste0("(",round(med,1), ", 0.5)"), pos=4, col=4)
  
  # add latency information to latency table:
  if (nrow(latencyTab)>1){
    latencyVec = paste0(100*round(1-cdf(latencies), 2), "%")
    latencyTab = as.data.frame(cbind(latencyTab, latencyVec))
    return(latencyTab)
  }
}

####### SLOPE CALC ###### 
# if slope is not yet calculated (or needs to be recalculated),
# loop through each river file and calculate slope:
for (i in 1:length(rivEndPtTabFpaths)){
  tab = foreign::read.dbf(rivEndPtTabFpaths[i])

  # clean up table:
  names(tab)[names(tab) == "RASTERVALU"] = "ELEV_M"
  names(tab)[names(tab) == "LENGTH_GEO"] = "LENGTH_KM"
  
  # check if there are any unpaired ORIG_FIDs:
  oddInd = seq(1, nrow(tab), 2)
  evenInd = seq(2, nrow(tab), 2)
  unpairedInd = which(tab$ORIG_FID[evenInd]-tab$ORIG_FID[oddInd] != 0)
  if (length(unpairedInd) > 0){
    message(paste(length(unpairedInd), "unpaired segment endpoints. Remove these?"))
  }
  
  # read in and set up connectivity table:
  j = grep(rivEndPtTabNames[i], rivNetworkConnectFnames)
  cTab = read.csv(rivNetworkConnectFpaths[j], header=F)
  names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
                  "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
                  "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
  
  # calculate slope:
  ptm <- proc.time()
  sTab = slopeCalc(tab, cTab, zeroSlope=1e-4)
  print(proc.time() - ptm)
  
  hist(sTab$SLOPE, 100, main=i)
  
  # add tab data to polyline shapefile data:
  # replace previously modified dbf file with new copy:
  if (file.exists(origPolylinesFpaths[i])){
    file.remove(slopeOutFpaths[i])
    file.copy(origPolylinesFpaths[i], slopeOutFpaths[i])
  }
  
  # transfer data to polyline vectors:
  pLineTab = foreign::read.dbf(slopeOutFpaths[i])
  m = match(pLineTab$ARCID, sTab$ARCID)
  pLineTab = sTab[m, ]
  foreign::write.dbf(pLineTab, slopeOutFpaths[i])
  
  print(paste(i, slopeOutFpaths[i], "slope calc done run"))
  
}
####### POI JOIN TO RIVER NETWORK #####

# read in city table:
cityTab = foreign::read.dbf(cityFpath);  
cityTab = cityTab[match(c('NAME','LATITUDE','LONGITUDE','POP_MAX','POP_MIN'), names(cityTab))]
cityXY = cbind(cityTab$LONGITUDE, cityTab$LATITUDE)

# read in dam table:
damTab = foreign::read.dbf(damFpath);  
damTab = damTab[match(c('DAM_NAME','LAT_DD','LONG_DD','YEAR','AREA_SKM','CAP_MCM'), names(damTab))]
damXY = cbind(damTab$LONG_DD, damTab$LAT_DD)

# read in gauge table:
gaugeTab = foreign::read.dbf(gaugeFpath);
#gaugeTab = gaugeTab[match(c('STAID','LAT_GAGE','LNG_GAGE', 'DRAIN_SQKM', 'HUC02'), names(gaugeTab))]
#gaugeXY = cbind(gaugeTab$LNG_GAGE, gaugeTab$LAT_GAGE)
# added lat lon fields using WGS84 in arc:
gaugeTab = gaugeTab[match(c('Site_NO','lat','lon', 'HUC8'), names(gaugeTab))]
gaugeXY = cbind(gaugeTab$lon, gaugeTab$lat)

# run through each continent and match up cities, dams, and gauges to the river networks: 
shpFpaths = sub('.dbf', '.shp', slopeOutFpaths)
for (i in c(4,6)){#1:length(shpFpaths)){
  print(i)
  # extract the verticies of the polylines # (very slow - took africa 3 hours to read in):
  # consider trying much faster function, sf::st_combine(st_read(shpFpaths[i])))
  ptm = proc.time()
  pLineTab = read.shp(shpFpaths[i])
  print(proc.time() - ptm)
  ptm = proc.time()
  pLineTab = convert.to.simple(pLineTab)
  print(proc.time() - ptm)
  
  # extract end point verticies and add a city field to tab:
  attrTab = foreign::read.dbf(slopeOutFpaths[i])
  
  # for each city, find the river segment that the city is most likely to be located on:
  ptm = proc.time()
  attrTab = POI2RivJoin(
    POItab = cityTab, 
    POIXY = cityXY, 
    pLineTab = pLineTab, 
    attrTab = attrTab, 
    inAttr = c('POP_MAX'), # POI table column name(s) to be transfered to the output attribute table. 
    # right now, only numerical attributes can be transferred between tables. 
    POIname = 'CITY', # prefix for what type of POI data this is. Used for output table column names
    scoreFunc = function(area, dist){ return(area/dist^2) }, # distance to score the relationship between POI and segments
    sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
    smallSrad = 1e4) # specify radius of secondary search radius (in meters)
  print(proc.time() - ptm)
  
  # for each dam, find the river segment that the dam is most likely to be located on:
  ptm = proc.time()
  attrTab = POI2RivJoin(
    POItab = damTab, 
    POIXY = damXY, 
    pLineTab = pLineTab, 
    attrTab = attrTab, 
    inAttr = c('AREA_SKM', 'CAP_MCM'), # POI table column name(s) to be transfered to the output attribute table. 
    # right now, only numerical attributes can be transferred between tables. 
    POIname = 'DAM', # prefix for what type of POI data this is. Used for output table column names
    scoreFunc = function(area, dist){ return(area/dist^3) }, # distance to score the relationship between POI and segments
    sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
    smallSrad = 3e3) # specify radius of secondary search radius (in meters)
  print(proc.time() - ptm)
  
  # for each gauge, find the river segment that the gauge is most likely to be located on:
  ptm = proc.time()
  attrTab = POI2RivJoin(
    POItab = gaugeTab, 
    POIXY = gaugeXY, 
    pLineTab = pLineTab, 
    attrTab = attrTab, 
    inAttr = c('Site_NO', 'HUC8'), # POI table column name(s) to be transfered to the output attribute table. 
    #inAttr = c('STAID', 'DRAIN_SQKM', 'HUC02'), # POI table column name(s) to be transfered to the output attribute table. 
    # right now, only numerical attributes can be transferred between tables. 
    POIname = 'GAUGE', # prefix for what type of POI data this is. Used for output table column names
    scoreFunc = function(area, dist){ return(area/dist^3) }, # distance to score the relationship between POI and segments
    sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
    smallSrad = 2e3) # specify radius of secondary search radius (in meters)
  print(proc.time() - ptm)
  
  # write out POI table:
  foreign::write.dbf(attrTab, cityDamGaugeFpaths[i])
}



####### CELERITY & TRAVEL TIME #####
# Calculate flow wave celerity and travel time for each segment:
for (i in 1:length(cityDamGaugeFpaths)){ # ){ # 
  # read in river attribute table:
  tab = foreign::read.dbf(cityDamGaugeFpaths[i])

  # TEMP TEMP TEMP : replace NA slope with 1e-4:
  tab$SLOPE[is.na(tab$SLOPE)] = 1e-4
  # TEMP TEMP TEMP 
  
  # read in and set up connectivity table:
  j = grep(rivEndPtTabNames[i], rivNetworkConnectFnames)
  cTab = read.csv(rivNetworkConnectFpaths[j], header=F)
  names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
                  "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
                  "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
  
  # calculate flow wave celerity for each segment:
  # manning rect: β = 5/3; Chezy rect: β = 3/2; manning tri: 4/3; chez tri: 5/4
  tab = celerityModel_mann_rect(tab, B=5/3)
  
  # Add new columns to tab:
  SWOT_TRAC_DEN = 
    hBASIN = # main basin ID from hydroBASINS
    GLCC = # Global Land Cover Classification
    FLOODHAZARD = # flood hazard index
    UPSTR_DIST_KM = # upstream distance in km
    UPSTR_TIME_DAY = # upstream flow wave travel time in days
    CITY_UPSTR_TIME_DAY = # upstream from nearest city flow wave travel time in days
    DAM_UPSTR_TIME_DAY = # upstream from nearest dam flow wave travel time in days
    rep(NA, nrow(tab))
  tab = cbind(tab, hBASIN, GLCC, FLOODHAZARD, SWOT_TRAC_DEN, 
              UPSTR_DIST_KM, UPSTR_TIME_DAY, CITY_UPSTR_TIME_DAY, DAM_UPSTR_TIME_DAY)
  
  # read in midPtTab:
  midPtTab = foreign::read.dbf(rivMidPtTabFpaths[midPtMatch[i]])
  
  # join data from midPoints tab:
  m = match(tab$ARCID, midPtTab$ARCID) # match UIDs (prob unnecessary)
  tab$hBASIN = midPtTab$MAIN_BAS[m] #  main basin UID for the hydroBASINS dataset (http://www.hydrosheds.org/page/hydrobasins)
  tab$GLCC = midPtTab$GLCC[m] # land cover dataset (https://lta.cr.usgs.gov/glcc/globdoc2_0) #DN: 20 <- desert
  tab$FLOODHAZARD = midPtTab$RASTERVALU[m] # flood hazard composite from the DFO (via NASA): #avhrr
  tab$SWOT_TRAC_DEN = midPtTab$COUNT_coun[m] # swot track density (N overpasses per cycle @ segment centroid)
  

  # calculate the cumulative length and time along
  # entire river network:
  tab = cumRivTime(tab, cTab)
  # consider adding reverse stream order field 
  # consider adding a continent UID field

  # add tab data to polyline shapefile data:
  # replace previously modified file with new copy:
  if (file.exists(origPolylinesFpaths[i])){
    file.remove(rivOutFpaths[i])
    file.copy(origPolylinesFpaths[i], rivOutFpaths[i])
  }
  
  # transfer data to polyline vectors:
  pLineTab = foreign::read.dbf(rivOutFpaths[i])
  m = match(pLineTab$ARCID, tab$ARCID)
  pLineTab = tab[m, ] #cbind(pLineTab, tab[m, ((ncol(pLineTab)+1):ncol(tab))])
  #print(cbind(1:ncol(pLineTab), names(pLineTab)))
  foreign::write.dbf(pLineTab, rivOutFpaths[i])

  print(paste(i, rivEndPtTabNames[i], "done run"))
}


###### VALIDATION - GAUGE CROSS CORRELATION ANALYSIS ######

# set up PDF device:
pdfOut = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/outputData/plots/crossCorrelations2.pdf"
pdf(pdfOut, width = 7, height=5)

# set the minimum days of overlap between the two gauge records for
# them to be considered in the lag correlation analysis:
minOvrlp = 5*365

ptm = proc.time()
# run through each continent and search for gauges:
for (i in c(4,6)){
  
  # read in segment attribute table and add empirical celerity & R columns:
  tab = foreign::read.dbf(rivOutFpaths[i])
  nCol1 = ncol(tab)
  tab = cbind(tab, "OBS_CEL_R"=0, "OBS_CEL_MPS"=0, "OBS_CEL_DIST_KM"=0, 
              "OBS_CEL_AREADIF_PER"=0, "OBS_CEL_LAG_D"=0, "Q_OVRLP_DAYS"=0, "CORRANGE"=0) 
  newColInd = (1+nCol1):ncol(tab)
  
  #### #### 
  # create a data frame containing , , ,
  # distance 
  celTab = tableMaker(colNames=c("ID", #segment ID
                                 'R', # max lag correlation
                                 "cel_mps", # celerity
                                 'gDist_km', # distance between two correlated gauges
                                 'gAreaDif_per', # area difference between correlated gauges
                                 'lag_day', # lag time of max correlation
                                 'Qoverlap',# N days of overlapping Q data between 2 gauges 
                                 'Rrange'),  # range of cross correlations
                      nrows=0, fill=NA)
  
  # read in and set up connectivity table:
  cityDamGaugeNames = substr(rivEndPtTabFnames, 1, 2)
  ii = grep(cityDamGaugeNames[i], rivNetworkConnectFnames)
  cTab = read.csv(rivNetworkConnectFpaths[ii], header=F)
  names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
                  "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
                  "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
  
  #### for each gauge, find downstream gauges, and their downstream distances and drainage areas: ######
  segIDtab = dnStrGaugeCrawler(tab, cTab)
  
  # get the downstream gauges, drainage areas, and cumulative downstream distance:
  segIDtab = matrix(segIDtab, nrow=nrow(segIDtab))
  gCol = grep('GAUGE_Site|GAUGE_STAI', names(tab))
  gTab = matrix(tab[segIDtab, gCol], nrow=nrow(segIDtab))
  areaTab = matrix(tab$AREA[segIDtab], nrow=nrow(segIDtab))
  distTab =  matrix(tab$LENGTH_KM[segIDtab], nrow=nrow(segIDtab))
  cumDistTab = t(apply(distTab, 1, cumsum))
  
  # match up gauges linked to river network to file list. First, 
  # get list of downstream gauge IDs and paths::
  gFmatch_raw = match(gTab[,1], gNames)
  notNA = !is.na(gFmatch_raw)
  gFmatch = gFmatch_raw[notNA]
  gRows = which(notNA)
  g1Paths = gFpaths[gFmatch]
  
  # set up empty tables for calculated gauge parameters:
  #lagTimeTab = corTab = modelTTtab = QoverlapTab = matrix(NA, nrow=nrow(gTab), ncol=ncol(gTab))
    
  # for each matching gauge, get the cross correlation of each downstream gauge
  # and assign that to segments between the two gauges:
  for (j in 1:length(gRows)){
    print(j)
    gRow = gRows[j]
    gVec = gTab[gRow, ]
    gDnStrmInd = gVec != 0 & !is.na(gVec) & gVec != gVec[1]
    if (length(which(gDnStrmInd))==0){next} #print('no downstream gauges'); 
    
    # get gauge file paths, downstream distances, and upstream areas of 
    # gauges located downstream:
    gDnStrmID = gVec[gDnStrmInd]
    gFmatch_dnStrm_raw = match(gDnStrmID, gNames)
    gFmatch_dnStrmInd = which(!is.na(gFmatch_dnStrm_raw))
    if (length(gFmatch_dnStrmInd)==0){next}#print('no downstream gauges on file'); 
    
    gFmatch_dnStrm = gFmatch_dnStrm_raw[gFmatch_dnStrmInd]
    g2paths = gFpaths[gFmatch_dnStrm]
    gDnStrmDist = cumDistTab[gRow, gDnStrmInd][gFmatch_dnStrmInd]
    gDnStrmArea = areaTab[gRow, gDnStrmInd][gFmatch_dnStrmInd]
    
    lagRange = lagRangeCalc(minCel=0.1, maxCel=10, dnStrDist=gDnStrmDist)
    
    # only consider gauges that are close together that do not have a 
    # large difference in drainage areas: 
    #keepers = which((gDnStrmDist<1000 & (gDnStrmArea-areaTab[gRow,1])/gDnStrmArea<0.5))
    #if (length(keepers)==0){next}
    
    #### lag correlation analysis:
    g1 = read.table(g1Paths[j], sep=",", fill=T, header=T)[-1,]
    # remove low flows from upper gauge discharge data:
    Q1_orig = QreadAndProc(qTab=g1, quantile=0.0)
    
    # skip if gauge doesn't contain more than 5 years of flow data:
    if (length(which(!is.na(Q1_orig)))<minOvrlp){next}
    
    # this could be sped up by reading in everything and put discharge data into a large matrix
    # so that correlations could be run all at once but this approach is much more complicated 
    # than individually analyzing correlations one by one:
    for (k in 1:length(gFmatch_dnStrm)){ #for (k in (1:length(gFmatch_dnStrm))[keepers]){
      g2 = read.table(g2paths[k], sep=",", fill=T, header=T)[-1,]
      # match up the dates between the two gauge records:
      dM = match(g1$datetime, g2$datetime)
      dInd = which(!is.na(dM))
      if (length(dInd)<minOvrlp){next}
      # remove low flows from discharge data:
      Q1 = Q1_orig
      Q2 = QreadAndProc(qTab=g2, quantile=0.0)
      
      # skip if gauge doesn't contain more than 5 years of flow data:
      if (length(which(!is.na(Q2)))<minOvrlp){next}
      
      # plot hydrograph timeseries over entire records:
      # xlim = range(c(as.Date(g1$datetime), as.Date(g2$datetime)), na.rm=T)
      # ylim = range(c(Q1, Q2)*cfs2mfsConv, na.rm=T)
      # plot(as.Date(g1$datetime), Q1*cfs2mfsConv, type='l',
      #     xlim=xlim, ylim=ylim, main=paste(j, k), ylab="Flow (cms)", col=1, lwd=0.4)
      #lines(as.Date(g2$datetime), Q2*cfs2mfsConv, type='l', col=4, lwd=0.4)
      
      # match up two discharge records:
      Q1_match = Q1[dInd]
      Q2_match = Q2[dM[dInd]]
      # remove blank discharge measurements:
      keep = Q1_match!='' & Q2_match!='' #& Q1_match!=0 & Q2_match!=0 
      dates = as.Date(g1$datetime[dInd], "%Y-%m-%d")
      Q1 = as.numeric(Q1_match)
      Q1[!keep] = NA
      Q2 = as.numeric(Q2_match)
      Q2[!keep] = NA
      # skip to next gauge pair if there is not enough overlapping gauge records:
      Qoverlap = length(which(!is.na(Q1) & !is.na(Q2)))
      if (Qoverlap<minOvrlp){next}
      
      
      ####
      # use lag correlation to find the optimum lag time:
      #corVec = lagCor(q1=Q1, q2=Q2, lagRange=lagRange, k=k)
      corVec_raw = ccf(Q1, Q2, lag.max=200, na.action=na.pass, plot=F)
      positiveLag = corVec_raw$lag > 0 #lagRange[k,1]
      corVec = corVec_raw$acf[positiveLag]
      # Pick the lag with the maximum correlation:
      lag_day = which.max(corVec)
      
      # skip gauge pair if the folling conditions are not met:
      if (!T %in% !is.na(corVec)){next}
      if (length(lag_day)==0){next}
      if (is.na(lag_day)){next}
      # range and maxmimum of cross correlations must be greater than specified values:
      corRange = range(corVec)
      if ((corRange[2]-corRange[1]) < 0.5){next}
      if (corRange[2] < 0.5){next}
      # If peak correlation occurs on the first or last day of the lag window, 
      # remove it from consideration because it the maximum lag might
      # just be truncated by the window, rendering it meaningless. 
      if (T %in% c(lag_day == c(1,200))){next}
      
      
      ####
      # assign celerity, correlation, and other parameters to segments between two gauges:
      # this could be accomplished using less memory by creating a table that does not
      # include each individual segIDs and only 1 row per celerity result...
      ID = segIDtab[gRows[j], 1:which(gVec == gDnStrmID[k])]
      
      # # calculate correlation coefficient and celerity:
      R = corVec[lag_day]
      Rrange = corRange[2]-corRange[1]
      cel_mps = (gDnStrmDist[k] / lag_day)*0.01157407 *1.3 # convert from m/s to days and km to m
      gDist_km = gDnStrmDist[k]
      gAreaDif_per = 100*(gDnStrmArea[k]-areaTab[gRow,1])/gDnStrmArea[k]
      celTabMat = cbind(ID, R, cel_mps, gDist_km, gAreaDif_per, lag_day, Qoverlap, Rrange)
      celTab = rbind(celTab, celTabMat)
      
      # fill in tabs:
      modelTT = sum(tab$LENGTH_KM[ID]/(tab$CELER_MPS[ID]))*0.01157407 *1.3# convert from m/s to days and km to m
      modelTTtab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = modelTT
      lagTimeTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = lag_day
      corTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = corVec[lag_day]
      QoverlapTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = Qoverlap
      
      ####
      # Plot hydrograph and lag cross correlation analysis:
      lagCorrPlot(tab, corVec, lag_day, gDnStrmDist, modelTT,  i, j, k, 
                  Q1, Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID)
      
    }
  }
  
  write.csv(celTab, celTabOutFpaths[i], row.names=F)
  
  celTab = read.csv(celTabOutFpaths[i], header=T)
  
  print(proc.time() - ptm)
 
  # join empirical celerity data onto flowlines, giving preference for higher correlation coeficients:
  uID = unique(celTab$ID)
  celTabColInd = 2:ncol(celTab)
  for (j in 1:length(uID)){
    uIDmatch = which(celTab$ID == uID[j])
    tabInd = tab$ARCID == uID[j]
    # use a correlation-weighted mean to caluclate celerity:
    tab[tabInd, newColInd] = apply(celTab[uIDmatch, celTabColInd], 2, function(x) weighted.mean(x, celTab[uIDmatch, 2]))

    # alternative: use the value with the maximum correlation:
    #celTabInd = uIDmatch[which.max(celTab$R[uIDmatch])]
    #tab[tabInd, newColInd] = celTab[celTabInd, celTabColInd]
  }
  

  foreign::write.dbf(tab, obsOutFpaths[i])

  keep = tab$OBS_CEL_MPS>0 # & tab$WIDTH>100# & tab$Q_OVRLP_DAYS>1e3 & tab$OBS_CEL_R<1 & tab$OBS_CEL_MPS > 0  & tab$OBS_CEL_MPS < 5
  x = tab$OBS_CEL_MPS[keep]
  x = rep(x, round(tab$LENGTH_KM[keep]), each=T)
  hist(x, 40, main="hist of celerity", xlab="celerity (m/s)")

  x = tab$OBS_CEL_R[keep]
  x = rep(x, round(tab$LENGTH_KM[keep]), each=T)
  hist(x, 50, main="hist of lag correlations with R > 0.5", xlab="R")
  
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)

# PLOT VALIDATION DISTRIBUTIONS (FIGURE 2):

# read in and append each continent:
for (i in c(4,6)){
  tab = foreign::read.dbf(obsOutFpaths[i])
  if (i == 4){ 
    gTab = tab 
  }else{ 
    gTab = rbind(gTab, tab) 
  }
  print(paste(i, rivEndPtTabNames[i]))
}

xlim = c(0,10)
ylim = c(0,5e3)
keep = gTab$OBS_CEL_R>0  & gTab$CELER_MPS<xlim[2] & gTab$OBS_CEL_MP<xlim[2] & gTab$OBS_CEL_R>.5 & gTab$CORRANGE > 0.5 #& gTab$WIDTH>100
x1 = gTab$CELER_MPS[keep]
x1 = rep(x1, round(gTab$LENGTH_KM[keep]), each=T)
x2 = gTab$OBS_CEL_MP[keep]
x2 = rep(x2, round(gTab$LENGTH_KM[keep]), each=T)
h1 = hist(x1, seq(0,100,.2), plot=F)
h2 = hist(x2, seq(0,100,.2), plot=F)
ylim = range(c(h1$counts, h2$counts))

hist(x1, seq(0,xlim[2],.2), main="", #Modeled vs Observed Celerity", 
     xlim=xlim,
     ylim=ylim,
     xlab="Flood Wave Velocity (m/s)", 
     ylab="River Network Length (km)",
     col="light blue",#rgb(0,.5,0,.5), 
     border=1)
par(new=T)
hist(x2, seq(0,xlim[2],.2), main='',#main='WIDTH>100',
     xlim=xlim, 
     ylim=ylim,
     xlab="", 
     ylab="",
     col=rgb(.5,.5,.5,.5),
     border=1)
legend("topright", legend=c("Observations", "Model"), 
       text.col=c(rgb(.5,.5,.5,.5),'light blue'), 
       text.font=2, cex=1.4, lwd=0, box.lwd=0)

# calcualte validation stats:

# RMSE:
rmse <- function(error){sqrt(mean(error^2))}
# Bias (mean error):
me <- function(error){mean(error)}
# Standard Error:
se <- function(error){ sqrt(mean((error-mean(error))^2)) }

# Calculate the error:
error <- gTab$OBS_CEL_MP[keep] - gTab$CELER_MPS[keep]

print(paste("RMSE:", round(rmse(error), 1)))
print(paste("Bias:", round(me(error), 1)))
print(paste("standard error:", round(se(error), 1)))



####### DISTRIBUTION PLOTS - TABLE 1, FIGURE S3 #########
# analyze the global distribution of flow wave travel time:

# read in and append each continent:
for (i in c(1:7)){
  tab = foreign::read.dbf(rivOutFpaths[i])
  #tab = foreign::read.dbf(obsOutFpaths[i])
  ## TEMP TEMP 
  # remove observational columns and gauge columns:
  obsCols = grep("OBS", names(tab))
  if (length(obsCols)>0){
    tab = tab[,-obsCols]
    names(tab)[names(tab) != names(gTab)] =  names(gTab)[names(tab) != names(gTab)]
  }
  obsCols = grep("GAUGE", names(tab))
  if (length(obsCols)>0){
    tab = tab[,-obsCols]
    #names(tab)[names(tab) != names(gTab)] =  names(gTab)[names(tab) != names(gTab)]
  }
  ## TEMP TEMP ^^^
  if (i == 1){ 
    gTab = tab 
  }else{ 
    gTab = rbind(gTab, tab) 
  }
  print(paste(i, rivEndPtTabNames[i]))
}

# thresholds and variables for plotting:
zeroSlope = 1e-4 
uprQuant = 0.99
wide = gTab$WIDTH > 100
realS = gTab$SLOPE != zeroSlope
notDesert = gTab$GLCC != 8
Sof60 = !(gTab$hBASIN %in% as.numeric(basin2rm))
keep = notDesert & Sof60 & realS & wide
breaks = seq(0, 50, 1)
latencies = c(1,2,5,10,45)

# values reported in the text:
# ALL RIVERS:
keep = notDesert & Sof60
# celerity:
x = gTab$CELER_MPS[keep]
print(paste0(round(quantile(x, .5),1), "+",
             round(quantile(x, .75)-quantile(x, .5),1), "-",
             round(quantile(x, .5)-quantile(x, .25),1)))
# travel time:
x = gTab$UPSTR_TIME[keep]
print(paste0(round(quantile(x, .5)), "+",
             round(quantile(x, .75)-quantile(x, .5)), "-",
             round(quantile(x, .5)-quantile(x, .25))))
x = gTab$CITY_UPSTR[gTab$CITY_UPSTR>0 & keep]
print(paste0(round(quantile(x, .5)), "+",
             round(quantile(x, .75)-quantile(x, .5)), "-",
             round(quantile(x, .5)-quantile(x, .25))))
x = gTab$DAM_UPSTR_[gTab$DAM_UPSTR_>0 & keep]
print(paste0(round(quantile(x, .5)), "+",
             round(quantile(x, .75)-quantile(x, .5)), "-",
             round(quantile(x, .5)-quantile(x, .25))))
# swot-observable rivers:
keep = notDesert & Sof60 & wide
x = gTab$UPSTR_TIME[keep]
print(paste0(round(quantile(x, .5)), "+",
       round(quantile(x, .75)-quantile(x, .5)), "-",
       round(quantile(x, .5)-quantile(x, .25))))
x = gTab$CITY_UPSTR[gTab$CITY_UPSTR>0 & keep]
print(paste0(round(quantile(x, .5)), "+",
       round(quantile(x, .75)-quantile(x, .5)), "-",
       round(quantile(x, .5)-quantile(x, .25))))
x = gTab$DAM_UPSTR_[gTab$DAM_UPSTR_>0 & keep]
print(paste0(round(quantile(x, .5)), "+",
       round(quantile(x, .75)-quantile(x, .5)), "-",
       round(quantile(x, .5)-quantile(x, .25))))

####### ALL RIVERS
pdfOut = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/figs/Fig2_distributions_woSWOTwSWOTRivers.pdf"
pdf(pdfOut, width=7, height=9)
layout(matrix(1:6, nrow=3, byrow=F))

# plot the entire dataset:
keep = notDesert & Sof60
latencyTab = as.data.frame(latencies)
names(latencyTab) = "Latency (days)"

# a. BASIN TRAVEL TIME:
# j = j+1
latencyTab = distPlot(x=gTab$UPSTR_TIME[keep], weight=gTab$LENGTH_KM[keep], 
         uprQuant=uprQuant, breaks=breaks, j=j,
         xLab="Travel Time to Basin Outlet(Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Observable global river network length"

# b. CITY TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$CITY_UPSTR>0
latencyTab = distPlot(x=gTab$CITY_UPSTR[keep & nonZeroTravelTime], weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream City (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream city"

# c. DAM TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$DAM_UPSTR_>0
latencyTab = distPlot(x=gTab$DAM_UPSTR_[keep & nonZeroTravelTime], weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream Dam (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream dam"

# write out first part of table 1 for all rivers:
write.csv(latencyTab, sub('.png', '.csv', pngOut), row.names=F)
names(latencyTab) = c("latency", "basins", "cities", "dams"); print(latencyTab)


####### SWOT-OBSERVABLE RIVERS

# plot the SWOT observable rivers:
keep = notDesert & Sof60 & wide
latencyTab = as.data.frame(latencies)
names(latencyTab) = "Latency (days)"

# d. BASIN TRAVEL TIME:
j = j+1
latencyTab = distPlot(x=gTab$UPSTR_TIME[keep], weight=gTab$LENGTH_KM[keep], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Basin Outlet(Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Observable global river network length"

# e. CITY TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$CITY_UPSTR>0
latencyTab = distPlot(x=gTab$CITY_UPSTR[keep & nonZeroTravelTime], weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream City (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream city"

# f. DAM TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$DAM_UPSTR_>0
latencyTab = distPlot(x=gTab$DAM_UPSTR_[keep & nonZeroTravelTime], weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream Dam (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream dam"

# write out second part of table 1 for swot-observable rivers:
write.csv(latencyTab, sub('.pdf', '.csv', pdfOut), row.names=F)
names(latencyTab) = c("latency", "basins", "cities", "dams"); print(latencyTab)

dev.off()
cmd = paste('open', pdfOut)
system(cmd)

