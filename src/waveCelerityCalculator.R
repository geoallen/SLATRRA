# waveCelerityCalculator.R
# George H. Allen, 2017-2018

### To run in RStudio:
# if (!"git2r" %in% rownames(installed.packages())){
#    install.packages("git2r")}; require(git2r)
# fPath2repo = "/Users/geoallen/Documents/research/2017_06_16_waveSpeed/R/SLATTRA/"
# urlPath2rep = "https://github.com/geoallen/SLATRRA"
# if (!file.exists(fPath2repo)){
#   dir.create(fPath2repo)
#   clone(urlPath2rep, fPath2repo)
#   file.edit(fPath2repo)
# }


################################################################################
# DESCRIPTION
################################################################################

# The following script contains the data analysis and production of figures and  
# tables in the Allen et al., "Global estimates of river flow wave travel times  
# and implications for low-latency satellite data"  
# 
# The most basic sections of the code include: 
# 1. joining points of intererests (POIs: cities, dams, and gauges) to the 
# river network
# 2. calculating  river slope 
# 3. modeling flow wave celerity
# 4. calculating travel time
# 5. validating the model with a gauge-based estimate of celerity
# 6. Generating figures and tables, and alculating the various statistics 
# 7. presented in the manuscript 

# For input data, see the associated dataset description section 

################################################################################
# LOAD LIBRARIES
################################################################################
if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"MASS" %in% rownames(installed.packages())){
  install.packages("MASS")}; require(MASS)
if (!"geosphere" %in% rownames(installed.packages())){
  install.packages("geosphere")}; require(geosphere)
if (!"shapefiles" %in% rownames(installed.packages())){
  install.packages("shapefiles")}; require(shapefiles)
if (!"robustbase" %in% rownames(installed.packages())){
  install.packages("robustbase")}; require(robustbase)
if (!"abind" %in% rownames(installed.packages())){
  install.packages("abind")}; require(abind)

################################################################################
# FUNCTIONS
################################################################################
tableMaker = function(colNames, nrows, fill=NA){
  
  # Creates a data frame with specified dimensions
  #
  # Input vars:
  # colNames - vector of column names
  # nrows - value specifying number of rows in output table
  # fill - what to fill the output table with (deafult = NA)
  #
  # Output vars:
  # tab - data frame table 
  
  tab = as.data.frame(array(fill, c(nrows, length(colNames))))
  names(tab) = colNames
  
  return(tab)
  
}

POI2RivJoin <- function(POItab, POIXY, pLineTab, attrTab, inAttr, POIname, 
                        scoreFunc, sRad, smallSrad, plot=F){
  
  # Description:
  # For each point of interest (POI), find the river segment that the 
  # POI is most likely to be located on. Point of interests include cities, 
  # dams and gauges.
  
  # Input vars:
  # POItab - table containing point of interest attributes to be joined to network
  # POIXY - 2 column table containing lat and lon of each point of interest 
  # pLineTab - polyline shapefile table providing lat & lon of flowline vertices
  # attrTab - table with flowline shapefile attributes
  # inAttr - POI table column name(s) to be transfered to output attribute table
  # Right now, only numerical attributes can be transferred between tables
  # POIname - prefix for type of POI data this is. Used for output table col names
  # scoreFunc - distance function to score relationship between POI and segments
  # sRad - value specifying dimensions of the initial square subset (in meters)
  # smallSrad - value specifying radius of secondary search radius (in meters)
  # plot - optional boolean to plot each join for examination (default=F)
  #
  # Output vars:
  # attrTab - original flowline shapefile attributes with join fields added
  
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
  
  # for each point of interest find best matching river segment:
  for (j in POIsubInd){ #nrow(POIXY)){
    # subset river polylines with an endpoints within a square distance 
    # from each POI:
    NESW = destPoint(POIXY[j,], c(0,90,180,270), sRad)
    closeRivs = which(tXY[,1] > NESW[4,1] & tXY[,1] < NESW[2,1] & 
                        tXY[,2] > NESW[3,2] & tXY[,2] < NESW[1,2])
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
    # score river segments by their drainage area and their proximity to 
    # the POI center:
    score = scoreFunc(closestA, closestD)
    #print(cbind(closestTabInd, closestD, closestA, score/max(score)))
    maxScoreInd = which.max(score)
    topScoreInd = closestTabInd[maxScoreInd]
    # add info to the attrTab:
    #print(attrTab$ARCID[topScoreInd])
    outTab[topScoreInd, ] = as.numeric(
      as.matrix(cbind(score[maxScoreInd], POItab[j, inAttrColInd])))
  }
  
  # add output table to attrTab:
  attrTab = cbind(attrTab, outTab)
  
  # Optional Plot:
  if (plot==T){
    plot(tXY[closeRivs, 1], tXY[closeRivs, 2])
    points(POIXY[j,1], POIXY[j,2], col=2, pch=15)
    lines(closeXY[, 1], closeXY[, 2])
    points(closeXY[, 1], closeXY[, 2], cex=0.4, pch=16)
    lines(closeXY[closestRivs, 1], closeXY[closestRivs, 2], col=4)
    points(tXY[closestTabInd, 1], tXY[closestTabInd, 2], pch=16)
    
    plot(tXY[closestTabInd, 1], tXY[closestTabInd, 2])
    line(closeXY[closestRivs, 1], closeXY[closestRivs, 2], col=4)
    points(POIXY[j,1], POIXY[j,2], col=2, pch=15)
    x = which(pLineTab$Id == closestTabInd[closestRivsOrd[which.max(score)]])
    line(pLineTab$X[x], pLineTab$Y[x], col=2)
  }
  
  return(attrTab)
  
}

gmMean <- function(x, na.rm=TRUE){
  
  # Computes the geometric mean
  #
  # Input vars:
  # x - a vector of values
  #
  # Output vars:
  # geomMean - the geometric mean
  
  geomMean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  
  return(geomMean)
  
}

dirCreator <- function(dirPath){
  
  # Create local directory(ies) in specified path locations
  #
  # Input vars:
  # dirPath - a single file path or a vector of paths
  #
  # Output vars:
  # none
  
  for (i in 1:length(dirPath)){
    if (!dir.exists(dirPath[i])){
      dir.create(dirPath[i], recursive=T)
      print(paste("created new directory:", dirPath[i]))
    }
  }
  
}

largerBasin <- function(cTab, mTab, UID, ncolCtab){
  
  # Description:
  # For a given river segment, returns the upstream segment with the greatest
  # drainage area. Takes in connectivity table and info about drainage area:
  #
  # Input vars:
  # cTab - connectivity tab that gives which flowline segment is upstream 
  # mTab - table containing flowline shapefile attributes
  # UID - unique flowline segment ID 
  # ncolCtab - number of non-zero columns in connectivity tab
  #
  # Output vars:
  # upUID - unique ID of the flowline segment ID with the largest drainage area
  
  upUIDs = cTab[UID, 4:ncolCtab]
  upUIDs = upUIDs[upUIDs != 0]
  
  # in the case of two basins being of equal size, the first listed seg is used:
  upUID = upUIDs[which.max(mTab$AREA[upUIDs])]
  
  return(upUID)
  
}

celerityCalc <- function(tab, width, depth, slope, N, B){
  
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
  # Description:
  # calculates flow wave celerity and travel time of each segment of 
  # flowline dataset based on Manning's Formula and the kinematic wave model
  #
  # Input vars:
  # tab - shapefile attribute tatble
  # width - bankful width (m): --> Andreadis et al
  # depth - bankful depth (m): --> Andreadis et al 
  # slope -  river slope (m/m): --> 15 arcsec HydroSHEDS Conditioned DEM
  # N - Manning's roughness coefficient (s/m^0.33)
  # B - kinematic wave coefficient 
  #
  # Output vars:
  # tab - return shapefile attribute table with new fields of celerity and 
  #  flow wave travel time attached
  
  # flow wave celerity model:
  # segment length (m): --> hydrosheds
  
  L = tab$LENGTH_KM*1e3 # convert to m
  
  # Calculate the hydraulic radius (m):
  R = width * depth / (2*depth + width)
  
  # manning's equation for flow velocity (m/s):
  v = N^(-1) * R^(2/3) * slope^(1/2)
  
  # Wave celerity equation from Lighthill and Whitham (1955) (m/s):
  celerity = v * B 
  
  # wave propagation time (days):
  propT = L / (celerity*86400) # convert from s to days
  
  # add information to tab:
  tab$CELER_MPS = celerity # (m/s)
  tab$PROPTIME_DAY = propT # (days)
  
  return(tab)
  
}

nextDnStrmPOIlength <- function(tab, cTab, mouthInd, POIbool, outField){
  
  # Description:
  # Calculate the travel time/ distance to the next downstream POI (city of dam).
  # Zero-out the cumulative time of river segments that are
  # not upstream from a city or dam (just upstream of the basin outlet)
  #
  # Input vars:
  # tab - shapefile attribute tatble
  # cTab - connectivity table
  # mouthInd - index of river mouth
  # POIbool - point of interest boolean. 
  # outField - name of column that is added to tab
  #
  # Output vars:
  # tab - return shapefile attribute table with new outField attached
  
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

cumRivTime_POI <- function(tab, cTab, SLATTRAdir){
  
  # Description:
  # Start with river segment at netowrk mouth and work upstream to calculate
  # cumulative river length and flow wave travel time. Also keeps track of the
  # distance to next downstream city and dam (this is computationally costly).
  # 
  # Input vars:
  # tab - shapefile attribute tatble
  # cTab - connectivity table
  # SLATTRAdir - path to working base directory
  #
  # Output vars:
  # tab - return shapefile attribute table with new outFields attached
   
  # load in other functions:
  source(paste0(SLATTRAdir, '/src/nextDnStrmPOIlength.R'))
  
  # identify which segments are at the bottom of networks:
  mouthInd = which(cTab$DNSTR == 0)
  thisUID = cTab$UID[mouthInd]
  
  # match current ID in cTab with ID in tab:
  tM = match(thisUID, tab$ARCID)
  tab$CITY_UPSTR_TIME_DAY[tM] = tab$DAM_UPSTR_TIME_DAY[tM] = 0
  
  # add current segment length & time to upstream length & time:
  tab$UPSTR_DIST_KM[tM] = tab$MC_LENGTH[tM]
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
  
  # match downstream IDs to current IDs:
  while (length(thisInd) > 0){
    
    thisUID = cTab$UID[thisInd]
    dnstrUID = cTab$UID[lastInd][thisM[thisInd]]
    
    tM = match(thisUID, tab$ARCID)
    dM = match(dnstrUID, tab$ARCID)
    
    # add downstream data to data of current segment length:
    tab$UPSTR_DIST_KM[tM] = tab$UPSTR_DIST_KM[dM] + tab$MC_LENGTH[tM]
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
  #print(proc.time() - ptm)
  
  return(tab)
  
}

MC_tabulator <- function(tab, tabdList, varList, h, i, binInt, maxBin, keep, 
                         tabOutFdir, rivEndPtTabNames, SWOT, nRun){
  
  # Description:
  # For each Monte Carlo ensemble run, tabulates distribution of data for given
  # fields. Concatenates this tabulated data and write out to a CSV file. 
  # 
  # Input vars:
  # tab - input shapefile attribute table
  # tabdList - tabulated distribution table to be filled in
  # varList - list of variables to analyze
  # h - run index
  # i - region index
  # binInt - histogram binning interval
  # maxBin - maximum bin break
  # keep - data filter to subset data
  # tabOutFdir - directory path to write output file
  # rivEndPtTabNames - abbreviated name of region
  # SWOT - boolean indicating whether to run in terms of SWOT path density
  # nRun - number of enseble runs
  #
  # Output vars:
  # tabdList - return tabulated distribution table
  
  # fill in tabulation tables with statistical distributions 
  # used to generate Fig. S3
  tabdNames = names(tabdList)
  
  # weight river lengths depending SWOT overpass frequency:
  if (SWOT == T){
    weight = tab$MC_LENGTH * tab$SWOT_TRAC_DEN
  } else {
    weight = tab$MC_LENGTH
  }
  
  # fill in each distrib. table with histogram counts:
  for (j in 1:length(varList)){
    breaks = seq(0, maxBin+binInt, binInt)
    if (length(grep('_c', tabdNames[j]))>0){ keep = keep & tab$CITY_UPSTR_TIME_DAY>0 }
    if (length(grep('_d', tabdNames[j]))>0){ keep = keep & tab$DAM_UPSTR_TIME_DAY>0 }
    keep = keep & varList[[j]]<=max(breaks)
    x = varList[[j]][keep]
    x = rep(x, round(weight[keep]), each=T)
    tabdList[[j]][h,] = hist(x, breaks, plot=F)$counts
  }
  
  # on last simulation run, write out distribution table(s) to CSV(s):
  if (h==nRun){
    tabdListNames = names(tabdList)
    outTabList = paste0(tabOutFdir, '/distributions/', rivEndPtTabNames[i], '/', tabdListNames, ".csv")
    print(paste("writing out distribution tables:", outTabList))
    for (j in 1:length(tabdList)){
      tabdList[[j]] = cbind(run=1:nrow(tabdList[[j]]), tabdList[[j]])
      write.csv(tabdList[[j]], outTabList[j], row.names=F)
    }
  }
  
  return(tabdList)
  
}

dnStrGaugeCrawler <- function(tab, cTab){

  # Description:
  # For each gauge, find the downstream gauge and add that index to a table. 
  # The first column of this table lists the UID of the segments with gauges. 
  # The next column is the corresponding downstream segment, etc. 
  # 
  # Input vars:
  # tab - shapefile attribute tatble
  # cTab - connectivity table
  #
  # Output vars:
  # segIDtab - return table containg the IDs of downstream segments with gauges.
  #   The first column of this table lists the UID of the segments with gauges. 
  #   The next column is the corresponding downstream segment, etc. 
  
  gCol = grep('GAUGE_Site|GAUGE_STAI', names(tab))
  gBoo = tab[ ,gCol] != 0
  gID = tab[gBoo ,gCol]
  segID = tab$ARCID[gBoo]
  
  # check to make sure that gauges are only assigned to one segment each:
  if (length(unique(gID)) != length(gID)){
    message("length of unique gauge list is not equal to length of non unique 
            gauge list, indicating that single gauges assigned to multiple 
            segments!")
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
  
  # Description:
  # For each pair of USGS river gauges, find the range of possible flow wave
  # travel times. 
  # 
  # Input vars:
  # minCel - minimum realistic celerity
  # minCel - maximum realistic celerity
  # dnStrDist - vector of distances between the pairs of gauges
  #
  # Output vars:
  # lagRange - table with min and max travel times for each pair of gauges
  
  # km to m, and sec to days conversion:
  lagRange = ceiling(cbind(dnStrDist/maxCel, dnStrDist/minCel)*kmpday2mpsConv) 
  lagRange[lagRange > 200] = 200
  
  return(lagRange)
  
}

QreadAndProc = function(qTab, quantile){

  # Description:
  # Find which column contains discharge records in the USGS flow records.
  #
  # Input vars:
  # qTab - USGS gauge discharge table 
  # quantile - optional: only consider discharges > a specified percentile:
  #
  # Output vars:
  # Q - vector of USGS gauge derived discharge 
  
  colNames = names(qTab)
  ind = grep("00060", colNames)
  Qcol = ind[grep("_cd", colNames[ind], inv=T)][1]
  Q = qTab[, Qcol]
  Q = suppressWarnings(as.numeric(levels(Q))[Q])
  highQ = quantile(Q, quantile, na.rm=T)[[1]]
  Q[Q<highQ] = NA
  
  return(Q)
  
}

lagCor = function(q1, q2, lagRange, k){
 
  # Description:
  # Find the time lag that corresponds to the maximum cross correlation between
  # a upstream/downstream pair of gauge flow records:
  #
  # Input vars:
  # q1 - discharge record of gauge 1
  # q2 - discharge record of gauge 1
  # lagRange - table of releastic lags to test the correlations over
  # k - record index of gauge pair
  #
  # Output vars:
  # corVec - vector of cross correlations
  
  lagWin = 1:lagRange[k,2]
  corVec = rep(NA, lagRange[k,2])
  for (l in lagWin){
    q1 = c(q1, NA)
    q2 = c(NA, q2)
    corVec[l] = cor(q1, q2, use="pairwise.complete.obs")
  }
  
  return(corVec)
  
}

lagCorrPlot <- function(tab, corVec, lag_day, gDnStrmDist, modelTT, i, j, k, Q1, 
                        Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID){
  
  # Description:
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Plot Figure S2 in the supplemental material showing the paired gauged 
  # discahrge records and the lag cross correlation between the two. 
  #
  # Input vars:
  # coming soon. maybe. 
  
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
  # convert from m/s to days and km to m
  mtext(paste("Area dif:", round(gAreaDif_per), "%   modCel:",
              round(sum(tab$LENGTH_KM[ID])/(modelTT)*kmpday2mpsConv,2)), line=-1) 
}

################################################################################
# HARD-CODED CONSTANTS
################################################################################
### local file paths:

# inputs:'
SLATTRAdir = '/Users/geoallen/Documents/research/2017_06_16_waveSpeed/git/SLATRRA'


# river network input:
origPolylinesFdir = paste0(SLATTRAdir, '/input/riv_lines') # from Andreadis etal
rivEndPtTabFdir = paste0(SLATTRAdir, '/input/riv_endPts') # made with ArcMap
rivMidPtTabFdir = paste0(SLATTRAdir, '/input/riv_midPts') # made with ArcMap
rivNetworkConnectFdir = paste0(SLATTRAdir, '/input/riv_conn') # made with RRR
# river netowrk
# From Landscan/Natural Earth:
cityFpath = paste0(SLATTRAdir, '/input/city/ne_10m_populated_places.dbf') 
damFpath = paste0(SLATTRAdir, '/input/dam/GRanD_dams_v1_1.dbf') # from GRanD
gaugeFpath = paste0(SLATTRAdir, '/input/gauge/gages_x010g.dbf') # from USGS WIS
gaugeRecordFdir = paste0(SLATTRAdir, '/input/gauge/records') # from USGS WIS

# outputs:
cityDamGaugeFdir = paste0(SLATTRAdir, '/output/riv_POIjoin') # mid-process 
# write out. Shapefiles with joined cities, dams, and gauges
obsOutFdir = paste0(SLATTRAdir, '/output/riv_val')
rivOutFdir = paste0(SLATTRAdir, '/output/riv_out')
figOutFdir = paste0(SLATTRAdir, '/output/figs')
tabOutFdir = paste0(SLATTRAdir, '/output/tabs')

### variables:
# hydroBasins with river networks that drain north of 60N lat. (SRTM data gap):
basin2rm = c(3000001840,3000004740,3000009130,2000028310,
             2000041840,2000043720,8000009560)

# cubic feet per second to cubic peters per second conversion:
cfs2mfsConv = 0.028316847
# km/day to m/s conversion:
kmpday2mpsConv = 0.01157407
# threshold slope (unncessary if using monte carlo simulation):
zeroSlope = 1e-5
# tabulation interval and maximum for output tabs and Fig S3:
celBinInt = 0.2
TTbinInt = 1
maxCel = 50
maxTT = 100
# number of simulation runs:
nRun = 100

# get list of input file paths:
rivEndPtTabFnames = list.files(rivEndPtTabFdir, 'dbf', recursive=T)
rivEndPtTabNames = substr(rivEndPtTabFnames, 1, 2)
rivEndPtTabFpaths = paste0(rivEndPtTabFdir, '/', rivEndPtTabFnames)

rivMidPtTabFnames = list.files(rivMidPtTabFdir, 'dbf', recursive=T)
rivMidPtTabNames = substr(rivMidPtTabFnames, 1, 2)
rivMidPtTabFpaths = paste0(rivMidPtTabFdir, '/', rivMidPtTabFnames)
midPtMatch = match(rivEndPtTabNames, rivMidPtTabNames)

rivNetworkConnectFnames = list.files(rivNetworkConnectFdir)
rivNetworkConnectFpaths = list.files(rivNetworkConnectFdir, 'csv', recursive=T, 
                                     full.names=T)

origPolylinesFpaths = paste0(origPolylinesFdir, '/', rivEndPtTabFnames)

cityDamGaugeFpaths = paste0(cityDamGaugeFdir, '/', rivEndPtTabFnames)
obsOutFpaths = paste0(obsOutFdir, '/', rivEndPtTabFnames)
rivOutFpaths = paste0(rivOutFdir, '/', rivEndPtTabFnames)
celTabOutFpaths = paste0(rivOutFdir, '/', rivEndPtTabNames, 'riv_rangeGT05.csv')
# list gauge files:
gFnames = list.files(gaugeRecordFdir)
gNames = as.numeric(sub('.csv', '', gFnames))
gFpaths = paste0(gaugeRecordFdir, '/', gFnames)

# outpaths for column mean and meadian tabs:
meanTabOutPath = paste0(tabOutFdir, '/run_averages/', rivEndPtTabNames, '_runMeans.csv')
medTabOutPath = paste0(tabOutFdir, '/run_averages/', rivEndPtTabNames, '_runMedians.csv')
aveTabPaths = cbind(meanTabOutPath, medTabOutPath)
# create directories if they don't exist:
dirCreator(dirPath = c(cityDamGaugeFdir, obsOutFdir, rivOutFdir, 
           paste0(tabOutFdir, '/distributions/', rivEndPtTabNames, '/'),
           paste0(tabOutFdir, '/run_averages/'), paste0(figOutFdir, '/validation/')))

################################################################################
# POI JOIN TO RIVER NETWORK
################################################################################

# read in city table:
cityTab = foreign::read.dbf(cityFpath);
cityTab = cityTab[match(c('NAME','LATITUDE','LONGITUDE','POP_MAX','POP_MIN'),
                        names(cityTab))]
cityXY = cbind(cityTab$LONGITUDE, cityTab$LATITUDE)

# read in dam table:
damTab = foreign::read.dbf(damFpath);
damTab = damTab[match(
  c('DAM_NAME','LAT_DD','LONG_DD','YEAR','AREA_SKM','CAP_MCM'), names(damTab))]
damXY = cbind(damTab$LONG_DD, damTab$LAT_DD)

# read in gauge table:
gaugeTab = foreign::read.dbf(gaugeFpath);
# added lat lon fields using WGS84 in arc:
gaugeTab = gaugeTab[match(c('Site_NO','lat','lon', 'HUC8'), names(gaugeTab))]
gaugeXY = cbind(gaugeTab$lon, gaugeTab$lat)

# run through each continent and match up cities, dams, and gauges to the
# river networks:
shpFpaths = sub('.dbf', '.shp', origPolylinesFpaths)
for (i in 1:length(shpFpaths)){ #c(4,6)){ #
  print(i)
  # extract the verticies of the polylines # (very slow - took africa 3 hours
  # to read in):
  #consider trying much faster function, sf::st_combine(st_read(shpFpaths[i])))
  ptm = proc.time()
  pLineTab = read.shp(shpFpaths[i])
  print(proc.time() - ptm)
  ptm = proc.time()
  pLineTab = convert.to.simple(pLineTab)
  print(proc.time() - ptm)

  # extract end point verticies and add a city field to tab:
  attrTab = foreign::read.dbf(origPolylinesFpaths[i])

  # for each city, find the river segment that the city is most likely
  # to be located on:
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
    smallSrad = 1e4, # specify radius of secondary search radius (in meters)
    plot=F) 
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
    smallSrad = 3e3, # specify radius of secondary search radius (in meters)
    plot=F)   
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
    smallSrad = 2e3, # specify radius of secondary search radius (in meters)
    plot=F) 
  print(proc.time() - ptm)

  # write out POI table:
  foreign::write.dbf(attrTab, cityDamGaugeFpaths[i])
}

system("say done run")






################################################################################
#CELERITY & TRAVEL TIME
################################################################################

ptm = proc.time()

# For each continent, calculate flow wave celerity and travel time:
for (i in 1:length(cityDamGaugeFpaths)){
  
  print(paste("Begin simulations in region:", rivEndPtTabNames[i]))
  
  # for sensitivity analysis, run model several times with varying input parameters:
  for (h in 1:nRun){
    
    # if first simulation, read in input data:
    if (h==1){
      print("  Loading input data")
      # read in river polyline attribute table:
      tab_raw = foreign::read.dbf(cityDamGaugeFpaths[i]) 
      
      # read in and process segment endpoint shapefile table:
      EPtab = foreign::read.dbf(rivEndPtTabFpaths[i])
      names(EPtab)[names(EPtab) == "RASTERVALU"] = "ELEV_M"
      names(EPtab)[names(EPtab) == "LENGTH_GEO"] = "LENGTH_KM"
      
      # read in midPtTab:
      midPtTab = foreign::read.dbf(rivMidPtTabFpaths[midPtMatch[i]])
      
      # read in and set up connectivity table:
      j = grep(rivEndPtTabNames[i], rivNetworkConnectFnames)
      cTab = read.csv(rivNetworkConnectFpaths[j], header=F)
      names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
                      "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
                      "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
      
    }
    tab=tab_raw
    
    print(paste("  Simulation run: ", h, "of", nRun))
    
    ####
    # SLOPE CALCULATION:
    
    # check if there are any unpaired ORIG_FIDs:
    oddInd = seq(1, nrow(EPtab), 2)
    evenInd = seq(2, nrow(EPtab), 2)
    unpairedInd = which(EPtab$ORIG_FID[evenInd]-EPtab$ORIG_FID[oddInd] != 0)
    if (length(unpairedInd) > 0){
      message(paste(length(unpairedInd), "unpaired segment endpoints."))
    }
    
    # simulate slope uncertainty through monte carlo error propogation:
    N = nrow(tab)
    MC_LENCOR = runif(n=1, min=1.0, max=1.5) # simulate uncertainty of river length (low res. DEM short circuits river meanders)
    MC_LENGTH = EPtab$LENGTH_KM[evenInd]*MC_LENCOR
    #MC_UPSTR_ELEV = runif(n=N, min=EPtab$ELEV_M[oddInd]-3.5, max=EPtab$ELEV_M[oddInd]+3.5) # include 10 m of uncertainty to elevation
    #MC_DNSTR_ELEV = runif(n=N, min=EPtab$ELEV_M[evenInd]-3.5, max=EPtab$ELEV_M[evenInd]+3.5) # include 10 m of uncertainty to elevation
    MC_ZSLOPE = runif(n=1, min=1e-5, max=1e-4)
    
    # calculate gradient:
    slope = (EPtab$ELEV_M[oddInd]-EPtab$ELEV_M[evenInd])/(1e3*MC_LENGTH) # convert length from km to m
    #MC_SLOPE = (MC_UPSTR_ELEV - MC_DNSTR_ELEV)/(1e3*MC_LENGTH) # convert length from km to m
    
    # check out table 3 -- for minimum slope, look at charlotte's paper (min slope = 1e-5), max slope=1e-2
    MC_smult = runif(n=N, min=0.8, max=1.2) # from http://onlinelibrary.wiley.com/doi/10.1002/joc.2028/full 
    MC_d = rnorm(n=N, mean=1, sd=0.2)
    MC_SLOPE = MC_smult*(slope^MC_d)
    
    # set zero or negative slopes to a minimum threshold:
    MC_SLOPE[MC_SLOPE <= 0 | is.na(MC_SLOPE) | !is.finite(MC_SLOPE)] = MC_ZSLOPE
    
    # add slope to polyline attribute table:
    m = match(tab$ARCID, EPtab$ARCID[oddInd])
    SLOPE = MC_SLOPE[m]
    tab = cbind(tab, SLOPE)
    
    # generate width uncertainty from the 5th, and 95th CI from Andreadis etal 2012:
    mc_widthQuant = pnorm(rnorm(1))
    if (mc_widthQuant < 0.5){
      MC_WIDTH = abs(qnorm(mc_widthQuant, mean=tab$WIDTH, sd=(tab$WIDTH-tab$WIDTH5)/2))
    }else{
      MC_WIDTH = qnorm(mc_widthQuant, mean=tab$WIDTH, sd=(tab$WIDTH95-tab$WIDTH)/2)
    }
    # generate depth uncertainty from the 5th, and 95th CI from Andreadis etal 2012:
    mc_depthQuant = pnorm(rnorm(1))
    if (mc_depthQuant < 0.5){
      MC_DEPTH = abs(qnorm(mc_depthQuant, mean=tab$DEPTH, sd=(tab$DEPTH-tab$DEPTH5)/2))
    }else{
      MC_DEPTH = qnorm(mc_depthQuant, mean=tab$DEPTH, sd=(tab$DEPTH95-tab$DEPTH)/2)
    }
    
    # generate roughness uncertainty from XYZ: 
    MC_N = runif(n=1, min=0.02, max=0.05) 
    
    ####
    # calculate flow wave celerity for each segment:
    # manning rect: β=5/3; Chezy rect: 3/2; manning tri: 4/3; chez tri: 5/4
    tab = celerityCalc(tab, MC_WIDTH, MC_DEPTH, MC_SLOPE, MC_N, B=5/3)
    
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
                UPSTR_DIST_KM, UPSTR_TIME_DAY, CITY_UPSTR_TIME_DAY, DAM_UPSTR_TIME_DAY,
                MC_WIDTH, MC_DEPTH, MC_LENCOR, MC_LENGTH, MC_SLOPE, MC_ZSLOPE, MC_N)
    
    # join data from midPoints tab:
    m = match(tab$ARCID, midPtTab$ARCID) # match UIDs (prob unnecessary)
    tab$hBASIN = midPtTab$MAIN_BAS[m] #  main basin UID for the hydroBASINS dataset (http://www.hydrosheds.org/page/hydrobasins)
    tab$GLCC = midPtTab$GLCC[m] # land cover dataset (https://lta.cr.usgs.gov/glcc/globdoc2_0) #DN: 20 <- desert
    tab$FLOODHAZARD = midPtTab$RASTERVALU[m] # flood hazard composite from the DFO (via NASA): #avhrr
    tab$SWOT_TRAC_DEN = midPtTab$COUNT_coun[m] # swot track density (N overpasses per cycle @ segment centroid)
    tab$CONTINENT = i # add continental UID field (1=af, 2=as, 3=au, 4=ca, 5=eu, 6=na, 7=sa)
      
    # calculate the cumulative length and time along
    # entire river network. (computationally costly):
    tab = cumRivTime_POI(tab, cTab, SLATTRAdir)
  
    # PROCESS OUTPUT DATA:
    # add tab data to polyline shapefile data:
    # replace previously modified file with new copy:
    if (h == 1){
      if (file.exists(origPolylinesFpaths[i])){
        file.remove(rivOutFpaths[i])
        file.copy(origPolylinesFpaths[i], rivOutFpaths[i])
      }
      # match up new data indices to original indices:
      origTab = foreign::read.dbf(rivOutFpaths[i])
      m = match(origTab$ARCID, tab$ARCID)
    }
    # reorder data to polyline vectors:
    tab = tab[m, ]
    
    
    
    

    
    # WEIGHTED AVERAGED SHAPEFILE DBF:
    # for each simulation, take an cumulative mean of each cell in the 
    # shapefile attribute table. This shape file is then used to make Fig. 1 & 3. 
    # take the mean value of all simulations for each segment in the 
    # global flowline network (used to create maps in Fig 1, 3):
    if (h == 1){
      weightMeanTab = tab
    } else {
      weights = c((h-1), 1)/h
      weightMeanTab = Reduce(`+`, Map(`*`, list(weightMeanTab, tab), weights))
    }
    

    # CONVERGENCE PLOT CUMULATIVE TABLES:
    # for each simulation run, take the mean and median of each column
    # and concatenate them to a mean and median table. Can be used for
    # convergence plots. 
    
    # set 0 dam and city travel times to NA to not bias column averaging:
    cumTab = tab
    cumTab$CITY_UPSTR_TIME_DAY[cumTab$CITY_UPSTR_TIME_DAY<=0] = NA
    cumTab$DAM_UPSTR_TIME_DAY[cumTab$DAM_UPSTR_TIME_DAY<=0] = NA
    
    if (h == 1){ 
      meanTab = medTab = tableMaker(names(cumTab), nRun, NA)
    }
    meanTab[h,] = colMeans(cumTab, na.rm=T)
    medTab[h,] = colMedians(as.matrix(cumTab), na.rm=T)
    if (h == nRun){
      meanTab = cbind(run=1:nRun, meanTab, TOT_LENGTH_KM=sum(cumTab$MC_LENGTH))
      medTab = cbind(run=1:nRun, medTab, TOT_LENGTH_KM=sum(cumTab$MC_LENGTH))
    }
    
    
    # HISTOGRAM TABULATED CELERITY AND TRAVEL TIME TABLES:
    # for each simulation tabulate the TT histograms seen in Fig. S3.
    # additionally, tabulate the distribution of celerity:
    
    # create data filters:
    wide = tab$MC_WIDTH > 100
    notDesert = tab$GLCC != 8
    Sof60 = !(tab$hBASIN %in% as.numeric(basin2rm))
    
    # create distribution output tables:
    if (h == 1){
      # celerity tables:
      tabdCel = tabdCel_swot = 
        tableMaker(paste0("cel_", seq(0, maxCel, by=celBinInt)), nRun, 0)
      tabdCelList = list(tabdCel=tabdCel)
      tabdCel_swotList = list(tabdCel_swot=tabdCel_swot)
      # travel time tables
      tabdTT_b = tabdTT_c = tabdTT_d = 
        tabdTT_b_swot = tabdTT_c_swot = tabdTT_d_swot = 
        tableMaker(paste0("day_", seq(0, maxTT, by=TTbinInt)), nRun, 0)
      tabdTTList = list(tabdTT_b=tabdTT_b, tabdTT_c=tabdTT_c, tabdTT_d=tabdTT_d)
      tabdTT_swotList = list(tabdTT_b_swot=tabdTT_b_swot, tabdTT_c_swot=tabdTT_c_swot, 
                             tabdTT_d_swot=tabdTT_d_swot)
    }
    
    # fill in and write out tabulation tables with statistical distributions 
    # used to generate Table 1 and Fig. S3:
    # celerities of all rivers:
    
    tabdCelList = MC_tabulator(tab, tabdList=tabdCelList, varList=list(tab$CELER_MPS), 
              h, i, binInt=celBinInt, maxBin=maxCel, keep=notDesert & Sof60, 
              tabOutFdir, rivEndPtTabNames, SWOT=F, nRun)
    # celerities of SWOT rivers only:
    tabdCel_swotList = MC_tabulator(tab, tabdList=tabdCel_swotList, varList=list(tab$CELER_MPS), 
              h, i, binInt=celBinInt, maxBin=maxCel, keep=notDesert & Sof60 & wide, 
              tabOutFdir, rivEndPtTabNames, SWOT=T, nRun)
    # travel times of all rivers:
    tabdTTList = MC_tabulator(tab, tabdList=tabdTTList, 
              varList=list(tab$UPSTR_TIME_DAY, tab$CITY_UPSTR_TIME_DAY, tab$DAM_UPSTR_TIME_DAY), 
              h, i, binInt=TTbinInt, maxBin=maxTT, keep=notDesert & Sof60, 
              tabOutFdir, rivEndPtTabNames, SWOT=F, nRun)
    # travel times of SWOT rivers only:
    tabdTT_swotList = MC_tabulator(tab, tabdList=tabdTT_swotList, 
              varList=list(tab$UPSTR_TIME_DAY, tab$CITY_UPSTR_TIME_DAY, tab$DAM_UPSTR_TIME_DAY), 
              h, i, binInt=TTbinInt, maxBin=maxTT, keep=notDesert & Sof60 & wide, 
              tabOutFdir, rivEndPtTabNames, SWOT=T, nRun)
  
    print(paste(h, '-', round(medTab$CELER_MPS[h], 3), 'm/s'))
  } # end sensitivity loop

  # write out mean shapefile dbf:
  print("Writing outputs")
  foreign::write.dbf(weightMeanTab, rivOutFpaths[i])
  
  # write out mean and median column tabs for monte carlo convergence plots:
  #print(paste("writing out column averaged table:", meanTabOutPath[i]))
  write.csv(meanTab, meanTabOutPath[i], row.names=F)
  write.csv(medTab, medTabOutPath[i], row.names=F)

} # end region loop

print(proc.time() - ptm)
#system("say travel time calculation done run!")






################################################################################
# VALIDATION - GAUGE CROSS CORRELATION ANALYSIS
################################################################################
# 
# # set up PDF device:
# pdfOut = paste0(figOutFdir, '/validation/crossCorrelations.pdf')
# pdf(pdfOut, width = 7, height=5)
# 
# # set the minimum days of overlap between the two gauge records for
# # them to be considered in the lag correlation analysis:
# minOvrlp = 5*365 # 5 years
# 
# ptm = proc.time()
# # run through each region with gauges (na & ca) and conduct validation:
# for (i in c(4,6)){
#   
#   # read in segment attribute table and add empirical celerity & R columns:
#   tab = foreign::read.dbf(rivOutFpaths[i])
#   nCol1 = ncol(tab)
#   tab = cbind(tab, "OBS_CEL_R"=0, "OBS_CEL_MPS"=0, "OBS_CEL_DIST_KM"=0, 
#               "OBS_CEL_AREADIF_PER"=0, "OBS_CEL_LAG_D"=0, "Q_OVRLP_DAYS"=0, "CORRANGE"=0) 
#   newColInd = (1+nCol1):ncol(tab)
#   
#   # create a data frame containing , , ,
#   # distance 
#   celTab = tableMaker(colNames=c("ID", #segment ID
#                                  'R', # max lag correlation
#                                  "cel_mps", # celerity
#                                  'gDist_km', # distance between two correlated gauges
#                                  'gAreaDif_per', # area difference between correlated gauges
#                                  'lag_day', # lag time of max correlation
#                                  'Qoverlap',# N days of overlapping Q data between 2 gauges 
#                                  'Rrange'),  # range of cross correlations
#                       nrows=0, fill=NA)
#   
#   # read in and set up connectivity table:
#   cityDamGaugeNames = substr(rivEndPtTabFnames, 1, 2)
#   ii = grep(cityDamGaugeNames[i], rivNetworkConnectFnames)
#   cTab = read.csv(rivNetworkConnectFpaths[ii], header=F)
#   names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
#                   "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
#                   "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
#   
#   # for each gauge, find downstream gauges, and their downstream distances and drainage areas: 
#   segIDtab = dnStrGaugeCrawler(tab, cTab)
#   
#   # get the downstream gauges, drainage areas, and cumulative downstream distance:
#   segIDtab = matrix(segIDtab, nrow=nrow(segIDtab))
#   gCol = grep('GAUGE_Site|GAUGE_STAI', names(tab))
#   gTab = matrix(tab[segIDtab, gCol], nrow=nrow(segIDtab))
#   areaTab = matrix(tab$AREA[segIDtab], nrow=nrow(segIDtab))
#   distTab =  matrix(tab$MC_LENGTH[segIDtab], nrow=nrow(segIDtab))
#   cumDistTab = t(apply(distTab, 1, cumsum))
#   
#   # match up gauges linked to river network to file list. First, 
#   # get list of downstream gauge IDs and paths::
#   gFmatch_raw = match(gTab[,1], gNames)
#   notNA = !is.na(gFmatch_raw)
#   gFmatch = gFmatch_raw[notNA]
#   gRows = which(notNA)
#   g1Paths = gFpaths[gFmatch]
#   
#   # set up empty tables for calculated gauge parameters:
#   #lagTimeTab = corTab = modelTTtab = QoverlapTab = matrix(NA, nrow=nrow(gTab), ncol=ncol(gTab))
#     
#   # for each matching gauge, get the cross correlation of each downstream gauge
#   # and assign that to segments between the two gauges:
#   for (j in 1:length(gRows)){
#     print(j)
#     gRow = gRows[j]
#     gVec = gTab[gRow, ]
#     gDnStrmInd = gVec != 0 & !is.na(gVec) & gVec != gVec[1]
#     if (length(which(gDnStrmInd))==0){next} #print('no downstream gauges'); 
#     
#     # get gauge file paths, downstream distances, and upstream areas of 
#     # gauges located downstream:
#     gDnStrmID = gVec[gDnStrmInd]
#     gFmatch_dnStrm_raw = match(gDnStrmID, gNames)
#     gFmatch_dnStrmInd = which(!is.na(gFmatch_dnStrm_raw))
#     if (length(gFmatch_dnStrmInd)==0){next}#print('no downstream gauges on file'); 
#     
#     gFmatch_dnStrm = gFmatch_dnStrm_raw[gFmatch_dnStrmInd]
#     g2paths = gFpaths[gFmatch_dnStrm]
#     gDnStrmDist = cumDistTab[gRow, gDnStrmInd][gFmatch_dnStrmInd]
#     gDnStrmArea = areaTab[gRow, gDnStrmInd][gFmatch_dnStrmInd]
#     
#     lagRange = lagRangeCalc(minCel=0.1, maxCel=10, dnStrDist=gDnStrmDist)
#     
#     # only consider gauges that are close together that do not have a 
#     # large difference in drainage areas: 
#     #keepers = which((gDnStrmDist<1000 & (gDnStrmArea-areaTab[gRow,1])/gDnStrmArea<0.5))
#     #if (length(keepers)==0){next}
#     
#     #### lag correlation analysis:
#     g1 = read.table(g1Paths[j], sep=",", fill=T, header=T)[-1,]
#     # remove low flows from upper gauge discharge data:
#     Q1_orig = QreadAndProc(qTab=g1, quantile=0.0)
#     
#     # skip if gauge doesn't contain more than 5 years of flow data:
#     if (length(which(!is.na(Q1_orig)))<minOvrlp){next}
#     
#     # this could be sped up by reading in everything and put discharge data into a large matrix
#     # so that correlations could be run all at once but this approach is much more complicated 
#     # than individually analyzing correlations one by one:
#     for (k in 1:length(gFmatch_dnStrm)){ #for (k in (1:length(gFmatch_dnStrm))[keepers]){
#       g2 = read.table(g2paths[k], sep=",", fill=T, header=T)[-1,]
#       # match up the dates between the two gauge records:
#       dM = match(g1$datetime, g2$datetime)
#       dInd = which(!is.na(dM))
#       if (length(dInd)<minOvrlp){next}
#       # remove low flows from discharge data:
#       Q1 = Q1_orig
#       Q2 = QreadAndProc(qTab=g2, quantile=0.0)
#       
#       # skip if gauge doesn't contain more than 5 years of flow data:
#       if (length(which(!is.na(Q2)))<minOvrlp){next}
#       
#       # plot hydrograph timeseries over entire records:
#       # xlim = range(c(as.Date(g1$datetime), as.Date(g2$datetime)), na.rm=T)
#       # ylim = range(c(Q1, Q2)*cfs2mfsConv, na.rm=T)
#       # plot(as.Date(g1$datetime), Q1*cfs2mfsConv, type='l',
#       #     xlim=xlim, ylim=ylim, main=paste(j, k), ylab="Flow (cms)", col=1, lwd=0.4)
#       #lines(as.Date(g2$datetime), Q2*cfs2mfsConv, type='l', col=4, lwd=0.4)
#       
#       # match up two discharge records:
#       Q1_match = Q1[dInd]
#       Q2_match = Q2[dM[dInd]]
#       # remove blank discharge measurements:
#       keep = Q1_match!='' & Q2_match!='' #& Q1_match!=0 & Q2_match!=0 
#       dates = as.Date(g1$datetime[dInd], "%Y-%m-%d")
#       Q1 = as.numeric(Q1_match)
#       Q1[!keep] = NA
#       Q2 = as.numeric(Q2_match)
#       Q2[!keep] = NA
#       # skip to next gauge pair if there is not enough overlapping gauge records:
#       Qoverlap = length(which(!is.na(Q1) & !is.na(Q2)))
#       if (Qoverlap<minOvrlp){next}
#       
#       
#       #
#       # use lag correlation to find the optimum lag time:
#       #corVec = lagCor(q1=Q1, q2=Q2, lagRange=lagRange, k=k)
#       corVec_raw = ccf(Q1, Q2, lag.max=200, na.action=na.pass, plot=F)
#       positiveLag = corVec_raw$lag > 0 #lagRange[k,1]
#       corVec = corVec_raw$acf[positiveLag]
#       # Pick the lag with the maximum correlation:
#       lag_day = which.max(corVec)
#       
#       # skip gauge pair if the folling conditions are not met:
#       if (!T %in% !is.na(corVec)){next}
#       if (length(lag_day)==0){next}
#       if (is.na(lag_day)){next}
#       # range and maxmimum of cross correlations must be greater than specified values:
#       corRange = range(corVec)
#       if ((corRange[2]-corRange[1]) < 0.5){next}
#       if (corRange[2] < 0.5){next}
#       # If peak correlation occurs on the first or last day of the lag window, 
#       # remove it from consideration because it the maximum lag might
#       # just be truncated by the window, rendering it meaningless. 
#       if (T %in% c(lag_day == c(1,200))){next}
#       
#       
#       #
#       # assign celerity, correlation, and other parameters to segments between two gauges:
#       # this could be accomplished using less memory by creating a table that does not
#       # include each individual segIDs and only 1 row per celerity result...
#       ID = segIDtab[gRows[j], 1:which(gVec == gDnStrmID[k])]
#       
#       # calculate correlation coefficient and celerity:
#       R = corVec[lag_day]
#       Rrange = corRange[2]-corRange[1]
#       cel_mps = (gDnStrmDist[k] / lag_day)*kmpday2mpsConv # convert from m/s to days and km to m
#       gDist_km = gDnStrmDist[k]
#       gAreaDif_per = 100*(gDnStrmArea[k]-areaTab[gRow,1])/gDnStrmArea[k]
#       celTabMat = cbind(ID, R, cel_mps, gDist_km, gAreaDif_per, lag_day, Qoverlap, Rrange)
#       celTab = rbind(celTab, celTabMat)
#       
#       # fill in tabs:
#       modelTT = sum(tab$MC_LENGTH[ID]/(tab$CELER_MPS[ID]))*kmpday2mpsConv # convert from m/s to days and km to m
#       modelTTtab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = modelTT
#       lagTimeTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = lag_day
#       corTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = corVec[lag_day]
#       QoverlapTab[gRow, gDnStrmInd][gFmatch_dnStrmInd[k]] = Qoverlap
#       
#       ###
#       # Plot hydrograph and lag cross correlation analysis:
#       lagCorrPlot(tab, corVec, lag_day, gDnStrmDist, modelTT,  i, j, k, 
#                   Q1, Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID)
#       
#     }
#   }
#   
#   write.csv(celTab, celTabOutFpaths[i], row.names=F)
#   
#   celTab = read.csv(celTabOutFpaths[i], header=T)
#   
#   print(proc.time() - ptm)
#  
#   # join empirical celerity data onto flowlines, giving preference for higher correlation coeficients:
#   uID = unique(celTab$ID)
#   celTabColInd = 2:ncol(celTab)
#   for (j in 1:length(uID)){
#     uIDmatch = which(celTab$ID == uID[j])
#     tabInd = tab$ARCID == uID[j]
#     # use a correlation-weighted mean to caluclate celerity:
#     tab[tabInd, newColInd] = apply(celTab[uIDmatch, celTabColInd], 2, function(x) weighted.mean(x, celTab[uIDmatch, 2]))
# 
#     # alternative: use the value with the maximum correlation:
#     #celTabInd = uIDmatch[which.max(celTab$R[uIDmatch])]
#     #tab[tabInd, newColInd] = celTab[celTabInd, celTabColInd]
#   }
#   
# 
#   foreign::write.dbf(tab, obsOutFpaths[i])
# 
#   keep = tab$OBS_CEL_MPS>0 # & tab$WIDTH>100# & tab$Q_OVRLP_DAYS>1e3 & tab$OBS_CEL_R<1 & tab$OBS_CEL_MPS > 0  & tab$OBS_CEL_MPS < 5
#   x = tab$OBS_CEL_MPS[keep]
#   x = rep(x, round(tab$LENGTH_KM[keep]), each=T)
#   hist(x, 40, main="hist of celerity", xlab="celerity (m/s)")
# 
#   x = tab$OBS_CEL_R[keep]
#   x = rep(x, round(tab$LENGTH_KM[keep]), each=T)
#   hist(x, 50, main="hist of lag correlations with R > 0.5", xlab="R")
#   
# }
# 
# dev.off()
# cmd = paste('open', pdfOut)
# system(cmd)
 
# PLOT VALIDATION DISTRIBUTIONS (FIGURE 2):
# join observations to model results:



# read in and append each validation tab:
for (i in c(4,6)){
  obsTab = foreign::read.dbf(obsOutFpaths[i])
  modTab = foreign::read.dbf(rivOutFpaths[i])
  if (i == 4){
    vTab = obsTab
    gTab = modTab
  }else{
    vTab = rbind(vTab, obsTab)
    gTab = rbind(gTab, modTab)
  }
  print(paste(i, rivEndPtTabNames[i]))
}
system("say austrolapithicus!")
obsCols = c(grep("OBS", names(vTab))[1]:ncol(vTab))
gTab = cbind(gTab, vTab[,obsCols])



# try basin-by-basin subsampling:
keep = gTab$CELER_MPS<xlim[2] & gTab$OBS_CEL_MP<xlim[2] & gTab$OBS_CEL_R>0.5 & gTab$CORRANGE > 0.5 #& gTab$WIDTH>100

basinUID = unique(gTab$hBASIN[keep])

obsBasinMeanCel = modBasinMeanCel = nSegs = vector()
for (i in 1:length(basinUID)){
  ind = which(gTab$hBASIN[keep] == basinUID[i])
  nSegs[i] = length(ind)
  obsBasinMeanCel[i] = mean(gTab$OBS_CEL_MP[keep][ind])
  modBasinMeanCel[i] = mean(gTab$CELER_MPS[keep][ind])
  
}

par(xpd=F)
plotMax = max(c(obsBasinMeanCel, modBasinMeanCel))
plot(c(0,plotMax), c(0,plotMax),
     main="Basin-by-basin comparison",
     xlab="Gauge-based celerity (m/s)",
     ylab="Modeled celerity (m/s)",
     type='n', asp=1)
text(obsBasinMeanCel, modBasinMeanCel, nSegs, cex=0.7)
abline(0, 1)



quantile(nSegs, (1:100)/100)
length(which(nSegs > 100))




# Plot up the side-by-side histograms:
pdfOut = paste0(figOutFdir, '/validationHistogram.pdf')
pdf(pdfOut, width = 3.4, height=6)
xlim = c(0,10)
ylim = c(0,5e3)
keep = gTab$OBS_CEL_R>0  & gTab$CELER_MPS<xlim[2] & gTab$OBS_CEL_MP<xlim[2] & gTab$OBS_CEL_R>.5 & gTab$CORRANGE > 0.5 #& gTab$WIDTH>100
x1 = gTab$CELER_MPS[keep]
x1 = rep(x1, round(gTab$LENGTH_KM[keep]), each=T)
x2 = gTab$OBS_CEL_MP[keep]
x2 = rep(x2, round(gTab$LENGTH_KM[keep]), each=T)
h1 = hist(x1, seq(0,100,0.2), plot=F)
h2 = hist(x2, seq(0,100,0.2), plot=F)
ylim = range(c(h1$counts, h2$counts))

par(lwd=0.5)
hist(x1, seq(0,xlim[2],.2), main="", #Modeled vs Observed Celerity",
     xlim=xlim,
     ylim=ylim,
     xlab="Celerity (m/s)",
     ylab="River Network Length (km)",
     axes=F,
     col="light blue",#rgb(0,.5,0,.5),
     border=1,
     lwd=0.5)
par(new=T)
hist(x2, seq(0,xlim[2],.2), main='',#main='WIDTH>100',
     xlim=xlim, ylim=ylim,
     xlab="", ylab="",
     axes=F, col=rgb(.5,.5,.5,.5),
     border=1, lwd=0.5)
axis(1, lwd=0.5)
axis(2, lwd=0.5)
#legend("topright", legend=c("Observations", "Model"),
#       text.col=c(rgb(.5,.5,.5,.5),'light blue'),
#       text.font=2, cex=1.4, lwd=0, box.lwd=0)

dev.off()
cmd = paste('open', pdfOut)
system(cmd)


# calcualte validation stats:

# RMSE:
rmse <- function(error){sqrt(mean(error^2))}
# Bias (mean error):
me <- function(error){mean(error)}
# Standard Error:
se <- function(error){ sqrt(mean((error-mean(error))^2)) }

# Calculate the error:
error = gTab$OBS_CEL_MP[keep] - gTab$CELER_MPS[keep]

# 25th, 50th, and 75th quartiles:
quartiles <- function(x){return(quantile(x, c(0.25, 0.5, 0.75)))}
uncertainty <- function(quarts){
  return(paste0(round(quarts[2],1), " +", 
                round(quarts[3]-quarts[2], 1), " -", 
                round(quarts[2]-quarts[1], 1)))
}


print(paste("RMSE:", round(rmse(error), 1)))
print(paste("Bias:", round(me(error), 1)))
print(paste("standard error:", round(se(error), 1)))
print(paste("Model spread:", uncertainty(quartiles(gTab$CELER_MPS[keep]))))
print(paste("Gauge spread:", uncertainty(quartiles(gTab$OBS_CEL_MP[keep]))))









################################################################################
# GRAPHS AND TABLES
################################################################################

## APPEND ALL REGIONS TOGETHER: ####
# read in and append each regional mean shapefile DBF:
for (i in c(1:7)){
  tab = foreign::read.dbf(rivOutFpaths[i])
  if (i == 1){
    gTab = tab
  }else{
    gTab = rbind(gTab, tab)
  }
  print(paste(i, rivEndPtTabNames[i]))
}
system("say austrolapithicus!")




## DISTRIBUTION GRAPHS WITH UNCERTAINTY: ####

# analyze the global distribution of flow wave travel time:
# concatenate each region distrib tab:
tabOrder = c('tabdCel', 'tabdTT_b', 'tabdTT_c', 'tabdTT_d',
             'tabdCel_swot', 'tabdTT_b_swot', 'tabdTT_c_swot', 'tabdTT_d_swot')
paramOrder = c('CELER_MPS', 'UPSTR_TIME', 'CITY_UPSTR', 'DAM_UPSTR_',
               'CELER_MPS', 'UPSTR_TIME', 'CITY_UPSTR', 'DAM_UPSTR_')
xLabOrder = c("Celerity (m/s)", "Travel Time (days)", "Travel Time (days)", "Travel Time (days)",
              "Celerity (m/s)", "Travel Time (days)", "Travel Time (days)", "Travel Time (days)")
keepOrder = c("keep","keep","keep_c","keep_d","keep_swot","keep_swot","keep_c_swot","keep_d_swot")
swotOrder = c(F,F,F,F,T,T,T,T)
paramBinInt = rep(celBinInt, length(paramOrder))
paramBinInt[paramOrder != 'CELER_MPS'] = TTbinInt

wide = gTab$MC_WIDTH > 100
notDesert = gTab$GLCC != 8
Sof60 = !(gTab$hBASIN %in% as.numeric(basin2rm))
keep = notDesert & Sof60
keep_c = keep & gTab$CITY_UPSTR>0
keep_d = keep & gTab$DAM_UPSTR_>0
keep_swot = keep & wide
keep_c_swot = keep_c & wide
keep_d_swot = keep_d & wide

# set up output table (Table 1):
dayVec = c(1:5,10,45)
latTab = as.data.frame(array(NA, c(length(dayVec), length(tabOrder))))
names(latTab) = tabOrder

probVec = c(0.25, 0.5, 0.75)
probTab = as.data.frame(array(NA, c(length(probVec), length(tabOrder))))
names(probTab) = tabOrder

pdfOut = paste0(figOutFdir, '/distributions.pdf')
pdf(pdfOut, width=6, height=10.5)
layout(matrix(1:8, nrow=4, byrow=F))
par(mar=c(5,4,2,5))


# travel times:
for (j in 1:length(tabOrder)){
  # sum each region together:
  for (i in 1:7){
    distribOutPath = paste0(tabOutFdir, '/distributions/', rivEndPtTabNames[i], '/', tabOrder[j], '.csv')
    if (i == 1){
      dTab = read.csv(distribOutPath, header=T)
    } else { dTab = dTab + read.csv(distribOutPath, header=T) }
  }

  # calc RMSE with units of minutes: 
  #rmse = sqrt(mean((colMeans(dTab) - t(dTab))^2))
  
  # divide table by max order of magnitude for pretty y axes:
  oOm = 10^floor(log10(max(dTab)))
  dTab = dTab[,-1]/oOm

  # set plotting x limit:
  celFlag = length(grep('cel', tabOrder[j], ignore.case=T)) > 0
  TTflag = length(grep('TT', tabOrder[j], ignore.case=T)) > 0
  if (celFlag){ xlim = c(0, 10) }
  if (TTflag){ xlim = c(0, 50) }

  # get histogram breaks:
  splitNames = strsplit(names(dTab), '_')
  breaks_raw = as.numeric(sapply(splitNames, '[[', 2))
  plotLimInd = which(breaks_raw<=xlim[2])
  breaks = breaks_raw[plotLimInd]
  lB = c(breaks[-length(breaks)])
  rB = breaks[-1]
  mid = colMeans(rbind(lB, rB))

  # plot median histogram:
  quarts = apply(dTab, 2, quantile, probs=c(0.25, 0.5, 0.75))[,plotLimInd[-length(breaks)]]
  ylim = range(c(0, quarts))
  plot(xlim, ylim, type='n',
       main='', xlab=xLabOrder[j], ylab='Global River Length (km)',
       bty='n', axes=F)
  box(lwd=0.5)
  axis(1, lwd=0.5)
  axis(2, lwd=0.5, las=1)
  mtext(formatC(oOm, format="g"), adj=0, padj=-1, outer=F, cex=0.7)
  polygon(x=rbind(lB, rB, rB, lB, NA),
          y=rbind(quarts[2,], quarts[2,], 0, 0, NA),
          bty="n", lwd=0.7, border=NA, col=rgb(0.7, 0.7, 0.7, 1))

  # add runs:
  # matlines(c(0,rB), t(dTab[,(plotLimInd)]),
  #          axt="n", yaxt ="n", xlab=NA, ylab=NA,
  #          bty="n", lty=1, lwd=0.3, col=rgb(1,0,0,0.05))

  # add 1st and 3rd quartile uncertainty bars:
  zC = T#quarts[1,]>0 & quarts[2,]>0 & quarts[3,]>0
  arrows(mid[zC], quarts[1,zC], mid[zC], quarts[3,zC],
         code=3, length=0.018, angle=90, lwd=0.5)

  # add CDFs:
  cumTab = apply(dTab, 1, cumsum)
  maxVec = apply(cumTab, 2, max)
  normCumTab = rbind(0, t(t(cumTab)/maxVec))
  par(new=T);
  matplot(c(lB, xlim[2]), normCumTab[1:(length(rB)+1),], type='l',
          xaxt="n", yaxt ="n", xlab=NA, ylab=NA,
          lty=1, bty="n", lwd=0.3, col=rgb(0,0,1,.2))
  axis(4, las=1, lwd=0.5, col=4, col.ticks=4, col.axis=4)
  #mtext("Probability", side=4, col=4)
  corns = par("usr"); par(xpd=T)
  text(x=corns[2]+(corns[2]-corns[1])/4, y=mean(corns[3:4]),
       'Probability', srt=270, col=4)
  
  
  # add median point derived from mean shapefile:
  # colInd = grep(paramOrder[j], names(gTab))
  # rowInd = mget(keepOrder[j])[[1]]
  # x = gTab[rowInd, colInd]
  # if (swotOrder[j]){ x = rep(x,  gTab$SWOT_TRAC_[rowInd], each=T) }
  # med = median(x)
  
  # fill in probability table:
  rowSumTab = rowSums(dTab)
  revDtab = dTab[,c(ncol(dTab):1)]
  revCumsum = t(apply(revDtab, 1, cumsum))[,c(ncol(dTab):1)]
  normRevCumsum = revCumsum/rowSumTab
  for (k in 1:length(probVec)){
    probInd = apply(abs(normRevCumsum-probVec[k]), 1, which.min)-1
    probDays = c(1:nrow(normRevCumsum))[probInd]*paramBinInt[j]
    probTab[(length(probVec)+1-k),j] =  mean(probDays)
  }
  
  # add median point derived from histogram:
  medDay = probTab[2,j]
  points(medDay, 0.5, pch=1, cex=1, col=1)
  text(medDay, 0.5, paste0("(",round(medDay,1), ", 0.5)"), pos=4, col=4)

  # fill in Table 1:
  perTab = normRevCumsum[,dayVec+1]
  obTab = (apply(perTab, 2, quantile, probs=c(0.25, 0.5, 0.75)))*100
  latTab[,j]  = paste0(round(obTab[2,],1), " +", round(obTab[3,]-obTab[2,],1),
                       " -", round(obTab[2,]-obTab[1,], 1))
  
  
  
}

# write out probability table:
probTab = rbind(probTab, paste0(round(probTab[2,],1), " +", round(probTab[3,]-probTab[2,],1),
                        " -", round(probTab[2,]-probTab[1,], 1)))
probTab = cbind(probVec, probTab)
write.csv(probTab, paste0(tabOutFdir, '/probabilityTab.csv'), row.names=F)
print(probTab)

# write out table 1:
latTab = cbind(dayVec, latTab[, -grep("cel", names(latTab), ignore.case=T)])
write.csv(latTab, paste0(tabOutFdir, '/tab1_latencies.csv'), row.names=F)
print(latTab)


dev.off()
cmd = paste('open', pdfOut)
system(cmd)







## MONTE CARLO CONVERGENCE PLOTS ####

# put file paths into a table:
aveTabPaths = cbind(meanTabOutPath, medTabOutPath)

# set up variables and lavels to be plotted:
inputColNames = c("MC_WIDTH", "MC_DEPTH", "MC_N", "MC_LENCOR",
                  "MC_SLOPE", "MC_ZSLOPE")
outputColNames = c("CELER_MPS", "UPSTR_TIME_DAY",
                   "CITY_UPSTR_TIME_DAY", "DAM_UPSTR_TIME_DAY")
inputLegNames = c("Width (m)", "Depth (m)", "Manning's n (m/s^1/3)", "Length Multiplier",
                  "Slope (m/m)", "Min. Slope Thresh. (m/m)")
outputLegNames = c("Celerity (m/s)", "Wave to Outlet Travel Time (days)",
                   "Wave to City Travel Time (days)", "Wave to Dam Travel Time (days)")
colNameList = list(inputColNames=inputColNames,
                   outputColNames=outputColNames)
legNameList = list(inputLegNames=inputLegNames,
                   outputLegNames=outputLegNames)
stat = c('Mean', "Median")
varType = c("Inputs", "Outputs")
cols = c("black", "red2", "blue2", "orange3", "cyan4", "darkgray", "dark green")
#plot(1:length(cols), 1:length(cols), col=cols, pch=15, cex=5)

pdfOut = paste0(figOutFdir, '/simConvergence.pdf')
pdf(pdfOut, width=8, height=7)
par(mfrow=c(2,1))
par(mar=c(4,15,2,2))




# read in each mean and median tab and combine by taking average, weighted by
# number of segments in each:
for (h in 1:ncol(aveTabPaths)){
  for (i in 1:nrow(aveTabPaths)){
    if (i==1){
      mTab = read.csv(aveTabPaths[i,h], header=T)
    }else{
      inTab = read.csv(aveTabPaths[i,h], header=T)

      totLength = mTab$TOT_LENGTH_KM
      weights = c(mTab$TOT_LENGTH_KM[1], inTab$TOT_LENGTH_KM[1])/
        (mTab$TOT_LENGTH_KM[1]+inTab$TOT_LENGTH_KM[1])
      mTab = Reduce(`+`, Map(`*`, list(mTab, inTab), weights))
      mTab$TOT_LENGTH_KM = totLength+inTab$TOT_LENGTH_KM

    }
  }
  
  
  # write out global column-averaged table:
  globAveTabPath = sub('af_', '_global_', aveTabPaths[1,h])
  write.csv(mTab, globAveTabPath, row.names=F)
  # mTab = read.csv(globAveTabPath, header=T)

  # for each input and output variable, plot simulation convergence:
  for (i in 1:length(colNameList)){
    colNames = colNameList[[i]]
    legNames = legNameList[[i]]
    plot(c(1,nRun), c(1, -1), type="n",
         main=paste0(letters[i], '. ', varType[i]),
         bty="n", xlab="N Ensembles", ylab="", yaxt="n")
    mtext(paste("Cumulative", stat[h]), side=2,
          line=length(colNames)*2.1+1)
    segments(-1, 0, nRun+1, 0, lty=2)



    for (j in 1:length(colNames)){
      colInd = grep(colNames[j], names(mTab))[1]
      y = mTab[,colInd]
      mY = mean(y[is.finite(y)])
      y[!is.finite(y)] = mY
      cumAve = cumsum(y)/(1:nRun)
      yRange = mY+c(-max(abs(cumAve-mY)), max(abs(cumAve-mY)))

      par(new=T)
      plot(1:nRun, cumAve, type='l', col=cols[j], ylim=yRange,
           bty="n", xlab='', ylab='', yaxt="n", xaxt="n")
      txt = c(yRange[1], mY, yRange[2])
      lab = formatC(signif(txt,digits=3))
      axis(2, at=txt, labels=lab,
           col=cols[j], col.ticks=cols[j], col.axis=cols[j],
           line=(length(colNames)-j)*2.1, cex.axis=1)
    }

    # add legend:
    legend("topright", legend=legNames,
           col=cols, cex=0.8, lwd = 1, box.lwd=NA)
  }

}
dev.off()
cmd = paste('open', pdfOut)
system(cmd)
print(proc.time() - ptm)












## SENSITIVITY ANALYSIS SCATTER PLOTS ####
h = 1
  
  
# read in global column-averaged table:
globAveTabPath = sub('af_', '_global_', aveTabPaths[1,h])
mTab = read.csv(globAveTabPath, header=T)

tabNames = names(mTab)

xFieldNames_raw = c("MC_WIDTH", "MC_DEPTH", "MC_LENCOR", "MC_SLOPE", "MC_ZSLOPE", "MC_N")
yFieldNames_raw = c("CELER_MPS","UPSTR_TIME_DAY")
xFieldNames = rep(xFieldNames_raw, length(yFieldNames_raw))
yFieldNames = rep(yFieldNames_raw, each=length(xFieldNames_raw))

xLabNames_raw = c("Mean simulated width (m)", "Mean simulated depth (m)", 
              "Length correction factor", "Mean simulated slope (m)",
              "Minimum slope threshold", "Simulated Manning's n") 
yLabNames_raw = c("Celerity (m/s)", "Mean wave travel time (days)")
xLabNames = rep(xLabNames_raw, length(yLabNames_raw))
yLabNames = rep(yLabNames_raw, each=length(xLabNames_raw))

xFieldMatch = match(xFieldNames, tabNames)
yFieldMatch = match(yFieldNames, tabNames)

# set up device: 
nPlots = length(xFieldNames)


pdfOut = paste0(figOutFdir, '/sensAnalysisXYplots.pdf')
pdf(pdfOut, width=6.5, height=9)
layout(matrix(1:nPlots, nrow=4, byrow=T))
par(mar=c(4,5,3,1))


# plot the sensitivity of each output to each input:
for (i in 1:nPlots){
  plot(mTab[,xFieldMatch[i]], mTab[,yFieldMatch[i]], 
       xlab = xLabNames[i],
       ylab = yLabNames[i],
       las=1, pch=16)
  mtext(letters[i], adj=0, padj=-1, outer=F, cex=1.1, font=2)
  
  # # add linear regression and statistics:
  # finites = is.finite(mTab[,xFieldMatch[i]]) & is.finite(mTab[,yFieldMatch[i]])
  # mod = lm(mTab[,yFieldMatch[i]][finites] ~ 
  #            mTab[,xFieldMatch[i]][finites])
  # abline(mod)
  # mtext(paste0("y = ", round(mod$coefficients[[2]], 2), "x + ", round(mod$coefficients[[1]], 2)))
  # p = round(summary(mod)$coefficients[2,4], 3)
  # if (p == 0){ p = "p < 0.001" }else{ p = paste("p =", p) }
  # mtext(p, adj=1, padj=-1, outer=F, cex=0.7)
  # r2 = paste("R2 =", round(summary(mod)$r.square, 2))
  # mtext(r2, adj=0, padj=-1, outer=F, cex=0.7)
}


dev.off()
cmd = paste('open', pdfOut)
system(cmd)

