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

# for a given river segment, returns the upstream segment with the greatest
# drainage area. Takes in connectivity table and info about drainage area:
tableMaker = function(colNames, nrows, fill=NA){
  tab = as.data.frame(array(fill, c(nrows, length(colNames))))
  names(tab) = colNames
  return(tab)
}
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
dirCreater <- function(dirPath){
  for (i in 1:length(dirPath)){
    if (!dir.exists(dirPath[i])){
      dir.create(dirPath[i], recursive=T)
      print(paste("created new directory:", dirPath[i]))
    }
  }
}
largerBasin = function(cTab, mTab, UID, ncolCtab){
  upUIDs = cTab[UID, 4:ncolCtab]
  upUIDs = upUIDs[upUIDs != 0]
  # in the case of two basins being of equal size, the first listed seg is used:
  upUID = upUIDs[which.max(mTab$AREA[upUIDs])]
  return(upUID)
}
POI2RivJoin <- function(POItab, POIXY, pLineTab, attrTab, inAttr, POIname, 
                        scoreFunc, sRad, smallSrad){
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
celerityModel_mann_rect <- function(tab, MC_WIDTH, MC_DEPTH, MC_SLOPE, MC_N, B=5/3){
  # flow wave celerity model:
  # segment length (m): --> hydrosheds
  L = tab$LENGTH_KM*1e3 # convert to m
  # Manning's roughness coefficient (s/m^0.33):
  n = MC_N
  # bankful width (m): --> Andreadis et al
  w = MC_WIDTH
  # bankful depth (m): --> Andreadis et al
  h = MC_DEPTH
  # river slope (m/m): --> hydrosheds
  S = MC_SLOPE # sinuosity correction
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
  #ptm = proc.time()
  
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
  #print(proc.time() - ptm)
  return(tab)
}
tabulator <- function(tab, tabdList, varList, h, i, binInt, maxBin, keep, 
                      tabOutFdir, rivEndPtTabNames, SWOT){
  # fill in tabulation tables with statistical distributions 
  # used to generate Fig. S3
  
  # weight river lengths depending SWOT overpass frequency:
  if (SWOT == T){
    weight = tab$MC_LENGTH *tab$SWOT_TRAC_DEN
  } else {
    weight = tab$MC_LENGTH
  }
  
  # fill in each distrib. table with histogram counts:
  for (j in 1:length(varList)){
    breaks = seq(0, maxBin+binInt, binInt)
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
  # for each gauge, find the downstream gauge and add that index to a table. 
  # The first column of this table lists the UID of the segments with gauges. 
  # the next column is the corresponding downstream segment, etc. 
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
  # km to m, and sec to days conversion:
  lagRange = ceiling(cbind(dnStrDist/maxCel, dnStrDist/minCel)*kmpday2mpsConv) 
  lagRange[lagRange > 200] = 200
  return(lagRange)
}
QreadAndProc = function(qTab, quantile){
  # find which column contains dsicahrge records:
  colNames = names(qTab)
  ind = grep("00060", colNames)
  Qcol = ind[grep("_cd", colNames[ind], inv=T)][1]
  # optional - only consider discharges > a specified percentile:
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
lagCorrPlot <- function(tab, corVec, lag_day, gDnStrmDist, modelTT, i, j, k, Q1, 
                        Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID){
  
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
maxCel = 20
maxTT = 100


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
# create directories if they don't exist:
dirCreater(dirPath = c(cityDamGaugeFdir, obsOutFdir, rivOutFdir, 
           paste0(tabOutFdir, '/distributions/', rivEndPtTabNames, '/'),
           paste0(tabOutFdir, '/run_averages/'), paste0(figOutFdir, '/validation/')))

################################################################################
# POI JOIN TO RIVER NETWORK
################################################################################
# 
# # read in city table:
# cityTab = foreign::read.dbf(cityFpath);  
# cityTab = cityTab[match(c('NAME','LATITUDE','LONGITUDE','POP_MAX','POP_MIN'), 
#                         names(cityTab))]
# cityXY = cbind(cityTab$LONGITUDE, cityTab$LATITUDE)
# 
# # read in dam table:
# damTab = foreign::read.dbf(damFpath);  
# damTab = damTab[match(
#   c('DAM_NAME','LAT_DD','LONG_DD','YEAR','AREA_SKM','CAP_MCM'), names(damTab))]
# damXY = cbind(damTab$LONG_DD, damTab$LAT_DD)
# 
# # read in gauge table:
# gaugeTab = foreign::read.dbf(gaugeFpath);
# # added lat lon fields using WGS84 in arc:
# gaugeTab = gaugeTab[match(c('Site_NO','lat','lon', 'HUC8'), names(gaugeTab))]
# gaugeXY = cbind(gaugeTab$lon, gaugeTab$lat)
# 
# # run through each continent and match up cities, dams, and gauges to the 
# # river networks: 
# shpFpaths = sub('.dbf', '.shp', origPolylinesFpaths)
# for (i in 1:length(shpFpaths)){ #c(4,6)){ #
#   print(i)
#   # extract the verticies of the polylines # (very slow - took africa 3 hours 
#   # to read in): 
#   #consider trying much faster function, sf::st_combine(st_read(shpFpaths[i])))
#   ptm = proc.time()
#   pLineTab = read.shp(shpFpaths[i])
#   print(proc.time() - ptm)
#   ptm = proc.time()
#   pLineTab = convert.to.simple(pLineTab)
#   print(proc.time() - ptm)
#   
#   # extract end point verticies and add a city field to tab:
#   attrTab = foreign::read.dbf(origPolylinesFpaths[i])
#   
#   # for each city, find the river segment that the city is most likely 
#   # to be located on:
#   ptm = proc.time()
#   attrTab = POI2RivJoin(
#     POItab = cityTab, 
#     POIXY = cityXY, 
#     pLineTab = pLineTab, 
#     attrTab = attrTab, 
#     inAttr = c('POP_MAX'), # POI table column name(s) to be transfered to the output attribute table. 
#     # right now, only numerical attributes can be transferred between tables. 
#     POIname = 'CITY', # prefix for what type of POI data this is. Used for output table column names
#     scoreFunc = function(area, dist){ return(area/dist^2) }, # distance to score the relationship between POI and segments
#     sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
#     smallSrad = 1e4) # specify radius of secondary search radius (in meters)
#   print(proc.time() - ptm)
#   
#   # for each dam, find the river segment that the dam is most likely to be located on:
#   ptm = proc.time()
#   attrTab = POI2RivJoin(
#     POItab = damTab, 
#     POIXY = damXY, 
#     pLineTab = pLineTab, 
#     attrTab = attrTab, 
#     inAttr = c('AREA_SKM', 'CAP_MCM'), # POI table column name(s) to be transfered to the output attribute table. 
#     # right now, only numerical attributes can be transferred between tables. 
#     POIname = 'DAM', # prefix for what type of POI data this is. Used for output table column names
#     scoreFunc = function(area, dist){ return(area/dist^3) }, # distance to score the relationship between POI and segments
#     sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
#     smallSrad = 3e3) # specify radius of secondary search radius (in meters)
#   print(proc.time() - ptm)
#   
#   # for each gauge, find the river segment that the gauge is most likely to be located on:
#   ptm = proc.time()
#   attrTab = POI2RivJoin(
#     POItab = gaugeTab, 
#     POIXY = gaugeXY, 
#     pLineTab = pLineTab, 
#     attrTab = attrTab, 
#     inAttr = c('Site_NO', 'HUC8'), # POI table column name(s) to be transfered to the output attribute table. 
#     #inAttr = c('STAID', 'DRAIN_SQKM', 'HUC02'), # POI table column name(s) to be transfered to the output attribute table. 
#     # right now, only numerical attributes can be transferred between tables. 
#     POIname = 'GAUGE', # prefix for what type of POI data this is. Used for output table column names
#     scoreFunc = function(area, dist){ return(area/dist^3) }, # distance to score the relationship between POI and segments
#     sRad = 3e4, # specify the dimensions of the initial square subset (in meters)
#     smallSrad = 2e3) # specify radius of secondary search radius (in meters)
#   print(proc.time() - ptm)
#   
#   # write out POI table:
#   foreign::write.dbf(attrTab, cityDamGaugeFpaths[i])
# }
# 
# system("say done run")
# 
# 
# 
# 
# 
#
################################################################################
#CELERITY & TRAVEL TIME
################################################################################

# Sources of uncertainty: 
# DONE: elevation: gaussian mu=slope, sd=(FIND REF)
# DONE: zero slopes: uniform 1e-5 to 1e3 (or 1e-2)
# DONE: river length: uniform: 1-1.5
# DONE: width: skewed normal (or uniform) 5%, 50%, 95% from Andreadis etal 
# DONE: depth: gaussian 5%, 50%, 95% from Andreadis etal 
# DONE: n: gaussian: 5%, 50%, 95% : 0.02, 0.03, 0.04 (or just assume 0.03)
# DONE: Channel shape -- argue it is a can of worms (add a sentence about how width and depth are so uncertain, beta and shape cancel out)
-- 

ptm = proc.time()
nRun = 3

# For each continent, calculate flow wave celerity and travel time:
for (i in 1:7){ # 1:length(cityDamGaugeFpaths)){ #
  
  print(paste("Begin simulations in region:", rivEndPtTabNames[i]))
  
  # for sensitivity analysis, run model several times with varying input parameters:
  for (h in 1:nRun){
    
    print(paste("Begin simulation run:", h))
    
    # read in river polyline attribute table:
    if (h==1){ tab_raw = foreign::read.dbf(cityDamGaugeFpaths[i]) }
    tab=tab_raw
    
    # read in and process segment endpoint shapefile table:
    EPtab = foreign::read.dbf(rivEndPtTabFpaths[i])
    names(EPtab)[names(EPtab) == "RASTERVALU"] = "ELEV_M"
    names(EPtab)[names(EPtab) == "LENGTH_GEO"] = "LENGTH_KM"
    
    # read in and set up connectivity table:
    j = grep(rivEndPtTabNames[i], rivNetworkConnectFnames)
    cTab = read.csv(rivNetworkConnectFpaths[j], header=F)
    names(cTab) = c("UID", "DNSTR", "N_UPSTRSEGS", 
                    "UPSTR1", "UPSTR2", "UPSTR3", "UPSTR4", "UPSTR5", "UPSTR6", 
                    "UPSTR7", "UPSTR8", "UPSTR9", "UPSTR10", "UPSTR11", "UPSTR12")
    
    
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
    MC_LENCOR = runif(n=N, min=1.1, max=1.5) # simulate uncertainty of river length (low res. DEM short circuits river meanders)
    MC_LENGTH = EPtab$LENGTH_KM[evenInd]*MC_LENCOR
    MC_UPSTR_ELEV = runif(n=N, min=EPtab$ELEV_M[oddInd]-5, max=EPtab$ELEV_M[oddInd]+5) # include 10 m of uncertainty to elevation
    MC_DNSTR_ELEV = runif(n=N, min=EPtab$ELEV_M[evenInd]-5, max=EPtab$ELEV_M[evenInd]+5) # include 10 m of uncertainty to elevation
    MC_ZSLOPE = runif(n=N, min=1e-5, max=1e-4)
    
    # calculate gradient and remove duplicates:
    MC_SLOPE = (MC_UPSTR_ELEV - MC_DNSTR_ELEV)/(1e3*MC_LENGTH) # convert length from km to m
    
    # set zero or negative slopes to a minimum threshold:
    MC_SLOPE[MC_SLOPE <= 0] = MC_ZSLOPE[MC_SLOPE <= 0]
    
    # add slope to polyline attribute table:
    m = match(tab$ARCID, EPtab$ARCID[oddInd])
    SLOPE = MC_SLOPE[m]
    tab = cbind(tab, SLOPE)
    
    
    # generate hydraulic geometry uncertainty from the 5th, and 95th CI from Andreadis etal 2012:
    MC_WIDTH = abs(rnorm(n=N, mean=tab$WIDTH, sd=gm_mean(c(tab$WIDTH-tab$WIDTH5, tab$WIDTH95-tab$WIDTH))/2))
    MC_DEPTH = abs(rnorm(n=N, mean=tab$DEPTH, sd=gm_mean(c(tab$DEPTH-tab$DEPTH5, tab$DEPTH95-tab$DEPTH))/2))
    # generate roughness uncertainty from XYZ: 
    MC_N = runif(n=N, min=0.02, max=0.05) # abs(rnorm(n=N, mean=0.03, sd=0.005)) # 
    
    ####
    # calculate flow wave celerity for each segment:
    # manning rect: β = 5/3; Chezy rect: β = 3/2; manning tri: 4/3; chez tri: 5/4
    tab = celerityModel_mann_rect(tab, MC_WIDTH, MC_DEPTH, MC_SLOPE, MC_N, B=5/3)
    
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
                MC_WIDTH, MC_DEPTH, MC_LENGTH, MC_SLOPE, MC_ZSLOPE, MC_N)
    
    # read in midPtTab:
    midPtTab = foreign::read.dbf(rivMidPtTabFpaths[midPtMatch[i]])
    
    # join data from midPoints tab:
    m = match(tab$ARCID, midPtTab$ARCID) # match UIDs (prob unnecessary)
    tab$hBASIN = midPtTab$MAIN_BAS[m] #  main basin UID for the hydroBASINS dataset (http://www.hydrosheds.org/page/hydrobasins)
    tab$GLCC = midPtTab$GLCC[m] # land cover dataset (https://lta.cr.usgs.gov/glcc/globdoc2_0) #DN: 20 <- desert
    tab$FLOODHAZARD = midPtTab$RASTERVALU[m] # flood hazard composite from the DFO (via NASA): #avhrr
    tab$SWOT_TRAC_DEN = midPtTab$COUNT_coun[m] # swot track density (N overpasses per cycle @ segment centroid)
    tab$CONTINENT = i # add continental UID field (1=af, 2=as, 3=au, 4=ca, 5=eu, 6=na, 7=sa)
      
    # calculate the cumulative length and time along
    # entire river network:
    tab = cumRivTime(tab, cTab)
  
    # WEIGHTED AVERAGED SHAPEFILE DBF:
    
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
    
    # take the mean value of all simulations for each segment in the 
    # global flowline network:
    if (h == 1){
      weightMeanTab = tab
    } else {
      weights = c((h-1), 1)/h
      weightMeanTab = Reduce(`+`, Map(`*`, list(weightMeanTab, tab), weights))
    }
    

    # CUMULATIVE TABLES:
    # for each simulation run, take the mean and median of each column
    # and concatenate them to a mean and median table. Can be used for
    # convergence plots: 
    if (h == 1){ 
      meanTab = medTab = tableMaker(names(tab), nRun, NA)
    }
    meanTab[h,] = colMeans(tab)
    medTab[h,] = colMedians(as.matrix(tab))
    if (h == nRun){
      meanTab = cbind(run=1:nRun, meanTab)
      medTab = cbind(run=1:nRun, medTab, TOT_LENGTH_KM=sum(tab$MC_LENGTH))
    }
    
    # TABULATED CELERITY AND TRAVEL TIME TABLES:
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
        tableMaker(paste0(seq(0, maxCel, by=celBinInt), "_mps"), nRun, 0)
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
    
    tabdCelList = tabulator(tab, tabdList=tabdCelList, varList=list(tab$CELER_MPS), 
              h, i, binInt=celBinInt, maxBin=maxCel, keep=notDesert & Sof60, 
              tabOutFdir, rivEndPtTabNames, SWOT=F)
    # celerities of SWOT rivers only:
    tabdCel_swotList = tabulator(tab, tabdList=tabdCel_swotList, varList=list(tab$CELER_MPS), 
              h, i, binInt=celBinInt, maxBin=maxCel, keep=notDesert & Sof60 & wide, 
              tabOutFdir, rivEndPtTabNames, SWOT=T)
    # travel times of all rivers:
    tabdTTList = tabulator(tab, tabdList=tabdTTList, 
              varList=list(tab$UPSTR_TIME_DAY, tab$CITY_UPSTR_TIME_DAY, tab$DAM_UPSTR_TIME_DAY), 
              h, i, binInt=TTbinInt, maxBin=maxTT, keep=notDesert & Sof60, 
              tabOutFdir, rivEndPtTabNames, SWOT=F)
    # travel times of SWOT rivers only:
    tabdTT_swotList = tabulator(tab, tabdList=tabdTT_swotList, 
              varList=list(tab$UPSTR_TIME_DAY, tab$CITY_UPSTR_TIME_DAY, tab$DAM_UPSTR_TIME_DAY), 
              h, i, binInt=TTbinInt, maxBin=maxTT, keep=notDesert & Sof60 & wide, 
              tabOutFdir, rivEndPtTabNames, SWOT=T)
  
    print(paste(h, '-', round(medTab$CELER_MPS[h], 3), 'm/s'))
  } # end sensitivity loop

  # write out mean shapefile dbf:
  print(paste("Writing out mean shapefile dbf:", rivOutFpaths[i]))
  foreign::write.dbf(weightMeanTab, rivOutFpaths[i])
  
  
  # write out mean and median column tabs:
  print(paste("writing out column averaged table:", meanTabOutPath[i]))
  write.csv(meanTab, meanTabOutPath[i], row.names=F)
  write.csv(medTab, medTabOutPath[i], row.names=F)

} # end region loop

print(proc.time() - ptm)
system("say travel time calculation done run!")

# TEMP:







################################################################################
# VALIDATION - GAUGE CROSS CORRELATION ANALYSIS
################################################################################

# set up PDF device:
pdfOut = paste0(figOutFdir, '/validation/crossCorrelations.pdf')
pdf(pdfOut, width = 7, height=5)

# set the minimum days of overlap between the two gauge records for
# them to be considered in the lag correlation analysis:
minOvrlp = 5*365 # 5 years

ptm = proc.time()
# run through each region with gauges (na & ca) and conduct validation:
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
  distTab =  matrix(tab$MC_LENGTH[segIDtab], nrow=nrow(segIDtab))
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
      
      # calculate correlation coefficient and celerity:
      R = corVec[lag_day]
      Rrange = corRange[2]-corRange[1]
      cel_mps = (gDnStrmDist[k] / lag_day)*kmpday2mpsConv # convert from m/s to days and km to m
      gDist_km = gDnStrmDist[k]
      gAreaDif_per = 100*(gDnStrmArea[k]-areaTab[gRow,1])/gDnStrmArea[k]
      celTabMat = cbind(ID, R, cel_mps, gDist_km, gAreaDif_per, lag_day, Qoverlap, Rrange)
      celTab = rbind(celTab, celTabMat)
      
      # fill in tabs:
      modelTT = sum(tab$MC_LENGTH[ID]/(tab$CELER_MPS[ID]))*kmpday2mpsConv # convert from m/s to days and km to m
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

















################################################################################
# GRAPHS AND TABLES, FIGURE S3
################################################################################



# MONTE CARLO CONVERGENCE PLOT:

aveTabPaths = cbind(meanTabOutPath, medTabOutPath)

paste0(tabOutFdir, '/run_averages/m')


# read in each mean and median tab and combine by taking average, weighted by
# number of segments in each:

for (h in 1:ncol(aveTabPaths)){
  for (i in 1:nrow(aveTabPaths)){
    if (i==1){ 
      mTab = read.csv(aveTabPaths[i,h], header=T)
      
      # TEMP TEMP TEMP
      mTab$TOT_LENGTH_KM = 1e4
      # TEMP TEMP TEMP
      
    }else{ 
      inTab = read.csv(aveTabPaths[i,h], header=T) 
      
      # TEMP TEMP TEMP
      inTab$TOT_LENGTH_KM = 1e6
      # TEMP TEMP TEMP
      
      totLength = mTab$TOT_LENGTH_KM
      weights = c(mTab$TOT_LENGTH_KM[1], inTab$TOT_LENGTH_KM[1])/(mTab$TOT_LENGTH_KM[1]+inTab$TOT_LENGTH_KM[1])
      mTab = Reduce(`+`, Map(`*`, list(mTab, inTab), weights))
      mTab$TOT_LENGTH_KM = totLength+inTab$TOT_LENGTH_KM
      
    }
  }
  
  # write out global column-averaged table:
  globAveTabPath = sub('af_', '_global_', aveTabPaths[1,h])
  write.csv(mTab, globAveTabPath, row.names=T)
  
  
  # plot simulation convergence:
  plot(c(1,nRun), c(1, -1), type="n",
       xlab="Simulation Runs",
       ylab="Cumulative Averages",
       main="Simulation convergence")
  segments(1, 0, nRun, 0)
  colNames = c("MC_WIDTH", "MC_DEPTH", "MC_LENGTH", "MC_SLOPE", "MC_ZSLOPE", "MC_N")
  colNames = c("CELER_MPS", "UPSTR_TIME_DAY", "CITY_UPSTR_TIME_DAY", "DAM_UPSTR_TIME_DAY")
  for (j in 1:length(colNames)){
    colInd = grep(colNames[j], names(mTab))[1]
    y = mTab[,colInd]
    
    cumAve = cumsum(y)/1:nRun
    #normCumAve = (cumAve-mean(y))/max(abs(cumAve-mean(y)))
    #lines(1:nRun, normCumAve, col=j)
    par(new=T)
    plot(1:nRun, cumAve, type='l', col=j)
  }
}
  









# read in and append each regional mean shapefile DBF:
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
uprQuant = 0.99
wide = gTab$WIDTH > 100
realS = gTab$SLOPE != zeroSlope
notDesert = gTab$GLCC != 8
Sof60 = !(gTab$hBASIN %in% as.numeric(basin2rm))
keep = notDesert & Sof60 & realS & wide
breaks = seq(0, 50, 1)
latencies = c(1,2,5,10,45)





# analyze the global distribution of flow wave travel time:
# concatenate each region distrib tab:
tabOrder = c('tabdCel', 'tabdTT_b', 'tabdTT_c', 'tabdTT_d',
             'tabdCel_swot', 'tabdTT_b_swot', 'tabdTT_c_swot', 'tabdTT_d_swot')

# travel times:
for (j in 1:length(tabOrder)){
  # concatenate each xyz
  for (i in 3:4){
    distribOutPath = paste0(tabOutFdir, '/distributions/', rivEndPtTabNames[i], '/', tabOrder[j], '.csv')
    
    if (i == 3){ 
      dTab = read.csv(distribOutPath, header=T)
    }else{
      dTab = rbind(dTab, read.csv(distribOutPath, header=T))
    }
  }
  

  # set plotting x limit:
  TTflag = grep('TT', tabOrder[j], ignore.case=T)
  celFlag = grep('cel', tabOrder[j], ignore.case=T)
  if (length(TTflag)>0){ xlim = c(0, 50) }
  if (length(celFlag)>0){ xlim = c(0, 10) }
  
  # get histogram breaks:
  splitNames = strsplit(names(dTab), '_')
  breaks = as.numeric(sapply(splitNames[-1], '[[', 2))
  plotLimInd = which(breaks<=xlim[2])
  nPlim = length(plotLimInd)
  breaks = breaks[plotLimInd]
  
  lB = breaks[-length(breaks)]
  rB = breaks[-1]
  mid = colMeans(rbind(lB, rB))
  mn = colMeans(dTab[-1])[plotLimInd[-nPlim]]
  std = apply(dTab[-1], 2, sd)[plotLimInd[-nPlim]]
  ylim = range(0, mn+std)
  
  # plot histogram:
  plot(xlim, ylim, type='n')
  polygon(x=rbind(lB, rB, rB, lB, NA), 
          y=rbind(mn, mn, 0, 0, NA),
          border=NA, col=rgb(.7, .7, .7))
  
  # add uncertainty bars:
  options(warn=-1)
  arrows(mid, mn, mid, mn+std, length=0.03, angle=90, lwd=0.5)
  arrows(mid, mn, mid, mn-std, length=0.03, angle=90, lwd=0.5)
  options(warn=0)
  
  # add CDF:
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
  
  
  
}




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

# FIXME: multiply length by N overpasses:
keep = notDesert & Sof60 & wide

x = gTab$UPSTR_TIME[keep]
x = rep(x, gTab$LENGTH_KM[keep])
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
pdfOut = paste0(figOutFdir, "distributions_woSWOTwSWOTRivers_minSlope", zeroSlope, ".pdf")

pdf(pdfOut, width=7, height=9)
layout(matrix(1:6, nrow=3, byrow=F))

# plot distributions of travel time over the entire dataset:
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
latencyTab = distPlot(x=gTab$CITY_UPSTR[keep & nonZeroTravelTime], 
         weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream City (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream city"

# c. DAM TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$DAM_UPSTR_>0
latencyTab = distPlot(x=gTab$DAM_UPSTR_[keep & nonZeroTravelTime], 
         weight=gTab$LENGTH_KM[keep & nonZeroTravelTime], 
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
latencyTab = distPlot(x=gTab$UPSTR_TIME[keep], 
         weight=gTab$LENGTH_KM[keep]*gTab$SWOT_TRAC_[keep], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Basin Outlet(Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Observable global river network length"

# e. CITY TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$CITY_UPSTR>0
latencyTab = distPlot(x=gTab$CITY_UPSTR[keep & nonZeroTravelTime], 
         weight=gTab$LENGTH_KM[keep & nonZeroTravelTime]*gTab$SWOT_TRAC_[keep & nonZeroTravelTime], 
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream City (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream city"

# f. DAM TRAVEL TIME:
j = j+1
nonZeroTravelTime = gTab$DAM_UPSTR_>0
latencyTab = distPlot(x=gTab$DAM_UPSTR_[keep & nonZeroTravelTime], 
         weight=gTab$LENGTH_KM[keep & nonZeroTravelTime]*gTab$SWOT_TRAC_[keep & nonZeroTravelTime],
         uprQuant=uprQuant, breaks=breaks, j,
         xLab="Travel Time to Next Downstream Dam (Days)",
         yLab="Global River Length (km)",
         latencyTab = latencyTab, latencies=latencies)
names(latencyTab)[j] = "Next downstream dam"

# write out second part of table 1 for swot-observable rivers:
#write.csv(latencyTab, sub('.pdf', '.csv', pdfOut), row.names=F)
#names(latencyTab) = c("latency", "basins", "cities", "dams"); print(latencyTab)

dev.off()
cmd = paste('open', pdfOut)
system(cmd)

