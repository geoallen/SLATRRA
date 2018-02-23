POI2RivJoin <- function(POItab, POIXY, pLineTab, attrTab, inAttr, POIname, 
                        scoreFunc, sRad, smallSrad, plot=F){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  # Description:
  # 
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # For each point of interest (POI), find the river segment that the 
  # POI is most likely to be located on. Point of interests include cities, 
  # dams and gauges. 
  #
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
  
  ################################################################################
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