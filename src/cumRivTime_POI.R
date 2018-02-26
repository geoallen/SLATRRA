cumRivTime_POI <- function(tab, cTab, SLATTRAdir){
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
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
  #
  ################################################################################
  
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
