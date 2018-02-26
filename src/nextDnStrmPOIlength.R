nextDnStrmPOIlength <- function(tab, cTab, mouthInd, POIbool, outField){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
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
  #
  ################################################################################
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