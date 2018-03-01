dnStrGaugeCrawler <- function(tab, cTab){
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
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
  #
  ################################################################################

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