largerBasin <- function(cTab, mTab, UID, ncolCtab){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
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
  ################################################################################
  
  upUIDs = cTab[UID, 4:ncolCtab]
  upUIDs = upUIDs[upUIDs != 0]
  
  # in the case of two basins being of equal size, the first listed seg is used:
  upUID = upUIDs[which.max(mTab$AREA[upUIDs])]
  
  return(upUID)
  
}