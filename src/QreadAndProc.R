QreadAndProc = function(qTab, quantile){
  ################################################################################
  # Author:
  # George H. Allen (20180316)
  #
  # Description:
  # Find which column contains discharge records in the USGS flow records.
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Input vars:
  # qTab - USGS gauge discharge table 
  # quantile - optional: only consider discharges > a specified percentile:
  #
  # Output vars:
  # Q - vector of USGS gauge derived discharge 
  #
  ################################################################################
  
  colNames = names(qTab)
  ind = grep("00060", colNames)
  Qcol = ind[grep("_cd", colNames[ind], inv=T)][1]
  Q = qTab[, Qcol]
  Q = suppressWarnings(as.numeric(levels(Q))[Q])
  highQ = quantile(Q, quantile, na.rm=T)[[1]]
  Q[Q<highQ] = NA
  
  return(Q)
  
}