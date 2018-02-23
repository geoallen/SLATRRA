gmMean <- function(x, na.rm=TRUE){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  # Description:
  # 
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Computes the geometric mean
  #
  # Input vars:
  # x - a vector of values
  #
  # Output vars:
  # geomMean - the geometric mean
  ################################################################################
  
  geomMean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  
  return(geomMean)
  
}