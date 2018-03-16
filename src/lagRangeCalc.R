lagRangeCalc = function(minCel, maxCel, dnStrDist, kmpday2mpsConv){
  
  ################################################################################
  # Author:
  # George H. Allen (20180316)
  #
  # Description:
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # For each pair of USGS river gauges, find the range of possible flow wave
  # travel times. 
  #
  # Input vars:
  # minCel - minimum realistic celerity
  # minCel - maximum realistic celerity
  # dnStrDist - vector of distances between the pairs of gauges
  # kmpday2mpsConv - # km to m, and sec to days conversion
  #
  # Output vars:
  # lagRange - table with min and max travel times for each pair of gauges
  #
  ################################################################################
  
  lagRange = ceiling(cbind(dnStrDist/maxCel, dnStrDist/minCel)*kmpday2mpsConv) 
  lagRange[lagRange > 200] = 200
  
  return(lagRange)
  
}