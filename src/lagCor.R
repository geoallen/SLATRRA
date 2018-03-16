lagCor = function(q1, q2, lagRange, k){
  ################################################################################
  # Author:
  # George H. Allen (20180316)
  #
  # Description:
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
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
  #
  ################################################################################
  
  lagWin = 1:lagRange[k,2]
  corVec = rep(NA, lagRange[k,2])
  for (l in lagWin){
    q1 = c(q1, NA)
    q2 = c(NA, q2)
    corVec[l] = cor(q1, q2, use="pairwise.complete.obs")
  }
  
  return(corVec)
  
}