dirCreator <- function(dirPath){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  # Description:
  # 
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Create local directory(ies) in specified path locations
  #
  # Input vars:
  # dirPath - a single file path or a vector of paths
  #
  # Output vars:
  # none
  ################################################################################
  
  for (i in 1:length(dirPath)){
    if (!dir.exists(dirPath[i])){
      dir.create(dirPath[i], recursive=T)
      print(paste("created new directory:", dirPath[i]))
    }
  }
  
}