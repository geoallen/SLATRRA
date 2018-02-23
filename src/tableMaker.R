tableMaker = function(colNames, nrows, fill=NA){
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  # Description:
  # 
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Creates a data frame with specified dimensions
  #
  # Input vars:
  # colNames - vector of column names
  # nrows - value specifying number of rows in output table
  # fill - what to fill the output table with (deafult = NA)
  #
  # Output vars:
  # tab - data frame table 
  ################################################################################
  # create a table:
  tab = as.data.frame(array(fill, c(nrows, length(colNames))))
  names(tab) = colNames
  
  return(tab)
  
}