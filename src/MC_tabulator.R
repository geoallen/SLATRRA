MC_tabulator <- function(tab, tabdList, varList, h, i, binInt, maxBin, keep, 
                      tabOutFdir, rivEndPtTabNames, SWOT, nRun){
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
  # Description:
  # For each Monte Carlo ensemble run, tabulates distribution of data for given
  # fields. Concatenates this tabulated data and write out to a CSV file. 
  # 
  # Input vars:
  # tab - input shapefile attribute table
  # tabdList - tabulated distribution table to be filled in
  # varList - list of variables to analyze
  # h - run index
  # i - region index
  # binInt - histogram binning interval
  # maxBin - maximum bin break
  # keep - data filter to subset data
  # tabOutFdir - directory path to write output file
  # rivEndPtTabNames - abbreviated name of region
  # SWOT - boolean indicating whether to run in terms of SWOT path density
  # nRun - number of enseble runs
  #
  # Output vars:
  # tabdList - return tabulated distribution table
  #
  ################################################################################
  
  # fill in tabulation tables with statistical distributions 
  # used to generate Fig. S3
  tabdNames = names(tabdList)
  
  # weight river lengths depending SWOT overpass frequency:
  if (SWOT == T){
    weight = tab$MC_LENGTH * tab$SWOT_TRAC_DEN
  } else {
    weight = tab$MC_LENGTH
  }
  
  # fill in each distrib. table with histogram counts:
  for (j in 1:length(varList)){
    breaks = seq(0, maxBin+binInt, binInt)
    if (length(grep('_c', tabdNames[j]))>0){ keep = keep & tab$CITY_UPSTR_TIME_DAY>0 }
    if (length(grep('_d', tabdNames[j]))>0){ keep = keep & tab$DAM_UPSTR_TIME_DAY>0 }
    keep = keep & varList[[j]]<=max(breaks)
    x = varList[[j]][keep]
    x = rep(x, round(weight[keep]), each=T)
    tabdList[[j]][h,] = hist(x, breaks, plot=F)$counts
  }
  
  # on last simulation run, write out distribution table(s) to CSV(s):
  if (h==nRun){
    tabdListNames = names(tabdList)
    outTabList = paste0(tabOutFdir, '/distributions/', rivEndPtTabNames[i], '/', tabdListNames, ".csv")
    print(paste("writing out distribution tables:", outTabList))
    for (j in 1:length(tabdList)){
      tabdList[[j]] = cbind(run=1:nrow(tabdList[[j]]), tabdList[[j]])
      write.csv(tabdList[[j]], outTabList[j], row.names=F)
    }
  }
  
  return(tabdList)
  
}