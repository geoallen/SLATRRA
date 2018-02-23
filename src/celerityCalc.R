celerityCalc <- function(tab, width, depth, slope, N, B){
  
  ################################################################################
  # Author:
  # George H. Allen (20180223)
  #
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  #
  # Description:
  # calculates flow wave celerity and travel time of each segment of 
  # flowline dataset based on Manning's Formula and the kinematic wave model
  #
  # Input vars:
  # tab - shapefile attribute tatble
  # width - bankful width (m): --> Andreadis et al
  # depth - bankful depth (m): --> Andreadis et al 
  # slope -  river slope (m/m): --> 15 arcsec HydroSHEDS Conditioned DEM
  # N - Manning's roughness coefficient (s/m^0.33)
  # B - kinematic wave coefficient 
  #
  # Output vars:
  # tab - return shapefile attribute table with new fields of celerity and 
  #  flow wave travel time attached
  ################################################################################
  # flow wave celerity model:
  # segment length (m): --> hydrosheds
  L = tab$LENGTH_KM*1e3 # convert to m
  
  # Calculate the hydraulic radius (m):
  R = width * depth / (2*depth + width)
  
  # manning's equation for flow velocity (m/s):
  v = N^(-1) * R^(2/3) * slope^(1/2)
  
  # Wave celerity equation from Lighthill and Whitham (1955) (m/s):
  celerity = v * B 
  
  # wave propagation time (days):
  propT = L / (celerity*86400) # convert from s to days
  
  # add information to tab:
  tab$CELER_MPS = celerity # (m/s)
  tab$PROPTIME_DAY = propT # (days)
  
  return(tab)
  
}