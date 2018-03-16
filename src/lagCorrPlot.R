lagCorrPlot <- function(tab, corVec, lag_day, gDnStrmDist, modelTT, i, j, k, Q1, 
                        Q2, cfs2mfsConv, dates, R, Qoverlap, gAreaDif_per, ID){
  
  ################################################################################
  # Author:
  # George H. Allen (20180316)
  #
  # Description:
  # This function is used in the analysis of Allen et al. "Global estimates of  
  # river flow wave travel times and implications for low-latency satellite data"
  # 
  # Plot Figure S2 in the supplemental material showing the paired gauged 
  # discahrge records and the lag cross correlation between the two. 
  #
  # Input vars:
  # coming soon. maybe. 
  ################################################################################
  
  par(mfrow=c(2,1))
  par(mar=c(4.1,5.1,2.1,7.1))
  # plot lag-correlation:
  plot(corVec, type='l',
       ylim = c(0, 1),
       xlab = "Lag (days)",
       ylab = "Correlation, R")
  points(lag_day, corVec[lag_day], cex=0.8)
  abline(v=lagRangeCalc(0.5, 2.5, gDnStrmDist[k]), lty=2, col='gray')
  abline(v=modelTT)
  mtext(paste("i:", i,"   j:", j, "    k:", k))
  
  # plot hydrograph time series and lag correlations
  Q1[Q1==0] = 1
  Q2[Q2==0] = 1
  ylim = range(c(Q1, Q2, na.rm=T)*cfs2mfsConv, na.rm=T)
  plot(dates, Q1*cfs2mfsConv,
       ylim=ylim,
       xlab="Year",
       ylab="Flow (cms)",
       #main=paste(i, j, k)
       #log='y',
       type='l', col=1, lwd=0.4, las=1)
  # par(new=T)
  # plot(c(rep(NA,lag_day),dates), c(rep(NA,lag_day), Q2)*cfs2mfsConv,
  #      ylim=ylim,
  #      yaxt = "n",
  #      ylab = NA, xlab = NA,
  #      #log='y',
  #      type='l', col='green', lwd=0.3,
  #      axes=F)
  par(new=T)
  plot(dates, Q2*cfs2mfsConv,
       ylim=ylim,
       yaxt = "n",
       ylab = NA, xlab = NA,
       #log='y',
       type='l', col='blue', lwd=0.3,
       axes=F)
  axis(4, las=1, col=4, col.ticks=4, col.axis=4)
  corns = par("usr"); par(xpd=T)
  text(x=corns[2]+40, y=mean(corns[3:4]),
       'Flow (cms)',
       srt=270,
       col=4)
  
  mtext(paste("dist km:", round(gDnStrmDist[k]),
              "   cel:", round(cel_mps, 2), "m/s",
              "   R:", round(R, 3),
              "   Q ovrlp:", Qoverlap))
  # convert from m/s to days and km to m
  mtext(paste("Area dif:", round(gAreaDif_per), "%   modCel:",
              round(sum(tab$LENGTH_KM[ID])/(modelTT)*kmpday2mpsConv,2)), line=-1) 
}