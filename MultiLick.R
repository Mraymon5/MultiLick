#Copyright 2018, Martin A. Raymond

#This file is part of LickFunction.

#LickFunction is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#LickFunction is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with LickFunction.  If not, see <http://www.gnu.org/licenses/>.

#' Conducts microstructural analysis of .wav lickometer files exported from Audacity
#' @param Folder String: filepath of .csv files(s) to be analyzed
#' @param KeyWord String: A search term identifying desired files
#' @param File String: The filepath of one file to be analyzed
#' @param Cutoff Numeric: minimum acceptable ILI in mS
#' @param ThresholdSD_L Numeric: Number of Standard Deviations from the mean to set the lick detection threshold for the left track
#' @param ThresholdSD_R Numeric: Number of Standard Deviations from the mean to set the lick detection threshold for the right track
#' @param Start Time to initiate analysis, in minutes
#' @param Stop Time to terminate analysis, in minutes
#' @param Output String: filepath of a .csv file created containing microstructure findings
#' @param ILIFile String: filepath of a .csv file created containging all ILIs from files analyzed
#' @param Title String: Title appended to data saved in the global environment
#' @param LeftLab String: Title of the left recording channel
#' @param RightLab String: Title of the Right recording channel
#' @param PlotSignal "R" plots right channel waveform, "L" plots left channel waveform, "B" plots both.
#' @param DSF Downsampling factor for PlotSignal. Warning: can hide threshold-crossing events.
#' @param PlotRawInterval "R" plots unfiltered ILI histogram for right channel, "L" plots unfiltered ILI histogram for left channel, "B" plots both.
#' @param PlotInterval "R" plots filtered ILI histogram for right channel, "L" plots filtered ILI histogram for left channel, "B" plots both.
#' @param PlotCumulative "R" plots cumulative licks over time for right channel, "L" plots cumulative licks over time for left channel, "B" plots both.
#' @param PlotZoom Plots a magnified selection of the raw signal, "R" for right channel, "L" for left channel, "B" for both.
#' @param ZoomMin Sets the lower bound of PlotZoom, in seconds.
#' @param ZoomMax Sets the upper bound of PlotZoom, in seconds.
#' @param Timer Reports program runtime.

require(tuneR)

MultiLick <- function(Folder,  #sets the folder containing the data
                    KeyWord = "", #sets a search term for identifying files
                    File = F, #sets a specific file to assess
                    Cutoff = 40, #minimum ILI, in milliseconds
                    ThresholdSD_L = 3, #threshold sds from the mean
                    ThresholdSD_R = 3, #threshold sds from the mean
                    Start = 0, #Time to initiate analysis of recording, in minutes
                    Stop = NULL, #Time to terminate analysis of recording, in minutes
                    Output = F, #sets a .csv file to save output to
                    ILIFile = F, #sets a .csv to save ILIs to
                    Title = "Lickout", #sets a title for data output
                    LeftLab = "Left", #title of left track
                    RightLab = "Right", #title of right track
                    PlotSignal = F, #choose tracks to plot
                    DSF = 2, #Downsampling factor for plotsignal. Warning: can hide some threshold-crossing events.
                    PlotRawInterval = F, #plot raw interval histogram
                    PlotInterval = F, #plot filtered interval hist.
                    PlotCumulative = F, #plot licks over time
                    PlotZoom = F, #plot a magnified section
                    ZoomMin = 0, #lower boundary of zoom, in seconds
                    ZoomMax = 0, #upper bound of zoom, in seconds
                    Timer = F ) { #reports runtime
  if (Timer != F) {
    ptm <- proc.time()
  } # This starts a timer
  
  FileNames <- list.files(path = Folder, pattern = KeyWord, full.names = T)
  Directory <- data.frame(list.files(path = Folder, pattern = KeyWord, full.names = F))
  colnames(Directory) <- "Files" # This constructs a directory of files assessed
  
  for (i in 1:length(FileNames)) {
    
    if (File != F) {
      D <- File
    } else {
      D <- FileNames[i]
    }
    
    # Data Import ==== 
    WavData <- readWave(filename = D, units = "seconds", from = Start*60, to = Stop*60) # Import raw signal

    if (WavData@stereo == T){RawData <- data.frame(Left = WavData@left, Right = WavData@right)}
    else {RawData <- data.frame(Left = WavData@left)}
      
    # Time Stamp ==== 
    SampleRate <- 1000 / WavData@samp.rate #Saves the number of milliseconds per sample
    RawData$mS <- c(1:nrow(RawData)) * SampleRate # create timestamp in milliseconds
    
    # Set Lick Threshold, Left ==== 
    L_Threshold <- (mean(RawData$Left) - (ThresholdSD_L * sd(RawData$Left))) # sets detection threshold 
    
    # Test Threshold, Left ==== 
    RawData$LF1 <- 0
    RawData$LF1[RawData$Left < L_Threshold] <- 1 #Creates a record of all threshold crossings
    
    
    # Begin Lick Coding, Left ==== 
    N <- 0
    RawData$LF2 <- c(0,RawData$LF1[1:(nrow(RawData)-1)])
    RawData$LF3 <- RawData$mS
    RawData$LF3[RawData$LF1 - RawData$LF2 != 1] <- 0
    if (sum(RawData$LF3) > 0){ #This chunk handles files with no licks
      LicksL <- data.frame(RawTime = RawData$LF3[RawData$LF3 != 0])
    } else {
      LicksL <- data.frame(RawTime = 0)
    }
    LicksL[, 2:7] <- 0
    colnames(LicksL)[2:7] <- c("FilteredTime", "BurstCount", "RawInterval", "FilteredInterval", "DoubleContact", "LickTotal")
    LicksL <- rbind(0, LicksL)
    for (X in 2:nrow(LicksL)) { 
      if ((LicksL$RawTime[X]-LicksL[X-1,2])>Cutoff) {
        LicksL[X,2] <- LicksL[X,1]
      } else {
        LicksL[X,2] <- LicksL[X-1,2]
      }}
    
    LL <- 0
    for (X in 2:nrow(LicksL)) {
      if ((LicksL[X, 2] - LicksL[X - 1, 2]) > Cutoff) {
        LL <- LL + 1
      }
      LicksL[X, 7] <- LL
    }
    
    # Burst Analysis, Left ==== 
    DC <- 0
    BurstL <- data.frame(Licks = 0, Duration = 0, Start = 0, Stop = 0)
    for(X in 2:nrow(LicksL)) {
      if ((LicksL[X, 2] - LicksL[X - 1, 2]) > Cutoff & (LicksL[X, 2] - LicksL[X - 1, 2]) < 1000) {
        LicksL[X, 3] <- (LicksL[X - 1, 3] + 1)
        LicksL[X - 1, 6] <- DC
      } else {
        if ((LicksL[X, 2] - LicksL[X - 1, 2]) < Cutoff) {
          LicksL[X, 3] <- LicksL[X - 1, 3]
          DC <- DC + 1; LicksL[X - 1, 6] <- DC
        } else {
          LicksL[X, 3] <- 1
          LicksL[X - 1, 6] <- DC
          DC <- 0
        }}}
    N <- 0
    for(X in 2:nrow(LicksL)) {
      if ((LicksL[X, 3] - LicksL[X - 1, 3]) < (-1)) {
        N <- N + 1
        BurstL[N, 1] <- LicksL[X - 1, 3]
      }}
    
    if (LicksL[nrow(LicksL), 3] > 2) {
      N <- N + 1
      BurstL[N, 1] <- LicksL[nrow(LicksL), 3]}#For bursts that end on the last lick
    
    LB <- N
    N <- 0
    for(X in 2:nrow(LicksL)){
      if ((LicksL[X, 3] - LicksL[X - 1, 3]) < (-1)) {
        N <- N + 1
        BurstL[N, 2] <- LicksL[X - 1, 2] - LicksL[X - (LicksL[X - 1, 3] + LicksL[X - 1, 6]), 2]
        BurstL[N, 4] <- LicksL[X - 1, 1]
        BurstL[N, 3] <- LicksL[(X - (LicksL[X - 1, 3] + LicksL[X - 1, 6])) ,1]
      }}
    
    if (LicksL[nrow(LicksL), 3] > 2) {
      N <- N + 1
      BurstL[N, 2] <- LicksL[X - 1, 2] - LicksL[X - (LicksL[X - 1, 3] + LicksL[X - 1, 6]), 2]
      BurstL[N, 3] <- LicksL[(nrow(LicksL)+1)-(LicksL[nrow(LicksL), 3] + LicksL[nrow(LicksL)-1, 6]), 1]
      BurstL[N, 4] <- LicksL[X , 1]
    }#For bursts that end on the last lick
    
    PauseL <- data.frame(c(BurstL$Start, (nrow(RawData) * SampleRate)), c(0, BurstL$Stop))
    colnames(PauseL)[1] <- "Start"
    colnames(PauseL)[2] <- "Stop"
    PauseL$Pause <- PauseL$Start - PauseL$Stop
    
    # Interval Analyses, Left ==== 
    LicksL[1, 4] <- 0
    for(X in 2:nrow(LicksL)) {
      LicksL[X, 4] <- (LicksL[X, 1] - LicksL[X - 1, 1])
    }
    LicksL[1, 5] <- 0
    for(X in 2:nrow(LicksL)) {
      LicksL[X, 5] <- (LicksL[X, 2] - LicksL[X - 1, 2])
    }
    
    # Plotting, Left ==== 
    RawData$LF4 <- 0
    RawData$LF4[(LicksL$RawTime[LicksL$FilteredInterval != 0]) / SampleRate] <- 1
    
    if (PlotSignal == "L"||PlotSignal == "B") {
      PlotMat <- data.frame(x = RawData$mS[(1:floor(nrow(RawData)/DSF))*DSF]/1000, y = RawData$Left[(1:floor(nrow(RawData)/DSF))*DSF], thresh = L_Threshold)
      matplot(x = PlotMat$x, y = PlotMat[,2:3], main = LeftLab, xlab = "Seconds", ylab = "Lin.", type = "l", lty = 1, col = c(1,2),xaxt = "n")
      axis(1, at = seq(0, ceiling(x = nrow(RawData)/100000)*100, by = 100), las = 2) #sets X axis to 100-sec increments
      if (sum(RawData$LF3) > 0) {
      PointPlot <- data.frame(x = LicksL$FilteredTime[LicksL$FilteredTime > 0]/1000, y = L_Threshold)
      points(PointPlot, col = "darkgreen")
    }}
    
    if (PlotRawInterval == "L"||PlotRawInterval == "B") {
      hist(LicksL$RawInterval[LicksL$RawInterval <= 500], main = LeftLab, xlab = "Interval (mS)", xlim = c(0, 500), breaks = seq(0, 500, by = 10))
    }
    
    if (PlotInterval == "L"||PlotInterval == "B") {
      hist(LicksL$RawInterval[LicksL$FilteredInterval <= 1000 & LicksL$FilteredInterval > 0], main = LeftLab, xlab = "Interval (mS)", xlim = c(0, 1000), breaks = seq(0, 1000, by = 10), col = "green4")
    }
    
    if (PlotCumulative == "L" || PlotCumulative == "B") {
      plot(x = c(0, (LicksL$RawTime[LicksL$FilteredInterval != 0]) / 1000, (nrow(RawData) * SampleRate)/ 1000), y = c(0, 1:LL, LL), type = "l", xlab = "Seconds", ylab = "Licks", main = LeftLab, xaxt = "n", panel.first= grid())
      axis(1, at = seq(0, ceiling(x = nrow(RawData)/100000)*100, by = 100), las = 2) #sets X axis to 100-sec increments
    }
    
    if (PlotZoom == "L"){
      MIN <- (ZoomMin * 1000)
      MAX <- (ZoomMax * 1000)
      PlotSmS <- data.frame((RawData$mS[RawData$mS > MIN & RawData$mS < MAX]) / 1000)
      PlotMatrixS <- data.frame(RawData$Left[RawData$mS > MIN & RawData$mS < MAX], L_Threshold)
      matplot(y = PlotMatrixS, x = PlotSmS, main = LeftLab, xlab = "Seconds", ylab = "Lin.", type = "l", lty = 1, col = c("black", "red"))
      if (length(LicksL$FilteredTime[LicksL$FilteredTime > MIN & LicksL$FilteredTime < MAX]) != 0) {
        PointPlot <- data.frame(x = LicksL$FilteredTime[LicksL$FilteredTime > 0 & LicksL$FilteredTime > MIN & LicksL$FilteredTime < MAX]/1000, y = L_Threshold)
        points(PointPlot, col = "darkgreen")
      }
    }
    
    # Stereo Operations ====
    if (WavData@stereo == T) {
    # Set Lick Threshold, Right ==== 
    R_Threshold <- (mean(RawData$Right) - (ThresholdSD_R * sd(RawData$Right)))
    
    # Test Threshold, Right ==== 
    RawData$RF1 <- 0
    RawData$RF1[RawData$Right < (R_Threshold)] <- 1
    
    # Begin Lick Coding, Right ==== 
    N <- 1
    RawData$RF2 <- c(0, RawData$RF1[1:(nrow(RawData)-1)])
    RawData$RF3 <- RawData$mS
    RawData$RF3[RawData$RF1 - RawData$RF2 != 1] <- 0
    if (sum(RawData$RF3) > 0){
      LicksR <- data.frame(RawTime = RawData$RF3[RawData$RF3 != 0])
    } else {
      LicksR <- data.frame(RawTime = 0)
    }
    LicksR[, 2:7] <- 0
    colnames(LicksR)[2:7] <- c("FilteredTime", "BurstCount", "RawInterval", "FilteredInterval", "DoubleContact", "LickTotal")
    LicksR <- rbind(0, LicksR)
    for(X in 2:nrow(LicksR)){
      if ((LicksR$RawTime[X] - LicksR[X - 1, 2]) > Cutoff) {
        LicksR[X, 2] <- LicksR[X, 1]
      } else {
        LicksR[X, 2] <- LicksR[X - 1, 2]}}
    RL <- 0
    for(X in 2:nrow(LicksR)) {
      if ((LicksR[X, 2] - LicksR[X - 1, 2]) > Cutoff) {
        RL <- RL + 1
      }
      LicksR[X, 7] <- RL
    }
    
    # Burst Analysis, Right ==== 
    DC <- 0
    BurstR <- data.frame(0, 0, 0, 0)
    colnames(BurstR) <- c("Licks", "Duration", "Start", "Stop")
    for(X in 2:nrow(LicksR)) {
      if ((LicksR[X, 2] - LicksR[X - 1, 2]) > Cutoff & (LicksR[X, 2] - LicksR[X - 1, 2]) < 1000) {
        LicksR[X, 3] <- (LicksR[X - 1, 3] + 1)
        LicksR[X - 1, 6] <- DC
      } else {
        if ((LicksR[X, 2] - LicksR[X - 1, 2]) < Cutoff) {
          LicksR[X, 3] <- LicksR[X - 1, 3]
          DC <- DC + 1
          LicksR[X - 1, 6] <- DC
        } else {
          LicksR[X, 3] <- 1
          LicksR[X - 1, 6] <- DC
          DC <- 0
        }}}
    N <- 0
    for(X in 2:nrow(LicksR)) {
      if ((LicksR[X, 3] - LicksR[X - 1, 3]) < (-1)) {
        N <- N + 1
        BurstR[N, 1] <- LicksR[X - 1, 3]
      }}
    
    if (LicksR[nrow(LicksR), 3] > 2) {
      N <- N + 1
      BurstR[N, 1] <- LicksR[nrow(LicksR), 3]}#For bursts that end on the last lick
    
    RB <- N
    N <- 0
    for(X in 2:nrow(LicksR)) {
      if ((LicksR[X, 3] - LicksR[X - 1, 3]) < (-1)) {
        N <- N + 1
        BurstR[N, 2] <- LicksR[X - 1, 2] - LicksR[X - (LicksR[X - 1, 3] + LicksR[X - 1, 6]), 2]
        BurstR[N, 4] <- LicksR[X - 1, 1]
        BurstR[N, 3] <- LicksR[(X - (LicksR[X - 1, 3] + LicksR[X - 1, 6])), 1]
      }}
    
    if (LicksR[nrow(LicksR), 3] > 2) {
      N <- N + 1
      BurstR[N, 2] <- LicksR[X - 1, 2] - LicksR[X - (LicksR[X - 1, 3] + LicksR[X - 1, 6]), 2]
      BurstR[N, 3] <- LicksR[(nrow(LicksR)+1)-(LicksR[nrow(LicksR), 3] + LicksR[nrow(LicksR)-1, 6]), 1]
      BurstR[N, 4] <- LicksR[X , 1]
    }#For bursts that end on the last lick
    
    PauseR <- data.frame(c(BurstR$Start, (nrow(RawData) * SampleRate)), c(0, BurstR$Stop))
    colnames(PauseR)[1] <- "Start"
    colnames(PauseR)[2] <- "Stop"
    PauseR$Pause <- PauseR$Start - PauseR$Stop
    
    # Interval Analyses, Right ==== 
    LicksR[1, 4] <- 0
    for(X in 2:nrow(LicksR)) {
      LicksR[X, 4] <- (LicksR[X, 1] - LicksR[X - 1, 1])
    }
    LicksR[1, 5] <- 0
    for(X in 2:nrow(LicksR)) {
      LicksR[X, 5] <- (LicksR[X, 2] - LicksR[X-1, 2])
    }
    
    # Pause Analysis ==== 
    PauseB <- data.frame(c(BurstR$Start, BurstL$Start), c(BurstR$Stop, BurstL$Stop))
    colnames(PauseB)[1] <- "Start"; colnames(PauseB)[2] <- "Stop"
    PauseB <- PauseB[order(PauseB$Start), ]
    PauseB <- data.frame(c(PauseB$Start, (nrow(RawData) * SampleRate)), c(0,PauseB$Stop))
    colnames(PauseB)[1] <- "Start"
    colnames(PauseB)[2] <- "Stop"
    PauseB$Pause <- PauseB$Start - PauseB$Stop
    
    # Plotting, Right ==== 
    DSF <- 4
    RawData$RF4 <- 0
    RawData$RF4[(LicksR$RawTime[LicksR$FilteredInterval != 0]) / SampleRate] <- 1 #Creates a record of counted events
    
    if (PlotSignal == "R" || PlotSignal == "B") {
      PlotMat <- data.frame(x = RawData$mS[(1:floor(nrow(RawData)/DSF))*DSF]/1000, y = RawData$Right[(1:floor(nrow(RawData)/DSF))*DSF], thresh = R_Threshold)
      matplot(x = PlotMat$x, y = PlotMat[,2:3], main = RightLab, xlab = "Seconds", ylab = "Lin.", type = "l", lty = 1, col = c(1,2),xaxt = "n")
      axis(1, at = seq(0, ceiling(x = nrow(RawData)/100000)*100, by = 100), las = 2) #sets X axis to 100-sec increments
      if (sum(RawData$RF3) > 0) {
        PointPlot <- data.frame(x = LicksR$FilteredTime[LicksR$FilteredTime > 0]/1000, y = R_Threshold)
        points(PointPlot, col = "darkgreen")
      }}

    if (PlotRawInterval == "R"||PlotRawInterval == "B") {
      hist(LicksR$RawInterval[LicksR$RawInterval <= 500], main = RightLab, xlab = "Interval (mS)", xlim = c(0, 500), breaks = seq(0, 500, by = 10))
    }
    
    if (PlotInterval == "R" || PlotInterval == "B") {
      hist(LicksR$RawInterval[LicksR$FilteredInterval <= 1000 & LicksR$FilteredInterval > 0], main = RightLab, xlab = "Interval (mS)", xlim = c(0, 1000), breaks = seq(0, 1000, by = 10), col = "green4")
    }
    
    if (PlotCumulative == "R" || PlotCumulative == "B") {
      plot(x = c(0, (LicksR$RawTime[LicksR$FilteredInterval != 0]) / 1000, (nrow(RawData) * SampleRate)/ 1000), y = c(0, 1:RL, RL), type = "l", xlab = "Seconds", ylab = "Licks", main = RightLab, xaxt = "n", panel.first= grid())
      axis(1, at = seq(0, ceiling(x = nrow(RawData)/100000)*100, by = 100), las = 2) #sets X axis to 100-sec increments
    }

    if (PlotZoom == "R"){
      MIN <- (ZoomMin * 1000)
      MAX <- (ZoomMax * 1000)
      PlotSmS <- data.frame((RawData$mS[RawData$mS > MIN & RawData$mS < MAX]) / 1000)
      PlotMatrixS <- data.frame(RawData$Right[RawData$mS > MIN & RawData$mS < MAX], R_Threshold)
      matplot(y = PlotMatrixS, x = PlotSmS, main = RightLab, xlab = "Seconds", ylab = "Lin.", type = "l", lty = 1, col = c("black", "red"))
      if (length(LicksR$FilteredTime[LicksR$FilteredTime > MIN & LicksR$FilteredTime < MAX]) != 0) {
        PointPlot <- data.frame(x = LicksR$FilteredTime[LicksR$FilteredTime > 0 & LicksR$FilteredTime > MIN & LicksR$FilteredTime < MAX]/1000, y = R_Threshold)
        points(PointPlot, col = "darkgreen")
      }
    }
    
    }
    if (WavData@stereo == F){
      #Data Export, Mono ====
      if (sum(RawData$LF3 > 0)) {
        Lickout <- data.frame(
          "LLickTotal" = LL,
          "LMin.Lick" = max(LicksL$LickTotal[LicksL$FilteredTime < 60000]),
          "LBurstTotal" = LB,
          "LBurstLicks" = mean(BurstL$Licks),
          "LBurstDur." = mean(BurstL$Duration),
          "LPauseTotal" = nrow(PauseL),
          "LPauseDur." = mean(PauseL$Pause),
          "LMPI" = mean(LicksL$FilteredInterval[LicksL$FilteredInterval > Cutoff & LicksL$FilteredInterval <= 160]),
          "LEfficiency" = (NROW(LicksL$FilteredInterval[LicksL$FilteredInterval > 0 & LicksL$FilteredInterval <= 160])) / LL,
          "LLatency" = BurstL[1, 3],
          #"LThresPcnt" = (L_Threshold / min(RawData$Left)),
          "Cut" = Cutoff,
          "SD_L" = ThresholdSD_L,
          "ID" = D,
          stringsAsFactors = F)
        ILItabX <- data.frame(ILI = LicksL$FilteredInterval[LicksL$FilteredInterval != 0], ID = i)
        
      } else {
        Lickout <- data.frame(
          "LLickTotal" = 0,
          "LMin.Lick" = 0,
          "LBurstTotal" = 0,
          "LBurstLicks" = 0,
          "LBurstDur." = 0,
          "LPauseTotal" = 1,
          "LPauseDur." = (nrow(RawData) * SampleRate),
          "LMPI" = NaN,
          "LEfficiency" = NaN,
          "LLatency" = (nrow(RawData) * SampleRate),
          #"LThresPcnt" = (L_Threshold / min(RawData$Left)),
          "Cut" = Cutoff,
          "SD_L" = ThresholdSD_L,
          "ID" = D,
          stringsAsFactors = F)
      }
    } else {
      # Data Export, Stereo ==== 
      if (sum(RawData$LF3) > 0) {
        LickoutL <- data.frame(
          "LLickTotal" = LL,
          "LMin.Lick" = max(LicksL$LickTotal[LicksL$FilteredTime < 60000]),
          "LBurstTotal" = LB,
          "LBurstLicks" = mean(BurstL$Licks),
          "LBurstDur." = mean(BurstL$Duration),
          "LPauseTotal" = nrow(PauseL),
          "LPauseDur." = mean(PauseL$Pause),
          "LMPI" = mean(LicksL$FilteredInterval[LicksL$FilteredInterval > Cutoff & LicksL$FilteredInterval <= 160]),
          "LEfficiency" = (NROW(LicksL$FilteredInterval[LicksL$FilteredInterval > 0 & LicksL$FilteredInterval <= 160])) / LL,
          "LLatency" = BurstL[1, 3],
          stringsAsFactors = F)
      } else {
        LickoutL <- data.frame(
          "LLickTotal" = 0,
          "LMin.Lick" = 0,
          "LBurstTotal" = 0,
          "LBurstLicks" = 0,
          "LBurstDur." = 0,
          "LPauseTotal" = 1,
          "LPauseDur." = (nrow(RawData) * SampleRate),
          "LMPI" = NaN,
          "LEfficiency" = NaN,
          "LLatency" = (nrow(RawData) * SampleRate),
          stringsAsFactors = F)
      }
      if (sum(RawData$RF3) > 0) {
        LickoutR <- data.frame(
          "RLickTotal" = RL,
          "RMin.Lick" = max(LicksR$LickTotal[LicksR$FilteredTime < 60000]),
          "RBurstTotal" = RB,
          "RBurstLicks" = mean(BurstR$Licks),
          "RBurstDur." = mean(BurstR$Duration),
          "RPauseTotal" = nrow(PauseR),
          "RPauseDur." = mean(PauseR$Pause),
          "RMPI" = mean(LicksR$FilteredInterval[LicksR$FilteredInterval > Cutoff & LicksR$FilteredInterval <= 160]),
          "REfficiency" = (NROW(LicksR$FilteredInterval[LicksR$FilteredInterval > 0 & LicksR$FilteredInterval <= 160])) / RL,
          "RLatency" = BurstR[1, 3],
          "PauseTotal" = nrow(PauseB),
          "PauseDur." = mean(PauseB$Pause),
          "PrefRvL" = (RL/(RL+LL)),
          #"LThresPcnt" = (L_Threshold / min(RawData$Left)),
          #"RThresPcnt" = (R_Threshold / min(RawData$Right)),
          "Cut" = Cutoff,
          "SD_L" = ThresholdSD_L,
          "SD_R" = ThresholdSD_R,
          "ID" = D,
          stringsAsFactors = F)
      } else {
        LickoutR <- data.frame(
          "RLickTotal" = 0,
          "RMin.Lick" = 0,
          "RBurstTotal" = 0,
          "RBurstLicks" = 0,
          "RBurstDur." = 0,
          "RPauseTotal" = 1,
          "RPauseDur." = (nrow(RawData) * SampleRate),
          "RMPI" = NaN,
          "REfficiency" = NaN,
          "RLatency" = (nrow(RawData) * SampleRate),
          "PauseTotal" = 1,
          "PauseDur." = (nrow(RawData) * SampleRate),
          "PrefRvL" = (RL/(RL+LL)),
          #"LThresPcnt" = (L_Threshold / min(RawData$Left)),
          #"RThresPcnt" = (R_Threshold / min(RawData$Right)),
          "Cut" = Cutoff,
          "SD_L" = ThresholdSD_L,
          "SD_R" = ThresholdSD_R,
          "ID" = D,
          stringsAsFactors = F)
      }
      if (sum(c(RawData$LF3, RawData$RF3)) > 0){
      ILItabX <- data.frame(ILI = c(LicksR$FilteredInterval[LicksR$FilteredInterval != 0], LicksL$FilteredInterval[LicksL$FilteredInterval != 0]), ID = i)
      }
      Lickout <- data.frame(cbind(LickoutL, LickoutR))
    }
    
    if (i > 1) { 
      LickExport[(nrow(LickExport) + 1), ] <- Lickout[1, ]
    } else {
      LickExport <- Lickout
    }
    
    if (Output != F) {
      write.csv(LickExport, file = Output, row.names = F)
    }
    
    if (WavData@stereo == T) {
    }

    if (ILIFile!= F & sum(c(RawData$LF3, RawData$RF3)) > 0) {
      if (i>1)
      {
        ILItab <- read.csv(ILIFile,header = T)
        ILItab <- data.frame(c(ILItab$ILI, ILItabX$ILI), c(ILItab$ID, ILItabX$ID))
        colnames(ILItab) <- c("ILI", "ID")
        write.csv(ILItab,file = ILIFile, row.names = F)
        #View(ILItab)
      } else {
        ILItab <- ILItabX
        colnames(ILItab) <- c("ILI","ID")
        write.csv(ILItab,file = ILIFile,row.names = F)
      }}
    if (File!= F) {
      break
    } else {
      print(i)
    }
  }
  if (File!= F) {
    assign(x = Title, value =  Lickout, envir = .GlobalEnv)
  } else {
    assign(x = Title, value = LickExport, envir = .GlobalEnv)
    #assign(x = paste(Title, "Directory", sep = ""), value = Directory, envir = .GlobalEnv)
  }
  if (Timer!= F) {
    proc.time() - ptm
  }
}