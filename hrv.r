#load and analyse emwave.emdb data

#see emwave.emdb.txt for details on the SQLite tables fields

#TODO
#plot "The Zone" lines
#split screen to get "banking" (45 degrees in curves)

#power spectrum VLF LF HF histogram IBI/60 Hz
#frequency band VLF 0-0.04 Hz, LF 0.04 - 0.15 Hz, HF 0.15 0.5 Hz
#######

library(RSQLite)

#sqlite db location is system dependent
user <- Sys.info()['user']

#define a local user in the local.r file to override the user name and path
if(file.exists('local.r')) {source("local.r")}

if( Sys.info()['sysname'] == "Windows") {
	emdb <- paste('C:/Documents and Settings/',user,'/My\ Documents/emWave/emwave.emdb',sep="")
} else {
    #assumed the linux OS has the same username
	emdb <- paste('/windows/D/Documents\ and\ Settings/',user,'/My\ Documents/emWave/emwave.emdb',sep="")
}

#if emwave directory cannot be found then assume we have a copy of the db in the working directory
if(!file.exists(emdb)) {
    cat(emdb,' not found, using a local copy\n')
    emdb <- 'emwave.emdb'
}
############# CONNECT & LOAD
m <- dbDriver("SQLite")
con <- dbConnect(m, dbname=emdb)

dbListTables(con)
#sessions may not be stored in chronological order as sessions can be carried out and uploaded only at a later time
rs <- dbSendQuery(con, "select * from PrimaryData order by IBIStartTime")
h <- fetch(rs, n=-1)
dbClearResult(rs)

dbDisconnect(con)
#############

h$PctLow         <- 100 - h$PctMedium - h$PctHigh
h$date           <- as.POSIXct(h$IBIStartTime,origin="1970-01-01")
h$end            <- as.POSIXct(h$IBIEndTime,origin="1970-01-01")
h$sessiontime    <- h$IBIEndTime - h$IBIStartTime
h$Level          <- h$ChallengeLevel
h$ChallengeLevel <- factor( h$ChallengeLevel
                            ,levels=1:4
                            ,labels=c("Low","Medium","High","Highest"))
h$Endian         <- factor( h$Endian
                            ,levels=0:1
                            ,labels=c("big","little"))
h$Weekday        <- factor( strftime((as.POSIXct(h$IBIStartTime,origin="1970-01-01")),format="%w")
                            ,levels=0:6
                            ,labels=c('Sun','Mon','Tue','Wed','Thu','Fri','Sat'))

#convert lists from 4 bytes hexadecimal format to integer
h$LiveIBI        <- lapply(h$LiveIBI       ,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
h$SampledIBI     <- lapply(h$SampledIBI    ,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
h$ArtifactFlag   <- lapply(h$ArtifactFlag  ,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
h$AccumZoneScore <- lapply(h$AccumZoneScore,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
h$ZoneScore      <- lapply(h$ZoneScore     ,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
h$EntrainmentParameter <- lapply(h$EntrainmentParameter,function(x) readBin(x,"int",size=4,endian=h$Endian,n=length(x)/4))
#convert interbeat intervals to beats per minute [exclude zeroes to avoid Inf]
h$BPM            <- lapply(h$LiveIBI,function(x) 60*1000/x[x>0] )
#cumulate each hex interbeat interval to decimal seconds [exclude zeroes to match BPM lists]
h$timeIBI        <- lapply(h$LiveIBI,function(x) 0.001*cumsum(x[x>0]) )
#longest singular duration spent in high coherence
#rle computes the lenghts of runs of equal values
#we are looking for the longest run of "2" i.e. high coherence
#0 = low coherence, 1 = medium coherence
#converts length of high coherence runs to seconds by multiplying by Entrainment sampling interval (ms * 0.001 = seconds)
h$maxhicoherence <- 0.001 * h$EntrainmentIntervalTime * sapply(h$ZoneScore,function(x) with(rle(x==2),max(lengths[!!values==TRUE])))

h$AverageBPM     <- sapply( h$BPM, mean )
h$FinalScore     <- sapply( h$AccumZoneScore, function(x) x[length(x)] )

########## PLOTS

hrvweekday <- function() {
    par(mfrow=c(7,1),mai=c(0.4,0.7,0.2,0.2),lab=c(10,10,7))
    #calculate week number for each entry so all timeseries by weekday are aligned
    h$Week <- as.integer(1 + (h$IBIStartTime - h$IBIStartTime[1]) / (86400*7))

    #build a matrix week by weekday and fill in NA when scores are missing
    hw <- data.frame(h$Week,h$Weekday,h$FinalScore)
    mflat <- rbind(hw,data.frame(expand.grid(h.Week=unique(h$Week),h.Weekday=unique(h$Weekday)),h.FinalScore=NA))
    
    #sort by descending score, lower scores and NA on the same week,weekday will be duplicates
    mflat <- mflat[order(mflat$h.FinalScore,decreasing=TRUE),]
    #remove duplicate (week,weekday), keep only the highest score of the (week,weekday)
    mflat <- mflat[!duplicated(cbind(mflat$h.Week,mflat$h.Weekday)),]
    
    #sort back to chronological order
    mflat <- mflat[order(mflat$h.Week,mflat$h.Weekday),]

    ScoreMatrix <- matrix(mflat$h.FinalScore,nrow=7,ncol=max(h$Week))
    #set threshold at 200 points
    horizonplot(ts(t(ScoreMatrix)),layout=c(7,1),origin=200)
}

hrvplot <- function(n=1) {
    pulse  <- unlist(h$BPM[n])
    pulset <- unlist(h$timeIBI[n])
    #score <- readBin(unlist(h$AccumZoneScore[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$AccumZoneScore[n]))/4)
    
    par(mfrow=c(3,1),mai=c(0.4,0.4,0.2,0.2),lab=c(10,10,7))
    plot(pulse ~ pulset,xlab="time",ylab="mean Heart Rate (BPM)",type ="l")
    #highlight low coherence sequences
    abline(v=5*(which(ts(unlist(h$ZoneScore[n]) == 0))-1),col="red")
    plot(ts(unlist(h$AccumZoneScore[n])),xlab="time",ylab="Accumulated Coherence Score",type ="l")
    plot(ts(unlist(h$EntrainmentParameter[n])),xlab="time",ylab="Entrainment Parameter",type ="l")
    #plot(unlist(h$BPM[n]) ~ unlist(h$timeIBI[n]),xlab="time",ylab="mean Heart Rate (BPM)",type ="l")
    #plot(ts(unlist(h$ZoneScore[n])),xlab="time",ylab="Accumulated Coherence Score",type ="l")
    
    #LEGEND
    cat('Start',strftime(h$date[n],format="%x %X"),'\n')
    cat('End  ',strftime(h$end[n],format="%x %X"),'\n')
    cat('session time',as.integer(h$sessiontime[n]/60),'min',h$sessiontime[n] %% 60,'s','\n')
    cat('mean HR:',as.integer(mean(pulse)),'BPM\n')
    cat('final score',h$FinalScore[n],'\n')
    cat('difficulty level',h$ChallengeLevel[n],'\n')
    cat('Coherence Ratio Low/Med/High%',as.integer(h$PctLow[n]),'/',as.integer(h$PctMedium[n]),'/',as.integer(h$PctHigh[n]),'\n')
    cat('longest duration spent in high coherence',h$maxhicoherence[n],'s\n')
}


hrvsummary <- function() {
    #start by displaying summary of all sessions
    par(mfrow=c(5,1),mai=c(0.4,0.7,0.2,0.2),lab=c(10,10,7))
    barplot(t(cbind(h$PctLow,h$PctMedium,h$PctHigh))
            ,col=c('red','blue','green')
            ,xlab=as.numeric(h$ChallengeLevel)
            ,ylab="%"
            ,main="coherence ratio by session"
    	    ,legend=c('Low','Medium','High')
            ,args.legend = list(x = "topleft")
           )
    
    #scores plot in black and difficulty levels in grey
    plot(ts(h$FinalScore),ylab="Accumulated Score",main="final accumulated score by session")
    lines(ts(h$Level*max(h$FinalScore)/4),ylab="Challenge Level",xlab="session",main="level by session",col="grey")
    
    plot(ts(h$maxhicoherence),ylab="time (seconds)",main="longest time spent in high coherence by session")
    plot(ts(h$AverageBPM),ylab="beats per minute",main="Average BPM by session")
    boxplot(h$FinalScore ~ h$Weekday,horizontal=TRUE,main="final score by day of week")
    
    #plot(h$PctHigh ~ h$date,type="l",col="green")
    #lines(h$PctMedium ~ h$date,type="l",col="blue")
}

hrvexport <- function(x1="") {
  #export session(s) as ascii RR data file, readable by Kubios and others
  x0 <- 1
	#include h$date in the filename
  if(x1=="") {
      x1 <- dim(h)[1] 
  } else {
      x0 <- x1
  }
  for (x in x0:x1) {
    	write(unlist(h$LiveIBI[x])
    			,paste("emwave.",strftime(h$date[x],format="%Y-%m-%dT%X"),'.dat',sep="")
    			,ncolumns=1)
  }
}


