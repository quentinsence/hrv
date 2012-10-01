#load and analyse emwave.emdb data
#sqlite emdb file is composed of 3 tables: Client PrimaryData VersionTable

#TODO
#export all samples in ascii RR data
#plot "The Zone" lines
#bioconduction PROcess peaks function
#split screen to get "banking" (45 degrees in curves)

#power spectrum VLF LF HF histogram IBI/60 Hz
#frequency band VLF 0-0.04 Hz, LF 0.04 - 0.15 Hz, HF 0.15 0.5 Hz
#######

library(RSQLite)

#sqlite db location is system dependent
user <- Sys.info()['user']

#define a local user in the local.r file to override the user name and path

if(file.exists('local.r')) {
	source("local.r")
}

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

#final scores
#as.numeric(unlist(h$AccumZoneScore[1])[length(unlist(h$AccumZoneScore[1]))-3])
#won't work if value > 255, must have hex values in group of 4
#h$FinalScore <- sapply( h$AccumZoneScore, FUN = function(x) as.numeric(unlist(x)[length(unlist(x))-3]) )
#lappy readBin(unlist(h$AccumZoneScore[10]),"int",size=4,endian="little",n=length(unlist(h$AccumZoneScore[10]))/4)

h$PctLow <- 100 - h$PctMedium - h$PctHigh
h$date <- as.POSIXct(h$IBIStartTime,origin="1970-01-01")
h$end  <- as.POSIXct(h$IBIEndTime,origin="1970-01-01")
h$sessiontime <- h$IBIEndTime - h$IBIStartTime
h$Level <- h$ChallengeLevel
h$ChallengeLevel <- factor(h$ChallengeLevel,levels=c(1,2,3,4),labels=c("Low","Medium","High","Highest"))
h$Endian <- factor(h$Endian,levels=c(0,1),labels=c("big","little"))


#foreach session
#I do not know how to stop unlist recycling all sessions in 1 so I'm using "for"
for (n in 1:dim(h)[1]) {
	  #convert hex interbeat intervals to decimal bpm
    #cumulate each hex interbeat interval to decimal seconds
    h$timeIBI[n] <- list(0.001 * cumsum(readBin(unlist(h$LiveIBI[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$LiveIBI[n]))/4)))
	  #convert back to integers
	  h$AccumZoneScore[n] <- list(readBin(unlist(h$AccumZoneScore[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$AccumZoneScore[n]))/4))
	  h$ZoneScore[n] <- list(readBin(unlist(h$ZoneScore[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$ZoneScore[n]))/4))

    h$LiveIBI[n] <- list(readBin(unlist(h$LiveIBI[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$LiveIBI[n]))/4))
    #convert interbeat intervals to beats per minute [exclude zeroes to avoid Inf]
    h$bpm[n] <- list(60*1000/unlist(h$LiveIBI[n])[unlist(h$LiveIBI[n])>0])
    #cumulate each hex interbeat interval to decimal seconds
    h$timeIBI[n] <- list(0.001 * cumsum(unlist(h$LiveIBI[n])[unlist(h$LiveIBI[n])>0]))
    #longest singular duration spent in high coherence
    #rle computes the lenghts of runs of equal values, we are looking for the longest run of "2"
    h$maxhicoherence[n] <- h$sessiontime[n] * with(rle(unlist(h$ZoneScore[n])==2),max(lengths[!!values==TRUE])) / length(unlist(h$ZoneScore[n]))
}

#recalc FinalScore as decimal
h$FinalScore <- sapply( h$AccumZoneScore, FUN = function(x) unlist(x)[length(unlist(x))] )

h$Weekday <- strftime((as.POSIXct(h$IBIStartTime,origin="1970-01-01")),format="%w")

hrvplot <- function(n=1) {
pulse  <- unlist(h$bpm[n])
pulset <- unlist(h$timeIBI[n])
#score <- readBin(unlist(h$AccumZoneScore[n]),"int",size=4,endian=h$Endian,n=length(unlist(h$AccumZoneScore[n]))/4)

par(mfrow=c(2,1),mai=c(0.4,0.4,0.2,0.2),lab=c(10,10,7))
plot(pulse ~ pulset,xlab="time",ylab="mean Heart Rate (BPM)",type ="l")
plot(ts(unlist(h$AccumZoneScore[n])),xlab="time",ylab="Accumulated Coherence Score",type ="l")
#plot(unlist(h$bpm[n]) ~ unlist(h$timeIBI[n]),xlab="time",ylab="mean Heart Rate (BPM)",type ="l")
#plot(ts(unlist(h$ZoneScore[n])),xlab="time",ylab="Accumulated Coherence Score",type ="l")

#LEGEND
cat('Start',strftime(h$date[n],format="%x %X"),'\n')
cat('End  ',strftime(h$end[n],format="%x %X"),'\n')
cat('session time',as.integer(h$sessiontime[n]/60),'min',h$sessiontime[n] %% 60,'sec','\n')
cat('mean HR:',as.integer(mean(pulse)),'bpm\n')
cat('final score',h$FinalScore[n],'\n')
cat('difficulty level',h$ChallengeLevel[n],'\n')
cat('Coherence Ratio Low/Med/High%',as.integer(h$PctLow[n]),'/',as.integer(h$PctMedium[n]),'/',as.integer(h$PctHigh[n]),'\n')
}

#start by displaying summary of all sessions
par(mfrow=c(4,1),mai=c(0.4,0.7,0.2,0.2),lab=c(10,10,7))
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
boxplot(h$FinalScore ~ h$Weekday,horizontal=TRUE,main="final score by day of week (0=Sunday)")

#plot(h$PctHigh ~ h$date,type="l",col="green")
#lines(h$PctMedium ~ h$date,type="l",col="blue")

#export session(s) as ascii RR data file, readable by kubios
hrvexport <- function(x1="") {
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


