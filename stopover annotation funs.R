

#to summarize spatial and temporal information that will allow me to characterize it as 
# a stopover/ not a stopover, and then type of stopover
# produces a short table and a plot with 2 days pre/post
# as an html page to look at a whole migration easily?
stop_summary_table <- function(migration_df){
  
  mig_df_sub <- migration_df[which(migration_df$stopover==1),]
  mig_df_sub$stop.id2 <- mig_df_sub$stop.id
  
  stops <- split(mig_df_sub, mig_df_sub$stop.id)
  
  s.ids <- names(stops)
  
  stops_table <- lapply(stops, function(x){
           s.dists <- dist_pts(x)
           n <- nrow(s.dists)
           m.id <- x$ID_mig[1]
           s.id <- x$stop.id2[1]
           s.dates <- range(x$study.local.timestamp)
           s.dur <- diff(s.dates)
           s.roost.prop <- sum(x$roost)/nrow(x)
           
           roost.am <- which(hour(s.dists$study.local.timestamp)==8)
           roost.pm <- which(hour(s.dists$study.local.timestamp)==17)
           ymax <- ifelse(max(s.dists$dist.km, na.rm=T)<20,
                          20,
                          max(s.dists$dist.km, na.rm=T))
           p <- plot(s.dists$dist.km ~s.dists$study.local.timestamp, 
                pch = 19, 
                ylim=c(0, ymax),
                ylab="Distance (km)", 
                xlab="Timestamp")+
             abline(v=s.dists$study.local.timestamp[roost.am], 
                    col="red")+
             abline(v=s.dists$study.local.timestamp[roost.pm], 
                    col="red")+
             abline(h=20, 
                    col="blue")
           print(p)
           data.frame(m.id,
                      s.id,
                      n, 
                      "start"= s.dates[1], 
                      "end"=s.dates[2],
                      s.dur, 
                      s.roost.prop)
           })

  return(stops_table)
 
  #with(migration_df,
  #     ScanTrack(study.local.timestamp, 
  #               X, 
  #               Y, 
  #               col = as.factor(stop.id)))
}

#if(s.roost.prop > 0.9){
  #stopover = 0
#}

# if(s.roost.prop > 0.6 & s.dur < 24){
# stop.cat = "overnight"
# })

# works for stops that we do know
###################################################
#left off here - need to pair the plots with the tables, 
#write code using mig id and stop id to easily change the stop code

# need to make something similar for stops that arent identified -
# just 1 day pre/post any day that day.dist < # km?
###################################################

fpt_nogap <- function(df.traj, r = 20000){
  #find r
  #i <- fpt(df.traj, seq(14000,30000, length=8))
  #i <- fpt(df.traj, seq(2500,6000, length=20))
  #plot(i, scale = 500, warn = FALSE,type="l")
  #n <- which(varlogfpt(i,graph = F) == min(varlogfpt(i,graph = F)))
  #r <- attr( varlogfpt(i,graph = F), "radii")[n]
  #use r
  df.fpt <- fpt(df.traj,r)

  return(df.fpt)
}

find_threshold <- function(myfpt, mig_nogap){
  # if(k==5){
  #  q<-0.41
  #   }else 
  if(sd(subset(myfpt,fpt>1)$fpt,na.rm=T)>10|
     mean(subset(myfpt,fpt>1)$fpt,na.rm=T)>20){
    q<-0.55 #high variance
  }else{
    q<-0.75 #low variances
  } #takes variance into account
  
  threshold<-as.numeric(ceiling(quantile(subset(myfpt,fpt>1)$fpt,q)))
  return(threshold)
}

#library(lutz)
library(suncalc)

roost_times <- function(DF, mig_nogap){
  tz(DF$timestamp)<-"UTC"
  sundat<-data.frame(date = as.Date(DF$study.local.timestamp),
                     datetime = DF$timestamp,
                     lat = DF$location.lat,
                     lon = DF$location.long,
                     local.time = DF$study.local.timestamp)
  tz(sundat$local.time)<-"UTC"
  
  sundat$date <-as.Date(sundat$date)
  sundat$hr <- hour(sundat$local.time)
  sundat <- sundat[which(!(is.na(sundat$date))),]
  sundat <- sundat[which(!(is.na(sundat$datetime))),]
  
  sundat$offset<-sundat$local.time - sundat$datetime
  if(abs(max(sundat$offset))>3000){
    sundat$offset<-sundat$offset / 3600
  }

  sundat2 <- getSunlightTimes(data=sundat, 
                              keep = c("sunrise", 
                                       "sunset"))
  
  sundat<-merge(sundat,sundat2,by=c("date","lat","lon"))
  sundat$sunrise<- hour(sundat$sunrise + sundat$offset)
  sundat$sunset<- hour(sundat$sunset + sundat$offset)
  sundat$sunrise<- ifelse(sundat$sunrise > 12, 
                          sundat$sunrise - 24, 
                          sundat$sunrise)
  sundat$sunset<- ifelse(sundat$sunset > 23, 
                         sundat$sunset - 24, 
                         ifelse(sundat$sunset < 12, 
                                sundat$sunset + 24, 
                                sundat$sunset))

  trj<- as.ltraj(id = as.character(DF$individual.local.identifier), 
                 xy = DF[,c("X", "Y")],
                 date = DF$study.local.timestamp,
                 burst = as.character(DF$individual.local.identifier))
  
  trjDF<-ld(trj)
  names(trjDF)[3]<-"study.local.timestamp"
  trjDF<-merge(trjDF,DF[,c('timestamp','study.local.timestamp','day')],
               by="study.local.timestamp")
  trjDF$speed<-(trjDF$dist/1000)/(trjDF$dt/3600) #km per hour    
  trjDF$hr<-lubridate::hour(trjDF$study.local.timestamp)
  
  r2 <- sort(trjDF[which(trjDF$speed > 5),]$hr)[6]
  r1 <- sort(trjDF[which(trjDF$speed > 5),]$hr)[length(trjDF[which(trjDF$speed > 5),]$hr)-6]
  
  #mig_nogap$hr<-lubridate::hour(mig_nogap$time)
  #mig_nogap$roost <- ifelse(mig_nogap$hr>r2 & mig_nogap$hr < r1, 0, 1)
  return(c(r1,r2))
}

set.ID.pop <- function(x, VAR){
  Ind_ID <- c(1,
              abs(
                diff(
                  as.numeric(
                    as.factor(
                      x$individual.local.identifier
                    )
                  ))))
  ID_seg <- c(1,
              abs(
                diff(
                  as.numeric(
                    as.factor(
                      VAR
                    )
                  ))))
  ID_seg <- Ind_ID + ID_seg
  ID_seg <- ifelse(ID_seg == 2, 1, ID_seg)
  ID_seg <- cumsum(
    ifelse(ID_seg !=0, 1, 0))
  return(ID_seg)
}


set.ID.ind <- function(x, VAR){
  ID_seg <- c(1,
              abs(
                diff(
                  as.numeric(
                    as.factor(
                      VAR
                    )
                  ))))
  ID_seg <- cumsum(
    ifelse(ID_seg !=0, 1, 0))
  return(ID_seg)
}

