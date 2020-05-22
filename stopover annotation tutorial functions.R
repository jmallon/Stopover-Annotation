# Julie Mallon
# May 22, 2020
# https://github.com/jmallon/Stopover-Annotation
# -------------------------------------------------------
# Functions required for 'stopover annotation tutorial.Rmd'
# -------------------------------------------------------


# calcualte distance in kms
dist_pts <- function(x){
  dx <- c(x$X[1:(nrow(x)-1)] - x$X[2:(nrow(x))], NA)
  dy <- c(x$Y[1:(nrow(x)-1)] - x$Y[2:(nrow(x))], NA)
  x$dist.km <- sqrt(dx^2+dy^2)/1000
  return(x)
}

# prevent issues due to gaps in track by interpolating missing points
interpolate_gappy_track <- function(x,y,time) {
  old_df <- data.frame(x=x,y=y,time=time)
  old_df[,'dt'] <- c(as.numeric(diff(old_df[,'time'])),NA)
  old_df <- old_df[which(old_df$dt>0),]
  new_length <- old_df$dt[1:(nrow(old_df)-1)] %>% sum(na.rm=T) +1
  
  new_df <- data.frame(time = seq(old_df[1,'time'], 
                                  old_df[nrow(old_df),'time'], 
                                  by = 3600),
                       dt = rep(1, new_length)) #empty dataframe to fill
  new_df <- merge(new_df, old_df[,1:3], by="time", all = T)
  
  i.fix <- which(is.na(new_df$x))
  
  for (i in i.fix) {
    xx <- new_df[i-1,'x']
    yy <- new_df[i-1,'y']
    timedif <- new_df[i,'dt'] #i.e. 2
    tt = c(1,timedif+1)
    df <- data.frame(tt,xx,yy)
    
    mod <- lm(xx ~ tt,df)
    xs <- predict(mod,data.frame(tt=1:(timedif) ) )
    
    mod <- lm(yy ~ tt,df)
    ys <- predict(mod,data.frame(tt=1:(timedif) ) )
    
    new_df[i, 3:4] <- data.frame(x=xs,y=ys)
  }
  new_df$OBJECTID <- 1:nrow(new_df)
  new_df
}

# marcher package's scantrack function
ScanTrack <- function(time,x,y=NULL, col=NULL, cex=NULL, ...)
{
  if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
  
  if(is.null(col)) col=rgb(0,0,0,.5) 
  if(is.null(cex)) cex = 0.5
  
  layout(rbind(c(1,2), c(1,3)))
  plot(x,y,asp=1, type="o", pch=19, col=col, cex=cex, ...)
  plot(time,x, type="o", pch=19, col=col, xaxt="n", xlab="", cex=cex, ...)
  plot(time,y, type="o", pch=19, col=col, cex=cex, ...)
}

# adds a 4th panel of FPT to marcher's ScanTrack function
ScanTrack4<- function(time,x,y=NULL, col=NULL, cex=NULL, ...)
{
  if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
  
  if(is.null(col)) col=rgb(0,0,0,.5) 
  if(is.null(cex)) cex = 0.5
  
  layout(rbind(c(1,2), c(1,3),c(1,4)))
  plot(x,y,asp=1, type="o", pch=19, col=col, cex=cex, ...)
  plot(time,x, type="o", pch=19, col=col, xaxt="n", xlab="", cex=cex, ...)
  plot(time,y, type="o", pch=19, col=col, cex=cex, ...)
}



#to summarize spatial and temporal information that will allow me to characterize it as 
# a stopover/ not a stopover, and then type of stopover
# produces a short table and a plot with 2 days pre/post
# as an html page to look at a whole migration easily?
stops_summary_table <- function(migration_df){
  
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
  
}


# when using this approach for several individuals, the variance in speed and stopver duration varies
# this attempts to take that into account so it doesn't over- or under-annotate stopovers
find_threshold <- function(myfpt, mig_nogap){
  if(sd(subset(myfpt,fpt>1)$fpt,na.rm=T)>10|
     mean(subset(myfpt,fpt>1)$fpt,na.rm=T)>20){
    q<-0.55 #high variance
  }else{
    q<-0.75 #low variances
  } #takes variance into account
  
  threshold<-as.numeric(ceiling(quantile(subset(myfpt,fpt>1)$fpt,q)))
  return(threshold)
}



# this calculates typical roost start and end times to annotate roost in dataframe
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
  
  sundat$offset <- sundat$local.time - sundat$datetime
  
  if(abs(max(sundat$offset))>3000){
    sundat$offset<-sundat$offset / 3600
  }
  
  sundat2 <- getSunlightTimes(data=sundat, 
                              keep = c("sunrise", 
                                       "sunset"))
  
  sundat <- merge(sundat,
                  sundat2,
                  by=c("date",
                       "lat",
                       "lon"))
  sundat$sunrise <- hour(sundat$sunrise + sundat$offset)
  sundat$sunset <- hour(sundat$sunset + sundat$offset)
  sundat$sunrise <- ifelse(sundat$sunrise > 12,       # ensure the offset didn't mess up times
                          sundat$sunrise - 24, 
                          sundat$sunrise)
  sundat$sunset <- ifelse(sundat$sunset > 23, 
                         sundat$sunset - 24, 
                         ifelse(sundat$sunset < 12, 
                                sundat$sunset + 24, 
                                sundat$sunset))
  
  trj <- as.ltraj(id = as.character(DF$individual.local.identifier), 
                 xy = DF[,c("X", "Y")],
                 date = DF$study.local.timestamp,
                 burst = as.character(DF$individual.local.identifier))
  
  trjDF <- ld(trj)
  names(trjDF)[3] <- "study.local.timestamp"
  trjDF <- merge(trjDF,
               DF[,c('timestamp',
                     'study.local.timestamp',
                     'day')],
               by="study.local.timestamp")
  trjDF$speed <- (trjDF$dist/1000) / (trjDF$dt/3600) #km per hour   
  
  trjDF$hr <- hour(trjDF$study.local.timestamp)
  
  r2 <- sort(trjDF[which(trjDF$speed > 5),]$hr)[6]
  r1 <- sort(trjDF[which(trjDF$speed > 5),]$hr)[length(trjDF[which(trjDF$speed > 5),]$hr)-6]
  
  return(c(r1,r2))
}


# give each stop a unique ID within an individual
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
