# migration annotation table code functions


rm_anoms <- function(mydata, speed = 80, plot=T, time.unit = "sec"){
  mydata <- as.data.frame(mydata)
  trm<-which(table(mydata$timestamp)>1) #%>% as.numeric()
  if(length(trm)>0){
    trm2<-which(as.character(mydata$timestamp) %in% names(trm))
    trm2<-trm2[which(!(trm2 %in% trm))]
    mydata<-mydata[-trm2,]
  }
  mydata$dT[1] <- NA
  dX <- diff(mydata$X)
  dY <- diff(mydata$Y)
  if(time.unit=="sec"){
    dT <- mydata$dT[-1]/3600
  }else if(time.unit=="hour"){
    dT <- mydata$dT[-1]
  }else{
    print("unexpected time unit")
    stop()
  }
  
  spd <- (sqrt(dX^2+dY^2)/1000) / dT
  idx <- NULL
  
  # find points that move to a distant point and immediately return to starting position
  idx <- which(abs(dX) > 80000)[which(diff(  which(abs(dX) > 80000))<3)]
  idy <- which(abs(dY) > 80000)[which(diff(  which(abs(dY) > 80000))<3)]
  idx <- idx[idx %in% idy]
  nrm <- idx + 1
  
  # identifies locations > 80 km/h and checks if the speed
  if(length(which(spd> speed))>0){
    idx <- which(abs(spd)>speed)+1
    mydata$seq <- 0
    mydata[unique(c(idx, idx-1)),]$seq <- 1:length(unique(c(idx, idx-1)))
    mydat2 <- mydata[-idx,]
    dX2 <- diff(mydat2$X)
    dY2 <- diff(mydat2$Y)
    dT2 <- c(NA, as.numeric(diff(mydat2$timestamp)))/3600
    spd2 <- (sqrt(dX2^2+dY2^2)/1000)/dT2[-1]
    good.pts <- which((spd2[which(mydat2$seq!=0)] - speed[idx[mydat2$seq[which(mydat2$seq!=0)]]-1]) > -5)
    if(length(good.pts) > 0 ){
      nrm <- c(nrm, idx[-good.pts])
    }else{
      nrm <- c(nrm, idx)
    }
    mydata <- mydata[,-which(names(mydata)=="seq")]
    
    
    if(length(nrm)>0){
      DF<-mydata[-nrm,]
      if(plot == T){
        mydata$points<-0
        mydata[nrm,]$points<-1
        DF<-mydata[-nrm,]
        with(mydata, ScanTrack(timestamp,X,Y,col=as.factor(points)))
        with(DF, ScanTrack(timestamp,X,Y,col=as.factor(points)))
        DF<-DF[,which(!(names(DF)=="points"))]
      }
    }else{
      DF<-mydata
      
    }
  }else{
    DF<-mydata
  }
  
  return(DF)
}


# manual filtering of individuals and years, returns list of individuals
library(data.table)
mig.filter.manual <- function(bird_clean){
  
  ind.l1 <- split(bird_clean,
                    bird_clean$individual.local.identifier)
  
  ind.l1 <- lapply(ind.l1, function(x){
    x$days.cont <- floor(as.numeric(x$study.local.timestamp - 
                                x$study.local.timestamp[1])/ 86400); x})
  
  bird_clean <- rbindlist(ind.l1)
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier %in% c("Car Talk", "David") & 
      bird_clean$year > 2017),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Boneyard1" &
      bird_clean$year > 2014),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Morongo" &
      (bird_clean$year==2008&
         bird_clean$day>200)|
      bird_clean$year==2009 |
      (bird_clean$year==2006 &
         bird_clean$day<155)),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Sarkis" &
      bird_clean$year > 2007),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Steamhouse 2" &
      bird_clean$year ==2013 & 
      bird_clean$day>250),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Rosalie" &
      bird_clean$year > 2008),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Tintin" &
      bird_clean$year ==2017 & 
      bird_clean$day<150),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Leo" &
      (bird_clean$year %in%
         c(2013:2015))),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier == "Whitey" &
      bird_clean$days.cont > 2660),]
  
  bird_clean <- bird_clean[-which(
    bird_clean$individual.local.identifier %in%
      c("Butterball", "Black Knight", "Cinco", "Irma",
        "Jennie", "Mandeb", "Mark", "Mary", "Prado",
        "Schaumboch", "Simona", "Young Luro", "Tesla", 
        "Castagnino", "Disney")),]
  
  # remove birds with only 1 year
  short.birds <- bird_clean %>% 
    group_by(individual.local.identifier, year) %>% 
    summarize(n()) %>% 
    group_by(individual.local.identifier) %>% 
    summarize(n = n()) %>% 
    as.data.frame()
  
  ind.list <- split(bird_clean,
                    bird_clean$individual.local.identifier)
  
  ind.list[which(short.birds$n ==1)] <- NULL #remove birds with 1 year

  return(ind.list)
}


# FUNCTIONS
nonmig.latlon <- function (x){
  x155.y <- range(x[which(x$day>155 & x$day<255),]$Y,na.rm = T)
  x155.x <- range(x[which(x$day>155 & x$day<255),]$X,na.rm = T)
  x350.y <- range(x[which(x$day<40|x$day>350),]$Y,na.rm = T)
  x350.x <- range(x[which(x$day<40|x$day>350),]$X,na.rm = T)
  
  x225.y <- range(x[which(x$day>155 & x$day<225),]$Y, na.rm = T)
  x225.x <- range(x[which(x$day>155 & x$day<225),]$X, na.rm = T)
  x15.y <- range(x[which(x$day<50 & x$day>15),]$Y,na.rm = T)
  x15.x <- range(x[which(x$day<50 & x$day>15),]$X,na.rm = T)
  
  x30.y <- range(x[which(x$day<20|x$day>340),]$Y,na.rm = T)
  x30.x <- range(x[which(x$day<20|x$day>340),]$X,na.rm = T)
  
  
  if(x$individual.local.identifier[1] %in% c("Morongo", "Andres")){
    lat.breed <- x225.y
    lat.winter <- x350.y
    lon.breed <- x225.x
    lon.winter <- x350.x
  }else if(x$individual.local.identifier[1] %in% c("Mac")){
    lat.breed <- x225.y
    lat.winter <- x15.y
    lon.breed <- x225.x
    lon.winter <- x15.x
  }else if(x$individual.local.identifier[1] %in% c("Tintin")){
    lat.breed <- x155.y
    lat.winter <- x30.y
    lon.breed <- x155.x
    lon.winter <- x30.x
  }else if(x$individual.local.identifier[1] %in% c("Argentina", "Bariloche2")){
    lat.breed <- x350.y
    lat.winter <- x225.y
    lon.breed <- x350.x
    lon.winter <- x225.x
  }else if(length(grep("N",x$utm.zone[1]))>0){
    lat.breed <- x155.y
    lat.winter <- x350.y
    lon.breed <- x155.x
    lon.winter <- x350.x
    
  }else{
    lat.breed <- x350.y
    lat.winter <- x155.y
    lon.breed <- x350.x
    lon.winter <- x155.x
  }
  
  return(data.frame( "lat.breed" = lat.breed, 
                     "lat.winter" = lat.winter, 
                     "lon.breed" = lon.breed, 
                     "lon.winter" = lon.winter))
}


# resolve the actual start and end of migration by the first day travelled over 50 km within 10 days of marked 'migration' days

set.ID <- function(x){
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
                      x$mig.state
                  )))))
  ID_seg <- Ind_ID + ID_seg
  ID_seg <- ifelse(ID_seg == 2, 1, ID_seg)
  ID_seg <- cumsum(
    ifelse(ID_seg !=0, 1, 0))
  return(ID_seg)
}

check_state_NC <- function(x){
  tofix <- table(x$ID_mig)[which(table(x$ID_mig)<40)]
  
  if(length(tofix)>0){
    x$mig.state[which(x$ID_mig %in% names(tofix))] <- "Non-mig" 
    x$ID_mig <- set.ID(x)
  }
  return(x)
}

annotate.mig <- function(x, print.dates = F, plot = T){
  # annotate all points > radius km from winter / breeding lats as migration
  
  tofix <- NULL
  x$mig.state <- "Mig"
  nm.utm <- nonmig.latlon(x)
  breed.50km <- which(x$Y > nm.utm$lat.breed[1]  &
                        x$Y < nm.utm$lat.breed[2]  &  
                        x$X > nm.utm$lon.breed[1]  &
                        x$X < nm.utm$lon.breed[2] )
  
  winter.50km <- which(x$Y > nm.utm$lat.winter[1] &
                         x$Y < nm.utm$lat.winter[2] &
                         x$X > nm.utm$lon.winter[1] &
                         x$X < nm.utm$lon.winter[2])
  x$mig.state[c(breed.50km, winter.50km)] <- "Non-mig"
  
  x$ID_mig <- set.ID(x)
  
  #with(x, ScanTrack(study.local.timestamp, 
  #                  X,Y, col=as.factor(mig.state))) 
  
  # filter out short duration groupings
  
  x <- check_state_NC(x) 
  
  # fine-tune annotation to include days birds traveled > 50 km near start and end of migration 
  
  x <- x %>% 
    group_by(ID_mig) %>% 
    mutate(mig.start = min(date, na.rm = T), 
           mig.end = max(date, na.rm = T)) %>%
    as.data.frame()
  
  
  nm45 <- which(x$day.dist>48 & 
                  x$mig.state=="Non-mig" )
  
  adj.start <- which(x$date>x$mig.start-3)
  adj.start <- adj.start[which(adj.start %in% nm45)]
  if(length(adj.start)>0){
    x[adj.start,]$mig.state <- "Mig" 
  }
  
  
  adj.end <- which(x$date<x$mig.end+3)
  adj.end <- adj.end[which(adj.end %in% nm45)]
  if(length(adj.end)>0){
    x[adj.end,]$mig.state <- "Mig"
  }
  
  x$ID_mig <- set.ID(x)
  x <- check_state_NC(x) 
  
  # return results
  
  if(plot==T){
    with(x, ScanTrack(study.local.timestamp, 
                      X,
                      Y, 
                      col=as.factor(mig.state))) 
  }
  
  if(print.dates==T){
    x %>% 
      group_by(mig.state, 
               ID_mig) %>% 
      summarize(min(date), 
                max(date))
  }
  
  
  return(x)
}
