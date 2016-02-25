setwd("U:/Actigraphy Raw Data (Marie McCarthy)/brushing teeth")

load("brushing.fingerprint.RDA")

require(stringr)
require(lubridate)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(zoo)
require(reshape2)
require(doParallel)
require(data.table)

test.df <- fread("./TAS1E35150309 (2016-01-21)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)
test.new <- fread("./TAS1E35150309 (2016-01-21)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)

#' format.df
#'
#' @param df - original data.frame to be formated 
#'
#' @return df - formated df
#' 
format.df <- function(df) {
  names(df) <- str_replace(names(df), " ", "")
  op <- options(digits.secs = 3)
  df$Timestamp <- dmy_hms(df$Timestamp, tz = "GMT")
  df$vector.dir <- NA
  df$vector.dir[df$AccelerometerX<0  & df$AccelerometerY<0  & df$AccelerometerZ<0]  <- 1 
  df$vector.dir[df$AccelerometerX<0  & df$AccelerometerY<0  & !df$AccelerometerZ<0] <- 2 
  df$vector.dir[df$AccelerometerX<0  & !df$AccelerometerY<0 & !df$AccelerometerZ<0] <- 3 
  df$vector.dir[!df$AccelerometerX<0 & !df$AccelerometerY<0 & !df$AccelerometerZ<0] <- 4 
  df$vector.dir[!df$AccelerometerX<0 & !df$AccelerometerY<0 & df$AccelerometerZ<0]  <- 5 
  df$vector.dir[!df$AccelerometerX<0 & df$AccelerometerY<0  & df$AccelerometerZ<0]  <- 6 
  df$vector.dir[df$AccelerometerX<0  & !df$AccelerometerY<0 & df$AccelerometerZ<0]  <- 7 
  df$vector.dir[!df$AccelerometerX<0 & df$AccelerometerY<0  & !df$AccelerometerZ<0] <- 8 
  df$vector.dir[df$AccelerometerX==0 & df$AccelerometerY==0 & df$AccelerometerZ==0] <- 9
  df$vector.dir <- as.factor(df$vector.dir)
  df <- df %>% mutate(vector.mag = sqrt(AccelerometerX^2+ AccelerometerY^2+ AccelerometerZ^2),
                      time_minute = floor_date(Timestamp, "minute"))  
  return(df)
}

test.df <- format.df(test.df)
test.new <- format.df(test.new)


#' Title
#'
#' @param v 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
get.peak.summary <- function(v, k = 10) {
  require(zoo)
  v.smooth <- rollapply(v, k, mean)
  switch.dir <- sapply(seq_along(v.smooth), 
                       function(x) ifelse(x<length(v.smooth),
                                          v.smooth[x]<v.smooth[x+1],
                                          F))
  peaks.per.sec <- sum(sapply(seq_along(v.smooth), 
                              function(x) ifelse(x<length(v.smooth), 
                                                 switch.dir[x] != switch.dir[x+1], 
                                                 F))) / (length(v.smooth)/100)
  if(peaks.per.sec == 0) return(F)
  
  period <- rollapply(which(sapply(seq_along(v.smooth), 
                                   function(x) ifelse(x<length(v.smooth), 
                                                      switch.dir[x] != switch.dir[x+1], 
                                                      F)) == T),
                      2, function(x) x[2] - x[1])
  
  avg.period <- mean(period)
  sd.period <- sd(period)
  return(c(peaks.per.sec = peaks.per.sec, avg.period = avg.period, sd.period = sd.period))
}

#' Title
#'
#' @param t 
#'
#' @return
#' @export
#'
#' @examples
peak.func <- function(t, test){
  tbl <- test %>% filter(time_minute == t) %>% 
    select(-Timestamp, -vector.dir, -time_minute) %>% 
    apply(.,2, get.peak.summary) %>% t() %>% data.frame %>% mutate(Timestamp = t) 
  if(nrow(tbl) == 4) {
    return(tbl)
  } else{
    return(data.frame(peaks.per.sec = rep(0,4), 
                      avg.period = rep(0,4), 
                      sd.period = rep(0,4), 
                      Timestamp  = rep(t,4)))
  }
}

#' Title
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
summary.df <- function(df){
  test <- df
  times <- unique(test$time_minute)
  
  test.summary <- select(test, -vector.dir, -time_minute) %>%
    melt(id.vars = "Timestamp") %>%
    mutate(Timestamp = floor_date(Timestamp, "minute")) %>%
    group_by(Timestamp, variable) %>%
    summarise(min = min(value), 
              Qu1 = quantile(value, .25),
              Median = median(value),
              Mean = mean(value),
              Qu2 = quantile(value, .75),
              Max = max(value)) %>%
    left_join(select(test, Timestamp, vector.dir, -time_minute) %>% 
                filter(vector.dir == 5) %>%
                mutate(Timestamp = floor_date(Timestamp, "minute")) %>% 
                group_by(Timestamp) %>% summarise(per.vec.dir5 = n()/6000),
              by = "Timestamp"
    )
  tmp <- do.call(rbind, lapply(times, function(t) peak.func(t, test)))
  test.summary <- cbind(test.summary, tmp %>% select(-Timestamp)) 
  return(test.summary)
}

test.df.summary <- summary.df(test.df)
test.new.summary <- summary.df(test.new)



#' Title
#'
#' @param brushing.fingerprint 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
similarity <- function(brushing.fingerprint, d){
  return(
    sqrt(sum(((brushing.fingerprint$Qu1st-d$Qu1)/brushing.fingerprint$Qu1st)^2)) +
      sqrt(sum(((brushing.fingerprint$Median-d$Median)/brushing.fingerprint$Median)^2)) +
      sqrt(sum(((brushing.fingerprint$Mean-d$Mean)/brushing.fingerprint$Mean)^2)) +
      sqrt(sum(((brushing.fingerprint$X3rd.Qu.-d$Qu2)/brushing.fingerprint$X3rd.Qu.)^2)) +
      sqrt(sum(((brushing.fingerprint$mean.dir5[4]-max(d$per.vec.dir5[4], 0, na.rm = T))/brushing.fingerprint$mean.dir5[4])^2)) +
      sqrt(sum(((brushing.fingerprint$peaks.per.sec - d$peaks.per.sec)/brushing.fingerprint$peaks.per.sec)^2)) +
      sqrt(sum(((brushing.fingerprint$avg.period - d$avg.period)/brushing.fingerprint$avg.period)^2)) +
      sqrt(sum(((brushing.fingerprint$sd.period - d$sd.period)/brushing.fingerprint$sd.period)^2)) 
  )
}


#' Title
#'
#' @param df 
#' @param brushing.fingerprint 
#'
#' @return
#' @export
#'
#' @examples
get.sim.results <- function(org.df, df, brushing.fingerprint, similarity){
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  times <- unique(df$Timestamp)
  
  sim.results <- foreach(t = seq_along(times), .combine = c) %dopar% {
    return(sim = similarity(as.data.frame(brushing.fingerprint), df[df$Timestamp == times[t], ]))
  }
  stopCluster(cl)

  sim.results <- data.frame(times = times, sim = sim.results)
  sim.results$event <- ifelse(sim.results$sim < 2, 1, 0) 

  counter <- 0
  for (i in 2:length(sim.results$event)) {
    if(sim.results$event[i] == 1 & sim.results$event[i-1] == 0) counter <- counter + 1
    if(sim.results$event[i] == 1 ) sim.results$event[i] <- counter
  }
  message(paste(counter, " events of brushing teeth have been identified"))

  result <- left_join(org.df, 
                      sim.results, 
                      by = c("time_minute" = "times"))

  return(result)
}

test.df <- get.sim.results(org.df = test.df, 
                df = test.df.summary, 
                brushing.fingerprint, 
                similarity)


#' Title
#'
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
plot.similarity <- function(df, scale.time = F){
  d <- filter(df, event > 0)
  
  if(scale.time){
    start = min(df$time_minute)
    end = max(df$time_minute)
  } else {
    start = min(d$time_minute)
    end = max(d$time_minute)
  }
  
  p0 <- ggplot(d, aes(Timestamp, as.factor(ifelse(sim<2,1,0)), group = event)) + 
    geom_point(col = "red", size = 2) + 
    labs(y="") + 
    scale_x_datetime(limits = c(start, end)) 
  p1 <- ggplot(d, aes(Timestamp, vector.mag, group = event)) + geom_point(aes(col = vector.dir)) + 
    geom_line() + theme(legend.position="none") + labs(title = "Overall Vector Mag") + 
    scale_x_datetime(limits = c(start, end))
  p2 <- ggplot(d, aes(Timestamp, as.numeric(as.character(vector.dir)), group = event)) + 
    geom_line() + labs(title = "Vector Dir") + 
    scale_x_datetime(limits = c(start, end))
  p3 <- ggplot(d, aes(Timestamp, AccelerometerX, group = event)) + geom_line() + 
    labs(title = "x Accel") + 
    scale_x_datetime(limits = c(start, end))
  p4 <- ggplot(d, aes(Timestamp, AccelerometerY, group = event)) + geom_line() + 
    labs(title = "Y Accel") + 
    scale_x_datetime(limits = c(start, end))
  p5 <- ggplot(d, aes(Timestamp, AccelerometerZ, group = event)) + geom_line() + 
    labs(title = "Z Accel") + 
    scale_x_datetime(limits = c(start, end))
  
  grid.arrange(p0, p1, p2, p3, p4, p5, ncol = 1)
}

plot.similarity(test.df, scale.time = T)


