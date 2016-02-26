setwd("U:/Actigraphy Raw Data (Marie McCarthy)/brushing teeth")

load("brushing.fingerprint.RDA")
load("brushing.fingerprint.sd.RDA")

# calculation of fingerprint didn't consider freq of device
brushing.fingerprint$avg.period <- brushing.fingerprint$avg.period / 100
brushing.fingerprint$sd.period <- brushing.fingerprint$sd.period / 100

require(stringr)
require(lubridate)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(zoo)
require(reshape2)
require(doParallel)
require(data.table)
require(caret)

org.df <- fread("TAS1E35150309 (2016-01-21)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)
# marie's data
test.df <- fread("c:/users/smithlou/desktop/TAS1E31150005 (2016-02-26)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)
# louis data
test.new <- fread("c:/users/smithlou/desktop/TAS1E31150003 (2016-02-26)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)

#' format.df
#'
#' @param df - raw actigraphy data to be formated, without metadata 
#'
#' @return df - formated df with formated time, vector.mag & vector.dir
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

# format the data
org.df <- format.df(org.df)
test.df <- format.df(test.df)
test.new <- format.df(test.new)


#' get.peak.summary
#'
#' @param v - vector to be smoothed / considered
#' @param k - rolling average
#' @param freq - frequency device was recording at
#'
#' @return vector summarising the peak for v
#' 
get.peak.summary <- function(v, k, freq) {
  require(zoo)
  v.smooth <- rollapply(v, k, mean)
  switch.dir <- sapply(seq_along(v.smooth), 
                       function(x) ifelse(x<length(v.smooth),
                                          v.smooth[x]<v.smooth[x+1],
                                          F))
  peaks.per.sec <- sum(sapply(seq_along(v.smooth), 
                              function(x) ifelse(x<length(v.smooth), 
                                                 switch.dir[x] != switch.dir[x+1], 
                                                 F))) / (length(v.smooth)/freq)
  if(peaks.per.sec == 0) return(F)
  
  period <- rollapply(which(sapply(seq_along(v.smooth), 
                                   function(x) ifelse(x<length(v.smooth), 
                                                      switch.dir[x] != switch.dir[x+1], 
                                                      F)) == T),
                      2, function(x) (x[2] - x[1])/freq)
  
  avg.period <- mean(period)
  sd.period <- sd(period)
  return(c(peaks.per.sec = peaks.per.sec, avg.period = avg.period, sd.period = sd.period))
}

#' peak.func
#'
#' @param t time period to calculate peak summary for
#' @param df dataset to pull vector from
#' @param freq frequency device recording at 
#' @param k window size for smoothing
#'
#' @return data.frame of peak summary for period t
#' 
peak.func <- function(t, df, freq, k){
  tbl <- df %>% filter(time_minute == t) %>% 
    select(-Timestamp, -vector.dir, -time_minute) %>% 
    apply(.,2, function(x) get.peak.summary(v = x, k = k, freq = freq)) %>% 
    t() %>% data.frame %>% mutate(Timestamp = t) 
  if(nrow(tbl) == 4) {
    return(tbl)
  } else{
    return(data.frame(peaks.per.sec = rep(0,4), 
                      avg.period = rep(0,4), 
                      sd.period = rep(0,4), 
                      Timestamp  = rep(t,4)))
  }
}

#' summary.df
#'
#' @param df data.frame to summarise
#' @param freq frequency device recording at
#' @param k window size of smoothing
#'
#' @export summary data.frame
#' 
summary.df <- function(df, freq, k=10){
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
                group_by(Timestamp) %>% summarise(per.vec.dir5 = n()/(60 * freq)),
              by = "Timestamp"
    )
  tmp <- do.call(rbind, lapply(times, function(t) peak.func(t, test, freq, k)))
  test.summary <- cbind(test.summary, tmp %>% select(-Timestamp)) 
  return(test.summary)
}

org.df.summary   <- summary.df(org.df, freq = 100, k=10)          # org data used to generate fingerprint
test.df.summary  <- summary.df(df = test.df, freq = 100, k=10)    # marie
test.new.summary <- summary.df(test.new, freq = 100, k=10)        # louis




#' similarity
#'
#' @param brushing.fingerprint 
#' @param d a one minute summary of actigraphy data for comparision with the fingerprint
#'
#' @return a measure of similarity - euclidean distance from fingerprint
#' 
similarity.euclidean <- function(brushing.fingerprint, d){
  return(
      # sqrt(sum(((brushing.fingerprint$Qu1st-d$Qu1)/brushing.fingerprint$Qu1st)^2)) +
      sqrt(sum(((brushing.fingerprint$Median-d$Median)/brushing.fingerprint$Median)^2)) +
      sqrt(sum(((brushing.fingerprint$Mean-d$Mean)/brushing.fingerprint$Mean)^2)) +
      # sqrt(sum(((brushing.fingerprint$Qu3rd-d$Qu2)/brushing.fingerprint$Qu3rd)^2)) +
      # sqrt(sum(((brushing.fingerprint$mean.dir5[4]-max(d$per.vec.dir5[4], 0, na.rm = T))/brushing.fingerprint$mean.dir5[4])^2)) +
      sqrt(sum(((brushing.fingerprint$avg.period - d$avg.period)/brushing.fingerprint$avg.period)^2)) +
      # sqrt(sum(((brushing.fingerprint$sd.period - d$sd.period)/brushing.fingerprint$sd.period)^2)) +
      sqrt(sum(((brushing.fingerprint$peaks.per.sec - d$peaks.per.sec)/brushing.fingerprint$peaks.per.sec)^2))
  )
}

#' similarity.boolean
#'
#' @param brushing.fingerprint 
#' @param brushing.fingerprint.sd 
#' @param sigma 
#' @param d 
#'
#' @return
#' 
similarity.boolean <- function(brushing.fingerprint, brushing.fingerprint.sd, sigma = 1, d){
  tmp <- brushing.fingerprint$Mean + (sigma * brushing.fingerprint.sd$Mean * c(-1,1))
  close.mean <- d$Mean >= tmp[1] & d$Mean <= tmp[2]
  tmp <- brushing.fingerprint$Median + (sigma * brushing.fingerprint.sd$Median * c(-1,1))
  close.median <- d$Median >= tmp[1] & d$Median <= tmp[2]
  # don't have sd for peaks per second of teeth brushing events
  tmp <- brushing.fingerprint$peaks.per.sec + 
    (sigma * brushing.fingerprint.sd$peak.per.sec.sd * c(-1,1))
  close.peak.rate <- d$peaks.per.sec >= tmp[1] & d$peaks.per.sec <= tmp[2]
  
  tmp <- brushing.fingerprint$avg.period + (sigma * brushing.fingerprint$sd.period * c(-1,1))
  close.peak.period <- d$avg.period >= tmp[1] & d$avg.period <= tmp[2]
  
  return(c(close.mean, close.median, close.peak.rate, close.peak.period))
}

#' get.sim.results
#'
#' @param raw.df -
#' @param summary.df 
#' @param brushing.fingerprint 
#' @param similarity 
#' @param close
#'
#' @return
#' 
get.sim.results <- function(raw.df, summary.df, 
                            brushing.fingerprint,
                            brushing.fingerprint.sd,
                            similarity.euclidean.m = similarity.euclidean, 
                            similarity.boolean.m = similarity.boolean,
                            sigma = 1,
                            close = 2){
    cores <- detectCores()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    times <- unique(summary.df$Timestamp)
    sim.results.e <- foreach(t = seq_along(times), .combine = c) %dopar% {
      return(sim = similarity.euclidean.m(as.data.frame(brushing.fingerprint)[4,], 
                                        summary.df[summary.df$Timestamp == times[t], ][4,]))
    }
    sim.results.b <- foreach(t = seq_along(times), .combine = rbind) %dopar% {
      return(sim = similarity.boolean.m(brushing.fingerprint = as.data.frame(brushing.fingerprint)[4,], 
                                      brushing.fingerprint.sd = as.data.frame(brushing.fingerprint.sd)[4,], 
                                      sigma = sigma, 
                                      d = summary.df[summary.df$Timestamp == times[t], ][4,]))
    }
    stopCluster(cl)
  
    sim.results <- data.frame(times = times, sim.e = sim.results.e, sim.b = sim.results.b)
    
    sim.results$event.e <- ifelse(sim.results$sim.e < close, 1, 0)
    sim.results$event.b <- ifelse(sim.results %>% select(sim.b.1:sim.b.4) %>% apply(., 1, sum) == 4, 1, 0)
 
    sim.results %>% filter(event.e == 1 )
    sim.results %>% filter(event.b == 1)
    sim.results %>% filter(event.e == 1 & event.b == 1)
    
    sim.results$event <- ifelse(sim.results$event.e == 1 &
                                  sim.results$event.b == 1, 
                                1,0)

  counter <- 0
  for (i in 2:length(sim.results$event)) {
    if(sim.results$event[i] == 1 & sim.results$event[i-1] == 0) counter <- counter + 1
    if(sim.results$event[i] == 1 ) sim.results$event[i] <- counter
  }
  message(paste(counter, " events of brushing teeth have been identified"))

  result <- left_join(raw.df, 
                      sim.results, 
                      by = c("time_minute" = "times"))

  return(result)
}




# save.image("test_brushing_pattern.RDATA")

org.result <- get.sim.results(raw.df = org.df,summary.df = org.df.summary, 
                              brushing.fingerprint = brushing.fingerprint, 
                              brushing.fingerprint.sd = brushing.fingerprint.sd, 
                              sigma = 3.5, close = 0.11)
org.result %>% select(time_minute, sim.e,  event) %>% filter(event > 0) %>% distinct()

# 100% on training dataset


test.result <- get.sim.results(raw.df = test.df, summary.df = test.df.summary, 
                               brushing.fingerprint = brushing.fingerprint, 
                               brushing.fingerprint.sd = brushing.fingerprint.sd, 
                               sigma = 3.5, close = 0.11)
test.result %>% select(time_minute, sim.e,  event) %>% filter(event > 0) %>% distinct()

# 66% and 6 FP's

new.result <- get.sim.results(raw.df = test.new, summary.df = test.new.summary, 
                              brushing.fingerprint = brushing.fingerprint, 
                              brushing.fingerprint.sd = brushing.fingerprint.sd, 
                              sigma = 3.5, close = 0.11)
new.result %>% select(time_minute, sim.e,  event) %>% filter(event > 0) %>% distinct()


###################################
#
# missing times in marie & louis data

ggplot(test.df %>% filter(Timestamp > ymd_hms("20160225 171835"),
                          Timestamp < ymd_hms("20160225 172200")),
       aes(Timestamp, vector.mag)) + geom_line()
test.result %>% filter(Timestamp >= ymd_hms("20160225 171800"),
                       Timestamp <= ymd_hms("20160225 172200")) %>%
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

ggplot(test.df %>% filter(Timestamp > ymd_hms("20160225 230115"),
                          Timestamp < ymd_hms("20160225 230500")),
       aes(Timestamp, vector.mag)) + geom_line()
test.result %>% filter(Timestamp >= ymd_hms("20160225 230100"),
                       Timestamp <= ymd_hms("20160225 230500")) %>%
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

ggplot(test.df %>% filter(Timestamp > ymd_hms("20160226 072930"),
                          Timestamp < ymd_hms("20160226 073230")),
       aes(Timestamp, vector.mag)) + geom_line()
test.result %>% filter(Timestamp >= ymd_hms("20160226 072900"),
                       Timestamp <= ymd_hms("20160226 073300")) %>%
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

test.minutes <- c(ymd_hms("20160225 171900"), ymd_hms("20160225 172000"), ymd_hms("20160225 172100"),
                  ymd_hms("20160225 230100"), ymd_hms("20160225 230200"), ymd_hms("20160225 230300"),ymd_hms("20160225 230400"),
                  ymd_hms("20160226 072900"), ymd_hms("20160226 073000"), ymd_hms("20160226 073100"), ymd_hms("20160226 073200"))

confusionMatrix(ifelse(unique(test.result$time_minute) %in% test.minutes, 1, 0),
ifelse(unlist(test.result %>% select(time_minute, event) %>% distinct() %>% select(event)) > 0, 1, 0)
)

# louis plots
ggplot(test.new %>% filter(Timestamp > ymd_hms("20160225 214200"),
                          Timestamp < ymd_hms("20160225 214500")),
       aes(Timestamp, vector.mag)) + geom_line()
test.result %>% filter(Timestamp >= ymd_hms("20160225 214200"),
                       Timestamp <= ymd_hms("20160225 214500")) %>%
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

ggplot(test.new %>% filter(Timestamp > ymd_hms("20160226 080800"),
                          Timestamp < ymd_hms("20160226 081100")),
       aes(Timestamp, vector.mag)) + geom_line()
test.result %>% filter(Timestamp >= ymd_hms("20160226 080800"),
                       Timestamp <= ymd_hms("20160226 081100")) %>%
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

test.result %>% filter(event.e ==1) %>% 
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

test.result %>% filter(event.b ==1) %>% 
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

test.result %>% filter(event ==1) %>% 
  select(time_minute, sim.e, sim.b.1, sim.b.2, sim.b.3, sim.b.4, event.e, event.b,  event) %>%
  distinct

new.minutes <- c(ymd_hms("20160225 214200"), ymd_hms("20160225 214300"), ymd_hms("20160225 214400"),
                  ymd_hms("20160226 080800"), ymd_hms("20160226 080900"), ymd_hms("20160226 081000"))

confusionMatrix(ifelse(unique(new.result$time_minute) %in% new.minutes, 1, 0),
                ifelse(unlist(new.result %>% select(time_minute, event) %>% distinct() %>% select(event)) > 0, 1, 0)
)


