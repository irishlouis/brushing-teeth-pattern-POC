setwd("U:/Actigraphy Raw Data (Marie McCarthy)/brushing teeth")

load("brushing summary data.RDATA")

#################################################################

require(data.table)
df <- fread("TAS1E35150309 (2016-01-21)RAW.csv", stringsAsFactors = F, skip = 10, header = T, data.table = F)
head(df)

require(stringr)
names(df) <- str_replace(names(df), " ", "")

require(lubridate)

op <- options(digits.secs = 3)
df$Timestamp <- dmy_hms(df$Timestamp, tz = "GMT")

class(df$Timestamp)

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

class(df)

summary(df)

require(dplyr)

df <- df %>% mutate(vector.mag = sqrt(AccelerometerX^2+ AccelerometerY^2+ AccelerometerZ^2),
                    time_minute = floor_date(Timestamp, "minute"))

require(ggplot2)
require(gridExtra)

# subsets for when teeth brushing was occuring
df1 <- df %>% 
  filter(Timestamp > dmy_hms("20/01/2016 170115", tz = "GMT") & Timestamp < dmy_hms("20/01/2016 170445", tz = "GMT") ) 

df2 <- df %>% 
  filter(Timestamp > dmy_hms("20/01/2016 180505", tz = "GMT") & Timestamp < dmy_hms("20/01/2016 180915", tz = "GMT") ) 

df3 <- df %>% 
  filter(Timestamp > dmy_hms("20/01/2016 182700", tz = "GMT") & Timestamp < dmy_hms("20/01/2016 183050", tz = "GMT") )

df4 <- df %>% 
  filter(Timestamp > dmy_hms("20/01/2016 215515", tz = "GMT") & Timestamp < dmy_hms("20/01/2016 215800", tz = "GMT") )

plot.data <- function(d){
  p1 <- ggplot(d, aes(Timestamp, vector.mag)) + geom_point(aes(col = vector.dir)) + 
    geom_line() + theme(legend.position="none") + labs(title = "Overall Vector Mag")
  p <- ggplot(d, aes(Timestamp, as.numeric(as.character(vector.dir)))) + geom_line() + labs(title = "Vector Dir")
  p2 <- ggplot(d, aes(Timestamp, AccelerometerX)) + geom_line() + labs(title = "X Accel Data")
  p3 <- ggplot(d, aes(Timestamp, AccelerometerY)) + geom_line() + labs(title = "Y Accel Data")
  p4 <- ggplot(d, aes(Timestamp, AccelerometerZ)) + geom_line() + labs(title = "Z Accel Data")
  
  grid.arrange(p1, p, p2, p3, p4, ncol = 1)
  }

plot.data(df1)
plot.data(df2)
plot.data(df3)
plot.data(df4)

## zoom in on 30 seconds of brushing activity
plot.data(df %>% 
            filter(Timestamp > dmy_hms("20/01/2016 180730", tz = "GMT") & 
                     Timestamp < dmy_hms("20/01/2016 180800", tz = "GMT") ))

## zoom in on 30 seconds of non brushing activity
plot.data(df %>% 
            filter(Timestamp > dmy_hms("20/01/2016 181530", tz = "GMT") & 
                     Timestamp < dmy_hms("20/01/2016 181600", tz = "GMT") ))

## zoom in on 30 seconds of non brushing activity
plot.data(df %>% 
            filter(Timestamp > dmy_hms("20/01/2016 201530", tz = "GMT") & 
                     Timestamp < dmy_hms("20/01/2016 201600", tz = "GMT") ))

## zoom in on 30 seconds of non brushing activity
plot.data(df %>% 
            filter(Timestamp > dmy_hms("20/01/2016 214530", tz = "GMT") & 
                     Timestamp < dmy_hms("20/01/2016 214600", tz = "GMT") ))

require(zoo)
summ.df.mean <- function(d) {
  splits <- rollapply(1:nrow(d), 6000, function(x) x )
  splits <- as.data.frame(splits)
  d <- d %>% select( -Timestamp, -time_minute, -vector.dir)
  tmp <- apply(splits, 1, function(x) {apply(d[x, ], 2, summary)})
  tmp <- as.data.frame(tmp)
  return(matrix(apply(tmp, 1, mean), 
                nrow = 6, 
                dimnames = list(rownames = c("Min.","Qu1st","Median", "Mean","Qu3rd","Max"),
                                colnames = c("AccelerometerX","AccelerometerY","AccelerometerZ","vector.mag"))))
}

summ.df.sd <- function(d) {
  splits <- rollapply(1:nrow(d), 6000, function(x) x  )
  splits <- as.data.frame(splits)
  d <- d %>% select( -Timestamp, -time_minute, -vector.dir)
  tmp <- apply(splits, 1, function(x) {apply(d[x, ], 2, summary)})
  tmp <- as.data.frame(tmp)
  return(matrix(apply(tmp, 1, sd), 
                nrow = 6, 
                dimnames = list(rownames = c("Min.","Qu1st","Median", "Mean","Qu3rd","Max"),
                                colnames = c("AccelerometerX","AccelerometerY","AccelerometerZ","vector.mag"))))
}

brushing.summary <- lapply(list(df1, df2, df3, df4), summ.df.mean)

brushing.summary

brushing.summary.dir <- lapply(list(df1, df2, df3, df4), function(x) round(table(x$vector.dir)/nrow(x), 2))
brushing.summary.dir

brushing.fingerprint <- data.frame( (brushing.summary[[1]] + 
                                       brushing.summary[[2]]+ 
                                       brushing.summary[[3]] +
                                       brushing.summary[[4]]) / 4)

brushing.fingerprint.sd <- lapply(list(df1, df2, df3, df4), summ.df.sd)
brushing.fingerprint.sd <- data.frame( (brushing.fingerprint.sd[[1]] + 
                                          brushing.fingerprint.sd[[2]]+ 
                                          brushing.fingerprint.sd[[3]] + 
                                          brushing.fingerprint.sd[[4]]) / 4)
brushing.fingerprint.sd

brushing.fingerprint <- cbind(t(brushing.fingerprint), 
                              mean.dir5 = c(0,0,0,mean(sapply(brushing.summary.dir, function(x) x[5]))),
                              sd.dir5 = c(0,0,0, sd(sapply(brushing.summary.dir, function(x) x[5]))))
brushing.fingerprint

brushing.fingerprint.sd  <- cbind(t(brushing.fingerprint.sd), 
                               sd.dir5 = c(0,0,0, sd(sapply(brushing.summary.dir, function(x) x[5]))))

brushing.fingerprint.sd

# sin wave fingerprint pattern
get.peak.summary <- function(v, k = 10) {
  require(zoo)
  v.smooth <- rollapply(v, k, mean)
  switch.dir <- sapply(seq_along(v.smooth), function(x) ifelse(x<length(v.smooth), 
                                                                     v.smooth[x]<v.smooth[x+1], 
                                                                     F))
  peaks.per.sec <- sum(sapply(seq_along(v.smooth), function(x) ifelse(x<length(v.smooth), 
                                                                            switch.dir[x] != switch.dir[x+1], 
                                                                            F))) / (length(v.smooth)/100)
  if(peaks.per.sec == 0) return(F)
  
  period <- rollapply(which(sapply(seq_along(v.smooth), function(x) ifelse(x<length(v.smooth), 
                                                                                 switch.dir[x] != switch.dir[x+1], 
                                                                                 F)) == T),
                      2, function(x) x[2] - x[1])
  
  avg.period <- mean(period)
  sd.period <- sd(period)
  return(c(peaks.per.sec = peaks.per.sec, avg.period = avg.period, sd.period = sd.period))
}

peak.summary <- lapply(list(df1, df2, df3, df4), function(x)
  (x %>% select(-Timestamp, -vector.dir, -time_minute) %>% apply(.,2, get.peak.summary) %>% t())
)

peak.summary
peak.summary.averages <- Reduce('+', peak.summary) / length(peak.summary)
peak.summary.averages

brushing.fingerprint <- cbind(brushing.fingerprint, peak.summary.averages)

save.image(file = "brushing summary data.RDATA")

########################################################

test <- df 

head(test)

require(reshape2)

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

head(test.summary)

times <- unique(test$time_minute)

require(doParallel)

tmp.func <- function(t){
  print(t)
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

tmp <- do.call(rbind, lapply(times, tmp.func ))

test.summary <- cbind(test.summary, tmp %>% select(-Timestamp)) 

head(test.summary)

####################


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

times <- unique(test.summary$Timestamp)
require(doParallel)

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

sim.results <- foreach(t = seq_along(times), .combine = c) %dopar% {
  return(sim = similarity(as.data.frame(brushing.fingerprint), test.summary[test.summary$Timestamp == times[t], ]))
}

stopCluster(cl)

sim.results <- data.frame(times = times, sim = sim.results)

head(sim.results)

summary(sim.results)

sim.results$event <- ifelse(sim.results$sim < 2, 1, 0) 

counter <- 0
for (i in 2:length(sim.results$event)) {
  if(sim.results$event[i] == 1 & sim.results$event[i-1] == 0) counter = counter +1
  if(sim.results$event[i] == 1 ) sim.results$event[i] <- counter
}

message(paste(counter, " events of brushing teeth have been identified"))

test.results <- left_join(test, 
                          sim.results, 
                          by = c("time_minute" = "times"))

test.results <- left_join(test.results, select(sim.results, time, event), by = c("ts" = "time"))

plot.similarity <- function(d){
  d <- filter(d, event > 0)
  p0 <- ggplot(d, aes(Timestamp, as.factor(ifelse(sim<2,1,0)), group = event)) + geom_point(col = "red", size = 2) + labs(y="")
  p1 <- ggplot(d, aes(Timestamp, vector.mag, group = event)) + geom_point(aes(col = vector.dir)) + 
    geom_line() + theme(legend.position="none") + labs(title = "Overall Vector Mag")
  p2 <- ggplot(d, aes(Timestamp, as.numeric(as.character(vector.dir)), group = event)) + geom_line() + labs(title = "Vector Dir")
  p3 <- ggplot(d, aes(Timestamp, AccelerometerX, group = event)) + geom_line() + labs(title = "x Accel")
  p4 <- ggplot(d, aes(Timestamp, AccelerometerY, group = event)) + geom_line() + labs(title = "Y Accel")
  p5 <- ggplot(d, aes(Timestamp, AccelerometerZ, group = event)) + geom_line() + labs(title = "Z Accel")
  
  grid.arrange(p0, p1, p2, p3, p4, p5, ncol = 1)
}

plot.similarity(test.results)


###################################################################

# direction of accel
tmp <- df %>%filter(Timestamp > dmy_hms("20/01/2016 180735", tz = "GMT") & Timestamp < dmy_hms("20/01/2016 180740", tz = "GMT") )
plot.data(tmp)
switch.dir <- sapply(seq_along(tmp$Timestamp), function(x) ifelse(x<nrow(tmp), tmp$vector.mag[x]<tmp$vector.mag[x+1], F))
# sapply(seq_along(tmp$Timestamp), function(x) ifelse(x<nrow(tmp), switch.dir[x] != switch.dir[x+1], F))
sum(sapply(seq_along(tmp$Timestamp), function(x) ifelse(x<nrow(tmp), switch.dir[x] != switch.dir[x+1], F)))
### finds nearly twice as many

## need to smooth curve
ggplot(data.frame(t = 1:490, v =rollmean(tmp$vector.mag, 10)), aes(t, v)) + geom_line()
tmp2 <- rollmean(tmp$vector.mag, 10)
switch.dir.smooth <- sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), tmp2[x]<tmp2[x+1], F))
# sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), switch.dir.smooth[x] != switch.dir.smooth[x+1], F))
sum(sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), switch.dir.smooth[x] != switch.dir.smooth[x+1], F)))
### finds correct number now

# find period for wave
rollapply(which(sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), switch.dir.smooth[x] != switch.dir.smooth[x+1], F)) == T),
          2, function(x) x[2] - x[1])
summary(rollapply(which(sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), switch.dir.smooth[x] != switch.dir.smooth[x+1], F)) == T),
                  2, function(x) x[2] - x[1]))
sd(rollapply(which(sapply(seq_along(tmp2), function(x) ifelse(x<length(tmp2), switch.dir.smooth[x] != switch.dir.smooth[x+1], F)) == T),
                  2, function(x) x[2] - x[1]))


##############################

test$minute <- floor_date(test$Timestamp, "minute")

get.peaks <- function(vec.mag, k = 10) {
  require(zoo)
  vec.mag.smooth <- rollapply(vec.mag, k, mean)
  switch.dir <- sapply(seq_along(vec.mag.smooth), function(x) ifelse(x<length(vec.mag.smooth), 
                                                                     vec.mag.smooth[x]<vec.mag.smooth[x+1], 
                                                                     F))
  peaks.per.sec <- sum(sapply(seq_along(vec.mag.smooth), function(x) ifelse(x<length(vec.mag.smooth), 
                                                           switch.dir[x] != switch.dir[x+1], 
                                                           F))) / 60
  if(peaks.per.sec == 0) return(F)
  
  period <- rollapply(which(sapply(seq_along(vec.mag.smooth), function(x) ifelse(x<length(vec.mag.smooth), 
                                                                       switch.dir[x] != switch.dir[x+1], 
                                                                       F)) == T),
                      2, function(x) x[2] - x[1])
  
  avg.period <- mean(period)
  sd.period <- sd(period)
  

  mdl <- lm(vec.mag ~ I(sin(pi*2.07*seq_along(vec.mag))))
  p.value <- 1 - pf(summary(mdl)$fstatistic[1], summary(mdl)$fstatistic[2], summary(mdl)$fstatistic[3]) 
  
  if(p.value < .000001 & peaks.per.sec > 7 & 
     peaks.per.sec < 10 &
     sd.period > 3 &
     sd.period < 5.7 &
     avg.period > 11 &
     avg.period < 17
     ) return(T) else return(F)
  }

require(doParallel)
times <- unique(test$minute)
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

sim.results2 <- foreach(t = seq_along(times), .combine = c) %dopar% {
  return(sim = get.peaks(vec.mag = test[test$minute == times[t], "vector.mag"], k = 10))
}

stopCluster(cl)

data.frame(times, sim.results2) %>% filter(sim.results2)














## fit sin wave
tmp$row <- 1:499
mdl <- lm(tmp$vector.mag ~ I(sin(pi*2.07*tmp$row)) + I(cos(pi*2.07*tmp$row)))
summary(mdl)
plot(tmp$Timestamp, tmp$vector.mag, type = "l" )
par(new=T)
plot(data.frame(Timestamp = tmp$Timestamp,
                vector.mag=mdl$coefficients[1] + mdl$coefficients[2]*sin(2.07*pi*1:499)), type = "l", col ="red")

plot(data.frame(Timestamp = tmp$Timestamp,
                vector.mag=mdl$coefficients[1] + mdl$coefficients[2]*sin(2.07*pi*1:499) + mdl$coefficients[3]*cos(2.07*pi*1:499)), type = "l", col ="red")

# need to smooth out curve, to remove minor bumps