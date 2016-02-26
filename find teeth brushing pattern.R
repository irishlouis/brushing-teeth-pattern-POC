setwd("U:/Actigraphy Raw Data (Marie McCarthy)/brushing teeth")

#################################################################
require(stringr)
require(lubridate)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(zoo)
require(reshape2)
require(doParallel)
require(data.table)
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

brushing.fingerprint <- Reduce('+', brushing.summary) / length(brushing.summary)
  
brushing.fingerprint.sd <- lapply(list(df1, df2, df3, df4), summ.df.sd)
brushing.fingerprint.sd <- Reduce('+', brushing.fingerprint.sd) / length(brushing.fingerprint.sd)
brushing.fingerprint.sd

brushing.fingerprint <- cbind(t(brushing.fingerprint), 
                              mean.dir5 = c(0,0,0,mean(sapply(brushing.summary.dir, function(x) x[5]))),
                              sd.dir5 = c(0,0,0, sd(sapply(brushing.summary.dir, function(x) x[5]))))
brushing.fingerprint

save(brushing.fingerprint, file = "brushing.fingerprint.RDA" )

brushing.fingerprint.sd  <- cbind(t(brushing.fingerprint.sd), 
                               sd.dir5 = c(0,0,0, sd(sapply(brushing.summary.dir, function(x) x[5]))))

brushing.fingerprint.sd

save(brushing.fingerprint.sd, file = "brushing.fingerprint.sd.RDA" )

# sin wave fingerprint pattern
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
  
  # mdl <- lm(v ~ I(sin(pi*2.07*seq_along(v))))
  # pvalue <- 1 - pf(summary(mdl)$fstatistic[1], summary(mdl)$fstatistic[2], summary(mdl)$fstatistic[3])
  
  return(c(peaks.per.sec = peaks.per.sec, avg.period = avg.period, sd.period = sd.period)) #, pvalue = pvalue))
}

peak.summary <- lapply(list(df1, df2, df3, df4), function(x)
  (x %>% select(-Timestamp, -vector.dir, -time_minute) %>% apply(.,2, get.peak.summary) %>% t())
)

peak.summary
peak.summary.averages <- Reduce('+', peak.summary) / length(peak.summary)
peak.summary.averages

peak.per.sec.sd <- sd(c(peak.summary[[1]][4,1], peak.summary[[2]][4,1], peak.summary[[3]][4,1], peak.summary[[4]][4,1]))

brushing.fingerprint <- cbind(brushing.fingerprint, peak.summary.averages)
brushing.fingerprint.sd <-cbind(brushing.fingerprint.sd, peak.per.sec.sd = peak.per.sec.sd)
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

plot.similarity <- function(d){
  d <- filter(d, event > 0)
  p0 <- ggplot(d, aes(Timestamp, as.factor(ifelse(sim<2,1,0)), group = event)) + 
    geom_point(col = "red", size = 2) + 
    labs(y="") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000"))) 
    
  p1 <- ggplot(d, aes(Timestamp, vector.mag, group = event)) + geom_point(aes(col = vector.dir)) + 
    geom_line() + theme(legend.position="none") + labs(title = "Overall Vector Mag") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000")))
  p2 <- ggplot(d, aes(Timestamp, as.numeric(as.character(vector.dir)), group = event)) + 
    geom_line() + labs(title = "Vector Dir") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000")))
  p3 <- ggplot(d, aes(Timestamp, AccelerometerX, group = event)) + geom_line() + 
    labs(title = "x Accel") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000")))
  p4 <- ggplot(d, aes(Timestamp, AccelerometerY, group = event)) + geom_line() + 
    labs(title = "Y Accel") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000")))
  p5 <- ggplot(d, aes(Timestamp, AccelerometerZ, group = event)) + geom_line() + 
    labs(title = "Z Accel") + 
    scale_x_datetime(limits = c(ymd_hms("20160120 155000"), ymd_hms("20160120 220000")))
  
  grid.arrange(p0, p1, p2, p3, p4, p5, ncol = 1)
}

plot.similarity(test.results)

###################################################################
#
# evaluation confusion matrix

require(caret)
brushing.minutes <- c(dmy_hms("20/01/2016 170100", tz = "GMT"),
                      dmy_hms("20/01/2016 170200", tz = "GMT"),
                      dmy_hms("20/01/2016 170300", tz = "GMT"),
                      dmy_hms("20/01/2016 180500", tz = "GMT"),
                      dmy_hms("20/01/2016 180600", tz = "GMT"),
                      dmy_hms("20/01/2016 180700", tz = "GMT"),
                      dmy_hms("20/01/2016 180800", tz = "GMT"),
                      dmy_hms("20/01/2016 182700", tz = "GMT"),
                      dmy_hms("20/01/2016 182800", tz = "GMT"),
                      dmy_hms("20/01/2016 182900", tz = "GMT"),
                      dmy_hms("20/01/2016 215500", tz = "GMT"),
                      dmy_hms("20/01/2016 215600", tz = "GMT"),
                      dmy_hms("20/01/2016 215700", tz = "GMT"))
times <- unique(test$time_minute)
confusionMatrix(ifelse(sim.results$event > 0, 1, 0), ifelse(times %in% brushing.minutes, 1, 0))

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

# trying dbscan clustering
require(fpc)
test.summary.wide <- melt(test.summary %>% rename(metric = variable), id.vars = c("Timestamp", "metric")) %>% 
                          dcast(Timestamp ~ metric + variable, fill = 0) 


times <- filter(test.results, event >0) %>% select(time_minute) %>% distinct
test.summary.wide <- test.summary.wide %>% filter(Timestamp %in% times$time_minute)

require(caret)
pca.mdl <- preProcess(test.summary.wide %>% select(-Timestamp), method = "pca", pcaComp = 2)
test.summary.wide.pca <- cbind(Timestamp = test.summary.wide$Timestamp, predict(pca.mdl, test.summary.wide %>% select(-Timestamp)))

ggplot(test.summary.wide.pca %>% mutate(flag = ifelse(Timestamp %in% times$time_minute, 1, 0)),
       aes(PC1, PC2)) + geom_point(aes(col = as.factor(flag)))

d.cluster <- dbscan(test.summary.wide.pca %>% select(-Timestamp), eps = .05, MinPts = 3, scale = TRUE)
ggplot(cbind(cluster = d.cluster$cluster, test.summary.wide.pca), aes(PC1, PC2)) + geom_point(aes(col = as.factor(cluster)))

cbind(cluster = d.cluster$cluster, Timestamp = as.character(test.summary.wide$Timestamp)) %>% 
  data.frame %>% filter(Timestamp %in% as.character(times$time_minute))

cbind(cluster = d.cluster$cluster, Timestamp = as.character(test.summary.wide$Timestamp)) %>% 
  data.frame %>% filter(cluster == 1)

# kmeans
k.cluster <- kmeans(test.summary.wide %>% select(-Timestamp), centers = 25)
table(k.cluster$cluster)
cbind(cluster = k.cluster$cluster, Timestamp = as.character(test.summary.wide$Timestamp)) %>% 
  data.frame %>% filter(Timestamp %in% as.character(times$time_minute))

cbind(cluster = k.cluster$cluster, Timestamp = as.character(test.summary.wide$Timestamp)) %>% 
  data.frame %>% filter(cluster == 21)

#################
# looking at vector mag rate of change

tmp.rate <- data.frame(vector.rate = rollapply(test$vector.mag, 2, function(x) (x[2] - x[1])/.01),
                       time = 1:(length(test$vector.mag)-1),
                       event = test.results$event[-1]) 

ggplot(tmp.rate %>% filter(event > 0), 
       aes(time, vector.rate)) + 
  geom_point(aes(col = as.factor(event)))

ggplot(tmp.rate, 
       aes(time, vector.rate)) + 
  geom_point(aes(col = as.factor(event)))



#######################################
# 
# what are most important parts of fingerprint?
brushing.fingerprint
brushing.fingerprint.sd

# peaks per second
ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = peaks.per.sec)) + 
                                        geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                                                                  as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                                                                  as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                                                                  as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
                                                                  ), col = "red") +
  geom_hline(yintercept = 9.3, col = "blue") +
  labs(title = "peaks per sec")

# avg.period 
ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = avg.period)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = 10.8, col = "blue") +
  labs(title = "avg period")

# sd.period
ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = sd.period)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = 4.5, col = "blue")

# mean
ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = Mean)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = 1.068, col = "blue")+
  labs(title = "mean")

# median
ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = Median)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = 1.03, col = "blue") +
  labs(title = "median")


ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = Qu1)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = .82, col = "blue") +
  labs(title = "1stQ")

ggplot(test.summary %>% filter(variable == "vector.mag",
                               Timestamp > dmy_hms("20/01/2016 170000", tz = "GMT"),
                               Timestamp < dmy_hms("20/01/2016 220000", tz = "GMT")), 
       aes(x =Timestamp)) + geom_line(aes(y = Qu2)) + 
  geom_vline(xintercept = c(as.numeric(dmy_hms("20/01/2016 170200", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 180700", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 182900", tz = "GMT")),
                            as.numeric(dmy_hms("20/01/2016 215600", tz = "GMT"))
  ), col = "red") +
  geom_hline(yintercept = 1.28, col = "blue") +
  labs(title = "3rdQ")