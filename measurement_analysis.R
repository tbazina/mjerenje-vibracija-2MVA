
##### Measurement data analysis and plot #####

## Source functions and objects file
source('functions.R')

################################################################################
################################################################################
# Input args from Rscript
# args = commandArgs(TRUE)
# Kratki_Spoj_1_mjerenje_napon_100
# rep1.path = args[1]
# rep1.path <- '/home/tbazina/mjerenje-vibracija-transformatora/mjerenja/Kratki_Spoj_1_mjerenje_napon_100/'
rep1.path <- './mjerenja/90_MVA/Kratki_Spoj_1_mjerenje_napon_100/'


# Load saved image
load(".RData")

################################################################################
################################################################################
# Data input Kratki_Spoj_1_mjerenje_napon_100
# Input first point Kratki_Spoj_1_mjerenje_napon_100
dat <- VibData(
  point.id = 'pt1',
  loc.x = 125,
  loc.y = 1345,
  rib = T,
  replication = 1,
  file.loc = paste0(rep1.path, 'Pt1')
)

# Save current workspace
save.image()

str(dat)
dat %>% select(time, velocity) %>%
  time.series.plot(
    x.name = "time", y.name = "velocity",
    x.min = c(0, 2, 2.1), x.max = c(13, 3, 2.2),
  )

dat %>% select(time, velocity, peak.frequency, sample.rate, d.time) %>%
  # slice(1:20000) %>%
  mutate(velocity_ts = ts(velocity, start = 0, frequency = sample.rate[1]), 
         ma.smooth = ma(velocity_ts, 1),
         sin.term = sin(2*pi*peak.frequency[1]*time),
         sin2.term = sin(4*pi*peak.frequency[1]*time),
         cos2.term = cos(4*pi*peak.frequency[1]*time),
         cos.term = cos(2*pi*peak.frequency[1]*time)) %>% drop_na(ma.smooth) %>%
         mutate(int = cumsum(ma.smooth) * d.time[1]) %>%
  time.series.plot(
    # x.name = "time", y.name = "velocity",
    x.name = "time", y.name = "ma.smooth",
    # x.name = "time", y.name = "int",
    x.min = c(0, 0, 0.0), x.max = c(13, 1, 0.1),
  )
  nest(data = everything()) %>%
  mutate(fit = map(data, ~ lm(.$velocity ~ .$sin.term + .$cos.term + .$sin2.term + .$cos2.term)),
         fit.ar = map(data, ~ ar(.$velocity_ts)),
         fit.arima = map(data, ~ auto.arima(.$velocity_ts)),
         tidied = map(fit, glance)) %>%
  unnest(data) %>%
  mutate(velocity.pred = predict(.$fit[[1]]),
         velocity.pred.ar = predict(.$fit.ar[[1]], n.ahead=20000, se.fit=F),
         velocity.pred.arima = predict(.$fit.arima[[1]], n.ahead=20000, se.fit=F)
         ) %>%
  select(velocity, velocity.pred, velocity.pred.ar, velocity.pred.arima, time) %>%
  time.series.plot(
    x.name = "time", y.name = "velocity.pred",
    x.min = c(0, 0, 0.0), x.max = c(13, 1, 0.1),
  )

dat %>% select(time, velocity, peak.frequency, sample.rate, d.time) %>%
  mutate(
    sin.term = sin(2*pi*peak.frequency[1]*time),
    cos.term = cos(2*pi*peak.frequency[1]*time),
  ) %>%
  nest(data = everything()) %>%
  mutate(
    fit.nls = map(data, ~ nls(velocity ~ I(A*cos(2*pi*peak.frequency[1]*time + fi)),
      data=., start = list(A=700, fi=pi/3), trace = F, model = T)),
    fit.lm = map(data, ~ lm(velocity ~ sin.term + cos.term,
      data=., model = T)),
    tidied = map(fit.lm, glance)
  ) %>%
  unnest(data) %>% 
  mutate(
    velocity.pred.nls = predict(fit.nls[[1]], se.fit = F, newdata=.),
    velocity.pred.lm = predict(fit.lm[[1]], se.fit = F, newdata=.),
    ) %>%
  select(time, velocity, velocity.pred.nls, velocity.pred.lm, fit.nls, fit.lm,
         d.time, peak.frequency) %>%
  mutate(
    vel.int = cumsum(velocity) * d.time[1],
    displacement = velocity / 2 / pi / peak.frequency[1],
    nls.int = cumsum(velocity.pred.nls) * d.time[1],
    nls.int = nls.int - mean(nls.int),
    lm.int = cumsum(velocity.pred.lm) * d.time[1],
    ) %>%
  # summarise(
  #   vel.rms = sqrt(mean(velocity^2)),
  #   nls.rms = sqrt(mean(velocity.pred.nls^2)),
  #   lm.rms = sqrt(mean(velocity.pred.lm^2)),
  #   vel.mean = mean(velocity),
  #   nls.mean = mean(velocity.pred.nls),
  #   lm.mean = mean(velocity.pred.lm),
  #   vel.min = min(velocity),
  #   nls.min = min(velocity.pred.nls),
  #   lm.min = min(velocity.pred.lm),
  #   vel.max = max(velocity),
  #   nls.max = max(velocity.pred.nls),
  #   lm.max = max(velocity.pred.lm)
  #   )
  slice(10000:11000) %>%
  ggplot(aes(x = time, y = displacement)) +
    geom_line() +
    geom_line(aes(y = nls.int), size=1, color="red")
    # geom_line(aes(y = lm.int), size=1, color="green") 
  ggplot(aes(x = time, y = velocity)) +
    geom_point() +
    geom_line(aes(y= velocity.pred.nls), size=1, color="red") +
    geom_line(aes(y= velocity.pred.lm), size=1, color="green") 
         
  mutate(
    peak.frequency = peak.frequency[1],
    ma.smooth = ma(velocity, 2),
    diff = c(0, diff(ma.smooth)),
    ) %>% drop_na(ma.smooth, diff) %>%
  filter(row_number() > min(row_number()[velocity >= 0 & diff > 0]+13)) %>%
  mutate(
    time = (row_number()-1)*d.time[1],
    sin.term = sin(2*pi*peak.frequency[1]*time + pi/2),
    ) %>%
  # time.series.plot(
  #   x.name = "time", y.name = "ma.smooth",
  #   x.min = c(0, 0, 0.0), x.max = c(13, 1, 0.1),
  # )
  nest(data = everything()) %>%
  mutate(
    fit = map(data, ~ lm(.$velocity ~ .$sin.term)),
    tidied = map(fit, glance),
  ) %>%
  pull(fit) %$% summary(.[[1]])
  time.series.plot(
    x.name = "time", y.name = "ma.smooth",
    x.min = c(0, 0, 0.0), x.max = c(13, 1, 0.1),
  )
  
  