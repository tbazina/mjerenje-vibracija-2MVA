
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

# Integral by lm fit, nls fit and dividing by 2*pi*f (100 Hz)
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
    geom_line(aes(y = nls.int), size=1, color="red") +
  scale_color_manual(values = c("2*pi*f (100 Hz)", "Fit")) +
  ggtitle("2*pi*f (100 Hz) - Fit") +
  labs(y = "Pomak [um]", x = "Vrijeme [s]")
    # geom_line(aes(y = lm.int), size=1, color="green") 
  ggplot(aes(x = time, y = velocity)) +
    geom_point() +
    geom_line(aes(y= velocity.pred.nls), size=1, color="red") +
    geom_line(aes(y= velocity.pred.lm), size=1, color="green") 
         
# Integral in frequency domain - bad result
test <- tibble(
  sample.freq = dat$sample.rate[1],
  k = seq(0, length(dat$velocity)),
  fft.pad = fft(c(dat$velocity, rep(0, length(dat$velocity))))[1:(length(dat$velocity)+1)],
  int.op = 1/(2i*pi*k*sample.freq/2/length(dat$velocity))
  ) %>%
  mutate(
    int.op = c(0, int.op[2:length(int.op)]),
    int.freq = fft.pad * int.op,
    int.freq.neg = rev(Conj(int.freq))
  ) %$%
  tibble(
    int.freq.dbl = c(.$int.freq, .$int.freq.neg[2:(length(.$int.freq.neg)-1)]),
    int.time = fft(int.freq.dbl, inverse = T)/length(int.freq.dbl)
  ) %$%
  tibble(
    int.time = detrend(Re(.$int.time[1:(length(.$int.time)/2)])),
    time = dat$time)
test %>% ggplot(aes(x = time, y = int.time)) + geom_line()
