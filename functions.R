
##### Functions and objects for time-series vibration measurement analysis #####

## Libraries
# Runs test
library(tseries)
# Data Table
library(data.table)
# Cleaning data
library(tidyr)
library(tibble)
library(dplyr)
library(magrittr)
library(purrr)
library(broom)
library(readxl)
# Plotting
library(lattice)
library(ggplot2)
library(latticeExtra)
library(RColorBrewer)
library(tactile)
library(viridis)
# Power Spectral Density Estimation
library(psd)
library(spectral)
# Autocorrelation plot
library(forecast)
# Importing data
library(readr)
library(xlsx)
# Saving images
# Finding Peaks
library(pracma)
# String operations
library(stringr)
# Control limits
library(qcc)


################################################################################
################################################################################
# Radiation ratio (efficiency) of baffled plate
baff.plate.rad.eff <- function(
  plate.width, plate.height, plate.thickness, freq
) {
  # Steel plate (EN10029 8 S235JRG2)
  steel.density <- 7850 # rho [kg/m^3]
  steel.young <- 210000 * 1e6 # E [N/m^2]
  steel.poisson <- 0.3 # ni [/]
  # Air
  air.density <- 1.2 # rho0 [kg/m^3] at 20°
  sound.speed <- 343 # c0 [m/s] at 20°
  plate.area <- plate.width * plate.height # S [m^2]
  plate.perimeter <- 2 * (plate.width + plate.height) # P [m]
  plate.bending.stiffness <- (steel.young * plate.thickness^3 / 12 /
    (1 - steel.poisson^2)) # D [Nm]
  fund.natural.freq<- (pi / 2 * (1/plate.width^2 + 1/plate.height^2) *
    sqrt(plate.bending.stiffness / steel.density / plate.thickness)) # f1,1 [Hz]
  corner.mode.freq <- 3 * sound.speed / plate.perimeter # fe [Hz]
  critical.freq <- (sound.speed^2 / 2 / pi *
                      sqrt(steel.density * plate.thickness /
                             plate.bending.stiffness)) # fc [Hz]
  alpha <- sqrt(freq / critical.freq)
  sigm1 <- 4 * plate.area * freq^2 / sound.speed^2
  sigm2 <- (4 * pi^2 * plate.bending.stiffness / sound.speed^2 / plate.area /
           steel.density / plate.thickness)
  sigm3 <- (plate.perimeter * sound.speed / 4 / pi^2 / plate.area /
             critical.freq *
             ((1 - alpha^2) * log((1+alpha) / (1-alpha)) + 2 * alpha) /
             (1 - alpha^2)^(3/2))
  sigm4 <- 0.45 * sqrt(plate.perimeter * critical.freq / sound.speed)
  sigm5 <- 1 / sqrt(1 - critical.freq / freq)
  if (freq < fund.natural.freq) {
    if (sigm1 > sigm2) return(sigm2)
    else return(sigm1)
  }
  else if (freq < corner.mode.freq) {
    return(sigm2)
  }
  else if (freq < 0.85 * critical.freq) {
    if (sigm3 > sigm2 & sigm3 < sigm4) {return(sigm3)}
    else if (sigm3 < sigm2) return(sigm2)
    else return(sigm4)
  }
  else if (freq >= 0.85 * critical.freq & freq < 1.2 * critical.freq) {
    return(sigm4)
  }
  else {
    if (sigm5 > sigm4) return(sigm4)
    else return(sigm5)
  }
}


################################################################################
################################################################################
## Contour plots with ribs

contour.plot.rib.freq <- function(
  matrix.data,
  title='Signal Time Domain', subtitle='',
  xlab='Length [mm]', ylab='Height [mm]',
  z.tick.num = 8,
  aspect = 25/36,
  save.fig = F,
  fig.filename = '',
  fig.width = 12,
  fig.height = 8,
  fig.dpi = 300
  ) {
  # Power transfomer plate height and width
  width.x <- 3610
  height.y <- 2600
  # Power transformer rib dimensions
  rib.flat.x <- 130
  rib.rad.x <- 30
  # Power transformer dimensions between ribs
  flat.x <- 380
  element.no <- 13
  
  # Geometric data from minitab for plotting
  plate.length <- c(
    125, 368.2142857,	611.4285714,	854.6428571,	1097.857143, 1341.071429,
    1584.285714,	1827.5,	2070.714286,	2313.928571, 2557.142857,	2800.357143,
    3043.571429,	3286.785714,	3530
    )
  plate.height <- c(
    122, 285.5, 449, 612.5, 776, 939.5, 1103, 1266.5, 1430, 1593.5, 1757, 1920.5,															
    2084, 2247.5, 2411															
  )
  
  # Plot ribs with radius
  rib.boundaries <- cumsum(
    c(0, rib.rad.x, rib.flat.x, rib.rad.x,
      rep(c(flat.x, rib.rad.x, rib.flat.x, rib.rad.x), times=6)))
  
  # Color palette
  custom.colors <- brewer.pal(7, 'Set1')
  # Add alpha
  custom.colors.settings <- list(
    plot.symbol = list(
      cex = 2.4,
      col = custom.colors[1],
      fill = '#ffffff',
      pch = 21),
    superpose.symbol = list(
      cex = rep(1.4, 7),
      col = custom.colors,
      pch = rep(20, 7)),
    add.text = list(
      cex = 0.9,
      col = '#000000'
    )
  )
  
  # Grids and ticks
  x.at <- seq(0, width.x, by = 200)
  y.at <- seq(0, height.y, by = 200)

  # Creating device for saving image
  if (save.fig) {
    png(
      filename = fig.filename,
      units = 'in',
      width = fig.width,
      height = fig.height,
      pointsize = 12 * fig.dpi / 72,
      res = fig.dpi,
      bg = 'white'
    )
  }
  
  # Plot
  plt <- levelplot(matrix.data, row.values = plate.length,
            column.values = plate.height, cuts = 10,
            main=title, sub = subtitle,
            ylab=ylab, xlab=xlab, aspect = aspect,
            region = T, pretty = T, labels = T, contour = T,
            col.regions = inferno(100, 0.9, 0.1),
            par.settings = custom.colors.settings,
            colorkey = list(
                raster = T, interpolate = T, tick.number = z.tick.num,
                width = 1, height = 0.8, space = "right"
              ),
            scales = list(
              x = list(at = x.at, rot = 45, limits=c(0, width.x)),
              y = list(at = y.at, limits = c(0, height.y-100))
            ),
            panel = function(...){
              panel.abline(v = x.at, h = y.at, col.line = 'lightgrey', ...)
              panel.2dsmoother(method = "loess", ..., args = list(span=0.1),
                               n = 300)
              panel.abline(v = rib.boundaries)
              }
            )
  
  # Turning off device
  if (save.fig) {
    trellis.par.set('fontsize', list(text = 12, points = 8))
    print(plt)
    dev.off()
  }
  
  # Showing plot
  print(plt)
}
################################################################################
################################################################################
# Frequency spectrum calculate and plot
freq.spec <- function(
  input.signal, sample.rate,
  plt = F, plt.title='Signal Frequency Spectrum', plt.subtitle='',
  plt.xlab='Frequency [Hz]', plt.ylab='Signal Amplitude',
  min.plot.freq = 0, max.plot.freq = 1000) {
  # Pad input.signal with zeros to nearest highly composite number for fast FFT
  extended.size <- nextn(length(input.signal))
  if (length(input.signal) < extended.size) {
    extended.signal <- replace(
      input.signal, (length(input.signal)+1):extended.size, 0
      )
  }
  else {
    extended.signal <- input.signal
  }
  
  # Compute fft of padded input signal and keep only first half
  # Second half are complex conjugates of first half
  signal.bin.num <- extended.size / 2
  signal.fft <- fft(extended.signal)[1:signal.bin.num]
  
  # Discrete signal Frequency bins
  signal.freq <- seq(
    from=0, by=sample.rate/extended.size, length.out=signal.bin.num
    )
  
  # Signal Power spectrum density for single-sided spectrum
  signal.psd <- (1 / (sample.rate * extended.size)) * Mod(signal.fft) ^ 2
  signal.psd[-1] <- 2 * signal.psd[-1]
  
  # Signal amplitude
  signal.amp <- signal.fft / extended.size
  signal.amp[-1] <- 2 * signal.amp[-1]
  signal.amp <- Mod(signal.amp)
  
  # Spectrum plot for signal
  if (plt) {
    plot.slice <- signal.freq >= min.plot.freq & signal.freq <= max.plot.freq
    print(xyplot(
      signal.amp[plot.slice] ~ signal.freq[plot.slice],
      type=c('l', 'g'), main=plt.title, sub=plt.subtitle,
      ylab=plt.ylab, xlab=plt.xlab
      )
    )
  }
  
  return(tibble(
    amp = signal.amp,
    psd = signal.psd,
    freq = signal.freq,
    fft = signal.fft
    ))
}
################################################################################

################################################################################
# Object for storing mesurement and analysis data
VibData <- function(
  file.loc,
  point.id,
  loc.x,
  loc.y,
  rib,
  replication
  ){
  # Data inputa as tibble
  df.orig <- read_tsv(
    file = file.loc,
    locale = locale(decimal_mark=c(".")),
    col_types = "nnnnnnnnn"
    )
  
  # Renaming Displacement -> Max Displacement
  names(df.orig)[6:8] <- paste('Peak', c('Vel (um/s)', 'Disp (nm)', 'Acc (g)'))
  
  # Calibration constants and velocity range in um/(s*V)
  cal_const <- c(4990, 98910)
  names(cal_const) <- c('LO', 'HI')
  
  
  if (df.orig[['Range(HI/LO)']][1] == 1) {
    vel.range <- 'LO'
  }
  else {
    #TODO: Check if 2 corresponds to high velocity range
    vel.range <- 'HI'
  }
  
  # Single value information
  info <- with(
    df.orig,
    tibble(
      point.id = point.id,
      loc.x = loc.x,
      loc.y = loc.y,
      rib = rib,
      replication = replication,
      subtitle = paste(
        point.id, ': (', loc.x, ', ', loc.y, ') rib: ', rib, ' rep: ',
        replication, sep=""),
      range = vel.range,
      d.time = `dt  (sec)`[1],
      d.freq = `df  (Hz)`[1],
      sample.rate = 1 / d.time,
      sample.num = nrow(df.orig),
      peak.freq.orig = `Peak Freq(Hz)`[1],
      peak.vel.orig = `Peak Vel (um/s)`[1],
      peak.disp.orig = `Peak Disp (nm)`[1],
      peak.acc.orig = `Peak Acc (g)`[1],
    )
  )
  
  # Time-domain tibble
  time.domain <- with(
    df.orig,
    tibble(
      time = (1:info$sample.num - 1) * info$d.time,
      voltage = `Voltage (V)`,
      # Removing DC component from signal before calculating velocity
      velocity = (voltage - mean(voltage)) * cal_const[info$range]
    )
  )
  
  # Sample acf and pacf
  time.domain <- time.domain %>%
    mutate(velocity_ts = ts(.$velocity, frequency = info$sample.rate[1])) %>%
    add_column(!!!tibble(
      acf = acf(.$velocity_ts, lag.max = length(.$velocity_ts),
                type = 'correlation', plot = F),
      pacf = pacf(.$velocity_ts, lag.max = length(.$velocity_ts),
                  type = 'correlation', plot = F)
      ) %$%
        tibble(
          acf = as.vector(.$acf$acf),
          pacf = c(1, as.vector(.$pacf$acf))
          )
      ) %>%
    select(-velocity_ts)
    
  # Peak-to-Peak Velocity
  info <- info %>%
    add_column(
      peak.to.peak.vel = with(
        time.domain, abs(max(velocity)) + abs(min(velocity))
        ),
      rms.vel = with(
        time.domain, sqrt(mean(velocity^2))
        ),
    )
  
  # Frequency-domain tibble
  spec.vel <- freq.spec(time.domain$velocity, info$sample.rate)
  freq.domain <- with(
    df.orig,
    tibble(
      frequency = spec.vel$freq,
      displacement.orig = NA,
      velocity.orig = NA,
      velocity.amp = spec.vel$amp,
      velocity.psd = spec.vel$psd,
      velocity.fft = spec.vel$fft
      )
  )
  
  # Rename observed value based on iput data (4th column)
  if (names(df.orig)[4] == 'Displacement (nm)') {
    rename.name <- 'displacement.orig'
  }
  else {
    rename.name <- 'velocity.orig'
  }
  freq.domain <- freq.domain %>%
    mutate(!!rename.name := df.orig[[4]][1:length(frequency)])
  
  # Top 5 Velocity Peaks
  peaks <- tibble(
    peak.velocity.ind = findpeaks(
      freq.domain$velocity.amp, nups=1, ndowns=1, npeaks=5, sortstr=T
      )[,2],
    peak.frequency = freq.domain$frequency[peak.velocity.ind],
    peak.velocity.amp = freq.domain$velocity.amp[peak.velocity.ind],
    peak.velocity.psd = freq.domain$velocity.psd[peak.velocity.ind]
  )
  
  
  # Cleaning data
  # df.orig <- df.orig %>%
  #   ## Replacing zeros with NA
  #   mutate(`dt  (sec)` = replace(`dt  (sec)`, `dt  (sec)` == 0, NA)) %>%
  #   mutate(`df  (Hz)` = replace(`df  (Hz)`, `df  (Hz)` == 0, NA)) %>%
  #   mutate(`Peak Freq(Hz)` = replace(
  #     `Peak Freq(Hz)`, `Peak Freq(Hz)` == 0, NA)) %>%
  #   mutate(`Peak Vel (um/s)` = replace(
  #     `Peak Vel (um/s)`, `Peak Vel (um/s)` == 0, NA)) %>%
  #   mutate(`Peak Disp (nm)` = replace(
  #     `Peak Disp (nm)`, `Peak Disp (nm)` == 0, NA)) %>%
  #   mutate(`Peak Acc (g)` = replace(
  #     `Peak Acc (g)`, `Peak Acc (g)` == 0, NA)) %>%
  #   mutate(`Range(HI/LO)` = replace(`Range(HI/LO)`, `Range(HI/LO)` == 0, NA)) %>%
  #   ## Adding time vector
  #   mutate('Time (sec)' = (1:length(`dt  (sec)`) - 1) * `dt  (sec)`[1]) %>%
  #   ## Adding frequency vector
  #   mutate('Frequency (Hz)' = (1:length(`df  (Hz)`) - 1) * `df  (Hz)`[1]) %>%
  #   ## Adding velocity vector
  #   mutate('Velocity (um/s)' = `Voltage (V)` * cal_const[vel.range])
  
  # Checking data
  # cat('\n________Glimpse________\n')
  # glimpse(info)
  # glimpse(time.domain)
  # glimpse(freq.domain)
  # glimpse(peaks)
  cat('\n________Head________\n')
  print(head(info))
  # print(head(time.domain))
  # print(head(freq.domain))
  print(head(peaks))
  # Descriptive statistics
  cat('\n________Descriptive Statistics________\n')
  # print(summary(info))
  print(summary(time.domain))
  print(summary(freq.domain))
  # print(summary(peaks))
  
  peaks <- peaks %>%
    add_row(peak.velocity.ind = rep(
      NA, length.out=(info$sample.num - nrow(peaks))))
  
  freq.domain <- freq.domain %>%
    add_row(frequency = rep(
      NA, length.out=(info$sample.num - nrow(freq.domain))))
  
  # Returning list with clean data
  return(tibble(!!!info,
              !!!time.domain,
              !!!freq.domain,
              !!!peaks))
}
################################################################################

################################################################################
# Time series plot and save data
time.series.plot <- function(
  data, x.name, y.name,
  title='Signal Time Domain', subtitle='',
  xlab='Time [s]', ylab='Signal Amplitude',
  x.min = c(0, 2, 2.1), x.max = c(Inf, 3, 2.2),
  panel.names = c("Entire Dataset", "1 s", " 100 ms"),
  x.tick.num = 100, y.tick.num = 6,
  aspect = 1/8,
  col.palette.name = "Blues",
  col.line.num = 5,
  plot.type = 'l',
  point.symbol = 21,
  point.cex = 0.8,
  save.fig = F,
  fig.filename = '',
  fig.width = 20,
  fig.height = 10,
  fig.dpi = 300
){
  # Color palette and padding
  custom.colors <- brewer.pal(6, col.palette.name)
  custom.plot.settings <- list(
    plot.line = list(
      col = custom.colors[col.line.num], lwd = 1, border = 'transparent'),
    plot.symbol = list(
      col = custom.colors[col.line.num], cex = point.cex, pch = point.symbol,
      fill = custom.colors[col.line.num]),
    layout.widths = list(left.padding = 0, right.padding=0)
  )
  # Plot slice limits
  x.slice <- data[[x.name]] >= min(x.min) & data[[x.name]] <= max(x.max)
  # Tick locations
  x.at <- pretty(floor(min(x.min)):ceiling(max(x.max)),
                 n=x.tick.num,
                 u5.bias = 100
  )
  y.at <- pretty(
    min(data[[y.name]][x.slice]):max(data[[y.name]][x.slice]),
    n=y.tick.num,
  )
  
  # Creating device for saving image
  if (save.fig) {
    png(
      filename = fig.filename,
      units = 'in',
      width = fig.width,
      height = fig.height,
      pointsize = 12 * fig.dpi / 72,
      res = fig.dpi,
      bg = 'white'
    )
  }
  
  # Plot
  # Data for plot and groups
  y.data <- rep(data[[y.name]][x.slice], length(x.max))
  x.data <- rep(data[[x.name]][x.slice], length(x.max))
  panel.groups <- rep(
    rev(panel.names), each=length(data[[y.name]][x.slice])
    )
  # List of x axis limits
  xlims <- cbind(rev(x.min), rev(x.max))
  xlims.list <- split(xlims, row(xlims))
  plt <- xyplot(
    y.data ~ x.data | panel.groups,
    type=c(plot.type), strip = F, strip.left = T,
    main=title, sub = subtitle,
    ylab=ylab, xlab=xlab, aspect = aspect,
    xlim = xlims.list,
    layout = c(1, length(x.min)),
    scales = list(
      relation = 'free',
      y = list(at = y.at, rot=0),
      x = list(at = x.at, rot=45),
      axs = "i"
    ),
    par.settings = custom.plot.settings,
    panel = function(...) {
      panel.abline(h = y.at, v = x.at, col.line = 'lightgrey', ...)
      panel.xyplot(...)
    }
  )
  
  # Turning off device
  if (save.fig) {
    trellis.par.set('fontsize', list(text = 12, points = 8))
    print(plt)
    dev.off()
  }
  
  # Showing plot
  print(plt)
}
################################################################################

################################################################################
# Frequency spectrum plot and save figure
freq.spec.plot <- function(
  data, x.name, y.name, condition='', groups=NULL,
  title='Signal Frequency Spectrum', subtitle='',
  xlab='Frequency [Hz]', ylab='Signal Amplitude',
  plot.type = c('h'),
  x.min = 1, x.max = 600,
  x.tick.num = 25, y.tick.num = 8,
  log.x = F, log.y = F,
  yscale.components.function = yscale.components.default,
  jitter.x = F, jitter.amount=8,
  aspect = 1/4,
  col.palette.name = "Blues",
  col.different.num = 6,
  col.alpha = 1,
  col.line.num = 5,
  col.polygon.num = 5,
  line.width = 1,
  table.panel.order = T,
  show.legend = F,
  save.fig = F,
  fig.filename = '',
  fig.width = 14,
  fig.height = 5,
  fig.dpi = 300
){
  # Color palette
  custom.colors <- brewer.pal(col.different.num, col.palette.name)
  # Add alpha
  custom.colors <- adjustcolor(custom.colors, alpha.f = col.alpha)
  custom.colors.settings <- list(
    plot.line = list(
      col = custom.colors[col.line.num], lwd = line.width, 
      border = 'transparent'),
    superpose.line = list(
      col = custom.colors[col.line.num:length(custom.colors)],
      lwd = line.width,  border = 'transparent'),
    superpose.symbol = list(
      col = custom.colors[col.line.num:length(custom.colors)],
      fill = custom.colors[col.line.num:length(custom.colors)],
      pch = rep(19, col.different.num)),
    plot.polygon = list(
      col = custom.colors[col.polygon.num], lwd = 1, border = 'transparent')
  )
  # Plot slice limits
  x.slice <- data[[x.name]] >= x.min & data[[x.name]] <= x.max
  # Tick locations
  x.at <- pretty(abs(x.min):x.max, n=x.tick.num, u5.bias=100)
  y.at <- pretty(
    seq(min(data[[y.name]][x.slice]), max(data[[y.name]][x.slice]),
        length.out = y.tick.num*100),
    n=y.tick.num
  )
  if (log.y) {
    y.at <- seq(
      floor(log10(min(data[[y.name]][x.slice]))),
      ceiling(log10(max(data[[y.name]][x.slice]))),
      by = 1
      )
  }
  
  # Creating device for saving image
  if (save.fig) {
    png(
      filename = fig.filename,
      units = 'in',
      width = fig.width,
      height = fig.height,
      pointsize = 12 * fig.dpi / 72,
      res = fig.dpi,
      bg = 'white'
    )
  }
  
  # Plot horizontal lines for each bin
  if (condition!='') {
    condition <- data[[condition]][x.slice]
  }
  # If groups exist
  if (!is.null(groups)) {
    x.v.abline.ind <- c(
      T, abs(diff(round(sort(data[[x.name]][x.slice], decreasing = T)))) >= 10
      )
    x.v.abline <- (data[[x.name]][x.slice])[x.v.abline.ind]
    # Jitter data if necessary
    if (jitter.x) {
      x.data <- data %>%
        mutate(groups_jitter = .[[groups]] - 
                 length(unique(.[[groups]])) %/% 2 - 1) %>%
        mutate(x.data = .[[x.name]] + groups_jitter * jitter.amount) %>%
        pull(x.data)
      x.data <- x.data[x.slice]
    } else {
      x.data <- data[[x.name]][x.slice]
    }
    groups <- data[[groups]][x.slice]
  } else {
    x.data <- data[[x.name]][x.slice]
  }
    plt <- xyplot(
      data[[y.name]][x.slice] ~ x.data | condition,
      groups = groups,
      type=plot.type, main=title, sub = subtitle,
      ylab=ylab, xlab=xlab, aspect = aspect,
      xlim = c(x.min, x.max),
      scales = list(
        y = list(at = y.at, log=log.y),
        x = list(at = x.at, rot=45, log=log.x)
      ),
      yscale.components = yscale.components.function,
      as.table = table.panel.order,
      auto.key = show.legend,
      par.settings = custom.colors.settings,
      panel = function(...) {
        panel.abline(h = y.at, v = x.at, col.line = 'lightgrey', ...)
        # panel.xyarea(origin=0, ...)
        panel.xyplot(...)
        if (!is.null(groups)) {
          panel.abline(v = x.v.abline, lty = "2565", lwd = 0.8, ...)
        }
      }
    )
  
  # Turning off device
  if (save.fig) {
    trellis.par.set('fontsize', list(text = 12, points = 8))
    print(plt)
    dev.off()
    
  }
  
  # Showing plot
  print(plt)
}
################################################################################

################################################################################
# AFC and PACF plot and save figure
acf.pacf.plot <- function(
  data, x.name, y.name, condition='', groups=NULL,
  title='Autocorrelation and partial autocorrelation', subtitle='',
  xlab='Time [s]', ylab='Correlation [ / ]',
  x.min = -0.01, x.max = 3.2,
  x.tick.num = 60, y.tick.num = 15,
  log.x = F, log.y = F,
  jitter.x = F, jitter.amount = 8,
  aspect = 1/4,
  col.palette.name = "Blues",
  col.different.num = 6,
  col.alpha = 1,
  col.line.num = 5,
  col.polygon.num = 5,
  line.width = 1.8,
  table.panel.order = T,
  strip.left = T,
  strip.top = F,
  show.legend = F,
  save.fig = F,
  fig.filename = '',
  fig.width = 20,
  fig.height = 10,
  fig.dpi = 300
){
  # Color palette
  custom.colors <- brewer.pal(col.different.num, col.palette.name)
  # Add alpha
  custom.colors <- adjustcolor(custom.colors, alpha.f = col.alpha)
  custom.colors.settings <- list(
    plot.line = list(
      col = custom.colors[col.line.num], lwd = line.width, 
      border = 'transparent'),
    superpose.line = list(
      col = custom.colors[col.line.num:length(custom.colors)],
      lwd = line.width,  border = 'transparent'),
    plot.polygon = list(
      col = custom.colors[col.polygon.num], lwd = 1, border = 'transparent')
  )
  
  # Plot slice limits
  x.slice <- data[[x.name]] >= x.min & data[[x.name]] <= x.max
  # Tick locations
  x.at <- pretty(
    floor(min(x.min)):ceiling(max(x.max)), n=x.tick.num, u5.bias = 100)
  y.at <- pretty(
    floor(min(data[[y.name]][x.slice])):ceiling(max(data[[y.name]][x.slice])),
    n=y.tick.num,
  )
  
  # 95% Confidence intervals
  sample.length <- length(data[[x.name]]) / 2
  upper <- 3 / sqrt(sample.length)
  lower <- -3 / sqrt(sample.length)
  
  # Creating device for saving image
  if (save.fig) {
    png(
      filename = fig.filename,
      units = 'in',
      width = fig.width,
      height = fig.height,
      pointsize = 12 * fig.dpi / 72,
      res = fig.dpi,
      bg = 'white'
    )
  }
  
  # Plot horizontal lines for each bin
  if (condition!='') {
    condition <- data[[condition]][x.slice]
  }
  # If groups exist
  if (!is.null(groups)) {
    # Jitter data if necessary
    if (jitter.x) {
      x.data <- data %>%
        mutate(groups_jitter = .[[groups]] - 
                 length(unique(.[[groups]])) %/% 2 - 1) %>%
        mutate(x.data = .[[x.name]] + groups_jitter * jitter.amount) %>%
        pull(x.data)
      x.data <- x.data[x.slice]
    } else {
      x.data <- data[[x.name]][x.slice]
    }
    groups <- data[[groups]][x.slice]
  } else {
    x.data <- data[[x.name]][x.slice]
  }
    plt <- xyplot(
      data[[y.name]][x.slice] ~ x.data | condition,
      groups = groups,
      type=c('h'), main=title, sub = subtitle,
      ylab=ylab, xlab=xlab, aspect = aspect,
      as.table = table.panel.order,
      upper = upper, lower = lower,
      strip = strip.top, strip.left = strip.left,
      xlim = c(x.min, x.max),
      scales = list(
        y = list(at = y.at, log=log.y),
        x = list(at = x.at, rot=45, log=log.x)
      ),
      auto.key = show.legend,
      par.settings = custom.colors.settings,
      panel = function(...) {
        panel.abline(h = y.at, v = x.at, col.line = 'lightgrey', ...)
        panel.xyplot(...)
        panel.abline(h = c(lower, upper), lty='2565', lwd=1.2,
                     col.line = 'red', ...)
      }
    )
  
  # Turning off device
  if (save.fig) {
    trellis.par.set('fontsize', list(text = 12, points = 8))
    print(plt)
    dev.off()
    
  }
  
  # Showing plot
  print(plt)
}
################################################################################
################################################################################
## Plot rib effect
rib.effect.plot <- function(
  data, x.name, y.name, condition='', groups=NULL,
  title= 'Velocity amplitude - Rib effect', subtitle='Distinct frequencies and RMS',
  xlab= 'Frequency',
  ylab= expression('Velocity [' * mu * m/s * ']'),
  plot.type = c('b'),
  y.tick.num = 8,
  log.x = F, log.y = T,
  yscale.components.function = yscale.components.log10ticks,
  aspect = 2/3,
  col.palette.name = "Blues",
  col.different.num = 6,
  col.alpha = 1,
  col.line.num = 5,
  col.polygon.num = 5,
  line.width = 1,
  table.panel.order = T,
  show.legend = F,
  save.fig = F,
  fig.filename = '',
  fig.width = 8,
  fig.height = 10,
  fig.dpi = 300
){
  # Color palette
  custom.colors <- brewer.pal(col.different.num, col.palette.name)
  # Add alpha
  custom.colors <- adjustcolor(custom.colors, alpha.f = col.alpha)
  custom.colors.settings <- list(
    plot.line = list(
      col = custom.colors[col.line.num], lwd = line.width, 
      border = 'transparent'),
    superpose.line = list(
      col = custom.colors[col.line.num:length(custom.colors)],
      lwd = line.width,  border = 'transparent'),
    superpose.symbol = list(
      col = custom.colors[col.line.num:length(custom.colors)],
      fill = custom.colors[col.line.num:length(custom.colors)],
      pch = rep(19, col.different.num)),
    plot.polygon = list(
      col = custom.colors[col.polygon.num], lwd = 1, border = 'transparent')
  )
  
  # Gridlines
  x.at = data[[x.name]]
  if (log.y) {
    y.at <- seq(
      floor(log10(min(data[[y.name]]))),
      ceiling(log10(max(data[[y.name]]))),
      by = 1
      )
  }
  # Conditions and groups
  groups <- data[[groups]]
  if (condition!='') {
    condition <- data[[condition]]
  }
  
  # Creating device for saving image
  if (save.fig) {
    png(
      filename = fig.filename,
      units = 'in',
      width = fig.width,
      height = fig.height,
      pointsize = 12 * fig.dpi / 72,
      res = fig.dpi,
      bg = 'white'
    )
  }
  
  # Plot data
  plt <- xyplot(
    data[[y.name]] ~ data[[x.name]] | condition,
    groups = groups,
    type=plot.type, main=title, sub = subtitle,
    ylab=ylab, xlab=xlab, aspect = aspect,
    scales = list(
      y = list(at = y.at, log=log.y),
      x = list(rot=45)
    ),
    yscale.components = yscale.components.function,
    as.table = table.panel.order,
    auto.key = show.legend,
    par.settings = custom.colors.settings,
    panel = function(...) {
      panel.abline(h = y.at, v = x.at, col.line = 'lightgrey', ...)
      panel.xyplot(...)
    }
  )
  
  # Turning off device
  if (save.fig) {
    trellis.par.set('fontsize', list(text = 12, points = 8))
    print(plt)
    dev.off()
    
  }
  
  # Showing plot
  print(plt)
}

################################################################################