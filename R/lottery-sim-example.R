
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # use if in RStudio

# setwd(getSrcDirectory()[1]) # otherwise

############################################################
###### example of a simulation, formatting, andvizualization

source("lottery-sim-utilities.R")
source("lottery-sim-main.R")

# defining parameters
NSP <- 3
TIME <- 100
SIDE <- 6

PARAM <- data.frame(b = 20, u = 0.95, s = 0.01, d = 0.8, sb = 0.8, sT = 0.8, rho = 0.5)
DRIVER <- data.frame(rt = 0.5, ru = 0.8, rs = 1, rsT = 1)

# runinng the model (with the above parameters, it should take 0.08 seconds)
system.time(out <- LotterySim(nsp = NSP, time = TIME, size = SIDE, param = PARAM, driver = DRIVER, seed = 210))

# data formating
data <- data.frame("time" = 1:TIME)
# cell coordinates
xx <- 1
yy <- 1
name <- paste("sp", 1:NSP, sep = "")
for (i in 1:NSP) {
  data[[name[i]]] <- out$sim_out[i, xx, yy, ]
}
TS <- melt(data, id = "time", value.name = "freq", variable.name = "species")

# plotting
figTS <- ggplot(data = TS, aes(x = time, y = freq)) +
  geom_line(aes(col = species)) +
  ylim(0, 1) +
  theme_bw()
figTS
