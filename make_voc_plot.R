E <- read.csv('/Users/wdaniels/Desktop/E.csv')
NW <- read.csv('/Users/wdaniels/Desktop/NW.csv')
S <- read.csv('/Users/wdaniels/Desktop/S.csv')
SW <- read.csv('/Users/wdaniels/Desktop/SW.csv')

library(lubridate)

E$Hour.Start..MST. <- mdy_hm(E$Hour.Start..MST., tz = "America/Denver")
NW$Hour.Start..MST. <- mdy_hm(NW$Hour.Start..MST., tz = "America/Denver")
S$Hour.Start..MST. <- mdy_hm(S$Hour.Start..MST., tz = "America/Denver")
SW$Hour.Start..MST. <- mdy_hm(SW$Hour.Start..MST., tz = "America/Denver")

# E$Hour.Start..MST.

E$VOC[E$VOC < 0] <- NA
NW$VOC[NW$VOC < 0] <- NA
S$VOC[S$VOC < 0] <- NA
SW$VOC[SW$VOC < 0] <- NA

this.mask <- day(E$Hour.Start..MST.) == 12

png('/Users/wdaniels/Desktop/voc.png',
    res = 100, pointsize = 24, width = 1920, height = 1080)

plot(E$Hour.Start..MST.[this.mask], E$VOC[this.mask], type = "l", ylim = c(0,1), ylab = "VOC [ppm]",
     lwd = 3,
     xlab = "")
lines(NW$Hour.Start..MST.[this.mask], NW$VOC[this.mask], col = "red", lwd = 3)
lines(S$Hour.Start..MST.[this.mask], S$VOC[this.mask], col = "blue", lwd= 3)
lines(SW$Hour.Start..MST.[this.mask], SW$VOC[this.mask], col = "green", lwd = 3)

legend("topright", c("E", "NW", "S", "SW"), col = c("black", "red", "blue", "green"), lwd = 5)


dev.off()
