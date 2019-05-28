#install.packages(c("fractal","DescTools","biwavelet"))
library(fractal) # contains surrogate function
library(ggplot2)
#library(plyr)
library(dplyr)
library(reshape2)
library(DescTools)
library(biwavelet)

sessionInfo()

# Example
MEI_raw <- read.table("C:/Users/siplem/Dropbox/Chapter3-SardineAnchovy/R files/resomewaveletresults/MEI.txt", header = TRUE, sep = "\t", fill = TRUE)
mei_raw <- tbl_df(MEI_raw)
mei_raw_m <- melt(mei_raw, id.vars="YEAR")
names(mei_raw_m)[2:3] <- c("MON", "VALUE")
ann_mei <- mei_raw_m %>%
  filter(!is.na(VALUE )) %>%
  group_by(YEAR) %>%
  summarize(mean_val = mean(VALUE))
names(ann_mei)[1] <- "year"

# Need to get Sardine.est and Anchovy.est **********

# Calculate annual mean
# ann_MEI <- scale(rowMeans(MEI_raw[,-1], na.rm = TRUE)) #scale(rowMeans(MEI_raw[,c(2,4,6,8,10,12)], na.rm = TRUE)) # mean of 0, sd of 1, 

# Make surrogate time series with same properties
sardine_phase <- surrogate(data.points$Sardine.est,method = "phase")
std_sardine <- data.points$Sardine.est - mean(data.points$Sardine.est)
#plot(1:length(sardine_phase),sardine_phase,type='l',ylim=range(c(sardine_phase,std_sardine)))
#lines(1:nrow(data.points),std_sardine,col='red')

anchovy_phase <- surrogate(data.points$Anchovy.est,method = 'phase')
std_anchovy <- data.points$Anchovy.est - mean(data.points$Anchovy.est)
#plot(1:length(anchovy_phase),anchovy_phase,type='l',ylim=range(c(anchovy_phase,std_anchovy)))
#lines(1:nrow(data.points),std_anchovy,col='red')

sa.col <- c("#ef8a62","#67a9cf")
plot(1:56,anchovy_phase,type='l',col=sa.col[1],lwd=2,
     ylim=range(c(anchovy_phase,sardine_phase)),
     ylab="Simulated Abundance",xlab="Year")
lines(1:56,sardine_phase,col=sa.col[2],lwd=2)

## create surrogate data sets using circulant 
## embedding method 
# surr <- surrogate(ann_mei$mean_val, method="ce")
# 
# surr_mei <- surr[1:length(surr)]
# 
# # plot the surrogate line with the MEI
# with(ann_mei, plot(year,mean_val, type="l", lwd=2), ylim=c(-2,2))
# lines(ann_mei$year, surr_mei, col = rgb(0,0,1,0.7), lwd=2)



## plot and compare various statistics of the 
## surrogate and original time series 
# plot(surr, show="both", type="time")
# plot(surr, show="both", type="sdf")
# plot(surr, show="both", type="lag")
# plot(surr, show="both", type="pdf")

# make surrogate MEI using each of the the "methods"
seed <- 42
set.seed(seed)
surr_ce <- surrogate(ann_mei$mean_val, method="ce")[1:65]
set.seed(seed)
surr_aaft <- surrogate(ann_mei$mean_val, method="aaft")[1:65]
set.seed(seed)
surr_phase <- surrogate(ann_mei$mean_val, method="phase")[1:65]
set.seed(seed)
surr_dh <- surrogate(ann_mei$mean_val, method="dh")[1:65]

Desc(ann_mei$mean_val)
Desc(surr_ce)
Desc(surr_aaft)
Desc(surr_phase)
Desc(surr_dh)

ann_mei_df <- as.data.frame(ann_mei)
plot(ann_mei_df$year,ann_mei_df$mean_val,type='l')
#lines(ann_mei$year, surr_phase, lwd=2, type="l", col=rgb(0,0,1,0.65), ylim=c(-2,2))




# Example ts
old<- par(mfrow=c(5,1), mar=c(1.5,2,1,1))
with(ann_mei, plot(year,mean_val, type="l", lwd=2, ylim=c(-2,2)))
abline(h=0, lty = 2)
plot(ann_mei$year, surr_ce, lwd=2, type="l", col=rgb(1,0,0,0.65), ylim=c(-2,2))
abline(h=0, lty = 2)
plot(ann_mei$year, surr_aaft, lwd=2, type="l", col=rgb(0,1,0,0.65), ylim=c(-2,2))
abline(h=0, lty = 2)
plot(ann_mei$year, surr_phase, lwd=2, type="l", col=rgb(0,0,1,0.65), ylim=c(-2,2))
abline(h=0, lty = 2)
plot(ann_mei$year, surr_dh, lwd=2, type="l", col=rgb(1,0,1,0.65), ylim=c(-2,2))
abline(h=0, lty = 2)
par(old)

# plot acf of each example ts
old<- par(mfrow=c(5,1), mar=c(1,2,1,1))
acf(ann_mei$mean_val, lag.max = 10)
acf(surr_ce, col=rgb(1,0,0,0.65), lag.max = 10)
acf(surr_aaft, col=rgb(0,1,0,0.65), lag.max = 10)
acf(surr_phase, col=rgb(0,0,1,0.65), lag.max = 10)
acf(surr_dh, col=rgb(1,0,1,0.65), lag.max = 10)
par(old)
# plot periodogram of each example ts

tmp <- ceiling(sqrt(nrow(ann_mei)))

if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp}

old<- par(mfrow=c(3,2), mar=c(1.2,3,2.2,1))
spec.pgram(ann_mei$mean_val, spans = c(m,m), main="Annual MEI")
spec.pgram(surr_ce, spans = c(m,m), col=rgb(1,0,0,0.65), main = "CE")
spec.pgram(surr_aaft, spans = c(m,m), col=rgb(0,1,0,0.65), main = "AAFT")
spec.pgram(surr_phase, spans = c(m,m), col=rgb(0,0,1,0.65), main = "PHASE")
spec.pgram(surr_dh, spans = c(m,m), col=rgb(1,0,1,0.65), main = "DH")
par(old)

# plot wavelet spectram

J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1
# # Torrence and Compo 1998 equation 11 constants for reconstruction of time series for
# # Morlet mother wavelets using a delta function 
# # Morlet wavelet constants
# C_delta <- 0.776 #  The factor C_delta comes from the reconstruction of a
# # delta function from its wavelet transform using the function Psi_0(eta)
# Psi_0 <- pi^(-0.25) # removes the energy scaling
# dj <- 0.01
# dt <- 1
# 
# reCon <- (dj*dt^(0.5))/(C_delta*Psi_0) # Constant for reconstruction of wavelets


mei_wt <- wt(cbind(ann_mei$year, ann_mei$mean_val), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)

ce_wt <- wt(cbind(ann_mei$year, surr_ce), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
aaft_wt <- wt(cbind(ann_mei$year, surr_aaft), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
phase_wt <- wt(cbind(ann_mei$year, surr_phase), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
dh_wt <- wt(cbind(ann_mei$year, surr_dh), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)

old<- par(mfrow=c(5,1), mar=c(1,2,1,1))
plot(mei_wt)
plot(ce_wt)
plot(aaft_wt)
plot(phase_wt)
plot(dh_wt)
par(old)

