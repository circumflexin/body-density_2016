### plot results of the selected model

setwd("/Users/pm29/Documents/Rmodels/BodyDensity/workshop/")
dataDir<-"/Users/pm29/Documents/Rmodels/BodyDensity/workshop/data/"
modelDir<-"/Users/pm29/Documents/Rmodels/BodyDensity/workshop/"

modelName<-c("model(12)")   ### insert the selected model
m<-1

# load & filter data -> fitName
load("all_whales.Rd")
fitName <- "_f_pitch30_depth100"
filterBool <- abs(tab$mean.pitch)>30 & tab$dive.max.depth>100 & abs(tab$phase)>0 & is.na(tab$acceleration)==F &
  is.na(tab$mean.speed)==F & tab$r.roll>0.9  & is.na(tab$DswGG)==F
#& tab$mean.depth>200 & tab$sonar==0 & tab$vertglide==1
tab <- tab[filterBool,]       

whales <- unique(tab$whale.id)

# load model & data & filter data
load(paste(modelDir, modelName[m], fitName, ".Rd", sep=""))

library(R2jags)
#library(RColorBrewer)

############################### COLOURS
whaleCol = seq(1:length(whales))

rbPal<-colorRampPalette(c("purple","blue", "deepskyblue", "cyan", "lawngreen",
                          "yellow", "orange", "orangered", "darkred"))

pitchBreaks <- seq(-90, 90, 10)
mycols <- rbPal(length(pitchBreaks))[as.numeric(cut(tab$mean.pitch,breaks = pitchBreaks))]

#############################

pdf(paste("plot_best_", modelName[m], ".pdf", sep=""), width=5,height=5)
options(graphics.record=TRUE)


# body density
test<-density(fit$BUGSoutput$sims.list$body.density[,1])
plot(1,1,xlim=c(1020, 1048),
     ylim=c(0,max(test$y)), col=NA,
     xlab="body density", ylab="posterior density")
grid()
for(w in 1:length(whales)) {
  lines(density(fit$BUGSoutput$sims.list$body.density[,w]),col=whaleCol[w])
}
numG<-c(1)
for(w in 1:length(whales)) {
  filterW<-tab$whale.num==w
  tabW<-tab[filterW,]
  numG[w]<-length(tabW$acceleration)
}
legend(x=1040, y=max(test$y), legend=paste(whales[1:length(whales)], "[", numG[1:length(numG)], "]", sep=""),
      lty=1, col=whaleCol, cex=0.5)

if(length(fit$BUGSoutput$sims.list$body.density.g)>0){
  lines(density(fit$BUGSoutput$sims.list$body.density.g), col=c("#00000050"), lwd=7)
  legend(x=1040, y=max(test$y)*0.1, "global mean", lty=1, lwd=7, col=c("#00000050"), cex=0.6)
}


# Vair.d (for indiv-specific)
if(ncol(fit$BUGSoutput$sims.list$Vair.d)==length(whales)){
  test<-density(fit$BUGSoutput$sims.list$Vair.d[,1])
  plot(1,1,xlim=c(2,50), ylim=c(0,max(test$y)*2), col=NA,
       xlab="Vair", ylab="posterior density")
  grid()
  for(w in 1:length(whales)) {
    lines(density(fit$BUGSoutput$sims.list$Vair.d[,w]), col=whaleCol[w])
  }
  lines(density(fit$BUGSoutput$sims.list$Vair), col=c("#80808030"), lwd=3)
  legend(x=35, y=max(test$y)*2, legend=whales, lty=1, col=w, cex=0.7)
  legend(x=35, y=max(test$y)*0.8, "global", lty=1, col=c("#80808030"), lwd=3, cex=0.7)
}

# Vair global
test<-density(fit$BUGSoutput$sims.list$Vair)
plot(1,1,xlim=c(2,50), ylim=c(0,max(test$y)*1.5), col=NA,
     xlab="Vair", ylab="posterior density",
     main=paste("mean global Vair = ",
                signif(mean(fit$BUGSoutput$sims.list$Vair), digits=3), " (sd = ",
                signif(sd(fit$BUGSoutput$sims.list$Vair), digits=3), ")", sep=""))
lines(density(fit$BUGSoutput$sims.list$Vair), col=c("#00000040"), lwd=5)
legend(x=35, y=max(test$y), "global Vair", lty=1, col=c("#00000040"), lwd=5, cex=0.7)


# CdAM
test<-density(fit$BUGSoutput$sims.list$CdAM[,1])
plot(1,1,xlim=c(0,80), ylim=c(0,max(test$y)), col=NA,
     xlab="CdAM", ylab="posterior density")
grid()
for(w in 1:length(whales)) {
  lines(density(fit$BUGSoutput$sims.list$CdAM[,w]),
        col=whaleCol[w])
}
legend(x=60, y=max(test$y), legend=whales, lty= 1, col=whaleCol, cex=0.5)
if(length(fit$BUGSoutput$sims.list$CdAM.g)>0){
  lines(density(fit$BUGSoutput$sims.list$CdAM.g), col=c("#00000040"), lwd=5)
  legend(x=60, y=max(test$y)*0.1, "global mean", lty=1, lwd=5, col=c("#00000040"), cex=0.5)
}

# compressibility
plot(density(fit$BUGSoutput$sims.list$compr), col="lightgrey", lwd=5,
     xlab="Compressibility", ylab="posterior density", xlim=c(0.1, 0.8), #ylim=c(0, 5),
     main=paste("mean compr = ",
                signif(mean(fit$BUGSoutput$sims.list$compr), digits=3), " (sd = ",
                signif(sd(fit$BUGSoutput$sims.list$compr), digits=3), ")", sep=""))


# body density.g with uniform prior
if(length(fit$BUGSoutput$sims.list$body.density.g)>0){
  test<-density(fit$BUGSoutput$sims.list$body.density.g)
  plot(1,1,xlim=c(800, 1400),
       ylim=c(0,max(test$y)), col=NA,
       xlab="body density", ylab="posterior density",
       main=paste("mean =", signif(fit$BUGSoutput$mean$body.density.g, digit=5), sep=""))
  lines(density(fit$BUGSoutput$sims.list$body.density.g), col=c("#00000040"), lwd=3)
  legend(x=1200, y=max(test$y)*1, "global mean", lty=1, lwd=3, col=c("#00000040"), cex=0.6)
  legend(x=1200, y=max(test$y)*0.2, "prior range", lty=1, lwd=5, col=c("#00008B50"), cex=0.6)
  segments(800, 0, 1400, 0, col=c("#00008B50"), lwd=10)
}


# CdAM.g with uniform prior
if(length(fit$BUGSoutput$sims.list$CdAM.g)>0){
  test<-density(fit$BUGSoutput$sims.list$CdAM.g)
  plot(1,1,xlim=c(0, 80),
       ylim=c(0,max(test$y)*1.2), col=NA,
       xlab="CdAM (E-06)", ylab="posterior density",
       main=paste("mean =", signif(fit$BUGSoutput$mean$CdAM.g, digit=3), sep=""))
  lines(density(fit$BUGSoutput$sims.list$CdAM.g), col=c("#00000040"), lwd=5)
  legend(x=30, y=max(test$y)*1.2, "global mean", lty=1, lwd=3, col=c("#00000040"), cex=0.6)
  legend(x=30, y=max(test$y)*0.2, "prior range", lty=1, lwd=5, col=c("#00008B50"), cex=0.6)
  xx<-seq(5, 20, 0.1)
  lines(xx, dnorm(xx, 10, 2), col=c("#00008B50"), lwd=5)
}


# gloal Vair with uniform prior
if(length(fit$BUGSoutput$sims.list$Vair)>0){
  test<-density(fit$BUGSoutput$sims.list$Vair)
  plot(1,1,xlim=c(2, 50),
       ylim=c(0,max(test$y)*1.5), col=NA,
       xlab="Vair", ylab="posterior density",
       main=paste("mean =", signif(fit$BUGSoutput$mean$Vair, digit=3), sep=""))
  lines(density(fit$BUGSoutput$sims.list$Vair), col=c("#00000040"), lwd=5)
  legend(x=35, y=max(test$y)*1, "global mean", lty=1, lwd=3, col=c("#00000040"), cex=0.6)
  legend(x=35, y=max(test$y)*0.2, "prior range", lty=1, lwd=5, col=c("#00008B50"), cex=0.6)
  segments(5, 0, 50, 0, col=c("#00008B50"), lwd=10)
}


# compr with uniform prior
# compressibility
test<-density(fit$BUGSoutput$sims.list$compr)
plot(density(fit$BUGSoutput$sims.list$compr), col=c("#00000040"), lwd=5,
     xlim=c(0, 0.8), ylim=c(0, max(test$y)), 
     xlab="Compressibility", ylab="posterior density",
     main=paste("mean =", signif(fit$BUGSoutput$mean$compr, digit=3), sep=""))
legend(x=0.6, y=max(test$y)*1, "global mean", lty=1, lwd=3, col=c("#00000040"), cex=0.6)
legend(x=0.6, y=max(test$y)*0.2, "prior range", lty=1, lwd=5, col=c("#00008B50"), cex=0.6)
segments(0.1, 0, 0.7, 0, col=c("#00008B50"), lwd=10)



par(mfrow=c(1,1))
# measured tau
test<-density(sp.data$tau)
plot(1,1,xlim=c(min(sp.data$tau),max(sp.data$tau)), ylim=c(0,max(test$y)*2), col=NA,
     xlab="TAU used in the model", ylab="density")
grid()
for(w in 1:length(whales)) {  
  wBool <- sp.data$whale.id==w  
  lines(density(sp.data$tau[wBool]), col=whaleCol[w], lty=1)
}
legend(max(sp.data$tau)*0.7, max(test$y)*8, legend=whales, lty=1, col=whaleCol, cex=0.5)



#plot(sp.data$a, fit$BUGSoutput$mean$a.mu, ylab="", xlab="", col=whaleCol[sp.data$whale.id], cex=0.8)
#mtext("observed acceleration", side=1, line=3, cex=1.2)
#mtext("posterior mean", side=2, line=3, cex=1.2)
#abline(0,1)
#grid()

#plot(sp.data$mean.speed, sp.data$a,
#     ylab="", xlab="", 
#     main="")
#points(sp.data$mean.speed, fit$BUGSoutput$mean$a.mu, col=whaleCol[sp.data$whale.id], cex=0.8)
#grid()
#mtext("speed (m/s)", side=1, line=3, cex=1.2)
#mtext("acceleration", side=2, line=3, cex=1.2)
#mtext("observed and fitted (in colour)", side=3, line=1, cex=1.2, font=2)




## observed acceleration vs glide depth
# colour coded by speed
rbPal2<-colorRampPalette(c("darkblue", "blue", "deepskyblue", "cyan", "lawngreen",
                           "yellow", "orange", "darkred"))

speedBreaks<-seq(0, 4, 0.01)      # speed is color coded by 0.01
spBreakLegend<-seq(0, 4, 0.5)     # data to make legend
SPcols <- rbPal2(length(speedBreaks))[as.numeric(cut(sp.data$mean.speed,breaks = speedBreaks))]

plot(sp.data$depth, sp.data$a, ylab="", xlab="",main="", col=SPcols, cex=0.5, pch=1)
grid()
mtext("depth (m)", side=1, line=3, cex=1.2)
mtext("obs. acceleration", side=2, line=3, cex=1.2)
mtext("obs. acceleration vs depth (coloured by speed)", side=3, line=1, cex=1.2, font=2)
legend(max(sp.data$depth)*0.9,max(sp.data$a), legend=spBreakLegend[1:length(spBreakLegend)], 
       fill=rbPal2(length(spBreakLegend))[1:length(spBreakLegend)][1:length(spBreakLegend)], cex=0.6)


#plot(sp.data$depth, sp.data$a, ylab="", xlab="",main="", col=pitchCols1, cex=0.8, pch=1)
#plot(sp.data$depth, sp.data$a, ylab="", xlab="",main="", col="grey", cex=0.5, pch=1)
#grid()
#points(sp.data$depth, fit$BUGSoutput$mean$a.mu, col=SPcols, pch=3, cex=0.5)
#mtext("depth (m)", side=1, line=3, cex=1.2)
#mtext("acceleration", side=2, line=3, cex=1.2)
#mtext("obs + fitted acceleration vs depth", side=3, line=1, cex=1.2, font=2)
#legend(max(sp.data$depth)*0.9,max(sp.data$a), legend=spBreakLegend[1:length(spBreakLegend)], 
#       fill=rbPal2(length(spBreakLegend))[1:length(spBreakLegend)][1:length(spBreakLegend)], cex=0.6)



## observed data alon, color coded by pitch
#rbPal<-colorRampPalette(c("purple","blue", "deepskyblue", "cyan", "lawngreen",
#                          "yellow", "orange", "orangered", "darkred"))
#pitchBreaks <- seq(-90, 90, 10)
#sp.data$mean.pitch<-asin(sp.data$sin.pitch)*180/pi      # add mean.pitch to sp.data
#mycols <- rbPal(length(pitchBreaks))[as.numeric(cut(sp.data$mean.pitch,breaks = pitchBreaks))]

#plot(sp.data$mean.speed, sp.data$a, ylab="", xlab="",main="", col=mycols, cex=0.5, pch=1, 
#     xlim=c(0, max(sp.data$mean.speed)*1.2))
#grid()
#mtext("speed (m/s)", side=1, line=3, cex=1.2)
#mtext("obs. acceleration", side=2, line=3, cex=1.2)
#mtext("obs a vs speed", side=3, line=1, cex=1.2, font=2)
#legend(max(sp.data$mean.speed)*1.1,max(sp.data$a), legend=pitchBreaks[-(8:12)], 
#       fill=rbPal(length(pitchBreaks))[(1:20)[-(8:12)]], cex=0.6)


rbPal<-colorRampPalette(c("purple","blue", "deepskyblue", "cyan", "lawngreen",
                          "yellow", "orange", "orangered", "darkred"))
pitchBreaks <- seq(-90, 90, 10)
sp.data$mean.pitch<-asin(sp.data$sin.pitch)*180/pi      # add mean.pitch to sp.data
mycols <- rbPal(length(pitchBreaks))[as.numeric(cut(sp.data$mean.pitch,breaks = pitchBreaks))]

## observed data & estimated a
plot(sp.data$mean.speed, sp.data$a, ylab="", xlab="",main="", 
     col="grey", cex=0.5, pch=1, xlim=c(0, max(sp.data$mean.speed)*1.2))
points(sp.data$mean.speed, fit$BUGSoutput$mean$a.mu, col=mycols, pch=3, cex=0.5)
grid()
mtext("speed (m/s)", side=1, line=3, cex=1.2)
mtext("acceleration", side=2, line=3, cex=1.2)
mtext("obs (grey) + fitted (colored by pitch) acc. vs speed", side=3, line=1, cex=1.2, font=2)
legend(max(sp.data$mean.speed)*1.1,max(sp.data$a), legend=pitchBreaks[-(8:12)], 
       fill=rbPal(length(pitchBreaks))[(1:20)[-(8:12)]], cex=0.6)

#for(w in 1:length(whales)) {

#  wBool <- sp.data$whale.id==w

## observed data alon
#  plot(sp.data$mean.speed[wBool], sp.data$a[wBool], ylab="", xlab="",main="", 
#       col="grey", cex=0.5, pch=1)
#  points(sp.data$mean.speed[wBool], fit$BUGSoutput$mean$a.mu[wBool], col=mycols, pch=3, cex=0.5)
#  grid()
#  mtext("speed (m/s)", side=1, line=3, cex=1.2)
#  mtext("acceleration", side=2, line=3, cex=1.2)
#  mtext(whales[w], side=3, line=1, cex=1.2, font=2)

#}



# plot observed A and estimated A with speed
rbPal<-colorRampPalette(c("purple","blue", "deepskyblue", "cyan", "lawngreen",
                          "yellow", "orange", "orangered", "darkred"))
pitchBreaks<-seq(-90, 90, 10)

sp.data$mean.pitch<-asin(sp.data$sin.pitch)*180/pi      # add mean.pitch to sp.data

for(w in 1:length(whales)) {
  
  wBool <- sp.data$whale.id==w
  
  mycols2<-rbPal(length(pitchBreaks))[as.numeric(cut(sp.data$mean.pitch[wBool],breaks = pitchBreaks))]   # new mycolor for filtered data (for each whale)
  
  ## observed data alon
  plot(sp.data$mean.speed[wBool], sp.data$a[wBool], ylab="", xlab="",main="", 
       col="darkgrey", cex=0.8, pch=19, xlim=c(min(sp.data$mean.speed), max(sp.data$mean.speed)+1), 
       ylim=c(min(sp.data$a), max(sp.data$a)))
  points(sp.data$mean.speed[wBool], fit$BUGSoutput$mean$a.mu[wBool], col="black", bg=mycols2, pch=21, cex=0.8)
  grid()
  mtext("speed (m/s)", side=1, line=3, cex=1.2)
  mtext("acceleration", side=2, line=3, cex=1.2)
  mtext(paste(whales[w], " [obs(grey), fitted(coloured by pitch)]", sep=""), side=3, line=1, cex=1.2, font=2)
  legend(max(sp.data$mean.speed),max(sp.data$a), legend=pitchBreaks[-(8:12)], 
         fill=rbPal(length(pitchBreaks))[(1:20)[-(8:12)]], cex=0.7)
}


options(graphics.record=FALSE)      
dev.off() 



