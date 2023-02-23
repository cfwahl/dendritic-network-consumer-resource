

# setup --------------------------------------------------------------------

# clean objects
rm(list = ls())

# load required libraries
library(vegan)
library(sciplot)
library(lmPerm)
library(gplots) 
library(plotrix)
library(tidyverse)

# data --------------------------------------------------------------------

# macroinvert community data, read in file Litterbag_bug_data_2013.csv 
litter.bag <- read.csv(file.choose())

# shredder only data
shredder.litter <- litter.bag[,1:25]

# leaf litter decomposition, read in shredder.decomp.csv
shredders.decomp <- read.csv(file.choose()) 

# shredder richness -------------------------------------------------------

# two way ANOVA source*leaf
summary(richness.aov<-aov(shredder.litter$Shredder.Richness~shredder.litter$Source*shredder.litter$Leaf))

# bargraph
bargraph.CI(x.factor=shredder.litter$Source, group=shredder.litter$Leaf, 
            response=shredder.litter$Shredder.Richness, legend=T, 
            ylab="Consumer Richness (Richness/Leaf Pack)", xlab="Network Location", 
            x.leg=1, y.leg=5, ylim=c(0,5), cex.leg=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.5)
# legend
text(5.5,4.55, "SOV   p-value", adj=1, cex=1)
text(5.5,4.25, "Location   < 0.001", adj=1)
text(5.5,4, "Resource   < 0.001", adj=1)
text(5.47,3.75,  "Location x Resource       0.13", adj=1)

#create subsets of treatments for richness and abundance analysizes 
HW.Ash <- subset(shredder.litter, Source.Leaf=="HW.Ash")
HW.Beech <- subset(shredder.litter, Source.Leaf=="HW.Beech")
MS.Ash <- subset(shredder.litter, Source.Leaf=="MS.Ash")
MS.Beech <- subset(shredder.litter, Source.Leaf=="MS.Beech")

#mean shredder richness for treatments
tapply(shredder.litter$Shredder.Richness, INDEX=list(shredder.litter$Source.Leaf), FUN=mean)

# standaard error of mean for shredder richness HW.Ash
std.error(HW.Ash$Shredder.Richness)

# standaard error of mean for shredder richness HW.Beech
std.error(HW.Beech$Shredder.Richness)

# standaard error of mean for shredder richness MS.Ash
std.error(MS.Ash$Shredder.Richness)

# standaard error of mean for shredder richness MS.Beech
std.error(MS.Beech$Shredder.Richness)

# shredder abundance ------------------------------------------------------

# two way ANOVA source*leaf
summary(richness.aov4 <- aov(shredder.litter$Shredder.Sum ~ shredder.litter$Source*shredder.litter$Leaf))

# bar graph
bargraph.CI(x.factor=shredder.litter$Source, group=shredder.litter$Leaf, 
            response=shredder.litter$Shredder.Sum, legend=T, 
            ylab="Consumer Abundance (# of Consumers/Leaf Pack)", xlab="Network Location", 
            x.leg=1, y.leg=20.5, ylim=c(0,20), cex.leg=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.5)

# legend
text(5.5,19, "SOV    p-value", adj=1, cex=1)
text(5.5,18, "Location   < 0.001", adj=1)
text(5.5,17, "Resource      0.03  ", adj=1)
text(5.5,16,  "Location x Resource      0.26  ", adj=1)

# mean shredder abundance for treatments
tapply(shredder.litter$Shredder.Sum, INDEX=list(shredder.litter$Source.Leaf), FUN=mean)

# standaard error of mean for shredder abundance HW.Ash
std.error(HW.Ash$Shredder.Sum)

# standaard error of mean for shredder abundance HW.Beech
std.error(HW.Beech$Shredder.Sum)

# standaard error of mean for shredder abundance MS.Ash
std.error(MS.Ash$Shredder.Sum)

# standaard error of mean for shredder abundance MS.Beech
std.error(MS.Beech$Shredder.Sum)

# resource loss -----------------------------------------------------------

# remove sites with zero shredders 
shredders.decomp.culled <- shredders.decomp[-c(40,60,82,94,96,101,106),]

#creates dataframe with species, shredder.sum, shredder.richness, and mass.loss
decomp.shredders1 <- data.frame(cbind(shredders.decomp.culled$Species,
                                    shredders.decomp.culled$Shredder.Sum,
                                    shredders.decomp.culled$Shredder.Richness,
                                    shredders.decomp.culled$Mass.loss.per.g.day,
                                    shredders.decomp.culled$HW.MS )) 
#renames columns 
names(decomp.shredders1) <- c("Ash.Beech","Shredder.Abun","Shredder.Richness","Decomp","HW.MS")

# make numeric for analyses
decomp.shredders1 <- decomp.shredders1 %>%
  mutate(Decomp = as.numeric(Decomp),
         Shredder.Abun = as.numeric(Shredder.Abun))

# makes subsets for leaf treatments
decomp.beech1 <- subset(decomp.shredders1, Ash.Beech=="F.grandifolia")
decomp.ash1 <- subset(decomp.shredders1, Ash.Beech=="F.pennsylvanica")

#makes a subset of beech leaves for reaches
decomp.HW.Beech <- subset(decomp.beech1, HW.MS=="Headwater")
decomp.MS.Beech <- subset(decomp.beech1, HW.MS=="Mainstem")

#makes a subset of ash leaves for reaches
decomp.HW.ash <- subset(decomp.ash1, HW.MS=="Headwater")
decomp.MS.ash <- subset(decomp.ash1, HW.MS=="Mainstem")

# Analysis for Ash.HW leaves
summary(lm(decomp.HW.ash$Decomp ~ log(decomp.HW.ash$Shredder.Abun)))

# make the plot 2x2
par(mfrow=c(2,2))

# plots ASH.HW breakdown against shredder abundance 
plot(decomp.HW.ash$Decomp ~ log(decomp.HW.ash$Shredder.Abun), 
     main="A: Headwater and F. pennsylvanica", cex.lab=1.2, cex.main=1.2, 
     cex=1.2, ylab="Resource Loss per Leaf Pack", xlab="",xlim=c(1,4), ylim=c(0,0.025))

# adds line of best fit to plot
abline(lm(decomp.HW.ash$Decomp ~ log(decomp.HW.ash$Shredder.Abun)), lwd=1.7)

# text
rsquarelm2 <- 0.1646
text(1.75,0.008, bquote( R^2 == .(rs), list(rs=round(rsquarelm2,4))), cex=1.1)
text(1.75,0.0045,"Slope = 0.00238", cex=1.1)
text(1.75,0.0019, "p-value = 0.04", cex=1.1)

# Analysis for BEECH.HW leaves
summary(lm(decomp.HW.Beech$Decomp ~ log(decomp.HW.Beech$Shredder.Abun)))

#plots BEECH.HW breakdown against logged shredder abundance, log transformed because data was skewed 
plot(decomp.HW.Beech $Decomp ~ log(decomp.HW.Beech $Shredder.Abun), 
     main="B: Headwater and F. grandifolia", cex.lab=1.2, cex.main=1.2, cex=1.2, 
     ylab="", xlab="",xlim=c(0,4), ylim=c(0,0.006))

#adds line of best fit to plot
abline(lm(decomp.HW.Beech$Decomp ~ log(decomp.HW.Beech$Shredder.Abun)), lwd=1.7)

# text
rsquarelm2 <- 0.1795
text(1,0.0019, bquote( R^2 == .(rs), list(rs=round(rsquarelm2,4))), cex=1.1)
text(1,0.0011,"Slope = 0.00025", cex=1.1)
text(1,0.0004, "p-value = 0.03", cex=1.1)

# Analysis for ASH.MS leaves
summary(lm(decomp.MS.ash$Decomp ~ log(decomp.MS.ash$Shredder.Abun)))

#plots ASH.MS breakdown against shredder abundance 
plot(decomp.MS.ash$Decomp ~ log(decomp.MS.ash$Shredder.Abun), 
     main="C: Mainstem and F. pennsylvanica", cex.main=1.2, cex.lab=1.2, cex=1.2, 
     ylab="Resource Loss per Leaf Pack", xlab="log(Shredder Abundance per Leaf Pack)",
     xlim=c(0,3), ylim=c(0,0.025))

#adds line of best fit to plot
abline(lm(decomp.MS.ash $Decomp~log(decomp.MS.ash $Shredder.Abun)), lwd=1.7)

rsquarelm2 <- 0.00144
text(0.75,0.008, bquote( R^2 == .(rs), list(rs=round(rsquarelm2,5))), cex=1.1)
text(0.8,0.0045,"Slope = -0.00016", cex=1.1)
text(0.75,0.0015, "p-value = 0.81", cex=1.1)

### Analysis for BEECH.MS leaves ###
summary(lm(decomp.MS.Beech$Decomp~log(decomp.MS.Beech$Shredder.Abun)))

#plots BEECH.MS breakdown against shredder abundance 
plot(decomp.MS.Beech$Decomp ~ log(decomp.MS.Beech$Shredder.Abun), 
     main="D: Mainstem and F. grandifolia", cex.lab=1.2, cex.main=1.2, cex=1.2, ylab="", 
     xlab="log(Shredder Abundance per Leaf Pack)",xlim=c(0,2.5), ylim=c(0,0.006))

#adds line of best fit to plot
abline(lm(decomp.MS.Beech$Decomp ~ log(decomp.MS.Beech$Shredder.Abun)), lwd=1.7)

rsquarelm2 <- 0.4649
text(0.7,0.0021, bquote( R^2 == .(rs), list(rs=round(rsquarelm2,4))), cex=1.1)
text(0.7,0.0013,"Slope = -0.00091", cex=1.1)
text(0.7,0.00066, "p-value < 0.001", cex=1.1)

# mass loss ---------------------------------------------------------------

# makes decomposition file for treatments
decomp<-shredders.decomp.culled[,1:6]

# attach file
attach(decomp)

# make into factors
species<-factor(Species)
HW.MS<-factor(decomp$HW.MS)
stream<-factor(decomp$Stream)

summary(aov(Mass.loss.per.g.day ~ HW.MS*Species, data=shredders.decomp))

# make the plot 1x1 again
par(mfrow=c(1,1))

# bar graph
bargraph.CI(x.factor=HW.MS, group=Species, response=Mass.loss.per.g.day, legend=T, 
            ylab="Resource Loss per Leaf Pack", xlab="Network Location", x.leg=1, 
            y.leg=0.025, ylim=c(0,0.025), cex.leg=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.5)

# legend
text(5.5,0.023, "SOV      p-value", adj=1, cex=1)
text(5.5,0.022, "Location     < 0.01 ", adj=1)
text(5.5,0.021, "Resource    < 0.001", adj=1)
text(5.5,0.020,  "Location x Resource     < 0.01 ", adj=1)

# mean mass loss for treatments
tapply(decomp$Mass.loss.per.g.day, INDEX=list(decomp$Source.Leaf), FUN=mean)

# standaard error of mean for breakdown HW.Ash
std.error(decomp.HW.ash$Decomp)

# standaard error of mean for breakdown HW.Beech
std.error(decomp.HW.Beech$Decomp)

# standaard error of mean for breakdown MS.Ash
std.error(decomp.MS.ash$Decomp)

# standaard error of mean for breakdown MS.Beech
std.error(decomp.MS.Beech$Decomp)

