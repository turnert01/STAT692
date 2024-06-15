#Load packages
#Some of these were not used
library(ape)
library(geiger)
library(treeplyr)
library(caper)
library(tidyverse)
library(readxl)
library(MASS)
library(stringr)
library(gridExtra)
library(ggpubr)
library(grDevices)
library(colorspace)
library(jtools)

#################################################
#################################################
#SHEET 1: EYE DISCRIMINABILITY DATA
#################################################
#################################################

##Load data
ed <- read_excel("Project 1.xlsx",
                 sheet=1)

##Load in phylogenic tree
tree = read.nexus("tree1.nex")
str(tree)

##Create a new column with species names
ed$tiplabel = gsub(" ", "_", ed$Species)

##Check that the names in the data and the tree match
name.check(phy = tree, data = ed, 
           data.names = ed$tiplabel)

##Combine tree and data
comb <- make.treedata(tree = tree,  data = ed, 
                      name_column = "tiplabel")


##Replace lost species names
comb$dat$tiplabel <- comb$phy$tip.label

##Complete data
completedata = as.data.frame(comb$dat)

##Final tree
finaltree = comb$phy

##Create comparative data for use in PGLS model
sheet1 <- comparative.data(phy = finaltree, data = completedata, 
                       names.col = tiplabel, vcv = TRUE, 
                       na.omit = FALSE, warn.dropped = TRUE)

##Plots of ab, bc, and cz distance  against
##conspecific aggression

plot1=ggplot(data=ed,
             aes(x=Conspecific_Aggression,
                 y=`ab distance`))+
  geom_point(color="maroon")+
  xlab("Conspecific Aggression")+
  ylab("AB Distance")+
  stat_cor(aes(label = ..r.label..),
          label.x = 12)

plot2=ggplot(data=ed,
             aes(x=Conspecific_Aggression,
                 y=`bc distance`))+
  geom_point(color="maroon")+
  xlab("Conspecific Aggression")+
  ylab("BC Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 12)

plot3 = ggplot(data=ed,
               aes(x=Conspecific_Aggression,
                   y=`cz distance`))+
  geom_point(color="maroon")+
  xlab("Conspecific Aggression")+
  ylab("CZ Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 12)

grid.arrange(plot1,
             plot2,
             plot3,
             nrow=1,
             top="Measures of Eye Contrast vs. Conspecific Aggression")

##Plots of ab, bc, and cz distance  against
##canine dimorphism 

plot1=ggplot(data=ed,
             aes(x=Canine_Dimorphism,
                 y=`ab distance`))+
  geom_point(color="maroon")+
  xlab("Canine Dimorphism")+
  ylab("AB Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 3.75)

plot2=ggplot(data=ed,
             aes(x=Canine_Dimorphism,
                 y=`bc distance`))+
  geom_point(color="maroon")+
  xlab("Canine Dimorphism")+
  ylab("BC Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 4)

plot3 = ggplot(data=ed,
               aes(x=Canine_Dimorphism,
                   y=`cz distance`))+
  geom_point(color="maroon")+
  xlab("Canine Dimorphism")+
  ylab("CZ Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 4)

grid.arrange(plot1,
             plot2,
             plot3,
             nrow=1,
             top="Measures of Eye Contrast vs. Canine Dimorphism")

##Plots of ab, bc, and cz distance  against
##group size 
plot1=ggplot(data=ed,
             aes(x=Social_Group_Size,
                 y=`ab distance`))+
  geom_point(color="maroon")+
  xlab("Social Group Size")+
  ylab("AB Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 70)

plot2=ggplot(data=ed,
             aes(x=Social_Group_Size,
                 y=`bc distance`))+
  geom_point(color="maroon")+
  xlab("Social Group Size")+
  ylab("BC Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 70)

plot3 = ggplot(data=ed,
               aes(x=Social_Group_Size,
                   y=`cz distance`))+
  geom_point(color="maroon")+
  xlab("Social Group Size")+
  ylab("CZ Distance")+
  stat_cor(aes(label = ..r.label..),
           label.x = 70)

grid.arrange(plot1,
             plot2,
             plot3,
             nrow=1,
             top="Measures of Eye Contrast vs. Social Group Size")


##Data analysis 

##There will be 3 models

###Rename variables
colnames(sheet1$data)[5]="abdistance"
colnames(sheet1$data)[7]="bcdistance"
colnames(sheet1$data)[9]="czdistance"

####AB Distance
abdistance <- pgls(sqrt(abdistance)~log(Conspecific_Aggression)+log(Canine_Dimorphism)+log(Social_Group_Size), 
             data = sheet1,
             lambda = "ML")

par(mfrow=c(2,2))
plot(abdistance)
shapiro.test(abdistance$residuals)
par(mfrow=c(1,1))

summary(abdistance)

####BC Distance
bcdistance <- pgls(sqrt(bcdistance)~log(Conspecific_Aggression)+log(Canine_Dimorphism)+log(Social_Group_Size), 
                   data = sheet1,
                   lambda = "ML")

par(mfrow=c(2,2))
plot(bcdistance)
shapiro.test(bcdistance$residuals)
par(mfrow=c(1,1))

summary(bcdistance)

####CZ Distance
czdistance <- pgls(sqrt(czdistance)~log(Conspecific_Aggression)+log(Canine_Dimorphism)+log(Social_Group_Size), 
                   data = sheet1,
                   lambda = "ML")

par(mfrow=c(2,2))
plot(czdistance)
shapiro.test(czdistance$residuals)
par(mfrow=c(1,1))

summary(czdistance)


#################################################
#################################################
#SHEET 2: ANALYSIS OF WHR AND BODY MASS
#################################################
#################################################

##Load data
ed <- read_excel("Project 1.xlsx",
                 sheet=2)

##Load in phylogenic tree
tree = read.nexus("tree1.nex")
str(tree)

##Create a new column with species names
ed$tiplabel = gsub(" ", "_", ed$Species)

##Check that the names in the data and the tree match
name.check(phy = tree, data = ed, 
           data.names = ed$tiplabel)

##Combine tree and data
comb <- make.treedata(tree = tree,  data = ed, 
                      name_column = "tiplabel")


##Replace lost species names
comb$dat$tiplabel <- comb$phy$tip.label

##Complete data
completedata = as.data.frame(comb$dat)

##Final tree
finaltree = comb$phy

#Create comparative data for use in PGLS model
sheet2 <- comparative.data(phy = finaltree, data = completedata, 
                           names.col = tiplabel, vcv = TRUE, 
                           na.omit = FALSE, warn.dropped = TRUE)

##Scatterplot of WHR vs. Body Mass
plot1=ggplot(data=ed,
             aes(x=Body_Mass,
                 y=WHR))+
  geom_point(color="maroon")+
  xlab("Body Mass")+
  ylab("WHR")+
  ggtitle("Width-to-Height Ratio vs. Body Mass")+
  stat_cor(aes(label = ..r.label..),
           label.x = 40000)

##Data analysis
whr <-pgls(log(WHR)~log(Body_Mass),
           data = sheet2,
           lambda = "ML")


shapiro.test(whr$residuals)
par(mfrow=c(2,2))
plot(whr)
par(mfrow=c(1,1))

summary(whr)

#################################################
#################################################
#SHEET 3: EYE COLOR and LATITUDE
#################################################
#################################################

##Load data
ed<- read_excel("Project 1.xlsx",
                 sheet=3)
a1= ed[ed$roi=="a1",]
b1= ed[ed$roi=="b1",]
c1= ed[ed$roi=="c1",]
z1= ed[ed$roi=="b1",]

##Load in phylogenic tree
tree = read.nexus("tree1.nex")
str(tree)

##Create a new column with species names
a1$tiplabel = gsub(" ", "_", a1$Species)
b1$tiplabel = gsub(" ", "_", b1$Species)
c1$tiplabel = gsub(" ", "_", c1$Species)
z1$tiplabel = gsub(" ", "_", z1$Species)

##Check that the names in the data and the tree match
name.check(phy = tree, data = a1, 
           data.names = a1$tiplabel)
name.check(phy = tree, data = b1, 
           data.names = b1$tiplabel)
name.check(phy = tree, data = c1, 
           data.names = c1$tiplabel)
name.check(phy = tree, data = z1, 
           data.names = z1$tiplabel)

##Combine tree and data
comba1 <- make.treedata(tree = tree,  data = a1, 
                      name_column = "tiplabel")
combb1 <- make.treedata(tree = tree,  data = b1, 
                        name_column = "tiplabel")
combc1 <- make.treedata(tree = tree,  data = c1, 
                        name_column = "tiplabel")
combz1 <- make.treedata(tree = tree,  data = z1, 
                        name_column = "tiplabel")


##Replace lost species names
comba1$dat$tiplabel <- comba1$phy$tip.label
combb1$dat$tiplabel <- combb1$phy$tip.label
combc1$dat$tiplabel <- combc1$phy$tip.label
combz1$dat$tiplabel <- combz1$phy$tip.label

##Complete data
completedataa1 = as.data.frame(comba1$dat)
completedatab1 = as.data.frame(combb1$dat)
completedatac1 = as.data.frame(combc1$dat)
completedataz1 = as.data.frame(combz1$dat)

##Final tree
finaltreea1 = comba1$phy
finaltreeb1 = combb1$phy
finaltreec1 = combc1$phy
finaltreez1 = combz1$phy

##Create comparative data for use in PGLS model
sheet3a1 <- comparative.data(phy = finaltree, data = completedataa1, 
                           names.col = tiplabel, vcv = TRUE, 
                           na.omit = FALSE, warn.dropped = TRUE)
sheet3b1 <- comparative.data(phy = finaltree, data = completedatab1, 
                             names.col = tiplabel, vcv = TRUE, 
                             na.omit = FALSE, warn.dropped = TRUE)
sheet3c1 <- comparative.data(phy = finaltree, data = completedatac1, 
                             names.col = tiplabel, vcv = TRUE, 
                             na.omit = FALSE, warn.dropped = TRUE)
sheet3z1 <- comparative.data(phy = finaltree, data = completedataz1, 
                             names.col = tiplabel, vcv = TRUE, 
                             na.omit = FALSE, warn.dropped = TRUE)


##Plots color for each region of interest
CEILa1 = LAB(L=a1$L,
             A=a1$a,
             B=a1$b)
plot(CEILa1)

CEILb1 = LAB(L=b1$L,
             A=b1$a,
             B=b1$b)
plot(CEILb1)

CEILc1 = LAB(L=c1$L,
             A=c1$a,
             B=c1$b)
plot(CEILc1)

CEILz1 = LAB(L=z1$L,
             A=z1$a,
             B=z1$b)
plot(CEILz1)


##Data analysis
##Regress Lightness on ABS Latitude
La1 <-pgls(log(L)~log(ABSLatitude),
           data = sheet3a1,
           lambda = "ML")


shapiro.test(La1$residuals)
par(mfrow=c(2,2))
plot(La1)
par(mfrow=c(1,1))

summary(La1)

Lb1 <-pgls(log(L)~log(ABSLatitude),
           data = sheet3b1,
           lambda = "ML")


shapiro.test(Lb1$residuals)
par(mfrow=c(2,2))
plot(Lb1)
par(mfrow=c(1,1))

summary(Lb1)

Lc1 <-pgls(log(L)~log(ABSLatitude),
           data = sheet3c1,
           lambda = "ML")


shapiro.test(Lc1$residuals)
par(mfrow=c(2,2))
plot(Lc1)
par(mfrow=c(1,1))

summary(Lc1)

Lz1 <-pgls(log(L)~log(ABSLatitude),
           data = sheet3z1,
           lambda = "ML")


shapiro.test(Lz1$residuals)
par(mfrow=c(2,2))
plot(Lz1)
par(mfrow=c(1,1))

summary(Lz1)

##Regress a chromacity on ABS Latitude
aa1 <-pgls(a~log(ABSLatitude),
           data = sheet3a1,
           lambda = "ML")


shapiro.test(aa1$residuals)
par(mfrow=c(2,2))
plot(aa1)
par(mfrow=c(1,1))

summary(aa1)

ab1 <-pgls(a~ABSLatitude,
           data = sheet3b1,
           lambda = "ML")


shapiro.test(ab1$residuals)
par(mfrow=c(2,2))
plot(ab1)
par(mfrow=c(1,1))

summary(ab1)

ac1 <-pgls(a~ABSLatitude,
           data = sheet3c1,
           lambda = "ML")


shapiro.test(ac1$residuals)
par(mfrow=c(2,2))
plot(ac1)
par(mfrow=c(1,1))

summary(ac1)

az1 <-pgls(a~ABSLatitude,
           data = sheet3z1,
           lambda = "ML")


shapiro.test(az1$residuals)
par(mfrow=c(2,2))
plot(az1)
par(mfrow=c(1,1))

summary(az1)


##Regress b chromacity on ABS Latitude
ba1 <-pgls(b~ABSLatitude,
           data = sheet3a1,
           lambda = "ML")


shapiro.test(ba1$residuals)
par(mfrow=c(2,2))
plot(ba1)
par(mfrow=c(1,1))

summary(ba1)

bb1 <-pgls(b~log(ABSLatitude),
           data = sheet3b1,
           lambda = "ML")


shapiro.test(bb1$residuals)
par(mfrow=c(2,2))
plot(bb1)
par(mfrow=c(1,1))

summary(bb1)

bc1 <-pgls(b~ABSLatitude,
           data = sheet3c1,
           lambda = "ML")


shapiro.test(bc1$residuals)
par(mfrow=c(2,2))
plot(ac1)
par(mfrow=c(1,1))

summary(bc1)

bz1 <-pgls(b~log(ABSLatitude),
           data = sheet3z1,
           lambda = "ML")


shapiro.test(bz1$residuals)
par(mfrow=c(2,2))
plot(bz1)
par(mfrow=c(1,1))

summary(bz1)
