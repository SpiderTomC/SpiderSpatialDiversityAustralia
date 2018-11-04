EastSpiders = read.csv( #filepath of Spider Records East coast# , stringasfactors=F)
#latitudinal and longitudinal bounds
min(EastSpiders$decimalLatitude)
max(EastSpiders$decimalLatitude)
min(EastSpiders$decimalLongitude)
max(EastSpiders$decimalLongitude)

#chao1 function#
Chao1<-function(n)        {
  
  return(length(n) + length(which(n == 1)) * (length(which(n == 1)) - 1) / (2 *  (length(which(n == 2)) + 1)))
  
}
#squares function for diversity estimation
squares<-function(n)	{
  S <- length(n)
  N <- sum(n)
  s1 <- length(which(n == 1))
  if (s1 == S)
    return(NA)
  return(S + s1^2 * sum(n^2) / (N^2 - S * s1))
}

#creation of cell latitude and longitude sequences, using the maximum and  minimum determined previously#
cellLats= seq(-12,-39,by=-0.25)
cellLongs= seq(138, 159.5,by= 0.25)

#creation of cell band arrays. Quarter degree cells balance accuracy with sample size. Both Chao1 and Squares#
#are run within the loop, giving observed richness and two estimates#

Longrichness <- array()
Longcounts <- array()
Cellrichness<-array()
Cellcounts<- array()
cellsquare<- array()
squarecell<-array()
squareslat<-array()
squareslong<-array()
chaocell<-array()
latcount=1
longcount=1
for (i in cellLongs)
{ eastspidlong<- EastSpiders[which(EastSpiders$decimalLongitude<=i & EastSpiders$decimalLongitude >i-.25),]
Longrichness[longcount]<- length(unique(eastspidlong$specificEpithet))
Longcounts[longcount]<- nrow(eastspidlong)
for (m in cellLats)
{eastcell<- eastspidlong[which(eastspidlong$decimalLatitude<=m & eastspidlong$decimalLatitude>m-.25),]
eastcell<- eastcell[which(eastcell$specificEpithet!=""),]
eastcelltable<- sort(table(droplevels(eastcell$specificEpithet,decreasing = T)))
eastcelltable<-eastcelltable[2:length(eastcelltable)]
Cellrichness[latcount]<-length(unique(eastcell$specificEpithet))
Cellcounts[latcount]<- nrow(eastcell)
squarecell[latcount] <- squares(eastcelltable)
chaocell[latcount] <- Chao1(eastcelltable)
squareslat[latcount]=m
squareslong[latcount]=i

latcount<-latcount+1

}
longcount<-longcount+1
}

#creation of richness matrices, flipping column and then transposing is necessary for correct orientation#
#of the contour map#
#raw richness#
spidercell1=matrix(Cellrichness,nrow=length(cellLongs),ncol=length(cellLats))
spidercell1[spidercell1==0]<-NA
spidercell1<-spidercell1[,ncol(spidercell1)[1]:1]
spidercell1=t(spidercell1)
#chao1 richness#
chao1cells= matrix(chao1cell,nrow=length(cellLats),ncol=length(cellLongs))
chao1cells[is.na(chao1cells)]<-0
chao1cells<-chao1cells[,ncol(chao1cells)[1]:1]
chao1cells=t(chao1cells)
#squares richness#
squaresmat= matrix(squarecell,nrow=length(cellLats),ncol=length(cellLongs))
squaresmat[squaresmat==0]<-NA
squaresmat=squaresmat[nrow(squaresmat)[1]:1,]
squaresmat=t(squaresmat)

#Worldclim data extraction and assembly#
#requires sp and raster package#
library(raster)
library(sp)
totalclim <-getData("worldclim",var="bio",res=10)
climcoord <- SpatialPoints(cbind(squareslong,squareslat))
eastclim <- totalclim[[c(1,12)]]
names(eastclim) <- c("Temp","Prec")
climvalues <- extract(eastclim,climcoord)
eastdry <- totalclim[[c(15)]]
dryvalues<- extract(eastdry,climcoord)
names(dryvalues)<- "Prec Seasonality"
eastalt <-getData('alt',country= 'AUS')
altvalues <-extract(eastalt,climcoord)
names(altvalues)<- "Elevation"
tempseason= totalclim [[c(4)]]
tempseasonvalues= extract(tempseason,climcoord)
names(tempseasonvalues)<- "Temp Seasonality"
dryquarter= totalclim[[c(17)]]
dryquartervalues= extract(dryquarter,climcoord)
names(dryquartervalues)<- "Prec Dry Quarter"

#alteration of environmental arrays to matrices
spiderprec=matrix(climvalues[,'Prec'],nrow=length(cellLats),ncol=length(cellLongs))
spiderprec<-spiderprec[nrow(spiderprec)[1]:1,]
spiderprec=t(spiderprec)
spidertemp=matrix(climvalues[,'Temp'],nrow=length(cellLats),ncol=length(cellLongs))
spidertemp<-spidertemp[nrow(spidertemp)[1]:1,]
spidertemp=t(spidertemp)
spideralt=matrix(altvalues,nrow=length(cellLats),ncol=length(cellLongs))
spideralt=spideralt[nrow(spideralt)[1]:1,]
spideralt=t(spideralt)
spiderdry= matrix (dryvalues,nrow=length(cellLats),ncol=length(cellLongs))
spiderdry=spiderdry[nrow(spiderdry)[1]:1,]
spiderdry=t(spiderdry)
tempseasonality= matrix (tempseasonvalues, nrow=length(cellLats),ncol=length(cellLongs))
tempseasonality= tempseasonality[nrow(tempseasonality)[1]:1,]
tempseasonality=t(tempseasonality)
tempseasonality= tempseasonality/100
dryquartermat= matrix(dryquartervalues,nrow=length(cellLats),ncol=length(cellLongs))
dryquartermat=dryquartermat[nrow(dryquartermat)[1]:1,]
dryquartermat=t(dryquartermat)

#Exploratory factor analysis#
squarefact= cbind(log(squarecell),climvalues,altvalues,dryvalues,dryquartervalues,tempseasonvalues,squareslat,Cellcounts)
squarefact= squarefact[which(squarefact[,9]>49),]
squarefact= na.omit(squarefact)
library(psych)
fa(scale(squarefact[,1:8]),nfactors=8)
fa(scale(squarefact[,1:8]),nfactors=5, rotate= 'varimax')
fa(scale(squarefact[,1:8]),nfactors=8,fm="wls",rotate= 'varimax')

#creation of contour plots#
png('Fig 01 Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spidercell1)),y=seq(-39,-12,length.out = ncol(spidercell1)),spidercell1,col = spectral,levels = c(0,10,40,70,100,130,160,190,220,250,280),key.axes=axis(4,c(0,10,40,70,100,130,160,190,220,250,280)))
dev.off()
png('Fig 01B Chao1 Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(chao1cells)),y=seq(-39,-12,length.out = ncol(chao1cells)),chao1cells,col = spectral,levels = c(0,10,40,70,100,130,160,190,220,250,280),key.axes=axis(4,c(0,10,40,70,100,130,160,190,220,250,280)))
dev.off()
png('Fig 01C Squares Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(squaresmat)),y=seq(-39,-12,length.out = ncol(squaresmat)),squaresmat,col = spectral,levels = c(0,10,40,70,100,130,160,190,220,250,280),key.axes=axis(4,c(0,10,40,70,100,130,160,190,220,250,280)))
dev.off()
#environmental
png('Fig S1 Temp.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spidertemp)),y=seq(-39,-12,length.out = ncol(spidertemp)),spidertemp,col = spectral,levels = c(0,3,6,9,12,15,18,21,24,27,30),key.axes=axis(4,seq(0,30,by=3)))
dev.off()
png('Fig S2 Prec.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spiderprec)),y=seq(-39,-12,length.out = ncol(spiderprec)),log(spiderprec),col = spectra,levels = log(10^(.25*6:15)),key.axes=axis(4,at=log(10^(.25*6:15)),labels=c(32,56,100,178,316,562,1000,1778,3162,5623)))
dev.off()
png('Fig S3 Alti.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spideralt)),y=seq(-39,-12,length.out = ncol(spideralt)),spideralt,col = spectral,levels = c(-1,200,400,600,800,1000,1200,1400,1600,1800,2000),key.axes=axis(4,seq(0,2000,by=200)))
dev.off()
png('Fig S4 Dry.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spiderdry)),y=seq(-39,-12,length.out = ncol(spiderdry)),spiderdry,col = spectral,levels = c(9,21,33,45,57,69,81,93,105,117,129),key.axes=axis(4,seq(9,129,by=12)))
dev.off()
png('Fig S5 DryQuarter.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(dryquartermat)),y=seq(-39,-12,length.out = ncol(dryquartermat)),dryquartermat,col = spectra,levels = c(0,28,56,84,112,140,168,196,224,252,280),key.axes=axis(4,seq(0,280,by=28)))
dev.off()
png('Fig S6 TempSeason.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(tempseasonality)),y=seq(-39,-12,length.out = ncol(tempseasonality)),tempseasonality,col = spectral,levels = c(12,20,25,30,35,40,45,50,55,60,64),key.axes=axis(4,c(12,20,25,30,35,40,45,50,55,60,64)))
dev.off()

#richness curves and latitudinal lines
png('Fig 02 Richness curves.png', width=1000, height=1200)
par( mgp=c(3,0.75,0), cex=1.3, cex.lab=1.5, cex.axis=1.2,las=1, mfcol=c(3,2), oma=c(0,2,0,0),mar=c(4.5,6,0.5,1))
plot(Cellcounts,Cellrichness,log='xy',pch=19,xlab='',ylab="Species richness",xaxt='n',yaxt='n',mar=c(0,6,0,2))
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200))
text(2,150,labels='a',cex=2)
plot(Cellcounts,chaocell,log='xy',pch=19,xlab="",ylab="Chao1 species richness",xaxt='n',yaxt='n',col="blue",mar=c(0,6,.5,2))
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200))
text(2,200,labels='b',cex=2)
plot(Cellcounts,squarecell,log='xy',pch=19,xlab="Data records per cell",ylab="Squares species richness",xaxt='n',yaxt='n',col="red",mar=c(5.5,6,.5,2))
axis(1,at=c(1,5,20,100,250,500,1000,2500),labels=c(1,5,20,100,250,500,1000,2500))
axis(2,at=c(1,5,20,100,200,500),labels=c(1,5,20,100,200,500))
text(2,400,labels='c',cex=2)
plot(Latitude,Cellrichness[Cellcounts>99],log='y',pch=19,xlab="",xaxt='n',ylab="Species richness",yaxt='n',xlim=c(-12,-39),mar=c(0,6,0,2))
axis(2,at=c(25,40,50,100,200),labels=c(25,40,50,100,200))
text(-14,200,labels='d',cex=2)
abline(coef=c(99.5,1.409),lwd=2,untf=T)
plot(Latitude,chaocell[Cellcounts>99],log='y',pch=19,xlab="",xaxt='n',ylab="Chao1 species richness",yaxt='n',col="blue",xlim=c(-12,-39),mar=c(0,6,0,2))
axis(2,at=c(25,40,50,100,200,300),labels=c(25,40,50,100,200,300))
text(-14,275,labels='e',cex=2)
abline(coef=c(141.517,1.257),untf=T,lwd=2,col="blue")
plot(Latitude,squarecell[Cellcounts>99],log='y',pch=19,xlab="Latitude",ylab="Squares species richness",yaxt='n',col="red",xlim=c(-12,-39),mar=c(5.5,6,.5,2))
axis(2,at=c(25,50,100,250,500,1000),labels=c(25,50,100,250,500,1000))
text(-14,500,labels='f',cex=2)
abline(coef=c(132.373,1.022),untf=T,lwd=2,col="red")
dev.off()
png('Fig 03 Richness vs Environment.png', width=800, height=800)
par(mar=c(5,4.2,0,2), mgp=c(3,0.75,0), cex=1.3, cex.lab=1.5,cex.axis=1.4, las=1,mfrow=c(2,1),oma=c(0,1,0.5,2))
plot(sqrt(squarestep[,3]),squarestep[,1],pch=19,log ='xy',xaxt='n',yaxt='n',xlab='Precipitation', ylab='Richness')
axis(1,at=c(sqrt(200),sqrt(500),sqrt(1000),sqrt(1500),sqrt(2000)), labels=c(200,500,1000,1500,2000))
axis(2,at=c(log(10),log(20),log(50),log(100),log(500)),labels=c(10,20,50,100,500))
text(sqrt(300),log(500),labels="a",cex=2)
abline(coef=c(1.55850,0.07874),untf=T,lwd=2)
plot(squarestep[,5],squarestep[,1],pch=19,log ='xy',xaxt='n',yaxt='n',xlab='Seasonality of Precipitation', ylab='Richness')
axis(2,at=c(log(10),log(20),log(50),log(100),log(500)),labels=c(10,20,50,100,500))
axis(1,at=c(15,30,45,60,75,90),labels=c(15,30,45,60,75,90))
text(15,log(500),labels="b",cex=2)
abline(coef=c(3.32202,0.01816),untf=T,lwd=2)
dev.off()