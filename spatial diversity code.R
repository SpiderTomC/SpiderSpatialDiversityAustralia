EastSpiders = read.csv( #filepath of Spider Records East coast# , stringAsFactors=F)
#latitudinal and longitudinal bounds
min(EastSpiders$Latitude, na.rm=T)
max(EastSpiders$Latitude, na.rm=T)
min(EastSpiders$Longitude,na.rm=T)
max(EastSpiders$Longitude,na.rm=T)

#chao1 function#
Chao1<-function(n)        {
  
  return(length(n) + length(which(n == 1)) * (length(which(n == 1)) - 1) / (2 *  (length(which(n == 2)) + 1)))
  
}


simpson<-function(n)	{
  N <- sum(n)
  if(N== length(n))
    return(NA)
  f <- n / N
  S <- N / (N - 1) * (1 - sum(f^2))
  return(1 / (1 - S))
}
#Fisher's alpha function

logSeriesParams<-function(ab)   {
  a = 10
  olda = 0
  n = sum(ab)
  z = 0
  while (abs(a - olda) > 0.0000001 && z < 1000)  {
    z = z + 1
    olda = a
    a = length(ab) / log(1 + n / a)
  }
  x = n/(n + a)
  return(c(a,x))
}
#creation of cell latitude and longitude sequences, using the maximum and  minimum determined previously#
cellLats= seq(-12,-39,by=-0.25)
cellLongs= seq(138, 159.5,by= 0.25)

#creation of cell band arrays. Quarter degree cells balance accuracy with sample size. Chao1, 1/D and Fisher's alpha#
#are run within the loop, giving observed richness and three estimates. also included within loop is the code to make an abundance matrix for all cells#

Longrichness <- array()
Longcounts <- array()
Cellrichness<-array()
Cellcounts<- array()
eastabund<- matrix(nrow=length(unique(EastSpiders$Species)),ncol = 9483,dimnames = list(c(unique(EastSpiders$Species))))
eastabund[is.na=T]<-0
eastspecieslist=c(unique(EastSpiders$Species))
simpcell<-array()
squareslat<-array()
squareslong<-array()
chaocell<-array()
fishercell<-array()
latcount=1
longcount=1
for (i in cellLongs)
{ eastspidlong<- EastSpiders[which(EastSpiders$Longitude<=i & EastSpiders$Longitude >i-.25),]
Longrichness[longcount]<- length(unique(eastspidlong$Species))
Longcounts[longcount]<- nrow(eastspidlong)
for (m in cellLats)
{eastcell<- eastspidlong[which(eastspidlong$Latitude<=m & eastspidlong$Latitude>m-.25),]
eastcell<- eastcell[which(eastcell$Species!=""),]
eastcelltable<- sort(table(eastcell$Species),decreasing = T)
if (nrow(eastcell)>0){
for (j in 1:length(eastspecieslist)){
  eastabund[j,latcount]<-nrow(eastcell[which(eastcell$Species==eastspecieslist[j]),])
}
}
Cellrichness[latcount]<-length(unique(eastcell$Species))
Cellcounts[latcount]<- nrow(eastcell)
simpcell[latcount] <- simpson(eastcelltable)
if (length(eastcelltable>0)){
  falpha<-logSeriesParams(eastcelltable)
  fishercell[latcount]<-falpha[1]
}
else {
fishercell[latcount]<-NA
}
chao1cell[latcount] <- Chao1(eastcelltable)
squareslat[latcount]=m
squareslong[latcount]=i

latcount<-latcount+1
}
longcount<-longcount+1
}
eastcelltable<-eastcelltable[2:length(eastcelltable)]
#creation of richness matrices, flipping column and then transposing is necessary for correct orientation#
#of the contour map#
#raw richness#
spidercell1=matrix(Cellrichness,nrow=length(cellLats),ncol=length(cellLongs))
spidercell1[spidercell1==0]<-NA
spidercell1=t(spidercell1)
spidercell1<-spidercell1[,ncol(spidercell1)[1]:1]
#chao1 richness#
chao1cell[simpcell==NA]<-NA
chao1cell[chao1cell==1]<-NA
chao1cell[chao1cell==0]<-NA
chao1cells= matrix(chao1cell,nrow=length(cellLats),ncol=length(cellLongs))
chao1cells=t(chao1cells)
chao1cells<-chao1cells[,ncol(chao1cells)[1]:1]
#simpsons 1/d#
simpmat= matrix(simpcell,nrow=length(cellLats),ncol=length(cellLongs))
simpmat=simpmat[nrow(simpmat)[1]:1,]
simpmat=t(simpmat)
#fisher's alpha#
fishercell[fishercell>500]<- NA
fishermat= matrix(fishercell, nrow=length(cellLats),ncol=length(cellLongs))
fishermat=fishermat[nrow(fishermat)[1]:1,]
fishermat=t(fishermat)
plot(Cellrichness[which(is.na(isrcell)!=T)],chao1cell[which(is.na(isrcell)!=T)],log="xy",xlab="Raw Richness", ylab="Estimated Richness", pch=16, cex=0.8)
points(Cellrichness[which(is.na(isrcell)!=T)],simpcell[which(is.na(isrcell)!=T)], cex=0.8, pch=16, col="red")
points(Cellrichness[which(is.na(isrcell)!=T)],isrcell[which(is.na(isrcell)!=T)], cex=0.8, pch=16, col="blue")
lines(c(1:350),c(1:350),lwd=2)
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
#latitude vs diversity
rawlatline= (lm(log(Cellrichness[Cellcounts>99])~squareslat[Cellcounts>99]))$coef
chaolatline= (lm(chao1cell[Cellcounts>99]~squareslat[Cellcounts>99]))$coef
simplatline= (lm(simpcell[Cellcounts>99]~squareslat[Cellcounts>99]))$coef
fisherlatline= (lm(fishercell[Cellcounts>99]~squareslat[Cellcounts>99]))$coef
#population density code
library(sf)
poplat=squareslat-0.125
poplong=squareslong+0.125
popgrid=cbind(poplong,poplat)
popgrid<- popgrid[order(popgrid[,2]),]
popdensity = t(read.table('c:\\Users\\Thoma\\Documents\\PhD Chapter one\\gpw_v4_population_density_2020_15_min.asc',sep='' ,na.strings='-9999'))
areas = read_sf("c:\\Users\\Thoma\\Documents\\PhD Chapter one\\ne_50m_urban_areas.shp")
coarse = st_make_grid(areas,offset=c(-180,-90),cellsize=c(0.25,0.25),n=c(360*4,180*4),what = "centers")
popcoords = st_coordinates(coarse)
popcoords2= popcoords[popcoords[,1]>100 & popcoords[,1]<162 & popcoords[,2]>=-50 & popcoords[,2]<=-10,]
popcoords3= popcoords[popcoords[,1]>=138 & popcoords[,1]<=159.625 & popcoords[,2]>=-39.125 & popcoords[,2]<=-12]
cellDensity = array(dim=nrow(popcoords),data=NA)
for (i in 1:nrow(popcoords))       {
  rounded_long = round((popcoords[i,1] + 180 + 0.25)*4)
  rounded_lat = ncol(popdensity) - round((popcoords[i,2] + 90 + 0.25)*4)
  
  if (rounded_lat > 0 && ncol(popdensity) - rounded_lat > 0)
    cellDensity[i] = popdensity[rounded_long,rounded_lat]

}
ausdensity<- cellDensity[popcoords[,1]>=138.125 & popcoords[,1]<=159.625 & popcoords[,2]>=-39.125 & popcoords[,2]<=-12.125]
ausgrid=cbind(popgrid,ausdensity)
ausgrid=ausgrid[order(ausgrid[,1],-ausgrid[,2]),]
#Exploratory factor analysis and linear regression#
envfact = cbind(log(Cellrichness), log(chao1cell),log(simpcell),log(fishercell),climvalues,altvalues,dryvalues,dryquartervalues,tempseasonvalues,squareslat, log(ausgrid[,3]), Cellcounts)
envfact= envfact[which(envfact[,13]>49),]
envfact[,6]<- sqrt(envfact[,6])
rownames(envfact)<-c(1:222)
envfact[envfact[,12]==-Inf]<-NA
envfact= na.omit(envfact)
restrictedabund=eastabund[,which(Cellcounts>49)]
restrictedabund<-restrictedabund[2:2323,]
restrictedabund <- restrictedabund[ order(row.names(restrictedabund)), ]
restrictedabund<-restrictedabund[,as.numeric(rownames(envfact))]

envscores= factanal(scale(envfact[,5:12]),factors=3, scores = 'regression')$scores
summary(lm(envfact[,1]~envscores[,1:3]))
summary(lm(envfact[,2]~envscores[,1:3]))
summary(lm(envfact[,3]~envscores[,1:3]))
summary(lm(envfact[,4]~envscores[,1:3]))
summary(lm(envfact[,2]~envscores[,1:3]+envfact[,13]))
summary(lm(envfact[,3]~envscores[,1:3]+envfact[,13]))
summary(lm(envfact[,4]~envscores[,1:3]+envfact[,13]))
summary(lm(envfact[,1]~envfact[,5]))
#linear models of environment vs richness one by one#
#observed richness
rawcoeff<- matrix(nrow=8, ncol = 5)
for (i in 5:12){
  model<-summary(lm(envfact[,1]~envfact[,i]))
  rawcoeff[i-4,1]<-model$adj.r.squared
  rawcoeff[i-4,2]<-model$sigma
  rawcoeff[i-4,3]<-model$coefficients[2,1]
  rawcoeff[i-4,4]<-model$coefficients[2,3]
  rawcoeff[i-4,5]<-model$coefficients[2,4]
}
#chao 1
chaocoeff<- matrix(nrow=8, ncol = 5)
for (i in 5:12){
  model<-summary(lm(envfact[,2]~envfact[,i]))
  chaocoeff[i-4,1]<-model$adj.r.squared
  chaocoeff[i-4,2]<-model$sigma
  chaocoeff[i-4,3]<-model$coefficients[2,1]
  chaocoeff[i-4,4]<-model$coefficients[2,3]
  chaocoeff[i-4,5]<-model$coefficients[2,4]
}
#inverse simpson's index#
simpcoeff<- matrix(nrow=8, ncol = 5)
for (i in 5:12){
  model<-summary(lm(envfact[,3]~envfact[,i]))
  simpcoeff[i-4,1]<-model$adj.r.squared
  simpcoeff[i-4,2]<-model$sigma
  simpcoeff[i-4,3]<-model$coefficients[2,1]
  simpcoeff[i-4,4]<-model$coefficients[2,3]
  simpcoeff[i-4,5]<-model$coefficients[2,4]
}
#fisher's alpha#
fishcoeff<- matrix(nrow=8, ncol = 5)
for (i in 5:12){
  model<-summary(lm(envfact[,4]~envfact[,i]))
  fishcoeff[i-4,1]<-model$adj.r.squared
  fishcoeff[i-4,2]<-model$sigma
  fishcoeff[i-4,3]<-model$coefficients[2,1]
  fishcoeff[i-4,4]<-model$coefficients[2,3]
  fishcoeff[i-4,5]<-model$coefficients[2,4]
}

#creation of contour plots#
library(car)
png('Fig 01 Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(spidercell1)),y=seq(-39,-12,length.out = ncol(spidercell1)),spidercell1,col = spectral,levels = c(0,10,40,70,100,130,160,190,220,250,290),key.axes=axis(4,c(0,10,40,70,100,130,160,190,220,250,290)))
dev.off()
png('Fig 01B Chao1 Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(chao1cells)),y=seq(-39,-12,length.out = ncol(chao1cells)),chao1cells,col = spectral,levels = c(0,10,30,70,120,170,220,270,320,370,420),key.axes=axis(4,c(0,10,30,70,120,170,220,270,320,370,420)))
dev.off()
png('Fig 01C Simpsons D.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(simpmat)),y=seq(-39,-12,length.out = ncol(simpmat)),simpmat,col = spectral,levels = c(0,5,10,20,40,60,80,100,120,140,160),key.axes=axis(4,c(0,5,10,20,40,60,80,100,120,140,160)))
dev.off()
png('Fig 01D ISR Rarefy.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(isrmat)),y=seq(-39,-12,length.out = ncol(isrmat)),isrmat,col = spectral,levels = c(0,10,20,40,80,120,160,200,240,280,320),key.axes=axis(4,c(0,10,20,40,80,120,160,200,240,280,320)))
dev.off()
png('Fig 01E Fisher.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,160,length.out = nrow(fishermat)),y=seq(-39,-12,length.out = ncol(fishermat)),fishermat,col = spectral,levels = c(0,10,20,30,45,60,75,90,105,120,135),key.axes=axis(4,c(0,10,20,30,45,60,75,90,105,120,135)))
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
png('Fig 02 Richness curves.png', width=1000, height=1000)
par( mgp=c(3,0.75,0), cex=1.3, cex.lab=2.5, cex.axis=1.8,las=1, mfcol=c(2,2), oma=c(0,2,0,0),mar=c(4.5,6,0.5,1))
plot(Cellcounts,Cellrichness,log='xy',pch=19,xlab='',ylab="Observed species richness",xaxt='n',yaxt='n',mar=c(0,6,0,2),cex.lab=2)
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200),cex= 2.2)
text(2,150,labels='a',cex=4)
plot(Cellcounts,chao1cell,log='xy',pch=19,xlab="Data records per cell",ylab="Chao1 species richness",xaxt='n',yaxt='n',mar=c(0,6,.5,2),cex.lab=2)
axis(1,at=c(1,5,20,100,250,500,1000,2500),labels=c(1,5,20,100,250,500,1000,2500),cex=2.5)
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200),cex=2.2)
text(2,200,labels='b',cex=4)
plot(Cellcounts,simpcell,log='xy',pch=19,xlab="",ylab="Simpson's 1/D",xaxt='n',yaxt='n',mar=c(5.5,6,.5,2),cex.lab=2)
axis(2,at=c(0,1,5,10,50,100),labels=c(0,1,5,10,50,100),cex=2.2)
text(2,90,labels='c',cex=4)
plot(Cellcounts,fishercell,log='xy',pch=19,xlab="Data records per cell",ylab="Fisher's alpha",xaxt='n',yaxt='n',col="black",mar=c(5.5,6,.5,2),cex.lab=2)
axis(2,at=c(0,1,5,10,50,100),labels=c(0,1,5,10,50,100),cex=2.2)
axis(1,at=c(1,5,20,100,250,500,1000,2500),labels=c(1,5,20,100,250,500,1000,2500),cex=2.5)
text(2,65,labels='d',cex=4)
dev.off()

png('Fig 02 Lat lines.png', width=1000, height=1000)
par( mgp=c(3,0.75,0), cex=1.3, cex.lab=2.5, cex.axis=1.8,las=1, mfcol=c(2,2), oma=c(0,2,0,0),mar=c(4.5,6.6,0.5,1))
plot(squareslat[Cellcounts>99],Cellrichness[Cellcounts>99],pch=19,log='y',xlab="",xaxt='n',ylab="Richness",yaxt='n',xlim=c(-12,-39),mar=c(0,6,0,2),cex.lab=2)
axis(2,at=c(50,100,150,200,250),labels=c(50,100,150,200,250),cex=2.2)
text(-14,210,labels='a',cex=4)
lines(c(0,-40),c(133.6779,52.5678),lwd=2)
plot(squareslat[Cellcounts>99],chao1cell[Cellcounts>99],pch=19,log='y', xlab="Latitude",ylab="Richness",yaxt='n',col=,xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2)
axis(2,at=c(50,100,150,200,250,300),labels=c(25,40,50,100,200,300),cex=2.2)
text(-14,275,labels='b',cex=4)
lines(c(0,-40),c(242.349794,87.717354),lwd=2)
plot(squareslat[Cellcounts>99],simpcell[Cellcounts>99],pch=19,log='y',xaxt='n',xlab="",ylab="",yaxt='n',xlim=c(-12,-39),mar=c(0,6,.5,2),cex.lab=2)
axis(2,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex=2.2)
text(-14,100,labels='c',cex=4)
lines(c(0,-40),c(111.197988,6.089548),lwd=2)
plot(squareslat[Cellcounts>99],fishercell[Cellcounts>99],pch=19,log='y',xlab="Latitude",ylab="",yaxt='n',xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2)
axis(2,at=c(1,5,10,50,100),labels=c(1,5,10,50,100),cex=2.2)
text(-14,80,labels='d',cex=4)
lines(c(0,-40),c(113.099647,14.74673),lwd=2)
dev.off()

png('Fig 03 Richness vs Environment.png', width=600, height=900)
par(mar=c(5,4.2,0,2), mgp=c(3,0.75,0), cex=1.3, cex.lab=1.5,cex.axis=1.4, las=1,mfrow=c(3,1),oma=c(0,1,0.5,2))
plot(envfact[,6],envfact[,3],pch=19,log ='y',yaxt='n',xaxt= 'n',xlab='Precipitation', ylab="Simpson's 1/D")
axis(1,at=c(sqrt(200),sqrt(500),sqrt(1000),sqrt(1500),sqrt(2000)), labels=c(200,500,1000,1500,2000))
axis(2,at=c(log(1),log(5),log(10),log(50),log(100)),labels=c(1,5,10,50,100))
text(sqrt(230),log(90),labels="a",cex=4)
lines(c(1,sqrt(3500)),c(1.72227,4.17862632594),lwd=2)
plot(envfact[,8],envfact[,3],pch=19,log ='y',xaxt='n',yaxt='n',xlab='Seasonality of Precipitation', ylab="Simpson's 1/D")
axis(2,at=c(log(5),log(10),log(20),log(50),log(100)),labels=c(5,10,20,50,100))
axis(1,at=c(15,30,45,60,75,90),labels=c(15,30,45,60,75,90))
text(15,log(90),labels="b",cex=4)
lines(c(0,150),c(2.3146,5.1946),lwd=2)
plot(envfact[,5],envfact[,3],pch=19,log ='y',xaxt='n',yaxt='n',xlab='Temperature', ylab="Simpson's 1/D")
axis(2,at=c(log(5),log(10),log(20),log(50),log(100)),labels=c(5,10,20,50,100))
axis(1,at=c(40,80,120,160,200,240),labels=c(4,8,12,16,20,24))
text(65,log(90),labels="c",cex=4)
lines(c(0,300),c(1.343502,4.304802),lwd=2)
dev.off()

png('Fig 04 Richness vs Environment.png', width=600, height=400)
par(mar=c(5,4.2,0,2), mgp=c(3,0.75,0), cex=1.3, cex.lab=1.5,cex.axis=1.4, las=1,oma=c(0,1,0.5,2))
plot(envfact[,6],envfact[,3],pch=19,log ='xy',xaxt='n',yaxt='n',xlab='Precipitation', ylab="Simpson's 1/D")
axis(1,at=c(sqrt(200),sqrt(500),sqrt(1000),sqrt(1500),sqrt(2000)), labels=c(200,500,1000,1500,2000))
axis(2,at=c(log(2),log(5),log(10),log(50),log(100)),labels=c(2,5,10,50,100))
lines(c(1,sqrt(3500)),c(1.681399,4.255653),lwd=2)
text(sqrt(300),log(80),labels="a",cex=2)
dev.off()
