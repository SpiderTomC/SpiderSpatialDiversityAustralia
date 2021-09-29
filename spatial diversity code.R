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
  #broken stick function
broken<-function(n)	{
  library(stats4)
  if (length(table(n)) < 3)
    return(list('richness' = NA, 'SAD.AICc' = NA, 'octave.AICc' = NA))
  R <- length(n)
  N <- sum(n[n <= 2^16])
  s <- array(dim=N,data=0)
  t <- table(n[n <= N])
  s[as.numeric(names(t))] <- t
  like<-function(S)	{
    p <- (1 - 1:N / N)^(S - 2) * (S - 1) / N
    p <- p[1:N] / sum(p[1:N])
    ll <- -sum(dbinom(s,R,p,log=T))
    if (is.infinite(ll) || is.nan(ll))
      ll <- 9999
    ll
  }
  S <- coef(mle(like,lower=list(S=length(n)),upper=list(S=1000),start=list(S=2 * length(n))))
  if (S == 1000)	{
    warning('richness estimate is out of range')
    return(list('richness' = NA, 'SAD.AICc' = NA, 'octave.AICc' = NA))
  }
  p <- (1 - 1:N / N)^(S - 2) * (S - 1) / N
  p[(N + 1):2^16] <- 0
  p <- p / sum(p)
  s[(N + 1):2^16] <- 0
  aicc <- -2 * sum(dbinom(s,R,p,log=T)) + 2 + 4 / (sum(s > 0) - 2)
  t <- ceiling(log(1:2^16) / log(2)) + 1
  s2 <- array()
  p2 <- array()
  for (i in 1:17)	{
    s2[i] <- sum(s[t == i])
    p2[i] <- sum(p[t == i])
  }
  o_aicc <- -2 * sum(dbinom(s2,sum(s2),p2,log=T)) + 2 + 4 / (sum(s2 > 0) - 2)
  return(list('richness' = as.numeric(S), 'SAD.AICc' = as.numeric(aicc), 'octave.AICc' = as.numeric(o_aicc), 'fitted.SAD' = p[1:max(max(n),1024)]))
}


#creation of cell latitude and longitude sequences, using the maximum and  minimum determined previously#
cellLats= seq(-12,-39,by=-0.25)
cellLongs= seq(138, 159.5,by= 0.25)

#creation of cell band arrays. Quarter degree cells balance accuracy with sample size. Chao1, 1/D and Fisher's alpha#
#are run within the loop, giving observed richness and three estimates. also included within loop is the code to make an abundance matrix for all cells#
library(sads)
library(vegan)
library(poilog)
Longrichness <- array()
Longcounts <- array()
Cellrichness<-array()
Cellcounts<- array()
eastabund<- matrix(nrow=length(unique(EastSpiders$Species)),ncol = 7303,dimnames = list(c(unique(EastSpiders$Species))))
eastabund[is.na=T]<-0
eastspecieslist=c(unique(EastSpiders$Species))
simpcell<-array()
squareslat<-array()
squareslong<-array()
chao1cell<-array()
acecell<- array()
fishercell<-array()
logcell<-array()
zsmcell<-array()
brokencell<-array()
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
  logcell[latcount]<- length(eastcelltable) / poilogMLE(eastcelltable)$p
  zsmcell[latcount] <- coef(fitmzsm(eastcelltable))[2]
}
else {
fishercell[latcount]<-NA
logcell[latcount]<-NA
zsmcell[latcount]<-NA
}
chao1cell[latcount] <- Chao1(eastcelltable)
acecell[latcount]<-estimateR(eastcelltable)[[4]]
brokencell[latcount] <- broken(eastcelltable)$richness
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
#log normal#
logmat= matrix(logcell,nrow=length(cellLats),ncol=length(cellLongs))
logmat=logmat[nrow(logmat)[1]:1,]
logmat=t(logmat)
#fisher's alpha#
fishercell[fishercell>500]<- NA
fishermat= matrix(fishercell, nrow=length(cellLats),ncol=length(cellLongs))
fishermat=fishermat[nrow(fishermat)[1]:1,]
fishermat=t(fishermat)
plot(Cellrichness[which(is.na(isrcell)!=T)],chao1cell[which(is.na(isrcell)!=T)],log="xy",xlab="Raw Richness", ylab="Estimated Richness", pch=16, cex=0.8)
points(Cellrichness[which(is.na(isrcell)!=T)],simpcell[which(is.na(isrcell)!=T)], cex=0.8, pch=16, col="red")
points(Cellrichness[which(is.na(isrcell)!=T)],isrcell[which(is.na(isrcell)!=T)], cex=0.8, pch=16, col="blue")
lines(c(1:350),c(1:350),lwd=2)
#ACE#
acemat= matrix(acecell,nrow=length(cellLats),ncol=length(cellLongs))
acemat=acemat[nrow(acemat)[1]:1,]
acemat=t(acemat)
#broken stick
brokenmat= matrix(brokencell,nrow=length(cellLats),ncol=length(cellLongs))
brokenmat=brokenmat[nrow(brokenmat)[1]:1,]
brokenmat=t(brokenmat)
#zero-sum multinomial
zsmmat= matrix(zsmcell,nrow=length(cellLats),ncol=length(cellLongs))
zsmmat=zsmmat[nrow(zsmmat)[1]:1,]
zsmmat=t(zsmmat)
#Worldclim data extraction and assembly#
#requires sp and raster package#
library(raster)
library(sp)
totalclim <-getData("worldclim",var="bio",res=10)
climcoord <- SpatialPoints(cbind(squareslong,squareslat))
climcoord2 <- SpatialPoints(cbind(squareslong-0.125,squareslat-0.125))
climcoord3 <- SpatialPoints(cbind(squareslong-0.25,squareslat-0.25))
climcoord4 <- SpatialPoints(cbind(squareslong,squareslat-0.25))
climcoord5 <- SpatialPoints(cbind(squareslong-0.25,squareslat))
#the five sets of coordinates correspond to the middle and corners of each cell#
eastclim <- totalclim[[c(1,12)]]
climvalues1 <- extract(eastclim,climcoord)
climvalues2 <- extract(eastclim,climcoord2)
climvalues3 <- extract(eastclim,climcoord3)
climvalues4 <- extract(eastclim,climcoord4)
climvalues5 <- extract(eastclim,climcoord5)
climvalues=matrix(ncol=2, nrow = nrow(climvalues1))
for (i in 1:nrow(climvalues1)) {
  climvalues[i,1]<-mean(c(climvalues1[i,1],climvalues2[i,1],climvalues3[i,1],climvalues4[i,1],climvalues5[i,1]),na.rm=T)
  climvalues[i,2]<-mean(c(climvalues1[i,2],climvalues2[i,2],climvalues3[i,2],climvalues4[i,2],climvalues5[i,2]),na.rm=T)
}
#this loop generates an average from the five sets of values, reducing missing cells within the matrix#
eastdry <- totalclim[[c(15)]]
dryvalues1<- extract(eastdry,climcoord)
dryvalues2<- extract(eastdry,climcoord2)
dryvalues3<- extract(eastdry,climcoord3)
dryvalues4<- extract(eastdry,climcoord4)
dryvalues5<- extract(eastdry,climcoord5)
dryvalues=array()
for (i in 1:length(dryvalues1)){
  dryvalues[i]<- mean(c(dryvalues1[i],dryvalues2[i],dryvalues3[i],dryvalues4[i],dryvalues5[i]),na.rm=T)
}
eastalt <-getData('alt',country= 'AUS')
altvalues1<- extract(eastalt,climcoord)
altvalues2<- extract(eastalt,climcoord2)
altvalues3<- extract(eastalt,climcoord3)
altvalues4<- extract(eastalt,climcoord4)
altvalues5<- extract(eastalt,climcoord5)
altvalues=array()
for (i in 1:length(altvalues1)){
  altvalues[i]<- mean(c(altvalues1[i],altvalues2[i],altvalues3[i],altvalues4[i],altvalues5[i]),na.rm=T)
}
tempseason= totalclim [[c(4)]]
tempseasonvalues1= extract(tempseason,climcoord)
tempseasonvalues2= extract(tempseason,climcoord2)
tempseasonvalues3= extract(tempseason,climcoord3)
tempseasonvalues4= extract(tempseason,climcoord4)
tempseasonvalues5= extract(tempseason,climcoord5)
tempseasonvalues=array()
for (i in 1:length(tempseasonvalues1)){
  tempseasonvalues[i]<-mean(c(tempseasonvalues1[i],tempseasonvalues2[i],tempseasonvalues3[i],tempseasonvalues4[i],tempseasonvalues5[i]),na.rm=T)
}
dryquarter= totalclim[[c(17)]]
dryquartervalues1<- extract(dryquarter,climcoord)
dryquartervalues2<- extract(dryquarter,climcoord2)
dryquartervalues3<- extract(dryquarter,climcoord3)
dryquartervalues4<- extract(dryquarter,climcoord4)
dryquartervalues5<- extract(dryquarter,climcoord5)
dryquartervalues=array()
for (i in 1:length(dryquartervalues1)){
  dryquartervalues[i]<- mean(c(dryquartervalues1[i],dryquartervalues2[i],dryquartervalues3[i],dryquartervalues4[i],dryquartervalues5[i]),na.rm=T)
}

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
envfact = cbind(log(Cellrichness),log(acecell),log(brokencell), log(chao1cell),log(fishercell),log(logcell),log(simpcell),log(zsmcell),climvalues,altvalues,dryvalues,dryquartervalues,tempseasonvalues,squareslat, squareslong, log(ausgrid[,3]), Cellcounts)
colnames(envfact)<-c("Observed","ACE","Broken stick","Chao1","Fisher's Alpha","Log normal","Simpson's 1/D","ZSM", "Temp","Prec","Elevation","Prec Seasonality","Prec Dry Quarter","Temp Seasonality","Lat","Long","Pop Density","Counts")
envfact= envfact[which(envfact[,18]>99),]
envfact[,10]<- sqrt(envfact[,10])
envfact[envfact[,11]==NA]<-1
envfact= na.omit(envfact)
rownames(envfact)<-c(1:85)

envscores= factanal(scale(envfact[,c(9:15,17)]),factors=3, scores = 'regression')$scores
summary(lm(envfact[,1]~envscores[,1:3]))
summary(lm(envfact[,2]~envscores[,1:3]))
summary(lm(envfact[,3]~envscores[,1:3]))
summary(lm(envfact[,4]~envscores[,1:3]))
summary(lm(envfact[,5]~envscores[,1:3]))
summary(lm(envfact[,6]~envscores[,1:3]))
summary(lm(envfact[,7]~envscores[,1:3]))
summary(lm(envfact[,8]~envscores[,1:3]))
summary(lm(envfact[,2]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,3]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,4]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,5]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,6]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,7]~envscores[,1:3]+log(envfact[,18])))
summary(lm(envfact[,8]~envscores[,1:3]+log(envfact[,18])))
#spatial autoregression
library(spdep)
library(spatialreg)
eastcoord= cbind(envfact[,16],envfact[,15])

eastreg = which(! duplicated(eastcoord))

eastnb = knn2nb(knearneigh(eastcoord[eastreg,],k=2,longlat=T))
eastnb2 = knn2nb(knearneigh(regcoord,k=2,longlat=T))
summary(spautolm(envfact[eastreg,1] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,2] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,3] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,4] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,5] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,6] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,7] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,8] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb)),Nagelkerke=T)
summary(spautolm(envfact[eastreg,1] ~ envfact[eastreg,7],listw=nb2listw(eastnb)),Nagelkerke=T)
#Moran's I for residuals
spautoresid1=(spautolm(envfact[eastreg,1] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid2=(spautolm(envfact[eastreg,2] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid3=(spautolm(envfact[eastreg,3] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid4=(spautolm(envfact[eastreg,4] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid5=(spautolm(envfact[eastreg,5] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid6=(spautolm(envfact[eastreg,6] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid7=(spautolm(envfact[eastreg,7] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
spautoresid8=(spautolm(envfact[eastreg,8] ~ envscores[eastreg,1:3],listw=nb2listw(eastnb))$fit$residuals)
eastdist=dist(eastcoord)
Moran.I(spautoresid1,weight=as.matrix(eastdist))
moran.test(spautoresid1,listw = nb2listw(eastnb))
moran.test(spautoresid2,listw = nb2listw(eastnb))
moran.test(spautoresid3,listw = nb2listw(eastnb))
moran.test(spautoresid4,listw = nb2listw(eastnb))
moran.test(spautoresid5,listw = nb2listw(eastnb))
moran.test(spautoresid6,listw = nb2listw(eastnb))
moran.test(spautoresid7,listw = nb2listw(eastnb))
moran.test(spautoresid8,listw = nb2listw(eastnb))
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
#ACE#
acecoeff<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(lm(envfact[,2]~envfact[,i]))
  acecoeff[rowcount,1]<-model$adj.r.squared
  acecoeff[rowcount,2]<-model$sigma
  acecoeff[rowcount,3]<-model$coefficients[2,1]
  acecoeff[rowcount,4]<-model$coefficients[2,3]
  acecoeff[rowcount,5]<-model$coefficients[2,4]
  rowcount<-rowcount+1
}
#broken stick
brokencoeff<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(lm(envfact[,3]~envfact[,i]))
  brokencoeff[rowcount,1]<-model$adj.r.squared
  brokencoeff[rowcount,2]<-model$sigma
  brokencoeff[rowcount,3]<-model$coefficients[2,1]
  brokencoeff[rowcount,4]<-model$coefficients[2,3]
  brokencoeff[rowcount,5]<-model$coefficients[2,4]
  rowcount<-rowcount+1
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
logcoeff<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(lm(envfact[,6]~envfact[,i]))
  logcoeff[rowcount,1]<-model$adj.r.squared
  logcoeff[rowcount,2]<-model$sigma
  logcoeff[rowcount,3]<-model$coefficients[2,1]
  logcoeff[rowcount,4]<-model$coefficients[2,3]
  logcoeff[rowcount,5]<-model$coefficients[2,4]
  rowcount<-rowcount+1
}
  #ZSM#
zsmcoeff<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(lm(envfact[,8]~envfact[,i]))
  zsmcoeff[rowcount,1]<-model$adj.r.squared
  zsmcoeff[rowcount,2]<-model$sigma
  zsmcoeff[rowcount,3]<-model$coefficients[2,1]
  zsmcoeff[rowcount,4]<-model$coefficients[2,3]
  zsmcoeff[rowcount,5]<-model$coefficients[2,4]
  rowcount<-rowcount+1
}
#linear models one by one using spatial autoregression
rawcoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,1] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  rawcoeff2[rowcount,1]<-model$NK
  rawcoeff2[rowcount,2]<-model$fit$s2
  rawcoeff2[rowcount,3]<-model$Coef[2,1]
  rawcoeff2[rowcount,4]<-model$Coef[2,3]
  rawcoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
acecoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,2] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  acecoeff2[rowcount,1]<-model$NK
  acecoeff2[rowcount,2]<-model$fit$s2
  acecoeff2[rowcount,3]<-model$Coef[2,1]
  acecoeff2[rowcount,4]<-model$Coef[2,3]
  acecoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#broken stick#
brokencoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,3] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  brokencoeff2[rowcount,1]<-model$NK
  brokencoeff2[rowcount,2]<-model$fit$s2
  brokencoeff2[rowcount,3]<-model$Coef[2,1]
  brokencoeff2[rowcount,4]<-model$Coef[2,3]
  brokencoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#Chao 1
chaocoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,4] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  chaocoeff2[rowcount,1]<-model$NK
  chaocoeff2[rowcount,2]<-model$fit$s2
  chaocoeff2[rowcount,3]<-model$Coef[2,1]
  chaocoeff2[rowcount,4]<-model$Coef[2,3]
  chaocoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#Fisher's alpha#
fishercoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,5] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  fishercoeff2[rowcount,1]<-model$NK
  fishercoeff2[rowcount,2]<-model$fit$s2
  fishercoeff2[rowcount,3]<-model$Coef[2,1]
  fishercoeff2[rowcount,4]<-model$Coef[2,3]
  fishercoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#log normal#
logcoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,6] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  logcoeff2[rowcount,1]<-model$NK
  logcoeff2[rowcount,2]<-model$fit$s2
  logcoeff2[rowcount,3]<-model$Coef[2,1]
  logcoeff2[rowcount,4]<-model$Coef[2,3]
  logcoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#Simpson's 1/D
simpcoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,7] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  simpcoeff2[rowcount,1]<-model$NK
  simpcoeff2[rowcount,2]<-model$fit$s2
  simpcoeff2[rowcount,3]<-model$Coef[2,1]
  simpcoeff2[rowcount,4]<-model$Coef[2,3]
  simpcoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
#ZSM#
zsmcoeff2<- matrix(nrow=8, ncol = 5)
rowcount=1
for (i in c(9:15,17)){
  model<-summary(spautolm(envfact[eastreg,8] ~ envfact[eastreg,i],listw=nb2listw(eastnb)),Nagelkerke=T)
  zsmcoeff2[rowcount,1]<-model$NK
  zsmcoeff2[rowcount,2]<-model$fit$s2
  zsmcoeff2[rowcount,3]<-model$Coef[2,1]
  zsmcoeff2[rowcount,4]<-model$Coef[2,3]
  zsmcoeff2[rowcount,5]<-model$Coef[2,4]
  rowcount<-rowcount+1
}
  #factor analysis plot comparing loadings of diversity estimates
allrich = cbind(Cellrichness,acecell,brokencell, chao1cell,fishercell,logcell ,simpcell,zsmcell)

colnames(allrich)=c("observed","ACE","broken stick","Chao 1","Fisher's alpha","Poisson log normal","1/D","ZSM")
richnesslab=c("Obs","Chao 1","1/D","alpha","PLN","ACE","BS","ZSM")
# number of factors should be varied to find a good cutoff

richfact = factanal(log(na.omit(allrich)),factors=3)$loadings
par(mfrow=c(2,1))
png('Fig 01 Diversity factor plot.png', width=700, height=700)


par(mar=c(5,5,1,2) + 0.1)

xoff = c(-0.04,-0.04,-0.04,0.04,0.04,0.04,-0.04,-0.04) * 1.5
yoff = c(0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015) * 1.5

plot(richfact[,1],richfact[,2],mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Factor 1',ylab='Factor 2',cex=0,xlim=c(0.2,0.95),ylim=c(0.2,0.95),cex.lab=1.6,cex.axis=1.4)
for (i in 1:nrow(richfact)){
  lines(c(richfact[i,1],richfact[i,1] + xoff[i]),c(richfact[i,2],richfact[i,2] + yoff[i]))
}
points(richfact[,1] + xoff,richfact[,2] + yoff,cex=2,pch=19,col='white')
points(richfact[5,1] + xoff[5],richfact[5,2] + yoff[5] + 0.005,cex=3,pch=19,col='white')
points(richfact[7,1] + xoff[7],richfact[7,2] + yoff[7],cex=3,pch=19,col='white')
points(richfact,cex=2,pch=21,bg=c('dodgerblue','blue','red','orange','purple','pink','black','green'))
text(richfact[,1] + 1.2 * xoff,richfact[,2] + 1.2 * yoff,labels=colnames(allrich),cex=1.3)
dev.off()

library(car)
png('Fig 01A Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.8, las=1)
filled.contour(x=seq(138,155,length.out = nrow(spidercell1)),y=seq(-39,-12,length.out = ncol(spidercell1)),spidercell1,col = spectral,levels = c(0,10,40,70,100,130,160,190,220,250,290),key.axes=axis(4,c(0,10,40,70,100,130,160,190,220,250,290)),cex=2)
points(151.21, 33.87, cex=3)
dev.off()
png('Fig 01B Chao1 Richness.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(chao1cells)),y=seq(-39,-12,length.out = ncol(chao1cells)),chao1cells,col = spectral,levels = c(0,10,30,70,120,170,220,270,320,370,420),key.axes=axis(4,c(0,10,30,70,120,170,220,270,320,370,420)))
dev.off()
png('Fig 01C ACE.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(acemat)),y=seq(-39,-12,length.out = ncol(acemat)),(acemat),col=spectral,levels = c(0,45,90,135,180,225,270,315,360,405,450),key.axes=axis(4,c(0,45,90,135,180,225,270,315,360,405,450)))
dev.off()
png('Fig 01D Simpsons D.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(simpmat)),y=seq(-39,-12,length.out = ncol(simpmat)),simpmat,col = spectral,levels = c(0,5,10,20,40,60,80,100,120,140,160),key.axes=axis(4,c(0,5,10,20,40,60,80,100,120,140,160)))
dev.off()
png('Fig 01E Fisher.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(fishermat)),y=seq(-39,-12,length.out = ncol(fishermat)),fishermat,col=spectral,levels = c(0,5,10,20,30,40,60,80,100,120,140),key.axes=axis(4,c(0,5,10,20,30,40,60,80,100,120,140),labels=c(0,5,10,20,30,40,60,80,100,120,140)))
dev.off()
png('Fig 01F Log normal.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(logmat)),y=seq(-39,-12,length.out = ncol(logmat)),log(logmat),col = spectral,levels = c(0,log(5),log(10),log(50),log(100),log(500),log(1000),log(2500),log(5000),log(10000),log(25000)),key.axes=axis(4,c(0,log(5),log(10),log(50),log(100),log(500),log(1000),log(2500),log(5000),log(10000),log(25000)),labels=c(0,5,10,50,100,500,1000,2500,5000,10000,25000)))
dev.off()
png('Fig S1 Broken stick.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(brokenmat)),y=seq(-39,-12,length.out = ncol(brokenmat)),brokenmat,col = spectral,levels = c(0,10,30,60,90,120,150,190,230,270,310),key.axes=axis(4,c(0,10,30,60,90,120,150,190,230,270,310)))
dev.off()
png('Fig S2 ZSM.png', width=550, height=700)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(zsmmat)),y=seq(-39,-12,length.out = ncol(zsmmat)),(zsmmat),col=spectral,levels = c(0,10,20,30,40,50,60,70,80,90,110),key.axes=axis(4,c(0,10,20,30,40,50,60,70,80,90,110)))
dev.off()
#environmental
png('Fig S1 Temp.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(spidertemp)),y=seq(-39,-12,length.out = ncol(spidertemp)),spidertemp,col = spectral,levels = c(0,3,6,9,12,15,18,21,24,27,30),key.axes=axis(4,seq(0,30,by=3)))
dev.off()
png('Fig S2 Prec.png', width=550, height=750)
par(mar=c(5.5,6,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
filled.contour(x=seq(138,155,length.out = nrow(spiderprec)),y=seq(-39,-12,length.out = ncol(spiderprec)),log(spiderprec),col = spectra,levels = log(10^(.25*6:15)),key.axes=axis(4,at=log(10^(.25*6:15)),labels=c(32,56,100,178,316,562,1000,1778,3162,5623)))
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
png('Fig 02 Richness curves.png', width=900, height=1200)
par( mgp=c(5,0.75,0), cex=1.3, cex.lab=3, cex.axis=1.8,las=1,mfrow=c(3,2), oma=c(0,4,0,0),mar=c(8,10,0,1))
plot(Cellcounts,Cellrichness,log='xy',pch=19,xlab='',ylab="Observed species richness",xaxt='n',yaxt='n',mar=c(4,7,0,1),cex.lab=2.5)
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200),cex.axis=1.5)
text(2,150,labels='a',cex=4)
plot(Cellcounts,acecell,log='xy',pch=19,xlab="",ylab="ACE species richness",xaxt='n',yaxt='n',mar=c(4,7,.5,1),cex.lab=2.5)
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200),cex.axis=1.5)
text(2,200,labels='b',cex=4)
plot(Cellcounts,chao1cell,log='xy',pch=19,xlab="",ylab="Chao1 species richness",xaxt='n',yaxt='n',mar=c(4,7,.5,1),cex.lab=2.5)
axis(2,at=c(1,5,20,50,100,200),labels=c(1,5,20,50,100,200),cex.axis=1.5)
text(2,200,labels='c',cex=4)
plot(Cellcounts,fishercell,log='xy',pch=19,xlab="",ylab="Fisher's alpha",xaxt='n',yaxt='n',col="black",mar=c(4,7,.5,1),cex.lab=2.5)
axis(2,at=c(0,1,5,10,50,100),labels=c(0,1,5,10,50,100),cex.axis=1.5)
text(2,70,labels='d',cex=4)
plot(Cellcounts,logcell,log='xy',pch=19,xlab="Data records per cell",ylab="Log normal",xaxt='n',yaxt='n',col="black",mar=c(4,7,.5,1),cex.lab=2.5, ylim=c(1,1200))
axis(1,at=c(1,5,20,100,250,500,1000,2500),labels=c(1,5,20,100,250,500,1000,2500),cex.axis=1.5)
axis(2,at=c(0,1,10,100,1000),labels=c(0,1,10,100,1000),cex.axis=1.5)
text(2,600,labels='e',cex=4)
plot(Cellcounts,simpcell,log='xy',pch=19,xlab="Data records per cell",ylab="Simpson's 1/D",xaxt='n',yaxt='n',mar=c(4,7,.5,1),cex.lab=2.5)
axis(2,at=c(0,1,5,10,50,100),labels=c(0,1,5,10,50,100),cex.axis=1.5)
axis(1,at=c(1,5,20,100,250,500,1000,2500),labels=c(1,5,20,100,250,500,1000,2500),cex.axis=1.5)
text(2,90,labels='f',cex=4)
dev.off()

lines(c(0,-40),exp(5.418491 + 0.072679 * c(0,-40)))
lines(c(1,40),exp(cf[1] + cf[2] * log(c(1,40))),col='red')


png('Fig 03 Lat lines.png', width=800, height=1200)
par( mgp=c(5,0.75,0), cex=1.3, cex.lab=2.5, cex.axis=1.8,las=1, mfrow=c(3,2), oma=c(0,2,0,0),mar=c(6,8.6,0.5,1))
plot(squareslat[Cellcounts>99],Cellrichness[Cellcounts>99],pch=19,log='y',xlab="",xaxt='n',ylab="Observed species richness",yaxt='n',xlim=c(-12,-39),mar=c(0,6,0,2),cex.lab=2.5)
axis(2,at=c(50,100,150,200,250),labels=c(50,100,150,200,250),cex=2.2)
text(-14,210,labels='a',cex=4)
text(-18,33,labels=expression(adjusted~R^2==0.1425),cex=1.5)
lines(c(0,-40),exp(c(rawlatline[1],-40*rawlatline[2]+rawlatline[1])),lwd=2)
plot(squareslat[Cellcounts>99],acecell[Cellcounts>99],pch=19,log='y',xaxt='n', xlab="",ylab="ACE species richness",yaxt='n',col=,xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2.5)
axis(2,at=c(50,100,150,200,250,300),labels=c(25,40,50,100,200,300),cex=2.2)
text(-14,275,labels='b',cex=4)
text(-18,38,labels=expression(adjusted~R^2==0.1515),cex=1.5)
lines(c(0,-40),exp(c(acelatline[1],-40*acelatline[2]+acelatline[1])),lwd=2)
plot(squareslat[Cellcounts>99],chao1cell[Cellcounts>99],pch=19,log='y', xlab="",ylab="Chao 1 species richness",xaxt='n', yaxt='n',col=,xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2.5)
axis(2,at=c(50,100,150,200,log(250),300),labels=c(25,40,50,100,200,300),cex=2.2)
text(-14,275,labels='c',cex=4)
text(-18,38,labels=expression(adjusted~R^2==0.134),cex=1.5)
lines(c(0,-40),exp(c(chaolatline[1],-40*chaolatline[2]+chaolatline[1])),lwd=2)
plot(squareslat[Cellcounts>99],fishercell[Cellcounts>99],pch=19,log='y',xaxt='n',xlab="",ylab="Fisher's alpha",yaxt='n',xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2.5)
axis(2,at=c(1,10,25,50,100),labels=c(1,10,25,50,100),cex=2.2)
text(-14,60,labels='d',cex=4)
text(-18,10,labels=expression(adjusted~R^2==0.3254),cex=1.5)
lines(c(0,-40),exp(c(fisherlatline[1],-40*fisherlatline[2]+fisherlatline[1])),lwd=2)
plot(squareslat[Cellcounts>99],logcell[Cellcounts>99],pch=19,log='y',xlab="Latitude",ylab="Log normal",yaxt='n',xlim=c(-12,-39),mar=c(5.5,6,0,2),cex.lab=2.5)
axis(2,at=c(10,50,100,250,500),labels=c(10,50,100,250,500),cex=2.2)
text(-14,500,labels='e',cex=4)
lines(c(0,-40),exp(c(loglatline[1],-40*loglatline[2]+loglatline[1])),lwd=2)
text(-18,45,labels=expression(adjusted~R^2==-0.0267),cex=1.5)
plot(squareslat[Cellcounts>99],simpcell[Cellcounts>99],pch=19,log='y',xlab="Latitude",ylab="Simpson's 1/D",yaxt='n',xlim=c(-12,-39),mar=c(0,6,.5,2),cex.lab=2.5)
axis(2,at=c(0,20,40,60,80,120),labels=c(0,20,40,60,80,120),cex=2.2)
text(-14,100,labels='f',cex=4)
text(-18,6,labels=expression(adjusted~R^2==0.4184),cex=1.5)
lines(c(0,-40),exp(c(simplatline[1],-40*simplatline[2]+simplatline[1])),lwd=2)
dev.off()
