library(CEMT);library(bit64);library(lme4);library(nlme);library(dplyr);library(rgdal);library(raster);library(MuMIn)

#getwd()
#setwd("C:/Lambert_Ye/UBC_Graduate/Greg's project/Analysis")
## data ====================
# ht data ---
csv <- fRead('masterfile7.csv');hd(csv)
ps <- aggregate(csv,by=list(Site=csv$Site,seedlot=csv$Seedlot),FUN='mean',na.rm=T);hd(ps)
ps2 <- ps[,c('Site','seedlot','HT10')];hd(ps2)

#site, seedlot climate ----
sclm <- fRead('original_data/Sx CC site ClimBC v5.4_Normal_1961_1990Y.csv',skip=8);hd(sclm)
pclm <- fRead('original_data/Sx CC seedlots ClimBC v5.4_Normal_1961_1990Y.csv',skip=8);hd(pclm)
names(pclm)[1] <- 'seedlot';hd(pclm)


#merge ht and climate ---
dat <- dtJoin(ps2,sclm, by = 'Site', type = "left");hd(dat)
dat2 <- dtJoin(dat,pclm, by = 'seedlot', type = "left");hd(dat2)
dat2 <- dat2[1:2140,] #omit NAs (9999s)
dat2 <- na.omit(dat2)

#dat2 removed NAs

## 1)	Top 10 most important climate variables will be identified and ranked=========
# 1 Using transfer function together 
vars <- c("MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH")
length(vars)
dat3 <- transform(dat2,td_MAT=MAT_S - MAT_P,td_MWMT=MWMT_S - MWMT_P,td_MCMT=MCMT_S - MCMT_P,
                  td_TD=TD_S - TD_P,td_MAP=MAP_S - MAP_P,td_MSP=MSP_S - MSP_P,
                  td_AHM=AHM_S - AHM_P,td_SHM=SHM_S - SHM_P,td_DD_0=DD_0_S - DD_0_P,
                  td_DD5=DD5_S - DD5_P,td_DD_18=DD_18_S - DD_18_P,td_DD18=DD18_S - DD18_P,
                  td_NFFD=NFFD_S - NFFD_P,td_bFFP=bFFP_S - bFFP_P,td_eFFP=eFFP_S - eFFP_P,
                  td_FFP=FFP_S - FFP_P,td_PAS=PAS_S - PAS_P,td_EMT=EMT_S - EMT_P,
                  td_EXT=EXT_S - EXT_P,td_Eref=Eref_S - Eref_P,td_CMD=CMD_S - CMD_P,
                  td_MAR=MAR_S - MAR_P,td_RH=RH_S - RH_P);hd(dat3) # add transferred distance

# 7-29 site, 37-59 seedlot, 60 - 82 transfer distance
Fit1 <- quadComb(dat3[,c(7:29,37:82)],dat3$HT,nComb=1,nVar=10) 
Fit1$list
Fit4 <- quadComb(dat3[,c(7:29)],dat3$HT,nComb=1,nVar=10)
Fit4$list


Fit3 <- quadComb(dat3[,c(37:59)],dat3$HT,nComb=1,nVar=10)
sglFit(dat3$TD_S,dat3$HT10)
Fit3$list
# Top 10 Variable using combined transfer function 
#> Fit1$list
#var1    r2_adj       P-value     sgma      aic
#1     td_TD 0.2880000 1.162585e-157 63.73402 23693.04 Transfer distance of temperature difference between MWMT and MCMT (°C)
#2   td_DD_0 0.2767000 2.067125e-150 64.23738 23726.47 Transfer distance of degree-days below 0°C, chilling degree-days
#3   td_MCMT 0.2663000 7.893562e-144 64.69781 23756.83 Transfer distance of mean coldest month temperature (°C)
#4    DD_0_S 0.2616000 6.511634e-141 64.90288 23770.28 Site degree-days below 0°C, chilling degree-days
#5   DD_18_S 0.2608000 2.100115e-140 64.93870 23772.62 Site degree-days below 18°C, heating degree-days
#6     MAT_S 0.2575000 2.258948e-138 65.08202 23781.99 Site mean annual temperature (°C)
#7      TD_S 0.2536000 5.781858e-136 65.25231 23793.10 Site temperature difference between MWMT and MCMT (°C)
#8    MCMT_S 0.2474000 4.047232e-132 65.52513 23810.83 Site mean coldest month temperature (°C)
#9      RH_S 0.2416000 1.363970e-128 65.77644 23827.10 Site mean annual relative humidity (%)
#10   eFFP_S 0.2320000 8.974819e-123 66.19302 23853.93 Site the day of the year on which FFP (frost-free period) ends

# Transfer distance only:
#var1   r2_adj       P-value     sgma      aic
#1     td_TD 0.288000 1.162585e-157 63.73402 23693.04 temperature difference between MWMT and MCMT (°C)
#2   td_DD_0 0.276700 2.067125e-150 64.23738 23726.47 degree-days below 0°C, chilling degree-days
#3   td_MCMT 0.266300 7.893562e-144 64.69781 23756.83 mean coldest month temperature (°C)
#4  td_DD_18 0.226500 1.666145e-119 66.42821 23869.00 degree-days below 18°C, heating degree-days
#5    td_MAT 0.222300 4.974577e-117 66.60686 23880.42 mean annual temperature (°C)
#6   td_eFFP 0.186100  4.654564e-96 68.13993 23977.13 the day of the year on which FFP ends
#7   td_NFFD 0.172500  2.209094e-88 68.70988 24012.53 the number of frost-free days
#8    td_MAP 0.170500  2.756224e-87 68.79165 24017.59 mean annual precipitation (mm)
#9     td_RH 0.157800  2.692983e-80 69.31540 24049.82 mean annual relative humidity (%)
#10   td_FFP 0.135300  3.840467e-68 70.23562 24105.87 frost-free period
## Provenance variable ranking
#> Fit3$list
#var1   r2_adj      P-value     sgma      aic
#1    MAR_P 0.022700 9.655811e-12 74.66831 24365.97
#2     TD_P 0.022340 1.434653e-11 74.68224 24366.77
#3    AHM_P 0.017590 2.445696e-09 74.86331 24377.06
#4     RH_P 0.014730 5.363283e-08 74.97232 24383.24
#5   Eref_P  0.014710 5.439945e-08 74.97283 24383.27
#6    MAP_P 0.013240 2.650070e-07 75.02879 24386.44
#7   MCMT_P 0.011250 2.248670e-06 75.10443 24390.72
#8    MSP_P 0.010710 4.009664e-06 75.12491 24391.88
#9    EXT_P 0.010190 7.055580e-06 75.14492 24393.02
#10   PAS_P 0.009723 1.158638e-05 75.16248 24394.01
par(mfrow=c(1,1))
sglFit(dat2$MCMT_P,dat2$HT10,xlab="MCMT_P")
sglFit(dat3$td_MCMT,dat2$HT10,xlab="td_MCMT")
sglFit(dat2$MCMT_S,dat2$HT10,xlab="MCMT_S")
mean(dat2$HT10)

sds <- 0
sdp <- 0
sList <- unique(dat2$Site)
pList <- unique(dat2$seedlot_name)

for(i in 1:length(sList)){
  s=sList[i]
  sds[i] <- sd(subset(dat2,Site==s)$HT10)
}
mean(sds)

for(i in 1:length(pList)){
  p=pList[i]
  sdp[i] <- sd(subset(dat2,seedlot_name==p)$HT10)
}
mean(sdp)

datt1 <- subset(dat2,Site=="Ftne"|Site=="Jord"|Site=="High"|Site=="Mayo"|Site=="Tete"|Site=="Whit")
Fitt1 <- quadComb(datt1[,c(37:59)],datt1$HT,nComb=1,nVar=10)
Fitt1$list
#> Fitt1$list
#var1     r2_adj      P-value     sgma      aic
#1     TD_P  0.0363700 4.413616e-07 93.20921 8804.688
#2   MCMT_P  0.0297800 5.432538e-06 93.52766 8809.729
#3    MAP_P  0.0258100 2.437183e-05 93.71860 8812.744
#4    PAS_P  0.0227700 7.664251e-05 93.86460 8815.044
#5    MAR_P  0.0223500 8.980373e-05 93.88482 8815.363

datt2 <- subset(dat2,Site=="Ftne"|Site=="Jord")
Fitt2 <- quadComb(datt2[,c(37:59)],datt2$HT,nComb=1,nVar=10)
Fitt2$list
#> Fitt2$list
#var1    r2_adj      P-value     sgma      aic
#1    MAR_P  0.095480 1.253654e-06 54.56326 2757.478
#2     TD_P  0.095170 1.308335e-06 54.57254 2757.565
#3   MCMT_P  0.084810 5.461817e-06 54.88413 2760.457
#4   Eref_P  0.079960 1.060326e-05 55.02937 2761.799
#5    AHM_P  0.047810 7.888981e-04 55.98234 2770.521


datt4 <- subset(dat2,Site=="Ftne")
Fitt4 <- quadComb(datt4[,c(37:59)],datt4$HT,nComb=1,nVar=10)
Fitt4$list

datt5 <- subset(dat2,Site=="Jord")
Fitt5 <- quadComb(datt5[,c(37:59)],datt5$HT,nComb=1,nVar=10)
Fitt5$list

par(mfrow=c(2,2))
sglFit(datt2$MCMT_P,datt2$HT10,xlim=c(-30,0),ylim=c(0,400))
sglFit(datt4$MCMT_P,datt4$HT10,xlim=c(-30,0),ylim=c(0,400))
sglFit(datt5$MCMT_P,datt5$HT10,xlim=c(-30,0),ylim=c(0,400))

datt3 <- subset(dat2,Site=="Ftne"|Site=="Whit")
Fitt3 <- quadComb(datt3[,c(37:59)],datt3$HT,nComb=1,nVar=10)
Fitt3$list
#> Fitt3$list
#var1   r2_adj      P-value     sgma      aic
#1   MCMT_P 0.486900 2.177783e-37 45.85570 2658.668
#2    EMT_P 0.454000 5.230507e-34 47.30592 2674.422
#3   DD_0_P 0.416000 2.314763e-30 48.92145 2691.414
#4     TD_P 0.414300 3.391107e-30 48.99623 2692.187
#5    MAP_P 0.314500 1.166992e-21 53.00411 2731.972

ps <- sample(1:107,17,replace=F)
pl <- pList[ps]
dat30 <- subset(dat3,seedlot_name %in% pl)

Fit32 <- quadComb(dat30[,c(37:59)],dat30$HT,nComb=1,nVar=10)
Fit32$list
Fit31 <- quadComb(dat30[,c(7:29)],dat30$HT,nComb=1,nVar=10)
Fit31$list


# adding a second variable
Fit2 <- quadComb(dat3[,c(7:29,37:82)],dat3$HT,nComb=3,nVar=10)

Fit2$list
# find other possible variable
urfFwd(dat3[,c(7:29,37:82)],dat3$HT10,preD=c(9,40),order=2)
head(dat3[,c(7:29,37:82)])
x <- dat3[,c(7:29,37:82)]
y <- dat3$HT10
preD=c(3,26)
urfFwd <- function(x,y,preD=c(1,2),order=2){
  xnm0 <- colnames(x)[preD];xnm0
  x1 <- x[,-c(preD)];head(x1) #remove the pre-determined var
  i=1
  for(i in 1:ncol(x1)){
    xnm <- c(xnm0,colnames(x1)[i]);xnm
    fmla <- as.formula(paste("y ~ .^2+", paste('I(',xnm,'^2)', collapse= "+")));fmla
    if(order==3){fmla <- as.formula(paste("y ~ .^2+",paste('I(',xnm,'^2)',collapse= "+"),'+',paste('I(',xnm,'^3)', collapse= "+")))};fmla
    lm1 <- lm(fmla,x[,xnm]);summary(lm1)
    r2p <- r2pFun(lm1);r2=r2p[1];pv=r2p[2];r2;pv
    nr2 <- data.frame(predictors=toString(xnm),r2,pv);nr2
    if(i==1){rr <- nr2}else{rr <- rbind(rr,nr2)};rr
  }
  rr2 <- facToNum(rr,xcol=1);#rr2
  rr2 <- sorT(rr2,byCols=2,rCols=-1);#hd(rr2)
  xnm2 <- strToList(rr2[1,1]);xnm2;length(xnm2)
  fmla2 <- as.formula(paste("y ~ .^2+", paste('I(',xnm2,'^2)', collapse= "+")));fmla2
  if(order==3){fmla2 <- as.formula(paste("y ~ .^2+", paste('I(',xnm2,'^2)', collapse= "+"),'+',paste('I(',xnm2,'^3)', collapse= "+")))}
  lm2 <- lm(fmla2,x[,xnm2]);#summary(lm2)
  r2p <- r2pFun(lm2);r2=r2p[1];pv=r2p[2];r2;pv
  p <- predict(lm2)
  length(y)
  plot(y~p);abline(lm(y~p),lwd=3,col='blue')
  legendTxt <- c(bquote(italic(r^2)==.(format(r2,digits=3))),
                 bquote(italic(p)==.(format((pv),5))));
  legend("topleft",legend=as.expression(legendTxt),bty='n',inset =- .015)
  lm2$list <- rr2;hd(rr2,10)
  return(lm2)
}
###
rr_r <- arrange(rr,pv)
#> rr_r
#predictors     r2            pv
#1     MCMT_S, MCMT_P, MAT_S 0.5407             0
#2    MCMT_S, MCMT_P, MWMT_S 0.5158             0
#3      MCMT_S, MCMT_P, TD_S 0.5203             0
#4     MCMT_S, MCMT_P, SHM_S 0.5052             0
#5    MCMT_S, MCMT_P, DD_0_S 0.5698             0
#6     MCMT_S, MCMT_P, DD5_S 0.6395             0
#7   MCMT_S, MCMT_P, DD_18_S 0.5291             0
#8    MCMT_S, MCMT_P, NFFD_S 0.6199             0
#9    MCMT_S, MCMT_P, bFFP_S 0.6154             0
#10   MCMT_S, MCMT_P, eFFP_S 0.5803             0
#11    MCMT_S, MCMT_P, FFP_S  0.612             0
#12    MCMT_S, MCMT_P, EMT_S 0.5915             0
#13    MCMT_S, MCMT_P, EXT_S 0.5567             0

## 2)	Transfer functions will be developed for both tree height and diameter at breast height; But We don't have DBH data.
#Transfer functions for single sites
#sList <- levels(dat3$Site);sList #list of site
sList <- unique(dat4$Site);sList
dat4 <- transform(dat2, tdt=MCMT_S - MCMT_P);hd(dat4)
pdfOut('transferFunction2.pdf',rowCol=c(3,2),paper='letter')
s=sList[8]
for(s in sList){
  sb <- subset(dat2,Site==s);head(sb)
  fun <- sglFit(sb$MCMT_P,sb$HT10,main=s,xlab='MCMT_P (°C)',ylab='HT10 (cm)')
  site= rbind('',site=s);site
  site2 = cat(site,summary(fun)$coefficients);site2
  if(s==sList[1]){s2=site}else{s2=rbind(s2,site)}
}
dev.off()
### outlierRM2
sList <- unique(dat4$Site);sList

pdfOut('transferFunction5.pdf',rowCol=c(3,2),paper='letter')
s=sList[8]
for(s in sList){
  sb <- subset(dat3,Site==s);head(sb)
  fun_pre <- sglFit(sb$td_MCMT,sb$HT10)
  fun_rm <- outlierRM2(sb$td_MCMT,sb$HT10,fun_pre,nst = 3)
  fun <- sglFit(fun_rm$x,fun_rm$y, main=s,xlab='td.MCMT (°C)',ylab='HT10 (cm)')
  site= rbind('',site=s);site
  site2 = cat(site,summary(fun)$coefficients);site2
  if(s==sList[1]){s2=site}else{s2=rbind(s2,site)}
}
dev.off()
sb3 <- subset(dat3,Site==s);head(sb)
sb13 <- subset(dat13,Site==s);head(sb)
sb12 <- subset(dat12,Site==s);head(sb)
par(mfrow=c(1,2))
sglFit100(sb3$td_MCMT,sb3$HT10,main="Mayo original",xlab='td.MCMT (°C)',ylab='10-year Height (cm)')
sglFit(sb13$td_MCMT,sb13$HT10,main="Mayo outlier removal",xlab='td.MCMT (°C)',ylab='10-year Height (cm)')
sglFit(sb12$td_MCMT,sb12$HT10,main="Mayo outlier removal",xlab='td.MCMT (°C)',ylab='10-year Height (cm)')


sList <- unique(dat4$Site);sList
dat4 <- dat3
pdfOut('transferFunction4.pdf',rowCol=c(3,2),paper='letter')
s=sList[1]
for(s in sList){
  sb <- subset(dat3,Site==s);head(sb)
  fun_pre <- sglFit(sb$td_MCMT,sb$HT10)
  fun_rm <- outlierRM2(sb$td_MCMT,sb$HT10,fun_pre,nst = 2)
  fun <- sglFit(fun_rm$x,fun_rm$y, main=s,xlab='td.MCMT (°C)',ylab='HT10 (cm)')
  site= rbind('',site=s);site
  site2 = cat(site,summary(fun)$coefficients);site2
  if(s==sList[1]){s2=site}else{s2=rbind(s2,site)}
}
dev.off()
i=38
sList <- unique(dat3$Site);sList
dat4 <- dat3
s=sList[1]
for(s in sList){
  sb <- subset(dat11,Site==s);head(sb)
  fun_pre <- sglFit(sb$td_MCMT,sb$HT10)
  fun_rm <- outlierRM2(sb$td_MCMT,sb$HT10,fun_pre,nst = 2)
  for (i in 1:length(sb$td_MCMT)){ 
    if((sb$td_MCMT[i] %in% fun_rm[,1]) & (sb$HT10[i] %in% fun_rm[,2])){
    }else{
      which(dat4$td_MCMT == sb$td_MCMT[i] & dat4$HT == sb$HT10[i])
      dat4 <- dat4[-which(dat4$td_MCMT == sb$td_MCMT[i] & dat4$HT == sb$HT10[i]),]
    }
  }
}
#parameter output========================================================
bk = '---------------------------------------------------------'
sink('test.txt')
for(s in sList){
  sb <- subset(dat4,Site==s);head(sb)
  fun <- sglFit(sb$tdt,sb$HT10,main=s,xlab='td.MCMT (°C)',ylab='HT10 (cm)',plot=F)
  print(paste0('Site: ',s))
  print(summary(fun)$coefficients)
  print(bk)
}
# Response functions for different provenance:
pList <- unique(dat2$seedlot);pList
dat8 <- transform(dat2, tdt=MCMT_S - MCMT_P);hd(dat8)
pdfOut('Responsefunction.pdf',rowCol=c(3,2),paper='letter')
p=pList[1]
for(p in pList){
  pb <- subset(dat4,seedlot_name==p);head(pb)
  fun <- sglFit(pb$MCMT_S,pb$HT10,main=p,xlab='MCMT_S (°C)',ylab='HT10 (cm)')
  prov= rbind('',prov=p);prov
  prov2 = cat(prov,summary(fun)$coefficients);prov2
  if(p==pList[1]){p2=prov}else{p2=rbind(p2,prov)}
}
dev.off()

# excluded weird provenances, T1
dat9 <- subset(dat3,dat3$seedlot!= 358&dat3$seedlot!= 368&dat3$seedlot!= 369&dat3$seedlot!= 364&dat3$seedlot!= 366&dat3$seedlot!= 379&
                 dat3$seedlot!= 387&dat3$seedlot!= 390&dat3$seedlot!= 393&dat3$seedlot!= 394&dat3$seedlot!= 395)
# excluded with more strict rules T2
dat10 <- subset(dat9,dat9$seedlot!=325&dat9$seedlot!=354&dat9$seedlot!=355&dat9$seedlot!=362&dat9$seedlot!=370&dat9$seedlot!=371&
                  dat9$seedlot!=372&dat9$seedlot!=374&dat9$seedlot!=377&dat9$seedlot!=378&dat9$seedlot!=381&dat9$seedlot!=382&
                  dat9$seedlot!=383&dat9$seedlot!=384&dat9$seedlot!=388&dat9$seedlot!=389&dat9$seedlot!=391&dat9$seedlot!=392&
                  dat9$seedlot!=396&dat9$seedlot!=400&dat9$seedlot!=401&dat9$seedlot!=402&dat9$seedlot!=403)


# 2) Mixed-effect transfer function ======================
#original
mix_MCMT1 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat12)
summary(mix_MCMT1) 
anova(mix_MCMT1)
r.squaredGLMM(mix_MCMT1)

URF_MCMT <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat12)
summary(URF_MCMT)
anova(URF_MCMT)

modFit(dat12,dat12$HT10,varList=c(which(names(dat12) == "td_MCMT"),which(names(dat12) == "MCMT_P")),zVar='10-year Height (cm)',IR=T,points=F)
r.squaredLR(mix_MCMT1)


# when exclude sites
mix_MCMT5 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat5)
summary(mix_MCMT5) 
anova(mix_MCMT5)
r.squaredGLMM(mix_MCMT5)

mix_MCMT6 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat6)
summary(mix_MCMT6) 
anova(mix_MCMT6)
r.squaredGLMM(mix_MCMT6)

mix_MCMT7 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat7)
summary(mix_MCMT7) 
anova(mix_MCMT7)
r.squaredGLMM(mix_MCMT7)


mix_MCMT2 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (td_MCMT|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT2) 
anova(mix_MCMT2)
r.squaredGLMM(mix_MCMT2)
r.squaredLR(mix_MCMT2)

# using two variable
mix_MCMT3 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                  + SHM_P + td_SHM + SHM_P * td_SHM + I(SHM_P^2) + I(td_SHM^2) +
                     +
                (td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat4)
summary(mix_MCMT3) 
anova(mix_MCMT3)
r.squaredGLMM(mix_MCMT3)
#> r.squaredGLMM(mix_MCMT3)
#R2m       R2c
#[1,] 0.409861 0.9519179
# SHM_S
mix_MCMT4 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                  + SHM_P + td_SHM + SHM_P * td_SHM + I(SHM_P^2) + I(td_SHM^2) +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT4) 
anova(mix_MCMT4)
r.squaredGLMM(mix_MCMT4)
#> r.squaredGLMM(mix_MCMT4)
#R2m       R2c
#[1,] 0.1808087 0.9294794

mix_MCMT11 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + td_SHM + I(td_SHM^2) + td_SHM*MCMT_P + td_SHM* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT11) 
anova(mix_MCMT11)
r.squaredGLMM(mix_MCMT11)
#> r.squaredGLMM(mix_MCMT11)
#R2m       R2c
#[1,] 0.2648185 0.9399157
mix_MCMT12 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + SHM_S + I(SHM_S^2) + SHM_S*MCMT_P + SHM_S* td_MCMT +(SHM_S|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT12) 
anova(mix_MCMT12)
r.squaredGLMM(mix_MCMT12)


mix_MCMT13 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + SHM_P +td_SHM + I(td_SHM^2)+ I(SHM_P^2) + td_SHM * MCMT_P+ td_SHM*td_MCMT+ SHM_P*MCMT_P + SHM_P* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat12)
summary(mix_MCMT13) 
anova(mix_MCMT13)
r.squaredGLMM(mix_MCMT13)
## 3)	Universal response functions will also be developed for both tree height and diameter at breast height; (DBH data is also missing)=========================
URF_MCMT <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat4)
summary(URF_MCMT)
anova(URF_MCMT)

modFit(dat4,dat4$HT10,varList=c(which(names(dat4) == "MCMT_S"),which(names(dat4) == "MCMT_P")),IR=T)
#Testing other variable
# ALL no better than R2 > 0.31

# testing for deleting not reasonable sites
# filter R2 >= 0.05 excluded: Terrace, Skimikin, Revelstoke, Nakusp, Duncanlake. (12 sites)
dat5 <- subset(dat3,dat3$Site != "Terr"& dat3$Site != "Reve" & dat3$Site != "Naku" & dat3$Site != "Skim" & dat3$Site != "Dunc" )

modFit(dat5,dat5$HT10,varList=c(which(names(dat5) == "RH_P"),which(names(dat5) == "MCMT_P")),IR=T)
# filter R2 > 0.10 excluded: Pinepass, Kalamalka (10 sites)
dat6 <- subset(dat5,dat5$Site != "Pine"& dat5$Site != "Kala")

modFit(dat6,dat6$HT10,varList=c(which(names(dat6) == "MCMT_S"),which(names(dat6) == "MCMT_P")),IR=T)

# filter R2 >= 0.15 excluded: Parsnip, Wells, Cranbrook (7 sites)
dat7 <- subset(dat6,dat6$Site != "Pars"& dat6$Site != "Well" & dat6$Site != "Cran")

modFit(dat7,dat7$HT10,varList=c(which(names(dat7) == "MCMT_S"),which(names(dat7) == "MCMT_P")),IR=T)
URF_MCMT7 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat7)
summary(URF_MCMT7) # R2 0.64 very good.
anova(URF_MCMT7)

modFit(dat9,dat9$HT10,varList=c(which(names(dat9) == "MCMT_S"),which(names(dat9) == "MCMT_P")),IR=T)
URF_MCMT9 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat9)
summary(URF_MCMT9) # R2=0.3047
anova(URF_MCMT9)

modFit(dat10,dat10$HT10,varList=c(which(names(dat10) == "MCMT_S"),which(names(dat10) == "MCMT_P")),IR=T)
URF_MCMT10 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat10)
summary(URF_MCMT10) 
anova(URF_MCMT10)

URF_MCMT2 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  MAT_S + MAT_P + I(MAT_S^2) +I(MAT_P^2) + MAT_P*MAT_S, data=dat3)
summary(URF_MCMT2)
anova(URF_MCMT2)

URF_MCMT4 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  SHM_S + SHM_P + I(SHM_S^2) +I(SHM_P^2) + SHM_P*SHM_S, data=dat3)
summary(URF_MCMT4)
anova(URF_MCMT4)



URF_MCMT5 <- lm(HT10~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  SHM_S + I(SHM_S^2)+SHM_S*MCMT_P+SHM_S*MCMT_S, data=dat3)
summary(URF_MCMT5) # 0.5073
anova(URF_MCMT5)

mix_MCMT11 <- lmer(HT10 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + td_SHM + I(td_SHM^2) + td_SHM*MCMT_P + td_SHM* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT11) 
anova(mix_MCMT11)
r.squaredGLMM(mix_MCMT11)
#> r.squaredGLMM(mix_MCMT11)
#R2m       R2c
#[1,] 0.2648185 0.9399157

#Using DD_0, AHM
URF_DD <- lm(HT10~DD_0_S +AHM_P+I(DD_0_S^2)+I(AHM_P^2)+AHM_P*DD_0_S, data=dat13)
summary(URF_DD)
modFit(dat13,dat13$HT10,varList=c(which(names(dat13) == "DD_0_S"),which(names(dat13) == "AHM_P")),zVar='10-year Height (cm)',IR=T,points=F,znj=0,zxj=3)
URF_DD <- lm(HT10~+AHM_P+I(DD_0_S^2)+I(AHM_P^2), data=dat4)
summary(URF_DD)

x <- dat3[,c(7:29,37:82)]
y <- dat3$HT10
preD=c(which(names(x)=="DD_0_S"),which(names(x)=="AHM_P"))
urfFwd <- function(x,y,preD=c(1,2),order=2)
  xnm0 <- colnames(x)[preD];xnm0
  x1 <- x[,-c(preD)];head(x1) #remove the pre-determined var
  i=1
  for(i in 1:ncol(x1)){
    xnm <- c(xnm0,colnames(x1)[i]);xnm
    fmla <- as.formula(paste("y ~ .^2+", paste('I(',xnm,'^2)', collapse= "+")));fmla
    if(order==3){fmla <- as.formula(paste("y ~ .^2+",paste('I(',xnm,'^2)',collapse= "+"),'+',paste('I(',xnm,'^3)', collapse= "+")))};fmla
    lm1 <- lm(fmla,x[,xnm]);summary(lm1)
    r2p <- r2pFun(lm1);r2=r2p[1];pv=r2p[2];r2;pv
    nr2 <- data.frame(predictors=toString(xnm),r2,pv);nr2
    if(i==1){rr <- nr2}else{rr <- rbind(rr,nr2)};rr
  }
rr#MAT_S
URF_DD2 <- lm(HT10~DD_0_S +AHM_P+I(DD_0_S^2)+I(AHM_P^2)+AHM_P*DD_0_S+MAT_S+I(MAT_S^2)+MAT_S*AHM_P+MAT_S*DD_0_S, data=dat13)
summary(URF_DD2)

URF_DD3 <- lm(HT10~DD_0_S +I(DD_0_S^2)+I(AHM_P^2)+MAT_S+I(MAT_S^2)+MAT_S*DD_0_S, data=dat13)
summary(URF_DD3)

mix_DD <- lmer(HT10~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name),REML = FALSE, dat= dat12)
summary(mix_DD)
r.squaredGLMM(mix_DD)

mix_DD2 <- lmer(HT10~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name)
                +td_MAT+I(td_MAT^2)+td_MAT*AHM_P+td_MAT*DD_0_S + (td_MAT|Site),REML = FALSE, dat= dat4)
summary(mix_DD2)
r.squaredGLMM(mix_DD2)


mix_DD1 <- lm(HT10~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0,dat= dat)
summary(mix_DD1)
r.squaredGLMM(mix_DD1)

mix_DD2 <- lmer(HT10~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name)
                +td_MAT+I(td_MAT^2)+td_MAT*AHM_P+td_MAT*DD_0_S + (td_MAT|Site),REML = FALSE, dat= dat4)
summary(mix_DD2)
r.squaredGLMM(mix_DD2)

## Using MAT_S MCMT_P
URF_F1 <- lm(HT10~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S, dat = dat2)
summary(URF_F1)

URF_F1c <- lm(HT10~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S, dat = dat4)
summary(URF_F1c)

modFit(dat4,dat4$HT10,varList=c(which(names(dat4) == "MAT_S"),which(names(dat4) == "MCMT_P")),IR=T,points=F,zVar="10-year tree height")

URF_F2 <- lm(HT10~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S , dat = dat2)
summary(URF_F2)

URF_F2c <- lm(HT10~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S , dat = dat4)
summary(URF_F2c)


#5) mapping ============
library(climatenaAPI)
varList <- c("MAT","DD5","MCMT")
rasterDownload(region='BC',res='4000m',period='Normal_1961_1990',varList=varList,destLoc="C:\\Lambert_Ye\\UBC_Graduate\\Greg's project\\Analysis")
wd <- "C:\\Lambert_Ye\\UBC_Graduate\\Greg's project\\Analysis\\BC\\4000m\\Normal_1961_1990SY"
stk_BC <- rasterStack(wd,varList,rType='tif',vConvert=T)
stk_WNA <- rasterStack(wd,varList,rType='tif',vConvert=T) 
stk_NA <- rasterStack(wd,varList,rType='tif',vConvert=T) 

### local seedlots
stk_BC2 <- stack(subset(stk_BC,1)/10,subset(stk_BC,2),subset(stk_BC,3)/10) 
#names      :    MAT,    DD5,   MCMT 
#min values :  -12.7,    7.0,  -25.3 
#max values :   10.4, 2525.0,    4.6 

which(names(dat4) == "MCMT_S")
dat21 <- dat4[,c(3,7,16,39)]
names(dat21) <- c("HT10","MAT","DD5","MCMT")
URF_F2M <- lm(HT10~MAT + MCMT + I(MAT^2) + I(MCMT^2) + MCMT * MAT + DD5 + I(DD5^2) + DD5 * MCMT + DD5 * MAT , dat = dat21)
summary(URF_F2M)

p1 <- predict(stk_BC2,URF_F2M)
p1[p1<0] <- 0
p1[p1>500] <- 500
p10 <- p1/100
plot(p10,main="Local seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                 font = 2, line = 2.5, cex = 0.8))

URF_F3 <- lm(HT10~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S +
             MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * DD5_S, dat = dat2)
summary(URF_F3)

dat22 <- dat4[,c(3,7,9,16,39)]
stk_BC3 <- stack(subset(stk_BC,1)/10,subset(stk_BC,2),subset(stk_BC,3)/10,subset(stk_BC,3)/10)
names(stk_BC3) <- c("MAT_S","DD5_S","MCMT_S","MCMT_P")

p3 <- predict(stk_BC3,URF_F3)
p3[p3<0] <- 0
p3[p3>500] <- 500
p30 <- p3/100
plot(p30,main="Local seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                 font = 2, line = 2.5, cex = 0.8))

b1=14.97;b2=-0.2537;b5=-6.365;b9=1.885e-2;b12=2.648
stk_BC4 <- stack(subset(stk_BC3,1),subset(stk_BC3,2),subset(stk_BC3,3),-(b1+b5*subset(stk_BC3,1)+b9*subset(stk_BC3,2)+b12*subset(stk_BC3,3))*0.5/b2)
names(stk_BC4) <- c("MAT_S","DD5_S","MCMT_S","MCMT_P")

p4 <- predict(stk_BC4,URF_F3)
p4[p4<0] <- 0
p4[p4>500] <- 500
p40 <- p4/100
plot(p40,main="Optimal seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

plot(subset(stk_BC4,4))
MCMTP <- subset(stk_BC4,4)
MCMTP[MCMTP>4.6] <- 4.6
MCMTP[MCMTP< -25.3] <- -25.3
plot(MCMTP,main="Optimal seedlot mean coldest month temperature",legend.args = list(text = 'seedlot MCMT (°C)', side = 4, 
                                                                          font = 2, line = 2.5, cex = 0.8))

stk_BC5 <- stack(subset(stk_BC3,1),subset(stk_BC3,2),subset(stk_BC3,3),MCMTP)
P4 <- predict(stk_BC5,URF_F3)

P4[P4<0] <- 0
P4[P4>500] <- 500
P40 <- P4/100
plot(P40,main="Optimal seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))


p50 <- p40 - p30
plot(p50,main="Height gain of using optimal seedlot",legend.args = list(text = '10-year tree height gain (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))

P50 <- P40 - p30
plot(P50,main="Height gain of using optimal seedlot",legend.args = list(text = '10-year tree height gain (m)', side = 4, 
                                                                        font = 2, line = 2.5, cex = 0.8))

par(mfrow=c(1,2))
plot(P40,main="Optimal seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))
plot(P50,main="Height gain of using optimal seedlot",legend.args = list(text = '10-year tree height gain (m)', side = 4, 
                                                                        font = 2, line = 2.5, cex = 0.8))

################################ ignore below
dat21 <- dat13[,c(3,7,15,43)]
names(dat21) <- c("HT10","MAT","DD_0","AHM")
URF_map1 <- lm(HT10~AHM+DD_0+I(DD_0^2)+I(AHM^2)+AHM*DD_0+MAT+I(MAT^2)+MAT*AHM+MAT*DD_0, data=dat21)
summary(URF_map1) # model with MAT
URF_map2 <- lm(HT10~AHM+DD_0+I(DD_0^2)+I(AHM^2)+AHM*DD_0, data=dat21)
summary(URF_map2) # original model


p1 <- predict(stk_BC2,URF_map1)
p1[p1<0] <- 0
p1[p1>500] <- 500
p10 <- p1/100
plot(p10,main="Local seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))
legend("10-year tree height")
### optimal seedlots
B1=1.474;B2=-5.159e-2;B5=1.090e-3;B9=0.2455
stk_BC3 <- stack(subset(stk_BC2,1),subset(stk_BC2,2),-(B1+B5*subset(stk_BC2,2)+B9*subset(stk_BC2,1))/(2*B2))
names(stk_BC3)[3] <- "AHM"
#names      :       MAT,      DD_0,       AHM 
#min values : -12.70000,  79.00000,  24.48546 
#max values :     10.40,   4469.00,     41.28 

P1 <- predict(stk_BC3,URF_map1)
P1[P1<0] <- 0
P1[P1>500] <- 500
P10 <- P1/100
plot(P1)

P2 <- (P1-p1)/100
par(mfrow=c(1,2),mai=c(1,0.5,1,1.5))
plot(P10,main="Optimal seedlot growth potential",legend.args = list(text = '10-year tree height (m)', side = 4, 
                                                                font = 2, line = 2.5, cex = 0.8))
plot(P2,main="Growth potential improvement of optimal seedlot",legend.args = list(text = 'Tree height gain (m)', side = 4, 
                                                                 font = 2, line = 2.5, cex = 0.8))

####################### TEST original model
p2 <- predict(stk_BC2,URF_map2)
p2[p2<0] <- 0
plot(p2)
b0=162.5;b1=3.464;b2=-4.936e-2;b3=-6.108e-3;b4=-1.271e-5;b5=-1.648e-6
stk_BC4 <- stack(subset(stk_BC2,1),subset(stk_BC2,2),-(b1+b5*subset(stk_BC2,2))/(2*b2))
names(stk_BC4)[3] <- "AHM"
names(stk_BC4)[2] <- "DD_0"

P2 <- predict(stk_BC4,URF_map2)
P2[P2<0] <- 0
plot(P2)
###########TEST
summary(URF_map1)
dat22 <- dat21
dat22$AHM <- -(1.474+1.090e-3*dat21$DD_0+0.2455*dat21$MAT)/(2*(-5.159e-2))


URF_MAP1 <- lm(HT10~AHM+DD_0+I(DD_0^2)+I(AHM^2)+AHM*DD_0+MAT+I(MAT^2)+MAT*AHM+MAT*DD_0, data=dat22)
summary(URF_MAP1)

P1 <- predict(stk_BC,URF_MAP1)
P1[P1<0] <- 0
plot(P1)
P4 <- predict(stk_WNA,URF_MAP1)
P4[P4<0] <- 0
plot(P4)

summary(URF_map2)
dat23 <- dat21
dat23$AHM <- -(3.464-1.648e-6*dat21$DD_0)/(2*-4.936e-2)
URF_MAP2 <- lm(HT10~AHM+DD_0+I(DD_0^2)+I(AHM^2)+AHM*DD_0, data=dat23)
summary(URF_MAP2)
URF_MAP2 <- lm(HT10~2*AHM+DD_0+I(DD_0^2)+I(AHM^2)+AHM*DD_0, data=dat23)

b0=162.5
b1=3.464
b2=-4.936e-2
b3=-6.108e-3
b4=-1.271e-5
b5=-1.648e-6


DD_0 <- c(100,200,300,400,500,600,700,1000,1500,2000,3000)
dat23 <- data.frame(DD_0)
dat23$HT <- B0+B1*dat23$DD_0+B2*dat23$DD_0*dat23$DD_0
URF_MAP2 <- lm(HT~DD_0+I(DD_0^2),dat=dat23)
summary(URF_MAP2)
P2 <- predict(stk_BC,URF_MAP2)
P2[P2<0] <- 0
plot(P2)
P2[P2<100] <- 100
plot(P2)
P5 <- predict(stk_WNA,URF_MAP2)
P5[P5<0] <- 0
plot(P5)





