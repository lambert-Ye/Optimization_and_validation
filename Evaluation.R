library(CEMT);library(bit64);library(lme4);library(nlme);library(dplyr);library(rgdal);library(raster);library(MuMIn);library(caret)
install.packages("recipes")
#getwd()
#setwd("C:/Lambert_Ye/UBC_Graduate/Paper_2_Evaluation/work_dir")
## data ====================
# ht data ---
csv <- fRead('masterfile16.csv');hd(csv)
csv16 <- csv[,c(3,8,16)];hd(csv16) #site,prov,height
csv16$HT16[csv16$HT16 == 0|csv16$HT16 == "."] <- NA
csv16[,3] <- as.numeric(csv16[,3])

ps <- aggregate(csv16,by=list(Site=csv16$Site,seedlot=csv16$Prov),FUN='mean',na.rm=T);hd(ps)
ps2 <- ps[,c('Site','seedlot','HT16')];hd(ps2)

#site, seedlot climate ----
sclm <- fRead('Sx CC site ClimBC v5.4_Normal_1961_1990Y.csv',skip=8);hd(sclm)
pclm <- fRead('Sx CC seedlots ClimBC v5.4_Normal_1961_1990Y.csv',skip=8);hd(pclm)
names(pclm)[1] <- 'seedlot';hd(pclm)


#merge ht and climate ---
dat <- dtJoin(ps2,sclm, by = 'Site', type = "left");hd(dat)
dat2 <- dtJoin(dat,pclm, by = 'seedlot', type = "left");hd(dat2)
dat2 <- dat2[1:2140,] #omit NAs (9999s)
dat2 <- na.omit(dat2)
#dat2 removed NAs

dat2

## Testing ========

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
Fit1 <- quadComb(dat3[,c(7:29,37:82)],dat3$HT16,nComb=1,nVar=10) 
Fit1$list

#> Fit1$list
#var1    r2_adj       P-value     sgma      aic
#1     td_TD 0.2523000 1.965060e-134 126.4934 26468.59
#2    DD_0_S 0.2346000 9.995924e-124 127.9792 26517.97
#3   td_DD_0 0.2297000 8.930437e-121 128.3918 26531.58
#4      TD_S 0.2194000 1.108766e-114 129.2481 26559.68
#5   td_MCMT 0.2156000 1.837532e-112 129.5614 26569.92
#6    MCMT_S 0.1842000  1.643858e-94 132.1233 26652.70
#7   DD_18_S 0.1712000  3.194744e-87 133.1779 26686.32
#8     MAT_S 0.1682000  1.462355e-85 133.4193 26693.98
#9  td_DD_18 0.1645000  1.502510e-83 133.7124 26703.25
#10   td_MAT 0.1592000  1.255406e-80 134.1392 26716.73

Fit4 <- quadComb(dat3[,c(7:29)],dat3$HT,nComb=1,nVar=10)
Fit4$list

#> Fit4$list
#var1    r2_adj       P-value     sgma      aic
#1   DD_0_S 0.2346000 9.995924e-124 127.9792 26517.97
#2     TD_S 0.2194000 1.108766e-114 129.2481 26559.68
#3   MCMT_S 0.1842000  1.643858e-94 132.1233 26652.70
#4  DD_18_S 0.1712000  3.194744e-87 133.1779 26686.32
#5    MAT_S 0.1682000  1.462355e-85 133.4193 26693.98
#6    MAR_S 0.1418000  2.946253e-71 135.5173 26759.94
#7     RH_S 0.1360000  3.485730e-68 135.9723 26774.11
#8   eFFP_S 0.0899000  2.421240e-44 139.5547 26884.07
#9   NFFD_S 0.0777400  2.926896e-38 140.4837 26912.12
#10   MSP_S 0.0738600  2.483204e-36 140.7795 26921.01


Fit3 <- quadComb(dat3[,c(37:59)],dat3$HT,nComb=1,nVar=10)
sglFit(dat3$TD_S,dat3$HT16)
Fit3$list

#> Fit3$list
#var1    r2_adj      P-value     sgma      aic
#1     TD_P 0.0330000 1.531052e-16 143.8514 27012.28
#2    MAR_P 0.0287700 1.535831e-14 144.1658 27021.50
#3    AHM_P 0.0240400 2.577764e-12 144.5161 27031.77
#4     RH_P 0.0232700 5.929122e-12 144.5731 27033.43
#5    MAP_P 0.0215600 3.765631e-11 144.6998 27037.14
#6   Eref_P 0.0197500 2.648593e-10 144.8335 27041.04
#7    EXT_P 0.0160100 1.469273e-08 145.1093 27049.09
#8    PAS_P 0.0152700 3.249935e-08 145.1639 27050.68
#9    MSP_P 0.0143700 8.500117e-08 145.2300 27052.60
#10  MCMT_P 0.0142400 9.791409e-08 145.2398 27052.89

## function plots======

sList <- unique(dat3$Site);sList
pList <- unique(dat3$seedlot);pList

pdfOut('responseFunction.pdf',rowCol=c(3,2),paper='letter')
p=pList[8]
for(p in pList){
  pb <- subset(dat3,seedlot==p);head(pb)
  fun <- sglFit(pb$MCMT_S,pb$HT16,main=p,xlab='MCMT_S (°C)',ylab='HT16 (cm)')
  prov= rbind('',prov=p);prov
  prov2 = cat(prov,summary(fun)$coefficients);prov2
  if(p==pList[1]){p2=prov}else{p2=rbind(p2,prov)}
}
dev.off()

pdfOut('transferFunction.pdf',rowCol=c(3,2),paper='letter')
s=sList[8]
for(s in sList){
  sb <- subset(dat3,Site==s);head(sb)
  fun <- sglFit(sb$td_MCMT,sb$HT16,main=s,xlab='MCMT_td (°C)',ylab='HT16 (cm)')
  site= rbind('',site=s);site
  site2 = cat(site,summary(fun)$coefficients);site2
  if(s==sList[1]){s2=site}else{s2=rbind(s2,site)}
}
dev.off()

#test WLS
s=sList[17]
sb <- subset(dat3,Site==s);head(sb)
sglFit(sb$td_MCMT,sb$HT16,main=s,xlab='MCMT_td (°C)',ylab='HT16 (cm)')
sglFit(sb$td_MCMT,sb$HT16,weights = 1/sb$HT16,main=s,xlab='MCMT_td (°C)',ylab='HT16 (cm)')


ols <- lm(HT16~td_MCMT + I(td_MCMT^2),sb)
summary(ols)
plot(ols)[1]

wls <- lm(HT16~td_MCMT + I(td_MCMT^2),weights = 1/HT16,sb)
summary(wls)
points(sb$td_MCMT,fitted(lm(HT16~td_MCMT + I(td_MCMT^2),weights = 1/HT16,sb)),pch=19)

nls <- nls(HT16 ~ A/(1+(td_MCMT-B)^2/C^2),data = sb, 
           start = list(A = 594, C = 15, B = sb$td_MCMT[which(sb$HT16 == max(sb$HT16))]))
summary(nls)
points(sb$td_MCMT,fitted(nls(HT16 ~ A/(1+(td_MCMT-B)^2/C^2),data = sb, 
                             start = list(A = 594, C = 15, B = sb$td_MCMT[which(sb$HT16 == max(sb$HT16))]))),pch=19,col = "red")




# test spatial correlation and variogram (seedlots 301, 302, 420, 427)
library(gstat)
#301
p=pList[1]
p
pb <- subset(dat3,seedlot==p);head(pb);dim(pb)
lmod <- lm(HT16~td_MCMT + I(td_MCMT^2),pb) #original OLS
sglFit(pb$td_MCMT,pb$HT16)
E <- rstandard(lmod)
(mydata <- data.frame(E,lat_S = pb[,4],long_S = pb[,5]))
coordinates(mydata) <- c("lat_S","long_S")

bubble(mydata,"E",col = c("black","grey"),main = "residuals") 
vario1 <- variogram(E ~ 1, mydata)
plot(vario1)
#302
p=pList[2]
#420
p=pList[119]
#427
p=pList[126]

#test weibull
x <- 1:300/100
a=2;b=2;c=10;d=-0.9;x1=(x-min(x))/(max(x)-min(x))
y <- c*(x1-d)^(a-1)*exp(-(x1-d)^a/b)
plot(x,y,pch=19,cex=0.2)


# Mixed-effect transfer function ======
#original
mix_MCMT1 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT1) 
anova(mix_MCMT1)
r.squaredGLMM(mix_MCMT1)
modFit(dat3,dat3$HT16,varList=c(which(names(dat3) == "td_MCMT"),which(names(dat3) == "MCMT_P")),zVar='16-year Height (cm)',IR=T,points=F)
r.squaredLR(mix_MCMT1)

mix_MCMT2 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = TRUE, data = dat3)
summary(mix_MCMT2) 
anova(mix_MCMT2)
r.squaredGLMM(mix_MCMT2)
modFit(dat3,dat3$HT16,varList=c(which(names(dat3) == "td_MCMT"),which(names(dat3) == "MCMT_P")),zVar='16-year Height (cm)',IR=T,points=F)
r.squaredLR(mix_MCMT2)

mix_MCMT3 <- gls(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), data = dat3)
summary(mix_MCMT3) 

# when exclude sites
mix_MCMT5 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat5)
summary(mix_MCMT5) 
anova(mix_MCMT5)
r.squaredGLMM(mix_MCMT5)

mix_MCMT6 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat6)
summary(mix_MCMT6) 
anova(mix_MCMT6)
r.squaredGLMM(mix_MCMT6)

mix_MCMT7 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat7)
summary(mix_MCMT7) 
anova(mix_MCMT7)
r.squaredGLMM(mix_MCMT7)


mix_MCMT2 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT + (td_MCMT|Site) + (td_MCMT|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT2) 
anova(mix_MCMT2)
r.squaredGLMM(mix_MCMT2)
r.squaredLR(mix_MCMT2)

# using two variable
mix_MCMT3 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                  + SHM_P + td_SHM + SHM_P * td_SHM + I(SHM_P^2) + I(td_SHM^2) +
                    +
                    (td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT3) 
anova(mix_MCMT3)
r.squaredGLMM(mix_MCMT3)
#> r.squaredGLMM(mix_MCMT3)
#R2m       R2c
#[1,] 0.409861 0.9519179
# SHM_S
mix_MCMT4 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                  + SHM_P + td_SHM + SHM_P * td_SHM + I(SHM_P^2) + I(td_SHM^2) +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT4) 
anova(mix_MCMT4)
r.squaredGLMM(mix_MCMT4)
#> r.squaredGLMM(mix_MCMT4)
#R2m       R2c
#[1,] 0.1808087 0.9294794

mix_MCMT11 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + td_SHM + I(td_SHM^2) + td_SHM*MCMT_P + td_SHM* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT11) 
anova(mix_MCMT11)
r.squaredGLMM(mix_MCMT11)
#> r.squaredGLMM(mix_MCMT11)
#R2m       R2c
#[1,] 0.2648185 0.9399157
mix_MCMT12 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + SHM_S + I(SHM_S^2) + SHM_S*MCMT_P + SHM_S* td_MCMT +(SHM_S|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT12) 
anova(mix_MCMT12)
r.squaredGLMM(mix_MCMT12)


mix_MCMT13 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + SHM_P +td_SHM + I(td_SHM^2)+ I(SHM_P^2) + td_SHM * MCMT_P+ td_SHM*td_MCMT+ SHM_P*MCMT_P + SHM_P* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat12)
summary(mix_MCMT13) 
anova(mix_MCMT13)
r.squaredGLMM(mix_MCMT13)
## 3)	Universal response functions will also be developed for both tree height and diameter at breast height; (DBH data is also missing) =====
URF_MCMT <- lm(HT16~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat3)
summary(URF_MCMT)
anova(URF_MCMT)

modFit(dat3,dat3$HT16,varList=c(which(names(dat3) == "MCMT_S"),which(names(dat3) == "MCMT_P")),IR=T)
#Testing other variable
# ALL no better than R2 > 0.31

# testing for deleting not reasonable sites
# filter R2 >= 0.05 excluded: Terrace, Skimikin, Revelstoke, Nakusp, Duncanlake. (15 sites)
dat4 <- subset(dat3,dat3$Site != "Naku" &  dat3$Site != "Dunc" )

# filter R2 >= 0.05 excluded: Terrace, Skimikin, Revelstoke, Nakusp, Duncanlake. (12 sites)
dat5 <- subset(dat3,dat3$Site != "Terr"& dat3$Site != "Reve" & dat3$Site != "Naku" & dat3$Site != "Skim" & dat3$Site != "Dunc" )

modFit(dat5,dat5$HT16,varList=c(which(names(dat5) == "MCMT_S"),which(names(dat5) == "MCMT_P")),IR=T)
# filter R2 > 0.10 excluded: Pinepass, Kalamalka (10 sites)
dat6 <- subset(dat5,dat5$Site != "Pine"& dat5$Site != "Kala")

modFit(dat6,dat6$HT16,varList=c(which(names(dat6) == "MCMT_S"),which(names(dat6) == "MCMT_P")),IR=T)

# filter R2 >= 0.15 excluded: Parsnip, Wells, Cranbrook (7 sites)
dat7 <- subset(dat6,dat6$Site != "Pars"& dat6$Site != "Well" & dat6$Site != "Cran")

modFit(dat7,dat7$HT16,varList=c(which(names(dat7) == "MCMT_S"),which(names(dat7) == "MCMT_P")),IR=T)
URF_MCMT7 <- lm(HT16~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S, data=dat7)
summary(URF_MCMT7) # R2 0.62 very good.
anova(URF_MCMT7)

modFit(dat3,dat3$HT16,varList=c(which(names(dat3) == "MAR_S"),which(names(dat3) == "MAR_P")),IR=T)
URF_MAR3 <- lm(HT16~MAR_P + MAR_S + I(MAR_P^2) + I(MAR_S^2) + MAR_P*MAR_S, data=dat3)
summary(URF_MAR3) 
anova(URF_MAR3)

modFit(dat7,dat7$HT16,varList=c(which(names(dat7) == "MAR_S"),which(names(dat7) == "MAR_P")),IR=T)
URF_MAR7 <- lm(HT16~MAR_P + MAR_S + I(MAR_P^2) + I(MAR_S^2) + MAR_P*MAR_S, data=dat7)
summary(URF_MAR7) 
anova(URF_MAR7)

URF_MCMT2 <- lm(HT16~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  MAT_S + MAT_P + I(MAT_S^2) +I(MAT_P^2) + MAT_P*MAT_S, data=dat3)
summary(URF_MCMT2)
anova(URF_MCMT2)

URF_MCMT4 <- lm(HT16~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  SHM_S + SHM_P + I(SHM_S^2) +I(SHM_P^2) + SHM_P*SHM_S, data=dat3)
summary(URF_MCMT4)
anova(URF_MCMT4)



URF_MCMT5 <- lm(HT16~MCMT_P + MCMT_S + I(MCMT_P^2) + I(MCMT_S^2) + MCMT_P*MCMT_S + 
                  SHM_S + I(SHM_S^2)+SHM_S*MCMT_P+SHM_S*MCMT_S, data=dat3)
summary(URF_MCMT5) # 0.5073
anova(URF_MCMT5)

mix_MCMT11 <- lmer(HT16 ~ MCMT_P + td_MCMT + I(MCMT_P^2) + I(td_MCMT^2) + MCMT_P * td_MCMT
                   + td_SHM + I(td_SHM^2) + td_SHM*MCMT_P + td_SHM* td_MCMT +(td_SHM|Site) + (td_MCMT|Site) + (1|seedlot_name), REML = FALSE, data = dat3)
summary(mix_MCMT11) 
anova(mix_MCMT11)
r.squaredGLMM(mix_MCMT11)
#> r.squaredGLMM(mix_MCMT11)
#R2m       R2c
#[1,] 0.2648185 0.9399157

#Using DD_0, AHM
URF_DD <- lm(HT16~DD_0_S +AHM_P+I(DD_0_S^2)+I(AHM_P^2)+AHM_P*DD_0_S, data=dat13)
summary(URF_DD)
modFit(dat13,dat13$HT16,varList=c(which(names(dat13) == "DD_0_S"),which(names(dat13) == "AHM_P")),zVar='16-year Height (cm)',IR=T,points=F,znj=0,zxj=3)
URF_DD <- lm(HT16~+AHM_P+I(DD_0_S^2)+I(AHM_P^2), data=dat3)
summary(URF_DD)

x <- dat3[,c(7:29,37:82)]
y <- dat3$HT16
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
URF_DD2 <- lm(HT16~DD_0_S +AHM_P+I(DD_0_S^2)+I(AHM_P^2)+AHM_P*DD_0_S+MAT_S+I(MAT_S^2)+MAT_S*AHM_P+MAT_S*DD_0_S, data=dat13)
summary(URF_DD2)

URF_DD3 <- lm(HT16~DD_0_S +I(DD_0_S^2)+I(AHM_P^2)+MAT_S+I(MAT_S^2)+MAT_S*DD_0_S, data=dat13)
summary(URF_DD3)

mix_DD <- lmer(HT16~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name),REML = FALSE, dat= dat12)
summary(mix_DD)
r.squaredGLMM(mix_DD)

mix_DD2 <- lmer(HT16~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name)
                +td_MAT+I(td_MAT^2)+td_MAT*AHM_P+td_MAT*DD_0_S + (td_MAT|Site),REML = FALSE, dat= dat3)
summary(mix_DD2)
r.squaredGLMM(mix_DD2)


mix_DD1 <- lm(HT16~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0,dat= dat)
summary(mix_DD1)
r.squaredGLMM(mix_DD1)

mix_DD2 <- lmer(HT16~td_DD_0 +AHM_P+I(td_DD_0^2)+I(AHM_P^2)+AHM_P*td_DD_0+ (td_DD_0|Site)+(1|seedlot_name)
                +td_MAT+I(td_MAT^2)+td_MAT*AHM_P+td_MAT*DD_0_S + (td_MAT|Site),REML = FALSE, dat= dat3)
summary(mix_DD2)
r.squaredGLMM(mix_DD2)

## Using MAT_S MCMT_P
URF_F1 <- lm(HT16~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S, dat = dat2)
summary(URF_F1)

URF_F1c <- lm(HT16~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S, dat = dat3)
summary(URF_F1c)

modFit(dat3,dat3$HT16,varList=c(which(names(dat3) == "MAT_S"),which(names(dat3) == "MCMT_P")),IR=T,points=F,zVar="16-year tree height")

URF_F2 <- lm(HT16~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S , dat = dat2)
summary(URF_F2)

URF_F2c <- lm(HT16~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S , dat = dat3)
summary(URF_F2c)

URF_16 <- lm(HT16~MAT_S + MCMT_P + I(MAT_S^2) + I(MCMT_P^2) + MCMT_P * MAT_S + 
               DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S + MCMT_S + I(MCMT_S^2) +
               MCMT_S * MCMT_P + MCMT_S * MAT_S + MCMT_S * DD5_S, dat = dat3)
summary(URF_16)

####################
#5) mapping
library(climatenaAPI)
varList <- c("MAT","DD5","MCMT")
rasterDownload(region='BC',res='4000m',period='Normal_1961_1990',varList=varList,destLoc="C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir")
wd <- "C:\\Lambert_Ye\\UBC_Graduate\\Greg's project\\Analysis\\BC\\4000m\\Normal_1961_1990SY"
stk_BC <- rasterStack(wd,varList,rType='tif',vConvert=T)
stk_WNA <- rasterStack(wd,varList,rType='tif',vConvert=T) 
stk_NA <- rasterStack(wd,varList,rType='tif',vConvert=T) 

raster_AHM <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\AHM.tif")
raster_bFFP <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\bFFP.tif")
raster_CMD <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\CMD.tif")
raster_CMI <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\CMI.tif")
raster_DD_0 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD_0.tif")
raster_DD_18 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD_18.tif")
raster_DD5 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5.tif")
raster_DD18 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD18.tif")
raster_DD1040 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD1040.tif")
raster_eFFP <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\eFFP.tif")
raster_EMT <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\EMT.tif")
raster_Eref <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\Eref.tif")
raster_FFP <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\FFP.tif")
raster_MAP <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP.tif")
raster_MAT <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT.tif")
raster_MCMT <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT.tif")
raster_MSP <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MSP.tif")
raster_MWMT <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MWMT.tif")
raster_NFFD <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\NFFD.tif")
raster_PAS <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\PAS.tif")
raster_RH <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\RH.tif")
raster_SHM <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\SHM.tif")
raster_TD <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\TD.tif")
plot(raster_TD)
raster <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\Site_Prod_Se1.tif")
plot(raster)
raster
newproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
test_raster <- projectRaster(raster,crs=newproj)
plot(test_raster)
plot(raster_TD)
### local seedlots
stk_BC2 <- stack(subset(stk_BC,1)/10,subset(stk_BC,2),subset(stk_BC,3)/10) 
#names      :    MAT,    DD5,   MCMT 
#min values :  -12.7,    7.0,  -25.3 
#max values :   10.4, 2525.0,    4.6 

which(names(dat3) == "MCMT_S")
dat21 <- dat3[,c(3,7,16,39)]
names(dat21) <- c("HT16","MAT","DD5","MCMT")
URF_F2M <- lm(HT16~MAT + MCMT + I(MAT^2) + I(MCMT^2) + MCMT * MAT + DD5 + I(DD5^2) + DD5 * MCMT + DD5 * MAT , dat = dat21)
summary(URF_F2M)

p1 <- predict(stk_BC2,URF_F2M)
p1[p1<0] <- 0
p1[p1>500] <- 500
p10 <- p1/100
plot(p10,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S +
               MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * DD5_S, dat = dat5)
summary(URF_16)

URF_16_2 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
               MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_2) #0.665

plot(stk_BC2,2)

dat22 <- dat3[,c(3,7,9,16,39)]
stk_BC3 <- stack(subset(stk_BC,1)/10,subset(stk_BC,2),subset(stk_BC,3)/10,subset(stk_BC,3)/10)
names(stk_BC3) <- c("MAT_S","DD5_S","MCMT_S","MCMT_P")

p3 <- predict(stk_BC3,URF_16)
p3[p3<0] <- 0
p3[p3>500] <- 500
p30 <- p3/100
plot(p30,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#optimal
b1=28.5;b2=-0.6003;b5=-13.35;b9=3.896e-2;b12=5.453
stk_BC4 <- stack(subset(stk_BC3,1),subset(stk_BC3,2),subset(stk_BC3,3),-(b1+b5*subset(stk_BC3,1)+b9*subset(stk_BC3,2)+b12*subset(stk_BC3,3))*0.5/b2)
names(stk_BC4) <- c("MAT_S","DD5_S","MCMT_S","MCMT_P")

p4 <- predict(stk_BC4,URF_16)
p4[p4<0] <- 0
p4[p4>500] <- 500
p40 <- p4/100
plot(p40,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))
p5 <- p4-p3
plot(p5,main="Optimal seedlot growth gain",legend.args = list(text = '16-year tree height (m)', side = 4, 
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
plot(P40,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))


p50 <- p40 - p30
plot(p50,main="Height gain of using optimal seedlot",legend.args = list(text = '16-year tree height gain (m)', side = 4, 
                                                                        font = 2, line = 2.5, cex = 0.8))

P50 <- P40 - p30
plot(P50,main="Height gain of using optimal seedlot",legend.args = list(text = '16-year tree height gain (m)', side = 4, 
                                                                        font = 2, line = 2.5, cex = 0.8))

par(mfrow=c(1,2))
plot(P40,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                    font = 2, line = 2.5, cex = 0.8))
plot(P50,main="Height gain of using optimal seedlot",legend.args = list(text = '16-year tree height gain (m)', side = 4, 
                                                                       font = 2, line = 2.5, cex = 0.8))
#====
# testing for deleting not reasonable sites
# filter R2 >= 0.05 excluded: Terrace, Skimikin, Revelstoke, Nakusp, Duncanlake. (12 sites)
dat5 <- subset(dat3,dat3$Site != "Terr"& dat3$Site != "Reve" & dat3$Site != "Naku" & dat3$Site != "Skim" & dat3$Site != "Dunc" )

modFit(dat5,dat5$HT16,varList=c(which(names(dat5) == "MCMT_S"),which(names(dat5) == "MCMT_P")),IR=T)
# filter R2 > 0.10 excluded: Pinepass, Kalamalka (10 sites)
dat6 <- subset(dat5,dat5$Site != "Pine"& dat5$Site != "Kala")

modFit(dat6,dat6$HT16,varList=c(which(names(dat6) == "MCMT_S"),which(names(dat6) == "MCMT_P")),IR=T)

# filter R2 >= 0.15 excluded: Parsnip, Wells, Cranbrook (7 sites)
dat7 <- subset(dat6,dat6$Site != "Pars"& dat6$Site != "Well" & dat6$Site != "Cran")

# assigning raster stacks
varList <- c("AHM","bFFP","CMD","CMI","DD_0", #5
             "DD_18","DD5","DD18","DD1040","eFFP", #10
             "EMT","Eref","FFP","MAP","MAT", #15
             "MCMT","MSP","MWMT","NFFD","PAS", #20
             "RH","SHM","TD")
wd <- "C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY"
stk_BC <- rasterStack(wd,varList,rType='tif',vConvert=T)
subset(stk_BC,1)#AHM

install.packages("cvTools")
library(cvTools)
detach("raster")
?cv
results_1 <- cvFit(URF_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_1
results_2 <- cvFit(URF_16_2, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_2
results_3 <- cvFit(URF_16_3, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_3
results_4 <- cvFit(URF_16_4, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4
results_5 <- cvFit(URF_16_5, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_5
results_6 <- cvFit(URF_16_6, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_6
results_7 <- cvFit(URF_16_7, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_7

#Candidate models (TD_P better)================
#1 MCMT_S,MAT_S,DD5_S - MCMT_P, R2=0.74 R2_5=0.79 (previous) CV_MAE = 74.96708
URF_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S +
               MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * DD5_S, dat = dat3)
summary(URF_1)
stk_BC2 <- stack(subset(stk_BC,7),subset(stk_BC,15)/10,subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("DD5_S","MAT_S","MCMT_S","MCMT_P")
p10 <- predict(stk_BC2,URF_1)
p10[p10<100] <- 100;p10[p10>800] <- 800
p100 <- p10/100
plot(p100,main="Local seedlot growth potential with URF model [1]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

aURF_1.5 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S +
              MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * DD5_S, dat = dat5)
summary(URF_1.5)
p15 <- predict(stk_BC2,URF_1.5)
p15[p15<100] <- 100;p15[p15>800] <- 800
p150 <- p15/100
plot(p150,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#2 MCMT_S,MAP_S,MAT_S - MCMT_P, R2=0.67 R2_5=0.75 (using MAP) CV_MAE = 84.90856 
URF_16_2 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_2)

stk_BC2 <- stack(subset(stk_BC,14),subset(stk_BC,15)/10,subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAP_S","MAT_S","MCMT_S","MCMT_P")
p20 <- predict(stk_BC2,URF_16_2)
p20[p20<100] <- 100;p20[p20>800] <- 800
p200 <- p20/100
plot(p200,main="Local seedlot growth potential with URF model [2]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_2.5 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat5)
summary(URF_16_2.5)
p25 <- predict(stk_BC2,URF_16_2.5)
p25[p25<100] <- 100;p25[p25>800] <- 800
p250 <- p25/100
plot(p250,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#3 MCMT_S,MAP_S,DD5_S - MCMT_P, R2=0.67 R2_5=0.78  CV_MAE = 84.65619
URF_16_3 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_3)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,7),subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
p30 <- predict(stk_BC2,URF_16_3)
p30[p30<100] <- 100;p30[p30>800] <- 800
p300 <- p30/100
plot(p300,main="Local seedlot growth potential  with URF model [3]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_3.4 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                   MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat4)
summary(URF_16_3.4)
p34 <- predict(stk_BC2,URF_16_3.4)
p34[p34<100] <- 100;p34[p34>800] <- 800
p340 <- p34/100
plot(p340,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))

URF_16_3.5 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                   MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat5)
summary(URF_16_3.5)
p35 <- predict(stk_BC2,URF_16_3.5)
p35[p35<100] <- 100;p35[p35>800] <- 800
p350 <- p35/100
plot(p350,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#4 TD_S,MAT_S,DD5_S - TD_P, R2=0.73 R2_5=0.77 (comparing #1, MCMT changed to TD) CV_MAE = 75.27863
URF_16_4 <- lm(HT16~ TD_P+ I(TD_P^2) + DD5_S + I(DD5_S^2)  + TD_P * DD5_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * DD5_S +
                 TD_S + I(TD_S^2) + TD_S * DD5_S + TD_S * TD_P + TD_S * MAT_S, dat = dat3)
summary(URF_16_4)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,7),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","DD5_S","TD_S","TD_P")
p40 <- predict(stk_BC2,URF_16_4)
p40[p40<100] <- 100;p40[p40>800] <- 800
p400 <- p40/100
plot(p400,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_4.5 <- lm(HT16~ TD_P+ I(TD_P^2) + DD5_S + I(DD5_S^2)  + TD_P * DD5_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * DD5_S +
                   TD_S + I(TD_S^2) + TD_S * DD5_S + TD_S * TD_P + TD_S * MAT_S, dat = dat5)
summary(URF_16_4.5)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,7),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","DD5_S","TD_S","TD_P")
p45 <- predict(stk_BC2,URF_16_4.5)
p45[p45<100] <- 100;p45[p45>800] <- 800
p450 <- p45/100
plot(p450,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#5 TD_S,MAT_S,MSP_S - TD_P, R2=0.59 R2_5=0.73 CV_MAE = 94.55284
URF_16_5 <- lm(HT16~ TD_P+ I(TD_P^2) + MSP_S + I(MSP_S^2)  + TD_P * MSP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MSP_S +
                 TD_S + I(TD_S^2) + TD_S * MSP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat3)
summary(URF_16_5)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,17),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","MSP_S","TD_S","TD_P")
p50 <- predict(stk_BC2,URF_16_5)
p50[p50<100] <- 100;p50[p50>800] <- 800
p500 <- p50/100
plot(p500,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_5.5 <- lm(HT16~ TD_P+ I(TD_P^2) + MSP_S + I(MSP_S^2)  + TD_P * MSP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MSP_S +
                   TD_S + I(TD_S^2) + TD_S * MSP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat5)
summary(URF_16_5.5)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,17),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","MSP_S","TD_S","TD_P")
p55 <- predict(stk_BC2,URF_16_5.5)
p55[p55<100] <- 100;p55[p55>800] <- 800
p550 <- p55/100
plot(p550,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

#6 TD_S,MAT_S,MAP_S - TD_P, R2=0.65 R2_5=0.768  CV_MAE = 86.88973
URF_16_6 <- lm(HT16~ TD_P+ I(TD_P^2) + MAP_S + I(MAP_S^2)  + TD_P * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MAP_S +
                 TD_S + I(TD_S^2) + TD_S * MAP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat3)
summary(URF_16_6)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","MAP_S","TD_S","TD_P")
p60 <- predict(stk_BC2,URF_16_6)
p60[p60<100] <- 100;p60[p60>800] <- 800
p600 <- p60/100
plot(p600,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_6.5 <- lm(HT16~ TD_P+ I(TD_P^2) + MAP_S + I(MAP_S^2)  + TD_P * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MAP_S +
                   TD_S + I(TD_S^2) + TD_S * MAP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat5)
summary(URF_16_6.5)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","MAP_S","TD_S","TD_P")
p65 <- predict(stk_BC2,URF_16_6.5)
p65[p65<100] <- 100;p65[p65>800] <- 800
p650 <- p65/100
plot(p650,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))
#7 MCMT_S,MSP_S,MAT_S - MCMT_P, R2=0.56 R2_5=0.74 (using MSP) CV_MAE = 96.85483
URF_16_7 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MSP_S + I(MSP_S^2) + MSP_S * MCMT_P + MSP_S * MAT_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MSP_S, dat = dat3)
summary(URF_16_7)

stk_BC2 <- stack(subset(stk_BC,17),subset(stk_BC,15)/10,subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MSP_S","MAT_S","MCMT_S","MCMT_P")
p70 <- predict(stk_BC2,URF_16_7)
p70[p70<100] <- 100;p70[p70>800] <- 800
p700 <- p70/100
plot(p700,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))

URF_16_7.5 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MSP_S + I(MSP_S^2) + MSP_S * MCMT_P + MSP_S * MAT_S +
                   MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MSP_S, dat = dat5)
summary(URF_16_7.5)
p75 <- predict(stk_BC2,URF_16_7.5)
p75[p75<100] <- 100;p75[p75>800] <- 800
p750 <- p75/100
plot(p750,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))
?par
par(mfrow=c(2,2),mai=c(0.5,0.5,0.5,1.2))
plot(p100,main="Local seedlot growth potential with URF [1]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                                       font = 2, line = 2.5, cex = 0.8))
plot(p300,main="Local seedlot growth potential with URF [3]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                                       font = 2, line = 2.5, cex = 0.8))
plot(p500,main="Local seedlot growth potential with URF [5]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                                       font = 2, line = 2.5, cex = 0.8))
plot(p600,main="Local seedlot growth potential with URF [6]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                                       font = 2, line = 2.5, cex = 0.8))
par(mfrow=c(1,1))
# can solve the problem that some locations with low cite densities.

#Spatial correlation ===============
#Raster calculation from GIS
#test
raster <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\SP_Se.tif")
plot(raster)
raster
newproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
test_raster <- projectRaster(raster,crs=newproj)
plot(test_raster)


bt <- test_raster
ndvi <- subset(stk_BC,14)
bt <- resample(bt, ndvi, method="bilinear") 
# Masking rasters
ndvi_m <- mask(ndvi, bt)
bt_m <- mask(bt, ndvi)
# Plot
plot(ndvi_m)
plot(bt_m) 

overlay <- stack(ndvi_m, bt_m)
overlay <- data.frame(na.omit(values(overlay))) 
names(overlay) <- c("ndvi", "bt")
# correlation test
cor.test(overlay[,1], overlay[,2])
#import rasters
Se_o <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\SP_Se.tif")
Sw_o <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\SP_Sw.tif")
Sx_o <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\SP_Sx.tif")
newproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
Se <- projectRaster(Se_o,crs=newproj)
Sw <- projectRaster(Sw_o,crs=newproj)
Sx <- projectRaster(Sx_o,crs=newproj)
Se[is.na(Se)] <- 0
Sw[is.na(Sw)] <- 0
Sx[is.na(Sx)] <- 0
Sxx <- Sxx_o <- max(Se,Sw,Sx)
Sxx_o[Sxx_o==0] <- NA
plot(Sxx_o)

Sxx_o <- resample(Sxx_o,p10,method="bilinear")
Sxx_o <- mask(Sxx_o,p10)
plot(Sxx_o)
oversxx_o <- stack(Sxx_o,p10,p15,p20,p25,p30,p34,p40,p45,p50,p55,p60,p65,p70,p75)
oversxx_o <- data.frame(na.omit(values(oversxx_o))) 
a <- cor.test(oversxx_o[,1], oversxx_o[,2]);a$estimate #0.40
a <- cor.test(oversxx_o[,1], oversxx_o[,3]);a$estimate #0.45
a <- cor.test(oversxx_o[,1], oversxx_o[,4]);a$estimate #0.46
a <- cor.test(oversxx_o[,1], oversxx_o[,5]);a$estimate #0.46
a <- cor.test(oversxx_o[,1], oversxx_o[,6]);a$estimate #0.61
a <- cor.test(oversxx_o[,1], oversxx_o[,7]);a$estimate #0.46
a <- cor.test(oversxx_o[,1], oversxx_o[,8]);a$estimate #0.41
a <- cor.test(oversxx_o[,1], oversxx_o[,9]);a$estimate #0.25
a <- cor.test(oversxx_o[,1], oversxx_o[,10]);a$estimate #0.47
a <- cor.test(oversxx_o[,1], oversxx_o[,11]);a$estimate #0.55
a <- cor.test(oversxx_o[,1], oversxx_o[,12]);a$estimate #0.50
a <- cor.test(oversxx_o[,1], oversxx_o[,13]);a$estimate #0.15
a <- cor.test(oversxx_o[,1], oversxx_o[,14]);a$estimate #0.53
a <- cor.test(oversxx_o[,1], oversxx_o[,15]);a$estimate #0.55

Sxx <- resample(Sxx,p10,method="bilinear")
Sxx <- mask(Sxx,p10)
plot(Sxx)
oversxx <- stack(Sxx,p10,p15,p20,p25,p30,p34,p40,p45,p50,p55,p60,p65,p70,p75)
oversxx <- data.frame(na.omit(values(oversxx))) 
a <- cor.test(oversxx[,1], oversxx[,2]);a$estimate #0.16
a <- cor.test(oversxx[,1], oversxx[,3]);a$estimate #0.41
a <- cor.test(oversxx[,1], oversxx[,4]);a$estimate #0.32
a <- cor.test(oversxx[,1], oversxx[,5]);a$estimate #0.33
a <- cor.test(oversxx[,1], oversxx[,6]);a$estimate #0.50
a <- cor.test(oversxx[,1], oversxx[,7]);a$estimate #0.37
a <- cor.test(oversxx[,1], oversxx[,8]);a$estimate #0.12
a <- cor.test(oversxx[,1], oversxx[,9]);a$estimate #-0.01
a <- cor.test(oversxx[,1], oversxx[,10]);a$estimate #-0.02
a <- cor.test(oversxx[,1], oversxx[,11]);a$estimate #0.45
a <- cor.test(oversxx[,1], oversxx[,12]);a$estimate #-0.08
a <- cor.test(oversxx[,1], oversxx[,13]);a$estimate #-0.01
a <- cor.test(oversxx[,1], oversxx[,14]);a$estimate #0.08
a <- cor.test(oversxx[,1], oversxx[,15]);a$estimate #0.50


# Reseach question 1: What is the best URF model for interior spruce? MCMT_S,MAP_S,DD5_S - MCMT_P, R2=0.67
URF_16_3 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_3)

stk_BC2 <- stack(subset(stk_BC,15),subset(stk_BC,7),subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
p30 <- predict(stk_BC2,URF_16_3)
p30[p30<100] <- 100;p30[p30>800] <- 800
p300 <- p30/100
plot(p300,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
# Reseach question 2: Why we use spatial validation?
#1):For models with similar R2, provide further comparison
#example
URF_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + DD5_S + I(DD5_S^2) + DD5_S * MCMT_P + DD5_S * MAT_S +
              MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * DD5_S, dat = dat3)
summary(URF_1)
stk_BC2 <- stack(subset(stk_BC,7),subset(stk_BC,15)/10,subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("DD5_S","MAT_S","MCMT_S","MCMT_P")
p10 <- predict(stk_BC2,URF_1)
p10[p10<100] <- 100;p10[p10>800] <- 800
p100 <- p10/100
plot(p100,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
#R2=0.73,even higher than the selected one.

#2): for locations with less site density, this might be super important criterion. For studies using provenance trials with a lot of site, this might be less problematic
URF_16_6 <- lm(HT16~ TD_P+ I(TD_P^2) + MAP_S + I(MAP_S^2)  + TD_P * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MAP_S +
                 TD_S + I(TD_S^2) + TD_S * MAP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat3)
summary(URF_16_6)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("MAT_S","MAP_S","TD_S","TD_P")
p60 <- predict(stk_BC2,URF_16_6)
p60[p60<100] <- 100;p60[p60>800] <- 800
p600 <- p60/100
plot(p600,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
#R2=0.65 for this one. But it expand to coastal area, which is bad.

#3): Application of 1):
URF_16_3 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_3)

stk_BC2 <- stack(subset(stk_BC,15),subset(stk_BC,7),subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
p30 <- predict(stk_BC2,URF_16_3)
p30[p30<100] <- 100;p30[p30>800] <- 800
p300 <- p30/100
plot(p300,main="Local seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))

#first-order derivative of MCMT_P
b1=-1.270e+01;b2=-5.091e-01;b5=-2.691e-03;b9=4.875e-03;b12=4.264e-01
stk_BC3 <- stack(subset(stk_BC,15),subset(stk_BC,7),subset(stk_BC,16)/10,-(b1+b5*subset(stk_BC,7)+b9*subset(stk_BC,15)+b12*subset(stk_BC,16)/10)*0.5/b2)
names(stk_BC3) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
p30 <- predict(stk_BC2,URF_16_3)
p31 <- predict(stk_BC3,URF_16_3)
p31[p31<100] <- 100;p31[p31>900] <- 900
p310 <- p31/100
plot(p310,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
p31 <- predict(stk_BC3,URF_16_3)
p31[p31<100] <- 100;p30[p30<100] <- 100
p32 <- p31-p30
p320 <- p32/100
p320 <- mask(p320,Sxx_o)
p320[is.na(p320)] <- 0
p320 <- mask(p320,Sxx)
plot(p320,main="Growth potential gain",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                     font = 2, line = 2.5, cex = 0.8))

p33 <- -(b1+b5*subset(stk_BC,7)+b9*subset(stk_BC,15)+b12*subset(stk_BC,16)/10)*0.5/b2

plot(p33,main="Optimal seed source",legend.args = list(text = 'Mean Coldest Month Temperature (°C)', side = 4, 
                                                          font = 2, line = 2.5, cex = 0.8))

#4): Another question: how will the model fit change with adding more variable? Both R2 and spatial correlation
#Ref:
Fit4 <- quadComb(dat3[,c(7:29)],dat3$HT,nComb=1,nVar=10)
Fit4$list

#> Fit4$list
#var1    r2_adj       P-value     sgma      aic
#1   DD_0_S 0.2346000 9.995924e-124 127.9792 26517.97
#2     TD_S 0.2194000 1.108766e-114 129.2481 26559.68
#3   MCMT_S 0.1842000  1.643858e-94 132.1233 26652.70
#4  DD_18_S 0.1712000  3.194744e-87 133.1779 26686.32
#5    MAT_S 0.1682000  1.462355e-85 133.4193 26693.98
#6    MAR_S 0.1418000  2.946253e-71 135.5173 26759.94
#7     RH_S 0.1360000  3.485730e-68 135.9723 26774.11
#8   eFFP_S 0.0899000  2.421240e-44 139.5547 26884.07
#9   NFFD_S 0.0777400  2.926896e-38 140.4837 26912.12
#10   MSP_S 0.0738600  2.483204e-36 140.7795 26921.01


Fit3 <- quadComb(dat3[,c(37:59)],dat3$HT,nComb=1,nVar=10)
sglFit(dat3$TD_S,dat3$HT16)
Fit3$list

#> Fit3$list
#var1    r2_adj      P-value     sgma      aic
#1     TD_P 0.0330000 1.531052e-16 143.8514 27012.28
#2    MAR_P 0.0287700 1.535831e-14 144.1658 27021.50
#3    AHM_P 0.0240400 2.577764e-12 144.5161 27031.77
#4     RH_P 0.0232700 5.929122e-12 144.5731 27033.43
#5    MAP_P 0.0215600 3.765631e-11 144.6998 27037.14
#6   Eref_P 0.0197500 2.648593e-10 144.8335 27041.04
#7    EXT_P 0.0160100 1.469273e-08 145.1093 27049.09
#8    PAS_P 0.0152700 3.249935e-08 145.1639 27050.68
#9    MSP_P 0.0143700 8.500117e-08 145.2300 27052.60
#10  MCMT_P 0.0142400 9.791409e-08 145.2398 27052.89

varList <- c("AHM","bFFP","CMD","CMI","DD_0", #5
             "DD_18","DD5","DD18","DD1040","eFFP", #10
             "EMT","Eref","FFP","MAP","MAT", #15
             "MCMT","MSP","MWMT","NFFD","PAS", #20
             "RH","SHM","TD")

## Selecting 10 candidate models for each number of parameters
nparameters <- data.frame(c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),c(1,2,3,4,5,6,7,8,9,0),
                          c(1,2,3,4,5,6,7,8,9,0))
names(nparameters) <- c("R2_1","PE_1","olc_1","whc_1","R2_2","PE_2","olc_2","whc_2","R2_3","PE_3","olc_3","whc_3",
                        "R2_4","PE_4","olc_4","whc_4","R2_5","PE_5","olc_5","whc_5","R2_6","PE_6","olc_6","whc_6",
                        "R2_7","PE_7","olc_7","whc_7") #R2,olc - overlapping correlation, whc - whole correlation

### 1] 2 variables - too few available model, just 5 
URF_2_1 <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P, dat = dat3)
s21 <- summary(URF_2_1) #R2=0.23
stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P")
pu21 <- predict(stk_BC2,URF_2_1)

URF_2_2 <- lm(HT16 ~ MAT_P + I(MAT_P^2) + MAT_S + I(MAT_S^2) + MAT_S * MAT_P, dat = dat3)
s22 <- summary(URF_2_2) #R2=0.19
stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,15)/10)
names(stk_BC2) <- c("MAT_S","MAT_P")
pu22 <- predict(stk_BC2,URF_2_2)

URF_2_3 <- lm(HT16 ~ TD_P + I(TD_P^2) + TD_S + I(TD_S^2) + TD_S * TD_P, dat = dat3)
s23 <- summary(URF_2_3) #R2=0.28
stk_BC2 <- stack(subset(stk_BC,23)/10,subset(stk_BC,23)/10)
names(stk_BC2) <- c("TD_S","TD_P")
pu23 <- predict(stk_BC2,URF_2_3)

URF_2_4 <- lm(HT16 ~ MAP_P + I(MAP_P^2) + MAP_S + I(MAP_S^2) + MAP_S * MAP_P, dat = dat3)
s24 <- summary(URF_2_4) #R2=0.10
stk_BC2 <- stack(subset(stk_BC,14)/10,subset(stk_BC,14)/10)
names(stk_BC2) <- c("MAP_S","MAP_P")
pu24 <- predict(stk_BC2,URF_2_4)

URF_2_5 <- lm(HT16 ~ DD_0_P + I(DD_0_P^2) + DD_0_S + I(DD_0_S^2) + DD_0_S * DD_0_P, dat = dat3)
s25 <- summary(URF_2_5) #R2=0.26
stk_BC2 <- stack(subset(stk_BC,5)/10,subset(stk_BC,5)/10)
names(stk_BC2) <- c("DD_0_S","DD_0_P")
pu25 <- predict(stk_BC2,URF_2_5)

oversxx <- stack(Sxx,pu21,pu22,pu23,pu24,pu25)
oversxx <- data.frame(na.omit(values(oversxx))) 
wh21 <- cor.test(oversxx[,1], oversxx[,2]);wh21$estimate 
wh22 <- cor.test(oversxx[,1], oversxx[,3]);wh22$estimate 
wh23 <- cor.test(oversxx[,1], oversxx[,4]);wh23$estimate 
wh24 <- cor.test(oversxx[,1], oversxx[,5]);wh24$estimate 
wh25 <- cor.test(oversxx[,1], oversxx[,6]);wh25$estimate

oversxx_o <- stack(Sxx_o,pu21,pu22,pu23,pu24,pu25)
oversxx_o <- data.frame(na.omit(values(oversxx_o))) 
ol21 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol21$estimate 
ol22 <- cor.test(oversxx_o[,1], oversxx_o[,3]);ol22$estimate 
ol23 <- cor.test(oversxx_o[,1], oversxx_o[,4]);ol23$estimate 
ol24 <- cor.test(oversxx_o[,1], oversxx_o[,5]);ol24$estimate 
ol25 <- cor.test(oversxx_o[,1], oversxx_o[,6]);ol25$estimate 

nparameters$R2_1 <- c(s21$adj.r.squared,s22$adj.r.squared,s23$adj.r.squared,s24$adj.r.squared,s25$adj.r.squared)
nparameters$PE_1 <- c(s21$sigma,s22$sigma,s23$sigma,s24$sigma,s25$sigma)
nparameters$olc_1 <- c(ol21$estimate,ol22$estimate,ol23$estimate,ol24$estimate,ol25$estimate)
nparameters$whc_1 <- c(wh21$estimate,wh22$estimate,wh23$estimate,wh24$estimate,wh25$estimate)

# 2] 3 variables
URF_3_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S +
                            MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P, dat = dat3)
summary(URF_3_1)

# 3] 4 variables
URF_4_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_4_1)

# 4] 5 variables
URF_5_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                 MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat3)
summary(URF_5_1)

# 5] 6 variables
URF_6_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P + MAT_P + I(MAT_P^2) + MAT_P * MCMT_P + MAT_P * MCMT_S + MAT_P * MAT_S +
                MAT_P * DD5_S + MAT_P * MAP_S, dat = dat3)
summary(URF_6_1)

# 6] 7 variables
URF_7_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P + AHM_P + I(AHM_P^2) + AHM_P * MCMT_P + AHM_P * MCMT_S + AHM_P * MAT_S +
                AHM_P * DD5_S + AHM_P * MAP_S + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MCMT_S + PAS_S * AHM_P + PAS_S * DD5_S +
                PAS_S * MAT_S + PAS_S * MAP_S + PAS_S * AHM_P, dat = dat3)
summary(URF_7_1)

# Test blockCV
install.packages("progress")
URF_16_3 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                 MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)
summary(URF_16_3)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,7),subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
p30 <- predict(stk_BC2,URF_16_3)
p30[p30<100] <- 100;p30[p30>800] <- 800
p300 <- p30/100
plot(p300,main="Local seedlot growth potential  with URF model [3]",legend.args = list(text = '16-year tree height (m)', side = 4, 
                                                                                       font = 2, line = 2.5, cex = 0.8))

library(raster);library(blockCV)
p300r <- aggregate(p300,fact=20,fun=mean)
po300 <- as.data.frame(rasterToPoints(p300r)) 
pc_data <- sf::st_as_sf(po300, coords = c("x", "y"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sb <- cv_spatial(x = pc_data,
                  column = "layer",
                  size = 500000,
                  k = 5,
                  selection = "random",
                  iteration = 50)
folds <- po300$fold <- sb$folds_ids

for (k in 1:5){
  train.data <- po300[which(fold != k),]
  test.data <- po300[which(fold ==k),]
  
}

#test for leave one out cv:
R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                            MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "DD5_S" | names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  R_sq[k] <-  b$adj.r.squared
  RMSE[k] <- b$sigma
}
mean(R_sq)
mean(RMSE)
plot(put44)
#unify legend
put41 <- pu41
put41[put41>max(values(put41)-2,na.rm = T)] <- 1000
ppt41 <- pp41
ppt41[ppt41>max(values(ppt41)-2,na.rm = T)] <- 1000
put42 <- pu42
put42[put42>max(values(put42)-20,na.rm = T)] <- 1000
put43 <- pu43
put43[put43>max(values(put43)-2,na.rm = T)] <- 1000
put44 <- pu44
put44[put44>max(values(put44)-25,na.rm = T)] <- 1000
plot(put44)
put45 <- pu45
put45[put45>max(values(put45)-2,na.rm = T)] <- 1000
pft_126_1 <- pf_126_1
pft_126_1[pft_126_1>max(values(pft_126_1)-5,na.rm = T)] <- 1000
pft_245_1 <- pf_245_1
pft_245_1[pft_245_1>max(values(pft_245_1)-5,na.rm = T)] <- 1000
pft_370_1 <- pf_370_1
pft_370_1[pft_370_1>max(values(pft_370_1)-5,na.rm = T)] <- 1000
pft_585_1 <- pf_585_1
pft_585_1[pft_585_1>max(values(pft_585_1)-5,na.rm = T)] <- 1000
pfot_126_1 <- pfo_126_1
pfot_126_1[pfot_126_1>max(values(pfot_126_1)-5,na.rm = T)] <- 1000
pfot_245_1 <- pfo_245_1
pfot_245_1[pfot_245_1>max(values(pfot_245_1)-5,na.rm = T)] <- 1000
pfot_370_1 <- pfo_370_1
pfot_370_1[pfot_370_1>max(values(pfot_370_1)-5,na.rm = T)] <- 1000
pfot_585_1 <- pfo_585_1
pfot_585_1[pfot_585_1>max(values(pfot_585_1)-5,na.rm = T)] <- 1000
pft_126_7 <- pf_126_7
pft_126_7[pft_126_7>max(values(pft_126_7)-5,na.rm = T)] <- 1000
pfot_126_7 <- pfo_126_7
pfot_126_7[pfot_126_7>max(values(pfot_126_7)-5,na.rm = T)] <- 1000
pft_126_4 <- pf_126_4
pft_126_4[pft_126_4>max(values(pft_126_4)-5,na.rm = T)] <- 1000
pfot_126_4 <- pfo_126_4
pfot_126_4[pfot_126_4>max(values(pfot_126_4)-5,na.rm = T)] <- 1000
pft_126_7 <- pf_126_7
pft_126_7[pft_126_7>max(values(pft_126_7)-5,na.rm = T)] <- 1000
pfot_126_7 <- pfo_126_7
pfot_126_7[pfot_126_7>max(values(pfot_126_7)-5,na.rm = T)] <- 1000
pft_245_4 <- pf_245_4
pft_245_4[pft_245_4>max(values(pft_245_4)-5,na.rm = T)] <- 1000
pfot_245_4 <- pfo_245_4
pfot_245_4[pfot_245_4>max(values(pfot_245_4)-5,na.rm = T)] <- 1000
pft_245_7 <- pf_245_7
pft_245_7[pft_245_7>max(values(pft_245_7)-5,na.rm = T)] <- 1000
pfot_245_7 <- pfo_245_7
pfot_245_7[pfot_245_7>max(values(pfot_245_7)-10,na.rm = T)] <- 1000
pft_370_4 <- pf_370_4
pft_370_4[pft_370_4>max(values(pft_370_4)-5,na.rm = T)] <- 1000
pfot_370_4 <- pfo_370_4
pfot_370_4[pfot_370_4>max(values(pfot_370_4)-10,na.rm = T)] <- 1000
pft_370_7 <- pf_370_7
pft_370_7[pft_370_7>max(values(pft_370_7)-5,na.rm = T)] <- 1000
pfot_370_7 <- pfo_370_7
pfot_370_7[pfot_370_7>max(values(pfot_370_7)-5,na.rm = T)] <- 1000
pft_585_4 <- pf_585_4
pft_585_4[pft_585_4>max(values(pft_585_4)-5,na.rm = T)] <- 1000
pfot_585_4 <- pfo_585_4
pfot_585_4[pfot_585_4>max(values(pfot_585_4)-10,na.rm = T)] <- 1000
pft_585_7 <- pf_585_7
pft_585_7[pft_585_7>max(values(pft_585_7)-5,na.rm = T)] <- 1000
pfot_585_7 <- pfo_585_7
pfot_585_7[pfot_585_7>max(values(pfot_585_7)-5,na.rm = T)] <- 1000

for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAT_S + I(MAT_S^2) + MAT_S * MCMT_P + MAT_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAT_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "DD5_S" | names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAT_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  R_sq[k] <-  b$adj.r.squared
  RMSE[k] <- b$sigma
}
mean(R_sq)
mean(RMSE)
RMSE_1 <- R_sq_1 <- c(1)
k <- 10
long_vector <- 1:2114
folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
for (k in 1:10){#k=1
  dat_test <- dat3[folds[[k]],]
  dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
  
  model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  R_sq_1[k] <-  b$adj.r.squared
  RMSE_1[k] <- b$sigma
}
mean(RMSE_1)
var(RMSE_1)
plot(DD5_126_4)
# future scenarios (13GCM):
DD5_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (1).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (2).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (3).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (4).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (5).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (6).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (7).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (8).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MCMT_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (1).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (2).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (3).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (4).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (5).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (6).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (7).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (8).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MAT_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (1).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (2).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (3).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (4).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (5).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (6).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (7).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (8).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MAP_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (1).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (2).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (3).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (4).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (5).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (6).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (7).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (8).tif") %>% resample(p10,method="bilinear") %>% mask(p10)


BC_stack_F <- stack(DD5_126_4,DD5_126_7,DD5_245_4,DD5_245_7,DD5_370_4,DD5_370_7,DD5_585_4,DD5_585_7,
                    MCMT_126_4,MCMT_126_7,MCMT_245_4,MCMT_245_7,MCMT_370_4,MCMT_370_7,MCMT_585_4,MCMT_585_7,
                    MAT_126_4,MAT_126_7,MAT_245_4,MAT_245_7,MAT_370_4,MAT_370_7,MAT_585_4,MAT_585_7,
                    MAP_126_4,MAP_126_7,MAP_245_4,MAP_245_7,MAP_370_4,MAP_370_7,MAP_585_4,MAP_585_7)

# 8 GCMs
DD5_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (10).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (11).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (12).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (13).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (14).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (15).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (16).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
DD5_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (17).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MCMT_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (10).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (11).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (12).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (13).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (14).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (15).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (16).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (17).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MAT_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (10).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (11).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (12).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (13).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (14).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (15).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (16).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (17).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

MAP_126_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (10).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_126_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (11).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_245_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (12).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_245_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (13).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_370_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (14).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_370_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (15).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_585_4 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (16).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_585_7 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (17).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

DD5_126_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (18).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_126_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (18).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_126_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (18).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_126_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (18).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

DD5_245_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (19).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_245_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (19).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_245_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (19).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_245_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (19).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

DD5_370_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (20).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_370_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (20).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_370_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (20).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_370_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (20).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

DD5_585_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\DD5 (21).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MCMT_585_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MCMT (21).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAT_585_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAT (21).tif") %>% resample(p10,method="bilinear") %>% mask(p10)
MAP_585_1 <- raster("C:\\Lambert_Ye\\UBC_Graduate\\Paper_2_Evaluation\\work_dir\\BC\\800m\\Normal_1961_1990SY\\MAP (21).tif") %>% resample(p10,method="bilinear") %>% mask(p10)

BC_stack_F <- stack(DD5_126_4,DD5_126_7,DD5_245_4,DD5_245_7,DD5_370_4,DD5_370_7,DD5_585_4,DD5_585_7,
                    MCMT_126_4,MCMT_126_7,MCMT_245_4,MCMT_245_7,MCMT_370_4,MCMT_370_7,MCMT_585_4,MCMT_585_7,
                    MAT_126_4,MAT_126_7,MAT_245_4,MAT_245_7,MAT_370_4,MAT_370_7,MAT_585_4,MAT_585_7,
                    MAP_126_4,MAP_126_7,MAP_245_4,MAP_245_7,MAP_370_4,MAP_370_7,MAP_585_4,MAP_585_7,
                    DD5_126_1,MCMT_126_1,MAT_126_1,MAP_126_1,DD5_245_1,MCMT_245_1,MAT_245_1,MAP_245_1,
                    DD5_370_1,MCMT_370_1,MAT_370_1,MAP_370_1,DD5_585_1,MCMT_585_1,MAT_585_1,MAP_585_1)
