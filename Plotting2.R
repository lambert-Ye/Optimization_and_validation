# Table 1 Candidate models ===================
## [1]
URF_2_1 <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P, dat = dat3)
s21 <- summary(URF_2_1) #R2=0.23
s21;AIC(URF_2_1)
results_2_1 <- cvFit(URF_2_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_2_1
RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
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
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),];length(dat_train$HT16)
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),];length(dat_test$HT16)
  
  model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P") ]
  test_result <- stats::predict(model,test)

  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)
stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P")
pu21 <- predict(stk_BC2,URF_2_1)
pu21[pu21<0] <- 0;pu21[pu21>1000] <- 1000
plot(pu21)
oversxx <- stack(Sxx,pu21)
oversxx <- data.frame(na.omit(values(oversxx)))
wh21 <- cor.test(oversxx[,1], oversxx[,2]);wh21$estimate 
oversxx_o <- stack(Sxx_o,pu21)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol21 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol21$estimate 

## [2]
URF_3_1 <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + 
                MAT_S * MCMT_P, dat = dat3)
s31 <- summary(URF_3_1) 
s31;AIC(URF_3_1)
results_3_1 <- cvFit(URF_3_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_3_1

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + 
                  MAT_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAT_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + 
                MAT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAT_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S")
pu31 <- predict(stk_BC2,URF_3_1)
pu31[pu31<0] <- 0;pu31[pu31>1000] <- 1000
plot(pu31)
oversxx <- stack(Sxx,pu31)
oversxx <- data.frame(na.omit(values(oversxx)))
wh31 <- cor.test(oversxx[,1], oversxx[,2]);wh31$estimate 
oversxx_o <- stack(Sxx_o,pu31)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol31 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol31$estimate 

## [3]
URF_3_2 <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_S + 
                MAP_S * MCMT_P, dat = dat3)
s32 <- summary(URF_3_2) 
s32;AIC(URF_3_2)
results_3_2 <- cvFit(URF_3_2, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_3_2

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_S + 
                  MAP_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MCMT_P + I(MCMT_P^2) + MCMT_S + I(MCMT_S^2) + MCMT_S * MCMT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_S + 
                MAP_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)
stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,14))
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAP_S")
pu32 <- predict(stk_BC2,URF_3_2)
pu32[pu32<0] <- 0;pu32[pu32>1000] <- 1000
plot(pu32)
oversxx <- stack(Sxx,pu32)
oversxx <- data.frame(na.omit(values(oversxx)))
wh32 <- cor.test(oversxx[,1], oversxx[,2]);wh32$estimate 
oversxx_o <- stack(Sxx_o,pu32)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol32 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol32$estimate 


## [4]
URF_4_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)   + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S+ MCMT_P * DD5_S + MAP_S, dat = dat3)

s41 <- summary(URF_4_1) 
s41;AIC(URF_4_1)
results_4_1 <- cvFit(URF_4_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4_1

URF_4_1 <- lm(HT16~    , dat = dat3)

mod <- lm(HT16~ I(DD5_S^2) + MAP_S  + MCMT_S * DD5_S + MCMT_S * MAP_S + MAP_S * DD5_S+ MCMT_S * MCMT_P +DD5_S+ MCMT_P * DD5_S + I(MAP_S^2)+ MAP_S * MCMT_P +I(MCMT_P^2) +MCMT_P+MCMT_S + I(MCMT_S^2) , dat = dat3)
s <- summary(mod)
s$adj.r.squared

for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "DD5_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "DD5_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)
stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,7),subset(stk_BC,15)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P","DD5_S","MAP_S")
pu41 <- predict(stk_BC2,URF_4_1)
pu41[pu41<0] <- 0;pu41[pu41>1000] <- 1000
plot(pu41)
oversxx <- stack(Sxx,Sxx_o,pu41)
oversxx <- data.frame(na.omit(values(oversxx))) 
wh41 <- cor.test(oversxx[,1], oversxx[,3]);wh41$estimate 
ol41 <- cor.test(oversxx[,2], oversxx[,3]);ol41$estimate

## [5]
URF_4_2 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)

s42 <- summary(URF_4_2) 
s42;AIC(URF_4_2)
results_4_2 <- cvFit(URF_4_2, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4_2

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "MAT_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MCMT_P+ I(MCMT_P^2) + MAT_S + I(MAT_S^2)  + MCMT_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "MAT_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)
stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14))
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","MAP_S")
pu42 <- predict(stk_BC2,URF_4_2)
pu42[pu42<0] <- 0;pu42[pu42>1000] <- 1000
plot(pu42)
oversxx <- stack(Sxx,pu42)
oversxx <- data.frame(na.omit(values(oversxx)))
wh42 <- cor.test(oversxx[,1], oversxx[,2]);wh42$estimate 
oversxx_o <- stack(Sxx_o,pu42)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol42 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol42$estimate 

## [6]
URF_4_3 <- lm(HT16~ MAT_P+ I(MAT_P^2) + DD5_S + I(DD5_S^2)  + MAT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MAT_P + MAP_S * DD5_S +
                MAT_S + I(MAT_S^2) + MAT_S * DD5_S + MAT_S * MAT_P + MAT_S * MAP_S, dat = dat3)

s43 <- summary(URF_4_3) 
s43;AIC(URF_4_3)
results_4_3 <- cvFit(URF_4_3, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4_3

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MAT_P+ I(MAT_P^2) + DD5_S + I(DD5_S^2)  + MAT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MAT_P + MAP_S * DD5_S +
                  MAT_S + I(MAT_S^2) + MAT_S * DD5_S + MAT_S * MAT_P + MAT_S * MAP_S, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MAT_S" | names(dat_test) == "MAT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "DD5_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MAT_P+ I(MAT_P^2) + DD5_S + I(DD5_S^2)  + MAT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MAT_P + MAP_S * DD5_S +
                MAT_S + I(MAT_S^2) + MAT_S * DD5_S + MAT_S * MAT_P + MAT_S * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MAT_S" | names(dat_test) == "MAT_P" | names(dat_test) == "MAP_S"| names(dat_test) == "DD5_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)
stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,15)/10,subset(stk_BC,7),subset(stk_BC,14))
names(stk_BC2) <- c("MAT_S","MAT_P","DD5_S","MAP_S")
pu43 <- predict(stk_BC2,URF_4_3)
pu43[pu43<0] <- 0;pu43[pu43>1000] <- 1000
plot(pu43)
oversxx <- stack(Sxx,pu43)
oversxx <- data.frame(na.omit(values(oversxx)))
wh43 <- cor.test(oversxx[,1], oversxx[,2]);wh43$estimate 
oversxx_o <- stack(Sxx_o,pu43)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol43 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol43$estimate 

## [7]
URF_4_4 <- lm(HT16~ MAT_P+ I(MAT_P^2) + MCMT_S + I(MCMT_S^2)  + MAT_P * MCMT_S + MCMT_P + I(MCMT_P^2) + MCMT_P * MAT_P + MCMT_P * MCMT_S +
                MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + MAT_S * MAT_P + MAT_S * MCMT_P, dat = dat3)

s44 <- summary(URF_4_4) 
s44;AIC(URF_4_4)
results_4_4 <- cvFit(URF_4_4, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4_4

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ MAT_P+ I(MAT_P^2) + MCMT_S + I(MCMT_S^2)  + MAT_P * MCMT_S + MCMT_P + I(MCMT_P^2) + MCMT_P * MAT_P + MCMT_P * MCMT_S +
                  MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + MAT_S * MAT_P + MAT_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MAT_S" | names(dat_test) == "MAT_P" | names(dat_test) == "MCMT_P"| names(dat_test) == "MCMT_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ MAT_P+ I(MAT_P^2) + MCMT_S + I(MCMT_S^2)  + MAT_P * MCMT_S + MCMT_P + I(MCMT_P^2) + MCMT_P * MAT_P + MCMT_P * MCMT_S +
                MAT_S + I(MAT_S^2) + MAT_S * MCMT_S + MAT_S * MAT_P + MAT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MAT_S" | names(dat_test) == "MAT_P" | names(dat_test) == "MCMT_P"| names(dat_test) == "MCMT_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,15)/10,subset(stk_BC,15)/10,subset(stk_BC,16)/10,subset(stk_BC,16)/10)
names(stk_BC2) <- c("MAT_S","MAT_P","MCMT_S","MCMT_P")
pu44 <- predict(stk_BC2,URF_4_4)
pu44[pu44<0] <- 0;pu44[pu44>1000] <- 1000
plot(pu44)
oversxx <- stack(Sxx,pu44)
oversxx <- data.frame(na.omit(values(oversxx)))
wh44 <- cor.test(oversxx[,1], oversxx[,2]);wh44$estimate 
oversxx_o <- stack(Sxx_o,pu44)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol44 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol44$estimate 


## [8]
URF_4_5 <- lm(HT16~ TD_P+ I(TD_P^2) + MAP_S + I(MAP_S^2)  + TD_P * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * TD_P + MAT_S * MAP_S +
                TD_S + I(TD_S^2) + TD_S * MAP_S + TD_S * TD_P + TD_S * MAT_S, dat = dat3)

s45 <- summary(URF_4_5) 
s45;AIC(URF_4_5)
results_4_5 <- cvFit(URF_4_5, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_4_5

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16 ~ TD_P+ I(TD_P^2) + MAT_S + I(MAT_S^2)  + TD_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * TD_P + MAP_S * MAT_S +
                  TD_S + I(TD_S^2) + TD_S * MAT_S + TD_S * TD_P + TD_S * MAP_S, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "TD_S" | names(dat_test) == "TD_P" | names(dat_test) == "MAP_S"| names(dat_test) == "MAT_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16 ~ TD_P+ I(TD_P^2) + MAT_S + I(MAT_S^2)  + TD_P * MAT_S + MAP_S + I(MAP_S^2) + MAP_S * TD_P + MAP_S * MAT_S +
                TD_S + I(TD_S^2) + TD_S * MAT_S + TD_S * TD_P + TD_S * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "TD_S" | names(dat_test) == "TD_P" | names(dat_test) == "MAP_S"| names(dat_test) == "MAT_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,23)/10,subset(stk_BC,23)/10,subset(stk_BC,15)/10,subset(stk_BC,14))
names(stk_BC2) <- c("TD_S","TD_P","MAT_S","MAP_S")
pu45 <- predict(stk_BC2,URF_4_5)
pu45[pu45<0] <- 0;pu45[pu45>1000] <- 1000
plot(pu45)

oversxx <- stack(Sxx,pu45)
oversxx <- data.frame(na.omit(values(oversxx)))
wh45 <- cor.test(oversxx[,1], oversxx[,2]);wh45$estimate 
oversxx_o <- stack(Sxx_o,pu45)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol45 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol45$estimate 

## [9]
URF_5_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat3)
s51 <- summary(URF_5_1)
s51;AIC(URF_5_1)
results_5_1 <- cvFit(URF_5_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_5_1

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                  MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                              names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                            names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,7))
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","MAP_S","DD5_S")
pu51 <- predict(stk_BC2,URF_5_1)
pu51[pu51<0] <- 0;pu51[pu51>1000] <- 1000
plot(pu51)

oversxx <- stack(Sxx,pu51)
oversxx <- data.frame(na.omit(values(oversxx)))
wh51 <- cor.test(oversxx[,1], oversxx[,2]);wh51$estimate 
oversxx_o <- stack(Sxx_o,pu51)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol51 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol51$estimate 

## [10]
URF_5_2 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat3)
s52 <- summary(URF_5_2)
s52;AIC(URF_5_2)
results_5_2 <- cvFit(URF_5_2, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_5_2

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                  MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                              names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                            names(dat_test) == "MAT_S"| names(dat_test) == "MAT_P") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,15)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","MAP_S","MAT_P")
pu52 <- predict(stk_BC2,URF_5_2)
pu52[pu52<0] <- 0;pu52[pu52>1000] <- 1000
plot(pu52)

oversxx <- stack(Sxx,pu52)
oversxx <- data.frame(na.omit(values(oversxx)))
wh52 <- cor.test(oversxx[,1], oversxx[,2]);wh52$estimate 
oversxx_o <- stack(Sxx_o,pu52)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol52 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol52$estimate 

## [11]
URF_6_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P + MAT_P + I(MAT_P^2) + MAT_P * MCMT_P + MAT_P * MCMT_S + MAT_P * MAT_S +
                MAT_P * DD5_S + MAT_P * MAP_S, dat = dat3)
s61 <- summary(URF_6_1)
s61;AIC(URF_6_1)
results_6_1 <- cvFit(URF_6_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_6_1

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                  MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P + MAT_P + I(MAT_P^2) + MAT_P * MCMT_P + MAT_P * MCMT_S + MAT_P * MAT_S +
                  MAT_P * DD5_S + MAT_P * MAP_S, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                              names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S" | names(dat_test) == "MAT_P") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P + MAT_P + I(MAT_P^2) + MAT_P * MCMT_P + MAT_P * MCMT_S + MAT_P * MAT_S +
                MAT_P * DD5_S + MAT_P * MAP_S, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                            names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S" | names(dat_test) == "MAT_P") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,7),subset(stk_BC,15)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","MAP_S","DD5_S","MAT_P")
pu61 <- predict(stk_BC2,URF_6_1)
pu61[pu61<0] <- 0;pu61[pu61>1000] <- 1000
plot(pu61)

oversxx <- stack(Sxx,pu61)
oversxx <- data.frame(na.omit(values(oversxx)))
wh61 <- cor.test(oversxx[,1], oversxx[,2]);wh61$estimate 
oversxx_o <- stack(Sxx_o,pu61)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol61 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol61$estimate 


## [12]
URF_7_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P + MAP_P + I(MAP_P^2) + MAP_P * MCMT_P + MAP_P * MCMT_S + MAP_P * MAT_S +
                MAP_P * MAT_P + MAP_P * MAP_S + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MCMT_S + PAS_S * MAP_P + PAS_S * MAT_P +
                PAS_S * MAT_S + PAS_S * MAP_S + PAS_S * MAP_P, dat = dat3)
s71 <- summary(URF_7_1)
s71;AIC(URF_7_1)
results_7_1 <- cvFit(URF_7_1, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_7_1

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_P +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                  MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P + MAP_P + I(MAP_P^2) + MAP_P * MCMT_P + MAP_P * MCMT_S + MAP_P * MAT_S +
                  MAP_P * MAT_P + MAP_P * MAP_S + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MCMT_S + PAS_S * MAP_P + PAS_S * MAT_P +
                  PAS_S * MAT_S + PAS_S * MAP_S + PAS_S * MAP_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                              names(dat_test) == "MAT_S"| names(dat_test) == "PAS_S" | names(dat_test) == "MAT_P"|
                              names(dat_test) == "MAP_P") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * MAP_S + MAT_S + I(MAT_S^2) + MAT_S * MAP_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P + MAP_P + I(MAP_P^2) + MAP_P * MCMT_P + MAP_P * MCMT_S + MAP_P * MAT_S +
                MAP_P * MAT_P + MAP_P * MAP_S + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MCMT_S + PAS_S * MAP_P + PAS_S * MAT_P +
                PAS_S * MAT_S + PAS_S * MAP_S + PAS_S * MAP_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "MAP_S"| 
                            names(dat_test) == "MAT_S"| names(dat_test) == "PAS_S" | names(dat_test) == "MAT_P"|
                            names(dat_test) == "MAP_P") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14),
                 subset(stk_BC,20),subset(stk_BC,15)/10,subset(stk_BC,14))
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","MAP_S","PAS_S","MAT_P","MAP_P")
pu71 <- predict(stk_BC2,URF_7_1)
pu71[pu71<0] <- 0;pu71[pu71>1000] <- 1000
plot(pu71)

oversxx <- stack(Sxx,pu71)
oversxx <- data.frame(na.omit(values(oversxx)))
wh71 <- cor.test(oversxx[,1], oversxx[,2]);wh71$estimate 
oversxx_o <- stack(Sxx_o,pu71)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol71 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol71$estimate 

# model [7] excluded, substitude model:

URF_5_3 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * PAS_S + MAT_S + I(MAT_S^2) + MAT_S * PAS_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat3)
s53 <- summary(URF_5_3)
s53;AIC(URF_5_3)
results_5_3 <- cvFit(URF_5_3, data = dat3 , y = dat3$HT16, K = 10, R = 10);results_5_3

RMSE_1 <- RMSE_re <- R_sq_1 <- c(1)
for(i in 1:10){
  k <- 10
  long_vector <- 1:2114
  folds <- split(long_vector, sample(rep(1:k, length.out = length(long_vector))))
  for (k in 1:10){#k=1
    dat_test <- dat3[folds[[k]],]
    dat_train <- dat3[!rownames(dat3) %in% rownames(dat_test),]
    
    model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * DD5_S +
                  MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * PAS_S + MAT_S + I(MAT_S^2) + MAT_S * PAS_S +
                  MAT_S * DD5_S + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
    
    test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "PAS_S"| 
                              names(dat_test) == "MAT_S"| names(dat_test) == "DD5_S") ]
    test_result <- stats::predict(model,test)
    
    a <- lm(dat_test$HT16 ~ test_result)
    b <- summary(a)
    R_sq_1[k] <-  b$adj.r.squared
    RMSE_1[k] <- b$sigma
  }
  mean(RMSE_1)
  RMSE_re[i] <- var(RMSE_1)
}
mean(RMSE_1)
mean(RMSE_re)

R_sq <- RMSE <- c(1)
for (k in 1:length(unique(dat3$Site))){ #k=1
  dat_train <- dat3[which(dat3$Site != unique(dat3$Site)[k]),]
  dat_test <- dat3[which(dat3$Site == unique(dat3$Site)[k]),]
  
  model <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + MAT_P + I(MAT_P^2)  + MCMT_P * MAT_P + PAS_S + I(PAS_S^2) + PAS_S * MCMT_P + PAS_S * MAT_P +
                MCMT_S + I(MCMT_S^2) + MCMT_S * MAT_P + MCMT_S * MCMT_P + MCMT_S * PAS_S + MAT_S + I(MAT_S^2) + MAT_S * PAS_S +
                MAT_S * MAT_P + MAT_S * MCMT_S + MAT_S * MCMT_P, dat = dat_train)
  
  test <- dat_test[,which(names(dat_test) == "MCMT_S" | names(dat_test) == "MCMT_P" | names(dat_test) == "PAS_S"| 
                            names(dat_test) == "MAT_S"| names(dat_test) == "MAT_P") ]
  test_result <- stats::predict(model,test)
  
  a <- lm(dat_test$HT16 ~ test_result)
  b <- summary(a)
  
  RMSE[k] <- sqrt(mean((test_result-dat_test$HT16)^2))
}
mean(R_sq)
mean(RMSE)
var(RMSE)

stk_BC2 <- stack(subset(stk_BC,16)/10,subset(stk_BC,16)/10,subset(stk_BC,15)/10,subset(stk_BC,14),subset(stk_BC,15)/10)
names(stk_BC2) <- c("MCMT_S","MCMT_P","MAT_S","PAS_S","MAT_P")
pu53 <- predict(stk_BC2,URF_5_3)
pu53[pu53<0] <- 0;pu53[pu53>1000] <- 1000
plot(pu53)

oversxx <- stack(Sxx,pu53)
oversxx <- data.frame(na.omit(values(oversxx)))
wh53 <- cor.test(oversxx[,1], oversxx[,2]);wh53$estimate 
oversxx_o <- stack(Sxx_o,pu53)
oversxx_o <- data.frame(na.omit(values(oversxx_o)))
ol53 <- cor.test(oversxx_o[,1], oversxx_o[,2]);ol53$estimate 


# Figure 2 Local_seedlot_growth_comparison.png=================

png(filename = "Figure_2_Local_seedlot_growth_comparison.png",height = 1750, width = 1700, res = 170)
par(mfrow=c(2,2),cex.lab=1.25,mai=c(0.55,0.55,0.35,0.1),oma=c(0,0,4,2))
plot(pu41,main="URF model [4]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                  font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pu42,main="URF model [5]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pu43,main="URF model [6]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pu45,main="URF model [8]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                   font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("   Local seedlot growth potential", side =3, line = 1, cex = 2, outer=TRUE)
dev.off()

Sxx_p <- Sxx/max(values(Sxx),na.rm = T)
png(filename = "Figure_2_Local_seedlot_growth_comparison.png",height = 1750, width = 1750, res = 170)
par(mfrow=c(2,2),cex.lab=1.25,mai=c(0.55,0.55,0.35,0.1),oma=c(0,0,4,2))
plot(Sxx_p,main="Interior spruce site productivity",cex.main = 2,legend.args = list(text = 'Site index', side = 4, 
                                                  font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(put41,main="URF model [4]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                  font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(put42,main="URF model [5]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                  font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)

plot(put45,main="URF model [8]",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                  font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("   Local seedlot growth potential", side =3, line = 1, cex = 2, outer=TRUE)
dev.off()

par(mfrow=c(1,1))

# Table 2 model [4] summary===================
URF_4_1 <- lm(HT16~ MCMT_P+ I(MCMT_P^2) + DD5_S + I(DD5_S^2)  + MCMT_P * DD5_S + MAP_S + I(MAP_S^2) + MAP_S * MCMT_P + MAP_S * DD5_S +
                MCMT_S + I(MCMT_S^2) + MCMT_S * DD5_S + MCMT_S * MCMT_P + MCMT_S * MAP_S, dat = dat3)

summary(URF_4_1) 
84.37/mean(dat3$HT16)

# Figure 3 Optimal seed source 1960-1961===============================
#first-order derivative of MCMT_P
b1=-1.270e+01;b2=-5.091e-01;b5=-2.691e-03;b9=4.875e-03;b12=4.264e-01
stk_BC3 <- stack(subset(stk_BC,15),subset(stk_BC,7),subset(stk_BC,16)/10,-(b1+b5*subset(stk_BC,7)+b9*subset(stk_BC,15)+b12*subset(stk_BC,16)/10)*0.5/b2)
names(stk_BC3) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pp41 <- predict(stk_BC3,URF_4_1)
pp41[pp41<0] <- 0;pp41[pp41>1000] <- 1000
plot(pp41,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                     font = 2, line = 2.5, cex = 0.8))
plot(pu41)
pd41 <- pp41-pu41
plot(pd41)
pd41 <- mask(pd41,Sxx_o)
pd41[is.na(pd41)] <- 0
pd41 <- mask(pd41,Sxx)
plot(pd41,main="Growth potential gain",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                          font = 2, line = 2.5, cex = 0.8))

pc41 <- -(b1+b5*subset(stk_BC,7)+b9*subset(stk_BC,15)+b12*subset(stk_BC,16)/10)*0.5/b2
pcd41 <- (pc41 - subset(stk_BC,16)/10)/10
pcc41 <- subset(stk_BC,16)/10+pcd41
plot(pcc41,main="Optimal seed source",legend.args = list(text = 'Mean Coldest Month Temperature (°C)', side = 4, 
                                                       font = 2, line = 2.5, cex = 0.8))
plot(pcd41,main="Difference between local and optimal",legend.args = list(text = 'Mean Coldest Month Temperature (°C)', side = 4, 
                                                                          font = 2, line = 2.5, cex = 0.8))

png(filename = "Figure_3_current_optimal_seed_source.png",height = 1750, width = 1700, res = 170)
par(mfrow=c(2,2),cex.lab=1.25,mai=c(0.55,0.55,0.35,0.1),oma=c(0,0,4,2))
my_palette <- colorRampPalette(c("deepskyblue2", "yellow", "coral")) # Example palette
col_1 <- my_palette(100)
plot(pcc41,main="Optimal seed source MCMT",col = col_1, legend.args = list(text = 'Mean Coldest Month Temperature (°C)', side = 4, 
                                                                           font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(ppt41,main="Optimal seedlot growth potential",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                     font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)

plot(pd41,main="Growth potential gain using optimal seedlot",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                          font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
my_palette <- colorRampPalette(c("deepskyblue2", "yellow")) # Example palette
col_6 <- my_palette(100)
plot(pcd41,main="Difference between local and optimal",col = col_6,legend.args = list(text = 'Mean Coldest Month Temperature (°C)', side = 4, 
                                                                          font = 2, line = 2.5, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
dev.off()

# Figure 4 and Appendix Optimal seed source , SSP 126, SSP245, SSP 370, SSP 585; ==============

# pf (local future), pg (change relatively to current), pfo (using optimal seed source), pfg(potential gain) 
## ssp126 - 4
stk_BC_F <- stack(subset(BC_stack_F,17),subset(BC_stack_F,1),subset(BC_stack_F,9)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_126_4 <- predict(stk_BC_F,URF_4_1)
pf_126_4[pf_126_4<0] <- 0;pf_126_4[pf_126_4>1000] <- 1000
plot(pf_126_4,main="Local seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                              font = 2, line = 2.5, cex = 0.8))
pg_126_4 <- pf_126_4-pu41
min <- min(values(pg_126_4),na.rm = T)
pg_126_4 <- mask(pg_126_4,Sxx_o)
pg_126_4 <- pg_126_4 - min(values(pg_126_4),na.rm = T)
pg_126_4[is.na(pg_126_4)] <- 0
pg_126_4 <- mask(pg_126_4,Sxx)
pg_126_4 <- pg_126_4 + min
plot(pg_126_4,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,17),subset(BC_stack_F,1),subset(BC_stack_F,9)/10,-(b1+b5*subset(BC_stack_F,1)+b9*subset(BC_stack_F,17)+b12*subset(BC_stack_F,9)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_126_4 <- predict(stk_BC_F,URF_4_1)
pfo_126_4[pfo_126_4<0] <- 0;pfo_126_4[pfo_126_4>1000] <- 1000
plot(pfo_126_4,main="Optimal seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                     font = 2, line = 2.5, cex = 0.8))
pfg_126_4 <- pfo_126_4 - pf_126_4
pfg_126_4 <- mask(pfg_126_4,Sxx_o)
pfg_126_4[is.na(pfg_126_4)] <- 0
pfg_126_4 <- mask(pfg_126_4,Sxx)
plot(pfg_126_4,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                  font = 2, line = 2.5, cex = 0.8))

## ssp126 - 7
stk_BC_F <- stack(subset(BC_stack_F,18),subset(BC_stack_F,2),subset(BC_stack_F,10)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_126_7 <- predict(stk_BC_F,URF_4_1)
pf_126_7[pf_126_7<0] <- 0;pf_126_7[pf_126_7>1000] <- 1000
plot(pf_126_7,main="Local seedlot growth potential in SSP1-2.6 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_126_7 <- pf_126_7-pu41
min <- min(values(pg_126_7),na.rm = T)
pg_126_7 <- mask(pg_126_7,Sxx_o)
pg_126_7 <- pg_126_7 - min(values(pg_126_7),na.rm = T)
pg_126_7[is.na(pg_126_7)] <- 0
pg_126_7 <- mask(pg_126_7,Sxx)
pg_126_7 <- pg_126_7 + min
plot(pg_126_7,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,18),subset(BC_stack_F,2),subset(BC_stack_F,10)/10,-(b1+b5*subset(BC_stack_F,2)+b9*subset(BC_stack_F,18)+b12*subset(BC_stack_F,10)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_126_7 <- predict(stk_BC_F,URF_4_1)
pfo_126_7[pfo_126_7<0] <- 0;pfo_126_7[pfo_126_7>1000] <- 1000
plot(pfo_126_7,main="Optimal seedlot growth potential in SSP1-2.6 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_126_7 <- pfo_126_7 - pf_126_7
pfg_126_7 <- mask(pfg_126_7,Sxx_o)
pfg_126_7[is.na(pfg_126_7)] <- 0
pfg_126_7 <- mask(pfg_126_7,Sxx)
plot(pfg_126_7,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2071-2100",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))
## ssp245 - 4
stk_BC_F <- stack(subset(BC_stack_F,19),subset(BC_stack_F,3),subset(BC_stack_F,11)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_245_4 <- predict(stk_BC_F,URF_4_1)
pf_245_4[pf_245_4<0] <- 0;pf_245_4[pf_245_4>1000] <- 1000
plot(pf_245_4,main="Local seedlot growth potential in SSP2-4.5 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_245_4 <- pf_245_4-pu41
min <- min(values(pg_245_4),na.rm = T)
pg_245_4 <- mask(pg_245_4,Sxx_o)
pg_245_4 <- pg_245_4 - min(values(pg_245_4),na.rm = T)
pg_245_4[is.na(pg_245_4)] <- 0
pg_245_4 <- mask(pg_245_4,Sxx)
pg_245_4 <- pg_245_4 + min
plot(pg_245_4,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,19),subset(BC_stack_F,3),subset(BC_stack_F,11)/10,-(b1+b5*subset(BC_stack_F,3)+b9*subset(BC_stack_F,19)+b12*subset(BC_stack_F,11)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_245_4 <- predict(stk_BC_F,URF_4_1)
pfo_245_4[pfo_245_4<0] <- 0;pfo_245_4[pfo_245_4>1000] <- 1000
plot(pfo_245_4,main="Optimal seedlot growth potential in SSP2-4.5 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_245_4 <- pfo_245_4 - pf_245_4
pfg_245_4 <- mask(pfg_245_4,Sxx_o)
pfg_245_4[is.na(pfg_245_4)] <- 0
pfg_245_4 <- mask(pfg_245_4,Sxx)
plot(pfg_245_4,main="Growth potential gain using optimal seedlot in SSP2-4.5 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

## ssp245 - 7
stk_BC_F <- stack(subset(BC_stack_F,20),subset(BC_stack_F,4),subset(BC_stack_F,12)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_245_7 <- predict(stk_BC_F,URF_4_1)
pf_245_7[pf_245_7<0] <- 0;pf_245_7[pf_245_7>1000] <- 1000
plot(pf_245_7,main="Local seedlot growth potential in SSP2-4.5 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_245_7 <- pf_245_7-pu41
min <- min(values(pg_245_7),na.rm = T)
pg_245_7 <- mask(pg_245_7,Sxx_o)
pg_245_7 <- pg_245_7 - min(values(pg_245_7),na.rm = T)
pg_245_7[is.na(pg_245_7)] <- 0
pg_245_7 <- mask(pg_245_7,Sxx)
pg_245_7 <- pg_245_7 + min
plot(pg_245_7,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,20),subset(BC_stack_F,4),subset(BC_stack_F,12)/10,-(b1+b5*subset(BC_stack_F,4)+b9*subset(BC_stack_F,20)+b12*subset(BC_stack_F,12)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_245_7 <- predict(stk_BC_F,URF_4_1)
pfo_245_7[pfo_245_7<0] <- 0;pfo_245_7[pfo_245_7>1000] <- 1000
plot(pfo_245_7,main="Optimal seedlot growth potential in SSP2-4.5 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_245_7 <- pfo_245_7 - pf_245_7
pfg_245_7 <- mask(pfg_245_7,Sxx_o)
pfg_245_7[is.na(pfg_245_7)] <- 0
pfg_245_7 <- mask(pfg_245_7,Sxx)
plot(pfg_245_7,main="Growth potential gain using optimal seedlot in SSP2-4.5 in 2071-2100",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

## ssp370 - 4
stk_BC_F <- stack(subset(BC_stack_F,21),subset(BC_stack_F,5),subset(BC_stack_F,13)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_370_4 <- predict(stk_BC_F,URF_4_1)
pf_370_4[pf_370_4<0] <- 0;pf_370_4[pf_370_4>1000] <- 1000
plot(pf_370_4,main="Local seedlot growth potential in SSP3-7.0 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_370_4 <- pf_370_4-pu41
min <- min(values(pg_370_4),na.rm = T)
pg_370_4 <- mask(pg_370_4,Sxx_o)
pg_370_4 <- pg_370_4 - min(values(pg_370_4),na.rm = T)
pg_370_4[is.na(pg_370_4)] <- 0
pg_370_4 <- mask(pg_370_4,Sxx)
pg_370_4 <- pg_370_4 + min
plot(pg_370_4,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,21),subset(BC_stack_F,5),subset(BC_stack_F,13)/10,-(b1+b5*subset(BC_stack_F,5)+b9*subset(BC_stack_F,21)+b12*subset(BC_stack_F,13)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_370_4 <- predict(stk_BC_F,URF_4_1)
pfo_370_4[pfo_370_4<0] <- 0;pfo_370_4[pfo_370_4>1000] <- 1000
plot(pfo_370_4,main="Optimal seedlot growth potential in SSP3-7.0 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_370_4 <- pfo_370_4 - pf_370_4
pfg_370_4 <- mask(pfg_370_4,Sxx_o)
pfg_370_4[is.na(pfg_370_4)] <- 0
pfg_370_4 <- mask(pfg_370_4,Sxx)
plot(pfg_370_4,main="Growth potential gain using optimal seedlot in SSP3-7.0 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

## ssp370 - 7
stk_BC_F <- stack(subset(BC_stack_F,22),subset(BC_stack_F,6),subset(BC_stack_F,14)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_370_7 <- predict(stk_BC_F,URF_4_1)
pf_370_7[pf_370_7<0] <- 0;pf_370_7[pf_370_7>1000] <- 1000
plot(pf_370_7,main="Local seedlot growth potential in SSP3-7.0 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_370_7 <- pf_370_7-pu41
min <- min(values(pg_370_7),na.rm = T)
pg_370_7 <- mask(pg_370_7,Sxx_o)
pg_370_7 <- pg_370_7 - min(values(pg_370_7),na.rm = T)
pg_370_7[is.na(pg_370_7)] <- 0
pg_370_7 <- mask(pg_370_7,Sxx)
pg_370_7 <- pg_370_7 + min
plot(pg_370_7,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,22),subset(BC_stack_F,6),subset(BC_stack_F,14)/10,-(b1+b5*subset(BC_stack_F,6)+b9*subset(BC_stack_F,22)+b12*subset(BC_stack_F,14)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_370_7 <- predict(stk_BC_F,URF_4_1)
pfo_370_7[pfo_370_7<0] <- 0;pfo_370_7[pfo_370_7>1000] <- 1000
plot(pfo_370_7,main="Optimal seedlot growth potential in SSP3-7.0 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_370_7 <- pfo_370_7 - pf_370_7
pfg_370_7 <- mask(pfg_370_7,Sxx_o)
pfg_370_7[is.na(pfg_370_7)] <- 0
pfg_370_7 <- mask(pfg_370_7,Sxx)
plot(pfg_370_7,main="Growth potential gain using optimal seedlot in SSP3-7.0 in 2071-2100",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

## ssp585 - 4
stk_BC_F <- stack(subset(BC_stack_F,23),subset(BC_stack_F,7),subset(BC_stack_F,15)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_585_4 <- predict(stk_BC_F,URF_4_1)
pf_585_4[pf_585_4<0] <- 0;pf_585_4[pf_585_4>1000] <- 1000
plot(pf_585_4,main="Local seedlot growth potential in SSP5-8.5 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_585_4 <- pf_585_4-pu41
min <- min(values(pg_585_4),na.rm = T)
pg_585_4 <- mask(pg_585_4,Sxx_o)
pg_585_4 <- pg_585_4 - min(values(pg_585_4),na.rm = T)
pg_585_4[is.na(pg_585_4)] <- 0
pg_585_4 <- mask(pg_585_4,Sxx)
pg_585_4 <- pg_585_4 + min
plot(pg_585_4,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,23),subset(BC_stack_F,7),subset(BC_stack_F,15)/10,-(b1+b5*subset(BC_stack_F,7)+b9*subset(BC_stack_F,23)+b12*subset(BC_stack_F,15)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_585_4 <- predict(stk_BC_F,URF_4_1)
pfo_585_4[pfo_585_4<0] <- 0;pfo_585_4[pfo_585_4>1000] <- 1000
plot(pfo_585_4,main="Optimal seedlot growth potential in SSP5-8.5 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_585_4 <- pfo_585_4 - pf_585_4
pfg_585_4 <- mask(pfg_585_4,Sxx_o)
pfg_585_4[is.na(pfg_585_4)] <- 0
pfg_585_4 <- mask(pfg_585_4,Sxx)
plot(pfg_585_4,main="Growth potential gain using optimal seedlot in SSP5-8.5 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

## ssp585 - 7
stk_BC_F <- stack(subset(BC_stack_F,24),subset(BC_stack_F,8),subset(BC_stack_F,16)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_585_7 <- predict(stk_BC_F,URF_4_1)
pf_585_7[pf_585_7<0] <- 0;pf_585_7[pf_585_7>1000] <- 1000
plot(pf_585_7,main="Local seedlot growth potential in SSP5-8.5 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_585_7 <- pf_585_7-pu41
min <- min(values(pg_585_7),na.rm = T)
pg_585_7 <- mask(pg_585_7,Sxx_o)
pg_585_7 <- pg_585_7 - min(values(pg_585_7),na.rm = T)
pg_585_7[is.na(pg_585_7)] <- 0
pg_585_7 <- mask(pg_585_7,Sxx)
pg_585_7 <- pg_585_7 + min
plot(pg_585_7,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                            font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,24),subset(BC_stack_F,8),subset(BC_stack_F,16)/10,-(b1+b5*subset(BC_stack_F,8)+b9*subset(BC_stack_F,24)+b12*subset(BC_stack_F,16)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_585_7 <- predict(stk_BC_F,URF_4_1)
pfo_585_7[pfo_585_7<0] <- 0;pfo_585_7[pfo_585_7>1000] <- 1000
plot(pfo_585_7,main="Optimal seedlot growth potential in SSP5-8.5 in 2071-2100",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_585_7 <- pfo_585_7 - pf_585_7
pfg_585_7 <- mask(pfg_585_7,Sxx_o)
pfg_585_7[is.na(pfg_585_7)] <- 0
pfg_585_7 <- mask(pfg_585_7,Sxx)
plot(pfg_585_7,main="Growth potential gain using optimal seedlot in SSP5-8.5 in 2071-2100",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

stk_BC_F <- stack(subset(BC_stack_F,35),subset(BC_stack_F,33),subset(BC_stack_F,34)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_126_1 <- predict(stk_BC_F,URF_4_1)
pf_126_1[pf_126_1<0] <- 0;pf_126_1[pf_126_1>1000] <- 1000
plot(pf_126_1,main="Local seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_126_1 <- pf_126_1-pu41
min <- min(values(pg_126_1),na.rm = T)
pg_126_1 <- mask(pg_126_1,Sxx_o)
pg_126_1 <- pg_126_1 - min(values(pg_126_1),na.rm = T)
pg_126_1[is.na(pg_126_1)] <- 0
pg_126_1 <- mask(pg_126_1,Sxx)
pg_126_1 <- pg_126_1 + min
plot(pg_126_1,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                      font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,35),subset(BC_stack_F,33),subset(BC_stack_F,34)/10,-(b1+b5*subset(BC_stack_F,33)+b9*subset(BC_stack_F,35)+b12*subset(BC_stack_F,34)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_126_1 <- predict(stk_BC_F,URF_4_1)
pfo_126_1[pfo_126_1<0] <- 0;pfo_126_1[pfo_126_1>1000] <- 1000
plot(pfo_126_1,main="Optimal seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_126_1 <- pfo_126_1 - pf_126_1
pfg_126_1 <- mask(pfg_126_1,Sxx_o)
pfg_126_1[is.na(pfg_126_1)] <- 0
pfg_126_1 <- mask(pfg_126_1,Sxx)
plot(pfg_126_1,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))


stk_BC_F <- stack(subset(BC_stack_F,39),subset(BC_stack_F,37),subset(BC_stack_F,38)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_245_1 <- predict(stk_BC_F,URF_4_1)
pf_245_1[pf_245_1<0] <- 0;pf_245_1[pf_245_1>1000] <- 1000
plot(pf_245_1,main="Local seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_245_1 <- pf_245_1-pu41
min <- min(values(pg_245_1),na.rm = T)
pg_245_1 <- mask(pg_245_1,Sxx_o)
pg_245_1 <- pg_245_1 - min(values(pg_245_1),na.rm = T)
pg_245_1[is.na(pg_245_1)] <- 0
pg_245_1 <- mask(pg_245_1,Sxx)
pg_245_1 <- pg_245_1 + min
plot(pg_245_1,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                      font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,39),subset(BC_stack_F,37),subset(BC_stack_F,38)/10,-(b1+b5*subset(BC_stack_F,37)+b9*subset(BC_stack_F,39)+b12*subset(BC_stack_F,38)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_245_1 <- predict(stk_BC_F,URF_4_1)
pfo_245_1[pfo_245_1<0] <- 0;pfo_245_1[pfo_245_1>1000] <- 1000
plot(pfo_245_1,main="Optimal seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_245_1 <- pfo_245_1 - pf_245_1
pfg_245_1 <- mask(pfg_245_1,Sxx_o)
pfg_245_1[is.na(pfg_245_1)] <- 0
pfg_245_1 <- mask(pfg_245_1,Sxx)
plot(pfg_245_1,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))


stk_BC_F <- stack(subset(BC_stack_F,43),subset(BC_stack_F,41),subset(BC_stack_F,42)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_370_1 <- predict(stk_BC_F,URF_4_1)
pf_370_1[pf_370_1<0] <- 0;pf_370_1[pf_370_1>1000] <- 1000
plot(pf_370_1,main="Local seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_370_1 <- pf_370_1-pu41
min <- min(values(pg_370_1),na.rm = T)
pg_370_1 <- mask(pg_370_1,Sxx_o)
pg_370_1 <- pg_370_1 - min(values(pg_370_1),na.rm = T)
pg_370_1[is.na(pg_370_1)] <- 0
pg_370_1 <- mask(pg_370_1,Sxx)
pg_370_1 <- pg_370_1 + min
plot(pg_370_1,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                      font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,43),subset(BC_stack_F,41),subset(BC_stack_F,42)/10,-(b1+b5*subset(BC_stack_F,41)+b9*subset(BC_stack_F,43)+b12*subset(BC_stack_F,42)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_370_1 <- predict(stk_BC_F,URF_4_1)
pfo_370_1[pfo_370_1<0] <- 0;pfo_370_1[pfo_370_1>1000] <- 1000
plot(pfo_370_1,main="Optimal seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_370_1 <- pfo_370_1 - pf_370_1
pfg_370_1 <- mask(pfg_370_1,Sxx_o)
pfg_370_1[is.na(pfg_370_1)] <- 0
pfg_370_1 <- mask(pfg_370_1,Sxx)
plot(pfg_370_1,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))


stk_BC_F <- stack(subset(BC_stack_F,47),subset(BC_stack_F,45),subset(BC_stack_F,46)/10,subset(stk_BC,16)/10)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pf_585_1 <- predict(stk_BC_F,URF_4_1)
pf_585_1[pf_585_1<0] <- 0;pf_585_1[pf_585_1>1000] <- 1000
plot(pf_585_1,main="Local seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                font = 2, line = 2.5, cex = 0.8))
pg_585_1 <- pf_585_1-pu41
min <- min(values(pg_585_1),na.rm = T)
pg_585_1 <- mask(pg_585_1,Sxx_o)
pg_585_1 <- pg_585_1 - min(values(pg_585_1),na.rm = T)
pg_585_1[is.na(pg_585_1)] <- 0
pg_585_1 <- mask(pg_585_1,Sxx)
pg_585_1 <- pg_585_1 + min
plot(pg_585_1,main="Height increase/decrease from climate change", legend.args = list(text = '16-year tree height gain(cm)', side = 4, 
                                                                                      font = 2, line = 2.5, cex = 0.8))
stk_BC_F <- stack(subset(BC_stack_F,47),subset(BC_stack_F,45),subset(BC_stack_F,46)/10,-(b1+b5*subset(BC_stack_F,45)+b9*subset(BC_stack_F,47)+b12*subset(BC_stack_F,46)/10)*0.5/b2)
names(stk_BC_F) <- c("MAP_S","DD5_S","MCMT_S","MCMT_P")
pfo_585_1 <- predict(stk_BC_F,URF_4_1)
pfo_585_1[pfo_585_1<0] <- 0;pfo_585_1[pfo_585_1>1000] <- 1000
plot(pfo_585_1,main="Optimal seedlot growth potential in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                                   font = 2, line = 2.5, cex = 0.8))
pfg_585_1 <- pfo_585_1 - pf_585_1
pfg_585_1 <- mask(pfg_585_1,Sxx_o)
pfg_585_1[is.na(pfg_585_1)] <- 0
pfg_585_1 <- mask(pfg_585_1,Sxx)
plot(pfg_585_1,main="Growth potential gain using optimal seedlot in SSP1-2.6 in 2041-2070",legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                              font = 2, line = 2.5, cex = 0.8))

col_2 <- c("#F2F2F2", colorRampPalette(c("coral","yellow"))(70),colorRampPalette(c("yellow","deepskyblue1"))(27))
col_3 <- c("#F2F2F2", colorRampPalette(c("coral","yellow"))(70),colorRampPalette(c("yellow","deepskyblue1"))(35))
col_4 <- c("#F2F2F2", colorRampPalette(c("coral","yellow"))(70),colorRampPalette(c("yellow","deepskyblue1"))(43))
col_5 <- c("#F2F2F2", colorRampPalette(c("coral","yellow"))(70),colorRampPalette(c("yellow","deepskyblue1"))(51))


png(filename = "Figure_4_future_scenarios.png",height = 1750, width = 3400, res = 180)
par(mfrow=c(2,4), cex.lab=1.25,mai=c(0.80,0.80,0.60,0.37),oma=c(0,0,0,2.4))
plot(pft_126_4,main="Local seedlot growth potential",cex.main = 1.8, legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                        font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_126_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_126_4,col=col_3,main="Height increase/decrease",cex.main = 1.8, legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                            font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_126_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                     font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP1-2.6 in 2041-2070", side =2, line = -2, adj = 0.85,cex = 1.7, outer=TRUE)
plot(pft_585_4,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                        font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[5], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_585_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[6], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_585_4,col=col_5,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                            font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[7], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_585_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                     font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[8], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP5-8.25 in 2041-2070", side =2, line = -2, adj = 0.15,cex = 1.7, outer=TRUE)
dev.off()


png(filename = "Appendix_Figure_1_future_scenarios_SSP126.png",height = 2650, width = 3400, res = 180)
par(mfrow=c(3,4), cex.lab=1.25,mai=c(0.80,0.80,0.60,0.37),oma=c(0,0,0,2.4))
plot(pft_126_1,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                        font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_126_1,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_126_1,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                            font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_126_1,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                     font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP1-2.6 in 2011-2040", side =2, line = -2, adj = 0.94,cex = 1.7, outer=TRUE)
plot(pft_126_4,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                        font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[5], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_126_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[6], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_126_4,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                            font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[7], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_126_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                     font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[8], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP1-2.6 in 2041-2070", side =2, line = -2, adj = 0.5,cex = 1.7, outer=TRUE)
plot(pft_126_7,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                        font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[9], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_126_7,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[10], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_126_7,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                            font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[11], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_126_7,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                     font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[12], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP1-2.6 in 2071-2100", side =2, line = -2, adj = 0.08,cex = 1.7, outer=TRUE)
dev.off()

png(filename = "Appendix_Figure_2_future_scenarios_SSP245.png",height = 2650, width = 3400, res = 180)
par(mfrow=c(3,4), cex.lab=1.25,mai=c(0.80,0.80,0.60,0.37),oma=c(0,0,0,2.4))
plot(pft_245_1,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_245_1,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_245_1,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_245_1,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP2-4.5 in 2011-2040", side =2, line = -2, adj = 0.94,cex = 1.7, outer=TRUE)
plot(pft_245_4,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[5], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_245_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[6], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_245_4,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[7], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_245_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[8], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP2-4.5 in 2041-2070", side =2, line = -2, adj = 0.5,cex = 1.7, outer=TRUE)
plot(pft_245_7,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[9], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_245_7,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[10], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_245_7,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[11], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_245_7,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[12], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP2-4.5 in 2071-2100", side =2, line = -2, adj = 0.08,cex = 1.7, outer=TRUE)
dev.off()

png(filename = "Appendix_Figure_3_future_scenarios_SSP370.png",height = 2650, width = 3400, res = 180)
par(mfrow=c(3,4), cex.lab=1.25,mai=c(0.80,0.80,0.60,0.37),oma=c(0,0,0,2.4))
plot(pft_370_1,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_370_1,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_370_1,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_370_1,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP3-7.0 in 2011-2040", side =2, line = -2, adj = 0.94,cex = 1.7, outer=TRUE)
plot(pft_370_4,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[5], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_370_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[6], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_370_4,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[7], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_370_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[8], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP3-7.0 in 2041-2070", side =2, line = -2, adj = 0.5,cex = 1.7, outer=TRUE)
plot(pft_370_7,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[9], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_370_7,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[10], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_370_7,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[11], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_370_7,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[12], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP3-7.0 in 2071-2100", side =2, line = -2, adj = 0.08,cex = 1.7, outer=TRUE)
dev.off()

png(filename = "Appendix_Figure_4_future_scenarios_SSP585.png",height = 2650, width = 3400, res = 180)
par(mfrow=c(3,4), cex.lab=1.25,mai=c(0.80,0.80,0.60,0.37),oma=c(0,0,0,2.4))
plot(pft_585_1,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[1], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_585_1,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[2], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_585_1,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[3], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_585_1,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[4], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP5-8.5 in 2011-2040", side =2, line = -2, adj = 0.94,cex = 1.7, outer=TRUE)
plot(pft_585_4,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[5], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_585_4,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[6], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_585_4,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[7], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_585_4,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[8], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP5-8.5 in 2041-2070", side =2, line = -2, adj = 0.5,cex = 1.7, outer=TRUE)
plot(pft_585_7,main="Local seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                       font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[9], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfot_585_7,main="Optimal seedlot growth potential",cex.main = 1.8,legend.args = list(text = '16-year tree height (cm)', side = 4, 
                                                                                          font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[10], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pg_585_7,col=col_3,main="Height increase/decrease", cex.main = 1.8,legend.args = list(text = '16-year tree height change(cm)', side = 4, 
                                                                                           font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[11], ")"), side = 1, line = -1, adj = 0, font = 2)
plot(pfg_585_7,main="Growth potential gain using optimal seedlot",cex.main = 1.8,legend.args = list(text = '16-year tree height gain (cm)', side = 4, 
                                                                                                    font = 2, line = 2.7, cex = 0.8))
mtext(paste0("(", letters[12], ")"), side = 1, line = -1, adj = 0, font = 2)
mtext("SSP5-8.5 in 2071-2100", side =2, line = -2, adj = 0.08,cex = 1.7, outer=TRUE)
dev.off()
#Figure added:===========

R_squared <- c(0.227,mean(0.473,0.439),mean(0.667,0.665,0.640,0.649),mean(0.662,0.763,0.703),0.803,0.748)
AIC <- c(26542,mean(25736,25868),mean(24768,24783,24938,24879),mean(24807,24052,24535),mean(23671),24202)
RMSE <- c(128.6,mean(106.2,109.6),mean(84.4,84.7,87.8,86.6),mean(85.0,71.1,79.7),64.9,73.4)
x <- c(2,3,4,5,6,7)

data <- data.frame(x, R_squared, AIC,RMSE)

fa <- ggplot(data, aes(x)) +
  geom_line(aes(y = R_squared, color = "R-squared")) +
  geom_line(aes(y = (AIC - min(AIC))/(max(AIC) - min(AIC)), color = "AIC")) +
  scale_y_continuous(sec.axis = sec_axis(~ .* (max(AIC) - min(AIC)) + min(AIC), name = "AIC"), 
                     limits = c(0, 1.0), name = "R-squared") +
  scale_color_manual(values = c("R-squared" = "blue", "AIC" = "red")) +
  labs(x = "Number of variables") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey95"),
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.left = element_text(color = "blue"),
        axis.text.y.left = element_text(color = "blue"),
        legend.position = "none")

ggsave("Figure_add_R2_AIC_No.png",plot = fa, dpi = 170, width = 6, height = 4, units = "in")
?ggsave

# Data

R_squared <- c(0.227, mean(0.473, 0.439), mean(0.667, 0.665, 0.640, 0.649), mean(0.662, 0.763, 0.703), 0.803, 0.748)
AIC <- c(26542, mean(25736, 25868), mean(24768, 24783, 24938, 24879), mean(24807, 24052, 24535), mean(23671), 24202)
RMSE <- c(128.6, mean(106.2, 109.6), mean(84.4, 84.7, 87.8, 86.6), mean(85.0, 71.1, 79.7), 64.9, 73.4)
x <- c(2, 3, 4, 5, 6, 7)


# Set up the plot
png(filename = "Figure_add_R2_AIC_No.png",height = 1080, width = 1560, res = 170)
par(mfrow=c(1,1),mar = c(5, 4, 4, 4) + 0.1)
plot(x, R_squared, xlim = c(1.5,7),type = "n", ylab = "", xlab = "Number of unique climate variables", ylim = c(0.2, 1.0))
axis(2, col.axis = 'red', col.ticks = 'red')
# Plot R_squared
lines(x, R_squared, col = "red",lwd=2)
mtext("R-squared", side = 2, line = 3,col = "red")
# Plot AIC on a second y-axis
par(new = TRUE)
plot(x, AIC, type = "n", ylab = "", xaxt = "n", yaxt = "n", xlab ="",ylim = c(23000, 27000),xlim = c(1.5,7))
lines(x, AIC, col = "slateblue3",lwd=2)
axis(side = 4, at = pretty(range(AIC)), labels = pretty(range(AIC)), col.axis = 'slateblue3', col.ticks = 'slateblue3')
mtext("AIC", side = 4, line = 3,col = "slateblue3")

# Plot RMSE on a third y-axis
par(new = TRUE)
plot(x, RMSE, type = "n", ylab = "", xaxt = "n", yaxt = "n", xlab ="",ylim = c(50, 130),xlim = c(1.5,7))
lines(x, RMSE, col = "orange",lwd=2)
axis(side = 2, at = c(50,70,90,110,130),NA, tcl = 0.5,col.ticks= "orange")
text(1.5,c(50,70,90,110,130),labels= c(50,70,90,110,130),srt=90, col = "orange")
mtext("RMSE", side = 2, line = -3.5, col="orange")

abline(h = seq(from = 0, to = 130, by = 10), col = "lightgray", lty = "dashed")
abline(v = seq(from = 0, to = max(x), by = 1), col = "lightgray", lty = "dashed")
dev.off()

# possible to be used====


par(mfrow=c(1,2))
plot(pf_370_7)


col_142 <- c("#F2F2F2", colorRampPalette(c("coral","yellow"))(70),colorRampPalette(c("yellow","#00A600"))(35))

summary(URF_4_1)
par(mfrow=c(3,3))
# Loop over the plots
for (i in 1:9) {
  # Generate random data and plot it
  x <- rnorm(10)
  y <- rnorm(10)
  plot(x, y, main = paste("Plot", i))
  
  # Add label to bottom left corner
  mtext(paste0("(", letters[i], ")"), side = 1, line = -1, adj = 1, font = 2)
}
