############################################################################################
###### Function ############################################################################
############################################################################################

 ExpBVSelGrOff<-function(Generation,SelNach,map,eff.file)
 {

verwen1<-st.marker.data.statistics(Generation) 
verwen2<-verwen1$genotypes

for (i in 1:length(verwen2[,1])) {
  verwen2$Reihenfolge<-which(verwen2$`Mar/Ind`[i]==mabstand$Name)
}
verwen2<-verwen2[order(verwen2$Reihenfolge),]
verwen2<-verwen2[,-length(verwen2[1,])]
rownames(verwen2)<-verwen2$`Mar/Ind`
verwen2<-verwen2[,-1]

####---
#calculate parental breeding values
gs.set.effects(eff=yld.eff,data.set=Generation)
yy <- gs.predict.genotypes(training.set=Generation, prediction.set=Generation)

Elter<-yy$i
AnzElter<-length(Elter)

####---
# selection intensity

N <- SelNach #number of selected offspring
G <- length(Elter) 
alpha <- N/G
i.alpha <- dnorm(qnorm(1-alpha,0,1),0,1)/alpha

if(i.alpha==0){
  ip <- 0
}
if(i.alpha!=0){
  ip <- i.alpha - ((G-N)/(2*N*(G+1)*i.alpha))
}

####---

verwendm <- row.names(verwen2)
mabstand<-map
mabstand$Reihenfolge <- rownames(mabstand)
rownames(mabstand) <- mabstand$Name
mabstand<- mabstand[mabstand$Name %in% verwendm,]
str(mabstand)
mabstand$Reihenfolge <- as.numeric(mabstand$Reihenfolge)
mabstand <-mabstand[order(mabstand$Reihenfolge, decreasing = F),]

####---

xe <-unique(mabstand$Chrom)
xxe <-length(xe)
dime <-length(mabstand$Name)
LDM5 <-matrix(rep(999,times= dime^2), nrow = dime)
t<-1

for (i in 1:xxe) {
  mab1 <-subset(mabstand,mabstand$Chrom==xe[i])
  länge <-length(mab1$Pos)
  
  for (ii in 1:länge) {
    LDM5[t:(länge+t-1),ii+t-1] <-abs(mab1$Pos-mab1$Pos[ii])
  }
  
  t<-t+länge
}

####---

xe <-unique(mabstand$Chrom)
xxe <-length(xe)
dime <-length(mabstand$Name)
LDM <-matrix(rep(0,times= dime^2), nrow = dime)
t<-1

for (i in 1:xxe) {
  mab1 <-subset(mabstand,mabstand$Chrom==xe[i])
  länge <-length(mab1$Pos)
  
  for (ii in 1:länge) {
    LDM[t:(länge+t-1),ii+t-1] <-abs(mab1$Pos-mab1$Pos[ii])
  }
  
  t<-t+länge
}

#calculate r 
LDMr <- (exp(-2*(LDM/100)))/4  

LDMr[which(LDM==0)] <-0
for (i in 1:length(LDMr[1,])) {
  LDMr[i,i] <-0.25
}

R<-LDMr

####---

eff <- read.table(eff.file, header=F, skip=1) # marker effects
eff <- eff[1:(nrow(eff)/2),]
colnames(eff) <- c("marker", "allele", "effect")
st.restrict.marker.data(data.set=Generation, mar.list=eff$marker) 

marp<-st.marker.data.statistics(Generation)

zz<- marp$genotypes
nn1<-length(marp$genotypes[,1])
nn2<-length(marp$genotypes[1,])-1

zzn1 <- zz[,-1]
u<-length(colnames(zzn1))
allelgeno2<-matrix(rep(0,times=u*length(zzn1[,1])),ncol = u)
colnames(allelgeno2)<- colnames(zzn1)
allelgeno1<-allelgeno2

for (i in 1:nn1) {
  gesuchtesAllel<-eff$allele[which(eff$marker==marp$genotypes[i,1])]
  
  for (ii in 1:nn2) {
    
    allel1<-substr(marp$genotypes[i,ii+1], 1, 1)==gesuchtesAllel
    allel2<-substr(marp$genotypes[i,ii+1], 3, 3)==gesuchtesAllel
    
    zz[i,ii+1]<-allel1-allel2 
    
    allelgeno1[i,ii]<-allel1*2-1 
    allelgeno2[i,ii]<-allel2*2-1 
  }
}  

row.names(zz)<-zz$`Mar/Ind`
zzn <- zz[,-1]

####---

mvec <- zzn

mvecn <- mvec
for (qq in 1:length(mvec[1,])) {
  mvecn[qq] <- as.numeric(mvec[,qq])*as.numeric(eff$effect)
}

FSn <- matrix(rep(0,times=length(mvec[1,])),nrow = length(mvec[1,]))
FS <- matrix(rep(0,times=(length(mvec[1,]))^2),nrow = length(mvec[1,]))

for (i in 1:length(mvec[1,])) { 
  FSn[i] <-t(as.numeric(mvecn[,i]))%*%R%*%as.numeric(mvecn[,i])
}

for (i in 1:length(mvec[1,])) { 
  for (ii in 1:length(mvec[1,])) {
    FS[ii,i] <-FSn[i]+FSn[ii]
  }
}

####---
#variance GO
x <-qnorm(1-alpha,0,1) 
if(ip==0){
  k <- 0
}
if(ip!=0){
  k<-ip*(ip-x)
}

####---
# Sigma
####---

VDH<-as.data.frame(matrix(rep(0,times=(nn2^2)),nrow = nn2))
colnames(VDH)<- colnames(zzn)
row.names(VDH)<- colnames(zzn)

c1 <- 0.5*(1-exp(-2*(LDM/100)))
c1[LDM5==999] <- 0.5
c2<- c1*2
cc.1 <- 1-c2

xe2 <-unique(mabstand$Chrom)
xxe2 <-length(xe2)

beta<-eff$effect

ElternVDH <- as.data.frame(matrix(rep(0,times=(nn2)),ncol = nn2))

colnames(ElternVDH)<- colnames(zzn)
#########################

for (Ekkv in 1:nn2) {  
    
    nm <- length(mvec[,1])
    X1 <- c(allelgeno1[,Ekkv])
    X2 <- c(allelgeno2[,Ekkv])
   
    VarDH<-0
    for (ChrNum in 1:xxe2) {
      
      D12 <-1/16 * (X1[which(mabstand$Chrom==xe2[ChrNum])] %*% t(-X2[which(mabstand$Chrom==xe2[ChrNum])]))
     
      D11 <-1/16 * (X1[which(mabstand$Chrom==xe2[ChrNum])] %*% t(X1[which(mabstand$Chrom==xe2[ChrNum])]))
      D22 <-1/16 * (X2[which(mabstand$Chrom==xe2[ChrNum])] %*% t(X2[which(mabstand$Chrom==xe2[ChrNum])]))
      
      ccc <- cc.1[which(mabstand$Chrom==xe2[ChrNum]),which(mabstand$Chrom==xe2[ChrNum])]
      beta2 <- beta[which(mabstand$Chrom==xe2[ChrNum])]
      
      Sigma12 <- ccc*D12
      Sigma12 <- ccc*Sigma12
      varDH12 <- t(beta2) %*% Sigma12 %*% beta2
 
      Sigma11 <- ccc*D11
      Sigma11 <- (2*D11)+Sigma11
      Sigma11 <- ccc*Sigma11
      varDH11 <- t(beta2) %*% Sigma11 %*% beta2
      Sigma22 <- ccc*D22
      Sigma22 <- (2*D22)+Sigma22
      Sigma22 <- ccc*Sigma22
      varDH22 <- t(beta2) %*% Sigma22 %*% beta2
      
      VarDH[ChrNum] <- (varDH12)*2+ varDH11+varDH22
      
    }
    sum(VarDH)
    var.DH<-sum(VarDH)
    
    ElternVDH[Ekkv] <- var.DH
    
}




#########################

for (kkv in 1:nn2) {
  for (mkv in 1:nn2) {
    combination <- c(kkv,mkv)  

    
    nm <- length(mvec[,1])
    X1 <- c(allelgeno1[,combination[1]]) 
    X2 <- c(allelgeno2[,combination[1]])
    X3 <- c(allelgeno1[,combination[2]])
    X4 <- c(allelgeno2[,combination[2]])
    
    VarDH<-0
    for (ChrNum in 1:xxe2) {
      
      D13 <-1/16 * (X1[which(mabstand$Chrom==xe2[ChrNum])] %*% t(-X3[which(mabstand$Chrom==xe2[ChrNum])]))
      D14 <-1/16 * (X1[which(mabstand$Chrom==xe2[ChrNum])] %*% t(-X4[which(mabstand$Chrom==xe2[ChrNum])]))
      D23 <-1/16 * (X2[which(mabstand$Chrom==xe2[ChrNum])] %*% t(-X3[which(mabstand$Chrom==xe2[ChrNum])]))
      D24 <-1/16 * (X2[which(mabstand$Chrom==xe2[ChrNum])] %*% t(-X4[which(mabstand$Chrom==xe2[ChrNum])]))
 
      ccc <- cc.1[which(mabstand$Chrom==xe2[ChrNum]),which(mabstand$Chrom==xe2[ChrNum])]
      beta2 <- beta[which(mabstand$Chrom==xe2[ChrNum])]
      
      Sigma13 <- ccc*D13
      varDH13 <- t(beta2) %*% Sigma13 %*% beta2
      
      Sigma14 <- ccc*D14
      varDH14 <- t(beta2) %*% Sigma14 %*% beta2
      
      Sigma23 <- ccc*D23
      varDH23 <- t(beta2) %*% Sigma23 %*% beta2
      
      Sigma24 <- ccc*D24
      varDH24 <- t(beta2) %*% Sigma24 %*% beta2

      VarDH[ChrNum] <- (varDH13+varDH14+varDH24+varDH23)*2
      
    }
    sum(VarDH)
    var.DH<-sum(VarDH) + ElternVDH[mkv] + ElternVDH[kkv]
    
    VDH[mkv,kkv] <- var.DH
    VDH[kkv,mkv] <- var.DH
    
  }
}

####---
  
  vgammsoff <- (VDH - FS)/4
  vgamseloff<- ((1-k)*FS)/4 +vgammsoff

  vBV <-var(yy$yhat)
  meanMSVind <- mean(FSn)

  vgamselindt<- ((1-k)*vBV)/4 + meanMSVind

  GO<- vgamseloff + vgamselindt

  zTerm <- ip*(0.5*(sqrt(FS))+sqrt(GO))

####---
  ### select parents

  Elter<-yy$i
  nn<-length(Elter)
  NN <- rep(0,times=nn)
  ee<-matrix(data=rep(NN,times=nn), nrow = nn,byrow = 1)

  rownames(ee)<- Elter
  colnames(ee)<- Elter

  #calculate ExpBVSelGrOff
  for (i in 1:nn) {
    for (ii in 1:nn) {
      BVe1<-yy$yhat[which(yy$i==Elter[i])]
      BVe2<-yy$yhat[which(yy$i==Elter[ii])]
    
      ####---
      ee[i,ii]<- 0.25*(BVe1+BVe2) + zTerm[i,ii]
    }
  }

####---
  ee1<-ee
  EE <-ee
  entf <- yy
  mm<- nn-SelNach
  
  if(mm!=0){
  
  for (l in 1:mm) {
  
    for (i in 1:nn) {
      EE[i,i]<- NA
    }
    
    EEt <- EE
    for (ko in 1:nn) {
      EEt[,ko] <-nn-as.numeric(as.factor(EE[,ko]))
    }
    
    for (i in 1:nn) {
      EEt[i,i]<- 1000
    }

    for (kk in 1:nn) {
      entf$yhat[kk] <-min(EEt[kk,])
    }
   
    ent<-entf$i[which(entf$yhat==max(entf$yhat))]
    ent<- ent[1]

    entff<- as.matrix(entf)
    entff<-entff[-which(entf$i==ent),]
    entf<-as.data.frame(entff)

    EE<- ee1[-which(colnames(EEt)==ent),-which(colnames(EEt)==ent)]
    ee1<-EE
    
    nn<-nn-1
  }
}

####---
  ## assign crosses ##
####---
  ee2 <- ee1
  pp <-length(ee2[1,])
  ppp<- pp/2

  for (i in 1:pp) {
    ee2[i,i]<- 0
  }

####---
  cross <- data.frame(elter1<- rep(0, times=SelNach/2),elter2<- rep(0, times=SelNach/2))
  colnames(cross) <-c("Elter1","Elter2")

  ee4 <- ee2
  ww <-length(ee2[1,])
  www<- ww/2

  for (i in 1:ww) {
    ee4[i,i]<- 0
  }

  for (v in 1:www) {
    samp <- sample(1:ww,1)
    mutt<-row.names(ee4)

    elt2<-order(ee4[samp,],decreasing = T)
    elte2<- elt2[1]

    Elter1<-mutt[samp]
    Elter2<-mutt[elte2]

    ee5 <-ee4[-c(samp,elte2),-c(samp,elte2)]
    ee4 <- NULL
    ee4<-ee5
    ww <- ww-2

    cross$Elter1[v] <- as.numeric(unlist(strsplit(Elter1, split='.', fixed=TRUE))[2])
    cross$Elter2[v] <- as.numeric(unlist(strsplit(Elter2, split='.', fixed=TRUE))[2])
  }
  return(cross)

 }
 