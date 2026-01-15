################################################################################
# RGS: wheat dataset
# Strategy MAXG2-TOP/UC1-SET
################################################################################

source("functions.R")

library ("SelectionTools")
library("sqldf")

################################################################################

st.input.dir  <- "input"
st.output.dir <- "output526"
dir.create(st.output.dir)
st.set.info.level (-2)
gs.set.num.threads(2)

NRUN <- 200  # how many sets of crosses are evaluated

################################################################################
# Wheat

eff.file <- "data/c001-yld-wheat.eff"

st.read.marker.data ("wheat.mpo",format="m",data.set="PP")
st.read.map         ("wheat.map",skip=1, format="mcp",data.set="PP")
st.read.performance.data ("wheat.dta",data.set="PP")

#processing data
st.restrict.marker.data (NoAll.MAX = 2, data.set = "PP")
st.restrict.marker.data (MaMis.MAX = 0.1, data.set = "PP")
st.restrict.marker.data (ExHet.MIN = 0.1, data.set = "PP")
st.restrict.marker.data (InMis.MAX= 0.1, data.set = "PP")

gs.esteff.rr("BLUP",data.set="PP") 
yld.eff <- gs.return.effects(data.set="PP")

st.set.simpop ( pop.name="PP", data.set="PP" ) 
load.effmap("yld", eff.file)


################################################################################

genotype.population("PP")
evaluate.population("PP", "yld")

# Select the best parental lines                        # SE-IL: GEGV

copy.population("PA","PP")
genotype.population("PA")
evaluate.population("PA", "yld")
population.sort("PA", decreasing=TRUE) 
population.divide("Psel", "PA", 144)

###########################################
# simulation
###########################################

# Cross selected parental lines                        # CR-IL: MAXG

st.get.simpop("Psel","Psel")
gs.set.effects(eff=yld.eff,data.set="Psel")
gs.cross.eval.ma(data.set="Psel")
gs.cross.eval.mu(data.set="Psel")
gs.cross.eval.va(pop.type="DH", data.set="Psel")
gs.cross.eval.es(alpha=1, data.set="Psel")
crs.Psel  <- gs.cross.info(data.set="Psel",sortby ="ma")
d.sorted <- data.frame(OTU1     =crs.Psel$P1No,
                       OTU2     =crs.Psel$P2No,
                       Measure  =crs.Psel$ma,stringsAsFactors=FALSE)
crs     <- sel.crs (d.sorted,ncp=1,nct=72) 

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

remove.population("F1")
for (ii in 1:72) {
  nme1 <- sprintf("p%03i",crs[ii,1])
  nme2 <- sprintf("p%03i",crs[ii,2])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
  append.population("F1",nme3)
}

# Cross F1 plants 
st.get.simpop("F1","F1")
gs.set.effects(eff=yld.eff,data.set="F1")
gs.cross.eval.ma(data.set="F1")
gs.cross.eval.mu(data.set="F1")
gs.cross.eval.va(pop.type="CRS", data.set="F1")
gs.cross.eval.es(alpha=0.2,data.set="F1")
crs  <- gs.cross.info(data.set="F1",sortby ="ma")      # CR-F1: MAXG
d.sorted <- data.frame(OTU1     =crs$P1No,
                       OTU2     =crs$P2No,
                       Measure  =crs$ma,stringsAsFactors=FALSE)
crs     <- sel.crs (d.sorted,ncp=1,nct=36)

for (ii in 1:36) {                
   nme1 <- sprintf("f1%02i",crs[ii,1])
   nme2 <- sprintf("f1%02i",crs[ii,2])
   nme3 <- sprintf("SYN1%02i",ii)   
   cross(nme3,nme1,nme2,10)         
 }

remove.population("SYN1")
for (ii in 1:36) {
  nme3 <- sprintf("SYN1%02i",ii)
  append.population("SYN1",nme3)
}

copy.population("SYN","SYN1")

C <- 2; {

    st.get.simpop("SYN","SYN")
    gs.set.effects(eff=yld.eff,data.set="SYN")
    gs.cross.eval.mu(data.set="SYN")
    gs.cross.eval.va(pop.type="CRS", data.set="SYN")
    gs.cross.eval.es(alpha=0.2,data.set="SYN")
    crs.Psel  <- gs.cross.info(data.set="SYN",sortby ="es")      # CR-C1: UC-OCS
    
    crs <- sel.set(ncp=1, nct=36, crs.list=crs.Psel)
    
    population.copy("tmp","SYN")
    for (ii in 1:360) {
        nme <- sprintf("sy2%03i",ii)
        population.divide (nme,"tmp",1)
    }
    
    for (ii in 1:36) {                
        nme1 <- sprintf("sy2%03i",crs[ii,1])
        nme2 <- sprintf("sy2%03i",crs[ii,2])
        nme3 <- sprintf("SYN2%02i",ii)   
        cross(nme3,nme1,nme2,10)         
    }
    
    remove.population("SYN2")
    for (ii in 1:36) {
        nme3 <- sprintf("SYN2%02i",ii)
        append.population("SYN2",nme3)
    }

    copy.population("SYN","SYN2")
    
} # C=2

C <- 3; {

    genotype.population("SYN")                 
    evaluate.population("SYN", "yld")
    population.sort("SYN", decreasing=TRUE)    
    population.divide("SYNsel", "SYN", 72)     
    
    nme <- sprintf("SYN%1isel",C-1)
    copy.population(nme,"SYNsel")
    
    population.copy("split","SYNsel")
    for (ii in 1:72) {
      nme <- sprintf("s%03i",ii)
      population.divide (nme,"split",1)
    }
    
    idx <- sample(1:72,72)                     
    
    for (ii in 1:36) {
      nme1 <- sprintf("s%03i",idx[ii])
      nme2 <- sprintf("s%03i",idx[36+ii])
      nme3 <- sprintf("SYNn%02i",ii)
      cross(nme3,nme1,nme2,10)
    }
    
    remove.population("SYN")
    for (ii in 1:36) {
      nme3 <- sprintf("SYNn%02i",ii)
      append.population("SYN",nme3)
    }
    
    nme <- sprintf("SYN%1i",C)
    copy.population (nme,"SYN") 

} # C=3


# Select 36 SYN3 plants

genotype.population("SYN")                 
evaluate.population("SYN", "yld")
population.sort("SYN", decreasing=TRUE)   # SE-3: GEGV 
population.divide("SYNsel", "SYN", 36)     
copy.population("SYN3sel","SYNsel")

# 6 DH Lines per selected SYN3 plant

remove.population("DH")
population.copy("split","SYNsel")
for (ii in 1:36) {
    population.divide ("ss","split",1)
    dh ("dd","ss",6)
    append.population("DH","dd")
}

eval <- c("SYN2","SYN3","DH")

v <- d <- dA <- NULL
for (i in 1:length(eval))
{
  genotype.population ( eval[i] )                 
  evaluate.population ( eval[i], "yld")
  P <- get.population.gvalue(eval[i])$gvalue
  d <- rbind(d,dA,data.frame (gen=eval[i],y=P))
  #
  D <- eval[i]
  
  {
    st.get.simpop(D,D)
    st.def.hblocks ( hap=5, hap.unit=2, data.set=D ) 
    st.recode.hil  (data.set=D)
    div <- mean(st.marker.data.statistics(data.set=D)$marker.list$ExHet)
  }
  
  {
    st.get.simpop(D,D)
    divA <- mean(st.marker.data.statistics(data.set=D)$marker.list$ExHet)
    
  }
  
  v <- rbind(v,data.frame (gen=eval[i],d=div, dA = divA))
}
m <- tapply(d$y,d$gen,mean);
e <- merge(data.frame(gen=names(m),y=m),v)
e
