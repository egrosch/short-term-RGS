################################################################################
# RGS: oat dataset
# Strategy SM-UC-1
################################################################################

#function
SimpleMatingHet <-function(Generation,alpha,NumCross,cullpairK,typk)
{
  
  #calculate relationship matrix G for K
  
  remove.population("testpop")
  population.copy("testpop",Generation)
  st.get.simpop("testpop","testpop")
  st.read.map(map,skip=1, format="mcp",data.set="testpop")
  genotype.population("testpop")                
  evaluate.population("testpop", "yld")
  
  # can only calculate ZZ for markers with at leest two alleles
  xx <- st.marker.data.statistics("testpop")
  if(length(which(xx$marker.list$NoAll<2))!=0){st.restrict.marker.data(data.set="testpop", ExHet.MIN = 0.0001) # entfernt alle Marker, bei denen es nur ein Allel gibt
  }
  ZZ <- gs.build.Z(data.set="testpop", out.filename="Z.matrix", auxfiles=T)
  AA <- as.matrix(ZZ)
  GG <- Gmatrix(SNPmatrix=AA, 
                maf=0.05, method="VanRaden")
  
  #need the original names for the G matrix
  ttt <- st.marker.data.statistics(Generation)
  colnames(GG) <- ttt$individual.list$Name
  row.names(GG) <- ttt$individual.list$Name
  relMat <- GG
  
  #calculate UC with SelectionTools #Alpha can now be 0
  gs.cross.eval.ma(data.set=Generation)
  gs.cross.eval.mu(data.set=Generation)
  gs.cross.eval.va(pop.type=typk, data.set=Generation)
  gs.cross.eval.es(alpha=alpha, data.set=Generation)
  crs.Psel  <- gs.cross.info(data.set=Generation,sortby ="ma")
  
  #combine Parents, criteria and relativeness
  for (u in 1:length(crs.Psel[,1])) {
    crs.Psel$K[u] <- relMat[which(colnames(relMat)==crs.Psel$P1Name[u]),which(colnames(relMat)==crs.Psel$P2Name[u])]
  }
  KrzAuswahl <- crs.Psel[,c(3,4,10,11)]
  colnames(KrzAuswahl) <- c("Parent1", "Parent2", "Y", "K")
  
  ## Running the optimization
  maxGain = selectCrosses(data = KrzAuswahl,
                          n.cross = NumCross,  #Number of Crosses
                          max.cross = 1,
                          min.cross = 0,
                          #K= relMat, #KrzAuswahl[,4],#relMat
                          culling.pairwise.k = cullpairK)  
  
  #save results in a data.frame
  cross <- data.frame(elter1<- rep(0, times=NumCross),elter2<- rep(0, times=NumCross))
  colnames(cross) <-c("Elter1","Elter2")
  for (v in 1:NumCross) {
    cross$Elter1[v] <- as.numeric(unlist(strsplit(maxGain[[2]]$Parent1[v], split='.', fixed=TRUE))[2])
    cross$Elter2[v] <- as.numeric(unlist(strsplit(maxGain[[2]]$Parent2[v], split='.', fixed=TRUE))[2])
  }
  return(cross)
}

library("SimpleMating")
library ("SelectionTools")
library("sqldf")
library("rrBLUP")
library(AGHmatrix)

einst <- 1.3 

################################################################################

st.input.dir  <- "input"
st.output.dir <- "output661"
dir.create(st.output.dir)
st.set.info.level (-2)
gs.set.num.threads(2)

dbfile <- "data/c003_661.sqldb" # database

NREP <- 300

################################################################################
# Oat

crop <- "oat"
strategy <- "SM-UC-1"
rfile <- "resultsoat" # file for the results in the database
eff.file <- "data/c001-yld-oat.eff"
map <- "oat.map"

st.read.marker.data ("oat.mpo",format="m",data.set="PP") 
st.read.map         ("oat.map",skip=1, format="mcp",data.set="PP")
st.read.performance.data ("oat.dta",data.set="PP")

#Markerdaten aufbereiten für Vergleichbarkeit der Diversität
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
# Loop for the simulation
###########################################

e    <- NULL

st.set.info.level(-2)

for (REP in 1:NREP) {

cat (sprintf("%05i\r",REP))

# Cross selected parental lines                        # CR-IL: SimpleMating-UC

st.get.simpop("Psel","Psel")
gs.set.effects(eff=yld.eff,data.set="Psel")
crs<-SimpleMatingHet("Psel",1,72,einst,"DH")

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
crs<-SimpleMatingHet("F1",0.2,36,einst,"CRS")

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
    crs<-SimpleMatingHet("SYN",0.2,36,einst,"CRS")
    
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

    st.get.simpop("SYN","SYN")
    gs.set.effects(eff=yld.eff,data.set="SYN")
    crs<-SimpleMatingHet("SYN",0.2,36,einst,"CRS")
    
    population.copy("tmp","SYN")
    for (ii in 1:360) {
        nme <- sprintf("sy3%03i",ii)
        population.divide (nme,"tmp",1)
    }
    
    for (ii in 1:36) {                
        nme1 <- sprintf("sy3%03i",crs[ii,1])
        nme2 <- sprintf("sy3%03i",crs[ii,2])
        nme3 <- sprintf("SYN3%02i",ii)   
        cross(nme3,nme1,nme2,10)         
    }
    
    remove.population("SYN3")
    for (ii in 1:36) {
        nme3 <- sprintf("SYN3%02i",ii)
        append.population("SYN3",nme3)
    }
    
    copy.population("SYN","SYN3")

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
    ssd.mating ("dd","ss",6, 5)
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
rownames(e) <- c()
e$Strategy <- strategy
e$Crop <- crop
e$Date <- date()

# Save simulation results in data base
conn <- dbConnect(RSQLite::SQLite(), dbfile)
dbWriteTable(conn, rfile, e, append=TRUE)
dbDisconnect(conn)

} # for (REP in 1:NREP)
