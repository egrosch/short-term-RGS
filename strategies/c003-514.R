################################################################################
# RGS: wheat dataset
# Strategy best
################################################################################

library ("SelectionTools")
library("sqldf")

################################################################################

st.input.dir  <- "input"
st.output.dir <- "output514"
dir.create(st.output.dir)
st.set.info.level (-2)
gs.set.num.threads(2)

dbfile <- "data/c003_514.sqldb" # database

NREP <- 300 # replications for the simulation

################################################################################
# Wheat

crop <- "wheat"
strategy <- "best"
rfile <- "resultswheat" # file for the results in the database
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
# Loop for the simulation
###########################################

e    <- NULL

st.set.info.level(-2)

for (REP in 1:NREP) {

cat (sprintf("%05i\r",REP))

# Cross selected parental lines                        # CR-IL: best

#copy.population("PA","PP")
genotype.population("PA")
evaluate.population("PA", "yld")
population.sort("PA", decreasing=TRUE) 

idx <- c(seq(1, 143, by = 2),seq(2, 144, by = 2))

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

remove.population("F1")
for (ii in 1:72) {
  nme1 <- sprintf("p%03i",idx[ii])
  nme2 <- sprintf("p%03i",idx[72+ii])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
  append.population("F1",nme3)
}
                                    

# Cross F1 plants 
st.get.simpop("F1","F1")
genotype.population("F1")
evaluate.population("F1", "yld")
population.sort("F1", decreasing=TRUE) 
    
idx <- c(seq(1, 71, by = 2),seq(2, 72, by = 2))         #CR-F1: best

for (ii in 1:36) {                
   nme1 <- sprintf("f1%02i",idx[ii])
   nme2 <- sprintf("f1%02i",idx[ii+36])
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
  genotype.population("SYN")
  evaluate.population("SYN", "yld")
  population.sort("SYN", decreasing=TRUE) 
    
  idx <- c(seq(1, 71, by = 2),seq(2, 72, by = 2))         #CR-C1: best
    
    population.copy("tmp","SYN")
    for (ii in 1:360) {
        nme <- sprintf("sy2%03i",ii)
        population.divide (nme,"tmp",1)
    }
    
    for (ii in 1:36) {                
        nme1 <- sprintf("sy2%03i",idx[ii])
        nme2 <- sprintf("sy2%03i",idx[ii+36])
        nme3 <- sprintf("SYN2%02i",ii)   
        cross(nme3,nme1,nme2,10)         
    }
    
    remove.population("SYN2")
    for (ii in 1:36) {
        nme3 <- sprintf("SYN2%02i",ii)
        append.population("SYN2",nme3)
    }
    remove.population("SYN")
    copy.population("SYN","SYN2")
    
} # C=2

C <- 3; {

    st.get.simpop("SYN","SYN")
  genotype.population("SYN")
  evaluate.population("SYN", "yld")
  population.sort("SYN", decreasing=TRUE) 
    
  idx <- c(seq(1, 71, by = 2),seq(2, 72, by = 2))       #CR-C2: best
    
    population.copy("tmp","SYN")
    for (ii in 1:360) {
        nme <- sprintf("sy3%03i",ii)
        population.divide (nme,"tmp",1)
    }
    
    for (ii in 1:36) {                
      nme1 <- sprintf("sy3%03i",idx[ii])
      nme2 <- sprintf("sy3%03i",idx[ii+36])
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
rownames(e) <- c()
e$Strategy <- strategy
e$Crop <- crop
e$Date <- date()

# Save simulation results in data base
conn <- dbConnect(RSQLite::SQLite(), dbfile)
dbWriteTable(conn, rfile, e, append=TRUE)
dbDisconnect(conn)

} # for (REP in 1:NREP)
