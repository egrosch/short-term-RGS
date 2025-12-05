################################################################################
# RGS: wheat dataset
# Strategy EMBV
################################################################################

calc.embv <- function(pop,N,G,REP){
    l <- list.populations()
    n <- l[ pop == l[,1],2]
    embv <- rep(0,n)
    gv   <- rep(0,REP)
    population.copy("EMBVsplit",pop)
    for (ii in 1:n) {
        population.divide ("EMBVSC","EMBVsplit",1)
        for (r in 1:REP) { 
            dh ("EMBVDH","EMBVSC",G)
            genotype.population("EMBVDH")
            evaluate.population("EMBVDH", "yld")
            population.sort("EMBVDH", decreasing=TRUE) 
            population.divide("EMBVDHsel", "EMBVDH", N)
            gval <-get.population.gvalue("EMBVDHsel")
            gv[r] <- mean(gval$gvalue)
        }
        embv[ii] <- mean(gv)
    }
    return(embv)
}

library ("SelectionTools")
library("sqldf")

################################################################################

st.input.dir  <- "input"
st.output.dir <- "output528"
dir.create(st.output.dir)
st.set.info.level (-2)
gs.set.num.threads(2)

dbfile <- "data/c003_528.sqldb" # database

NREP <- 300 # replications for the simulation

################################################################################
# Wheat

crop <- "wheat"
strategy <- "EMBV"
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

# Select the best 144 P lines / GEGV

copy.population("PA","PP")
genotype.population("PA")
evaluate.population("PA", "yld")
population.sort("PA", decreasing=TRUE) 
population.divide("Psel", "PA", 144)

###########################################
# Loop for the simulation
###########################################

e    <- NULL

for (REP in 1:NREP) {

cat (sprintf("%05i\r",REP))
# Random intermating of the selected P 

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

idx <- sample(1:144,144)

for (ii in 1:72) {
  nme1 <- sprintf("p%03i",idx[ii])
  nme2 <- sprintf("p%03i",idx[72+ii])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
}


for (ii in 1:36) {
  nme1 <- sprintf("f1%02i",ii)
  nme2 <- sprintf("f1%02i",36+ii)
  nme3 <- sprintf("SYN1%02i",ii)
  cross(nme3,nme1,nme2,10)
}

remove.population("SYN1")
for (ii in 1:36) {
  nme3 <- sprintf("SYN1%02i",ii)
  append.population("SYN1",nme3)
}

# Select the best 72 SYN1 plants / GEGV

copy.population("SYN","SYN1")

for (C in 2:3 ) {

#genotype.population("SYN")                 # Instead of selection accoriding
#evaluate.population("SYN", "yld")          # to the GEBV

embv <- calc.embv (pop="SYN", N=2, G=10, REP=20)  # We use the EMBV
set.population.gvalue("SYN",embv) 

population.sort("SYN", decreasing=TRUE)    # Selection accoriding to
population.divide("SYNsel", "SYN", 72)     # GEGV

nme <- sprintf("SYN%1isel",C-1)
copy.population(nme,"SYNsel")

population.copy("split","SYNsel")
for (ii in 1:72) {
  nme <- sprintf("s%03i",ii)
  population.divide (nme,"split",1)
}

idx <- sample(1:72,72)                     # Random mating of the selected
                                           # fraction
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
copy.population (nme,"SYN")                # Generates the new SYN population

} # for (C in 2:3 )

# Select 36 SYN3 plants

#genotype.population("SYN")                 
#evaluate.population("SYN", "yld")

embv <- calc.embv (pop="SYN", N=2, G=6, REP=20)  
set.population.gvalue("SYN",embv) 

population.sort("SYN", decreasing=TRUE)    
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
