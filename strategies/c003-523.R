################################################################################
# RGS: wheat dataset
# Strategy MSD-1
################################################################################

sel.crs <- function (d.sorted,ncp,nct,old.crs=NULL)
                                        # ncp Number of crosses per parent  
                                        # nct Total number of crosses       
    {
        indx <- rep(FALSE,nrow(d.sorted))           

        if (!is.null(old.crs))
          for (i in 1:nrow(old.crs))
            for (j in 1:nrow(d.sorted))
              if ((old.crs[i,1]==d.sorted[j,1])&&(old.crs[i,2]==d.sorted[j,2]))
                indx[j] <- TRUE
        
        new.crs <- 0
        for (i in 1:nrow(d.sorted))  {
            crosses.so.far <- c(d.sorted$OTU1[indx],d.sorted$OTU2[indx])

            search1 <- paste("\\b",d.sorted$OTU1[i],"\\b",sep="")
            search2 <- paste("\\b",d.sorted$OTU2[i],"\\b",sep="")
            
            if ( (new.crs<nct)                               &&
                (ncp > length(grep(search1,crosses.so.far))) &&
                (ncp > length(grep(search2,crosses.so.far))) &&
                (FALSE == indx[i])                                 )
                {
                    indx[i] <- TRUE
                    new.crs = 1 + new.crs
                }
        }
        crosses <- d.sorted[indx,]

        return (crosses)
    }

library ("SelectionTools")
library("sqldf")

################################################################################

st.input.dir  <- "input"
st.output.dir <- "output523"
dir.create(st.output.dir)
st.set.info.level (-2)
gs.set.num.threads(2)

dbfile <- "data/c003_523.sqldb" # database

NREP <- 300

################################################################################
# Wheat

crop <- "wheat"
strategy <- "MSD-1"
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
# Loop for the simulations
###########################################

e    <- NULL

st.set.info.level(-2)

for (REP in 1:NREP) {

cat (sprintf("%05i\r",REP))

# Cross selected partental lines                        # CR-IL: MAX

st.get.simpop("Psel","Psel")
gs.set.effects(eff=yld.eff,data.set="Psel")
gs.cross.eval.ma(data.set="Psel")
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

    st.get.simpop("SYN1","SYN1")
    gs.set.effects(eff=yld.eff,data.set="SYN1")
    gs.cross.eval.mu(data.set="SYN1")
    gs.cross.eval.va(pop.type="CRS", data.set="SYN1")
    gs.cross.eval.es(alpha=0.2,data.set="SYN1")
    crs  <- gs.cross.info(data.set="SYN1",sortby ="es")      # CR-C1: MSD
    d.sorted <- data.frame(OTU1     =crs$P1No,
                           OTU2     =crs$P2No,
                           Measure  =crs$es,stringsAsFactors=FALSE)
    crs     <- sel.crs (d.sorted,ncp=1,nct=36)  
    
    population.copy("tmp","SYN1")
    for (ii in 1:360) {
        nme <- sprintf("sy1%03i",ii)
        population.divide (nme,"tmp",1)
    }
    
    for (ii in 1:36) {                
        nme1 <- sprintf("sy1%03i",crs[ii,1])
        nme2 <- sprintf("sy1%03i",crs[ii,2])
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
rownames(e) <- c()
e$Strategy <- strategy
e$Crop <- crop
e$Date <- date()

# Save simulation results in data base
conn <- dbConnect(RSQLite::SQLite(), dbfile)
dbWriteTable(conn, rfile, e, append=TRUE)
dbDisconnect(conn)

} # for (REP in 1:NREP)
