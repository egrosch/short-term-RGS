################################################################################
# wheat data (MultiResistGS)
# Estimation of genetic effects
################################################################################

library ("sqldf")
library ("SelectionTools")

st.input.dir  <- "input"
st.output.dir <- "output"
gs.set.num.threads(2)

################################################################################
# Effect estimation

st.read.marker.data ("wheat.mpo",format="m",data.set="P")
st.read.map         ("wheat.map",skip=1, format="mcp",data.set="P")
st.read.performance.data ("wheat.dta",data.set="P")

#processing data
st.restrict.marker.data (NoAll.MAX = 2, data.set = "P")
st.restrict.marker.data (MaMis.MAX = 0.1, data.set = "P")
st.restrict.marker.data (ExHet.MIN = 0.1, data.set = "P")
st.restrict.marker.data (InMis.MAX= 0.1, data.set = "P")

st.copy.marker.data("PBLUP","P")

gs.esteff.rr("BLUP", data.set="PBLUP") 

eff <- gs.return.effects ( data.set="PBLUP" ) 
gs.write.pseff ( eff, file="data/c001-yld-wheat.eff" ) 


   
