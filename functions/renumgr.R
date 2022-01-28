#########################################################################
# Function to renumber the splits/groups in mvpart
#########################################################################
renumgr <- function(clusgr){
  aa <- 1
  renum <- rep(1, length(clusgr))
  for(i in 2:length(clusgr)) {
    if(clusgr[i] != clusgr[i-1]) aa <- aa+1
    renum[i] <- aa }
  return(renum)    }