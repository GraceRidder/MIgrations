

library(ETBDsim)
library(parallel)
Jmaxvalues <- c(800,1200,1600,2000)
extinctions <- c(-10,-8, -6, -4, -2)

for (v in extinctions) {
  for (i in Jmaxvalues) {
    # #run the simulation multiple times in parallel
    trials = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    repsim <- function(trials) {
      reslist <- list()
      for (u in 1:length(trials)) {
        res1 = ETBDspaceSYM(
          t = 2000,
          psymp = .2,
          SRS = F,
          GEO = F,
          LOGNRM = T,
          watchgrow = F,
          SADmarg = .1,
          siteN = 1,
          JmaxV = c(i),
          exparm = v,
          ExpSpParm = 2,
          ExpSp = T,
          SPgrow = .0,
          splitparm = .25
        )
        reslist[[u]] <- res1
      }
      results = list(reslist = reslist)
      return(results)
    }


    rs= mcmapply(repsim, trials, SIMPLIFY = F, mc.cores = 1)

    filename <- paste("~/Desktop/spDepRES/SYM_LG_SP_",i,".", v, ".RData", sep = "")

    save(rs, file = filename)

  }

}


