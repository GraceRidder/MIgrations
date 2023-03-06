##rs is a list of simulation replicates

analyzeSim <- function(rs, t){
  sites <- length(rs[[1]]$reslist[[1]]$matrix_list)
  reps <- length(rs)


  SPsitelist <- vector(mode='list', length=sites)
  ABsitelist <- vector(mode='list', length=sites)

  for (site in 1:sites){
    X <- c()
    Xx <- c()
    for (rep in 1:reps){
      res <- rs[[rep]]$reslist[[1]]
      tpac <- c()
      apac <- c()
      for (time in 1:length(res$mig)) {
        tpac <- append(tpac, length(unlist(res$mig[time][[1]][site])))
        apac <- append(apac,  sum(unlist(res$mig[time][[1]][site])))
      }
      X <- cbind(X, tpac)
      Xx <- cbind(Xx, apac)
    }
    SPsitelist[[site]] <- X
    ABsitelist[[site]] <- Xx
  }

  #Turn simulations into MCMC chains

  ###this is set to autoburn which means only the second half of the simulations/chains are used for this diagnostic
  library(coda)

  #species richness
  gelmans_S <- c()
  for (site in 1:sites) {
    sg <- as.mcmc.list(lapply(as.data.frame(SPsitelist[[site]]), mcmc))
    c <- gelman.diag(sg)$psrf[1, 1]
    gelmans_S <- rbind(gelmans_S, c)
    nam <- paste("site ", site, sep = "")
    row.names(gelmans_S)[site] <- nam
  }

  gelmans_A <- c()
  for (site in 1:sites) {
    #abundance
    sg <- as.mcmc.list(lapply(as.data.frame(ABsitelist[[site]]), mcmc))
    c <- gelman.diag(sg)$psrf[1, 1]
    gelmans_A <- rbind(gelmans_A, c)
    nam <- paste("site ", site, sep = "")
    row.names(gelmans_A)[site] <- nam
  }

  colnames(gelmans_A) = "AB_PSRF"
  colnames(gelmans_S) = "SP_PSRF"


  #visualilze the averaged runs
  #plot(rowMeans(X), typ= "l", main = "Species per site", xlab = "Jmax", ylab = "Species")
  #plot(rowMeans(Xx), typ= "l", main = "Abundance through time in site 12", xlab = "Time-step", ylab = "Species")


  x <- 1:t
  y <- rowMeans(SPsitelist[[2]])
  startpoint = t/2

  #plot(x[startpoint:length(X)] ,y[startpoint:length(X)], typ = "l")
  reg = lm(y[startpoint:length(X)]~x[startpoint:length(X)])
  abline(reg,col="red")

  #### check
  #reg

  ## consolidate results
  S_EQ <- c()
  S_B <- c()
  for (site in 1:sites) {
    y <- rowMeans(SPsitelist[[site]])
    reg = lm(y[startpoint:length(X)] ~ x[startpoint:length(X)])
    S_EQ <- rbind(S_EQ, reg[[1]][[1]])
    S_B <- rbind(S_B, reg[[1]][[2]])
  }
  colnames(S_EQ) = "SP_Equilibrium"
  colnames(S_B) = "SP_Beta_Coef"
  SP_results <- cbind(S_EQ, S_B, gelmans_S)

  A_EQ <- c()
  A_B <- c()
  for (site in 1:sites) {
    y <- rowMeans(ABsitelist[[site]])
    reg = lm(y[startpoint:length(X)] ~ x[startpoint:length(X)])
    A_EQ <- rbind(A_EQ, reg[[1]][[1]])
    A_B <- rbind(A_B, reg[[1]][[2]])
  }

  colnames(A_EQ) = "A_Equilibrium"
  colnames(A_B) = "A_Beta_Coef"
  AB_results <- cbind(A_EQ, A_B, gelmans_A)

  RESULTSH <- list(AB_results, SP_results)


  results=list(diagnostics = RESULTSH,
               SPsitelist = SPsitelist,
               ABsitelist = ABsitelist
  )


  return (results)

}




ss <- analyzeSim(rs, 400)

plot(rowMeans(ss$SPsitelist[[2]]), typ = "l")


