#SAD generators

NoForce(res1$matrix_list, Jmax)

NoForce = function(matlist, Jmax){

  SV <- c()
  for (i in 1:length(matlist)) {
    SV <- append(SV, length(matlist[[i]]))
  }

  S <- SV

  J <- c()
  for ( i in 1:length(Jmax)){
    J <- append(J, Jmax[i] * (S[i] / (100 + S[i])))
  }


oldJ <- c()
for ( i in 1:length(matlist)){
  old <- sum(matlist[[i]])
  oldJ <- append(oldJ, old)
}

newamt <- matlist

for (i in 1:length(matlist)) {
  if (oldJ[[i]]  > J[[i]]) {
    newamt[[i]] <- newamt[[i]] - (oldJ[[i]] - J[[i]]) / S[i]
  } else {
    newamt[[i]] <- newamt[[i]] + (oldJ[[i]] + J[[i]]) / S[i]
  }

}

return(newamt)

}


