



load('~/Desktop/twosite/DRS_SRS_3000_10_2.RData')


DR30002 <- DRs



DRR3000_3 <- cbind(DR3000,DR30001, DR30002)

dim(DRR3000_2)

plot(rev(rowMeans(DRR3000_3)), typ = 'l', col = 'red')
lines(rev(rowMeans(DRR2000_3)), col = 'blue')

maxx2 <- c()
for (i in 1:length(DRR3000_2[,1])){
maxx2 <-append(maxx2, max(DRR3000_2[i,]))
}

minn2 <- c()
for (i in 1:length(DRR3000_2[,1])){
  minn2 <-append(minn2, min(DRR3000_2[i,]))
}

two <- rowMeans(DRR3000_2)

library(ggplot2)

df <-as.data.frame(cbind(c(1:39)), two)




g <- ggplot(data=df, aes(x=c(1:39), y=three, ymin=rev(minn), ymax=rev(maxx))) +
   geom_line() +
   geom_ribbon(alpha=0.5) +
  scale_x_log10() +
  scale_y_log10()



g2 <- ggplot(data=df, aes(x=c(1:39), y=rev(two), ymin=rev(minn2), ymax=rev(maxx2))) +
  geom_line() +
  geom_ribbon(alpha=0.5) +
  scale_x_log10() +
  scale_y_log10()
g2
g2 + geom_line(aes(y=f), colour="red")

t <- rev(two)

s <- c(NA, NA, NA)
f <- append(s, t[1:36])


x = c()
D = 10
for ( l in 1:39){
  x <- append(x, D)
  D = D+10
}

library(epm)
library(geiger)
library(paleotree)

x = 1:3
DRs <- c()

for (o in 1:length(rs)) {
  tree <- rs[[o]]$reslist[[1]]$tree
  DRstats <- c()
  for (i in 1:length(x)) {
    troll <- timeSliceTree(tree, x[i], plot = F, drop.extinct = T)
    LL <- DRstat(troll)
    DRstats <- append(DRstats, mean(LL))
  }

  DRs <- cbind(DRs, DRstats)

}









fulltrops <- c()
fulltemps <- c()

sum_trop <- c()
sum_temp <- c()

per_temp <- c()
per_trop <- c()

for (b in 1:length(treelist)){



  xtree <- drop.extinct(treelist[[b]])
  tree <- treelist[[b]]
  mig <- miglist[[b]]


  ##will not error
  sites <- c()
  for (j in 2:length(tree$node.label)) {
    t = 1
    while (tree$node.label[j] %!in% c(unlist(row.names(mig[[t]][[1]])), "t001", "t002",
                                      unlist(row.names(mig[[t]][[2]])))) {
      t = t + 1
    }
    sites <- append(sites, t)
  }

  ###NODE LABELS
  #### indexing the site they exist in
  R <- c(.5)
  for (i in 2:length(tree$node.label)) {
    tt <- sites[i-1]
    if (tree$node.label[i] %in% unlist(row.names(mig[[tt]][[1]]))) {
      R <- append(R, 0)
    }
    if (tree$node.label[i] %in% unlist(row.names(mig[[tt]][[2]]))) {
      R <- append(R, 1)
    }

  }


  ###NODE LABELS
  #### indexing the time steps they exist in
  Tsites <- c()
  for (j in 1:length(tree$tip.label)) {
    t = 1
    while (tree$tip.label[j] %!in% c(unlist(row.names(mig[[t]][[1]])), unlist(row.names(mig[[t]][[2]])))) {
      t = t + 1
    }
    Tsites <- append(Tsites, t)
  }


  ###NODE LABELS
  #### indexing the site they exist in
  TR <- c()
  for (i in 1:length(tree$tip.label)) {
    tt <- Tsites[i]
    if (tree$tip.label[i] %in% unlist(row.names(mig[[tt]][[1]]))) {
      TR <- append(TR, 0)
    }
    if (tree$tip.label[i] %in% unlist(row.names(mig[[tt]][[2]]))) {
      TR <- append(TR, 1)
    }

  }


  ####################### try to track in tree

  #library(phytools)
  #library(ggtree)
  #library(dplyr)

  # x <- tree
  # d <- data.frame(label=x$tip.label, var1 = TR)
  # treex <- full_join(x, d, by='label')
  #
  # var1 = R
  # var2 = +(!R)
  # var2[1] <- 0.5
  # ancstats <- as.data.frame(var1)
  # ancstats <- cbind(ancstats, var2)
  # ancstats$node <- 1:tree$Nnode+Ntip(treex)



  nam <- tree$tip.label
  xx <- TR
  names(xx) <- nam


  library(epm)
  #library(phytools)
  #library(geiger)

  LL<- DRstat(tree)
  trait.bi1 <- ifelse(LL > .2, 3, 2)
  df <- data.frame(factor(xx),factor(trait.bi1))
  q <-match(xtree$tip.label,row.names(df))
  trimed <- df[q,]

  #p12 <- ggtree(xtree)
  #p12 <- gheatmap(p12, trimed, offset=.00, width=0.2, colnames = F, font.size=2) +
  #  scale_fill_manual(values=c("0" = "red","1" = "blue", "2" = "blue","3" = "red"))
  #p12

  #### in tropical species, percent of high DR estimates (red things)
  tropic <- c(trimed$factor.xx.=="0")
  tropic <- trimed[tropic,]
  sum(as.numeric(tropic$factor.trait.bi1.)-1)
  length(tropic[,1])
  OPIC <- sum(as.numeric(tropic$factor.trait.bi1.)-1)/length(tropic[,1])

  #### in temlerate species, percent of high DR estimates (red things)
  temp <- c(trimed$factor.xx.=="1")
  temper <- trimed[temp,]
  sum(as.numeric(temper$factor.trait.bi1.)-1)
  EMP <- sum(as.numeric(temper$factor.trait.bi1.)-1)/length(temper[,1])


  per_trop <- append(per_trop, OPIC)
  per_temp <- append(per_temp, EMP)


  df <- data.frame(factor(xx),LL)
  q <-match(xtree$tip.label,row.names(df))
  trimed <- df[q,]

  #### in tropical species, percent of high DR estimates (red things)
  tropic <- c(trimed$factor.xx.=="0")
  tropic <- trimed[tropic,]
  SUMOPIC <-   sum(tropic$LL)/length(tropic$LL)


  temp <- c(trimed$factor.xx.=="1")
  temper <- trimed[temp,]
  SUMEMP <-   sum(temper$LL)/length(temper$LL)

  sum_trop <- append(sum_trop, SUMOPIC)
  sum_temp <- append(sum_temp, SUMEMP)


}

boxplot(per_trop, per_temp, names =c("tropics", "temperate"), main = "threshold > 0.2")

boxplot(sum_trop, sum_temp, names =c("tropics", "temperate"), main = "average DR")


fulltrops <- append(fulltrops, per_trop)
fulltemps <- append(fulltemps, per_temp)


boxplot(fulltrops, fulltemps, names =c("tropics", "temperate"))
t.test(fulltrops, fulltemps)


library(phytools)
library(ggtree)
library(dplyr)
library(geiger)
install.packages('ETBDsim')

library("ETBDsim")
library("profvis")

treelist <- list()
miglist <- list()
for (i in 1:10){
 profvis({res1 = ETBD_migrateSYM(
    t =30,
    DIST = "SRS",     ### NO, GEO, SRS, NORM
    watchgrow = F,
    SADmarg = .1,
    siteN = 2,
    JmaxV = c(1500, 2500),
    NegExpEx = T,   ###dependent extinction
    exparm = -.7,
    psymp = .2,
    ExpSp = F,      ###dependent speciation
    ExpSpParm = 2,
    constantEX = 0,
    SPgrow = 0,
    splitparm = .3, ### splitting
    bud = T,
    split = F,
    migprob = .0,
    exparm2 = -.3
  )}
)


  treelist[[i]] <- res1$tree
  miglist[[i]] <- res1$mig
}


library(geiger)




plot(res1$tree, show.tip.label = F)
axis(1)
xtree <- drop.extinct(res1$tree)
plot(xtree, cex = .01)


f <- c()
for ( w in 1:length(treelist)){
  tree <- treelist[[w]]
g = c()
b = c()
x = 50:90
for ( i in 1:length(x)){

 troll <- try(timeSliceTree(tree, x[i], plot = F, drop.extinct = T))
 g <- append(g, gammaStat(troll))
 b <- append(b, beta_statistic(troll))

}
f <- cbind(f, g)
}



plot(rev(rowMeans(f)), typ = "l", col = "red", main = "gamma (10 reps")
abline(h=0)

plot(rev(b), typ = "l", col = "blue", main = "beta")



gammaStat(drop.extinct(res1$tree))

install.packages('treestats')
library('treestats')

beta_statistic(drop.extinct(res1$tree))



x <- c()
b <- c()
for (i in 1:length(treelist)){
  x <- append (x, (gammaStat(drop.extinct(treelist[[i]]))))
  b <- append(b, beta_statistic(drop.extinct(treelist[[i]])))
}

 hist(x)
hist(b)
mean(x)
mean(b)



migs3000 <- miglist
trees3000 <- treelist


X1 <- c()
X2 <- c()
for (j in 1:length(miglist)){
  mig <- miglist[[j]]
  X1 <- c()
for (i in 1:length(mig)){
  X1 <- append(X1,length(res1$mig[i][[1]][[1]]))
}
  X2 <- cbind(X2,X1)
}


plot(rowMeans(X), typ = 'l', col = "blue")
lines(rowMeans(X2), typ = 'l', col = "red")



library(epm)

LL <- c()
for ( j in 1:length(trees2000)){
  tree <- trees2000[[j]]
  LL<- append(LL, DRstat(tree))
}

length(LL)



LL2<- DRstat(res1$tree)

plot(density(LL3))
lines(density(LL2))





system.time({

  res1 = ETBD_migrateSYM.NE(
    t =30,
    DIST = "SRS",     ### NO, GEO, SRS, NORM
    watchgrow = F,
    SADmarg = .1,
    siteN = 2,
    JmaxV = c(2000, 2000),
    NegExpEx = T,   ###dependent extinction
    exparm = c(-.7,-.7) ,
    exparm2 = c(-.35,-.3),

    psymp = c(.3, .3),  ### sympatric speciation
    ExpSp = F,      ### dependent speciation
    ExpSpParm = 2,
    constantEX = 0,
    SPgrow = 0,
    splitparm = .3, ### splitting
    bud = T,
    split = F,
    migprob1 = 0,  ## probability of site one leaving
    migprob2 = 0 ## probability of site two leaving


  )

})


myTree <- ape::read.tree(text = res1$tree)
myTree$Nnode
plot(myTree, cex = .2)

length(res1$matrix_list[[1]])
length(res1$matrix_list[[2]])


plot(dropExtinct(myTree), cex= .2)



system.time({

  res1 = ETBD_migrateSYM.NE(
    t =80,
    DIST = "SRS",     ### NO, GEO, SRS, NORM
    watchgrow = F,
    SADmarg = .1,
    siteN = 2,
    JmaxV = c(1500, 1500),
    NegExpEx = T,   ###dependent extinction
    exparm = c(-.7,-.7) ,
    exparm2 = c(-.35,-.3),

    psymp = c(.3, .3),  ### sympatric speciation
    ExpSp = F,      ### dependent speciation
    ExpSpParm = 2,
    constantEX = 0,
    SPgrow = 0,
    splitparm = .3, ### splitting
    bud = T,
    split = F,
    migprob1 = 0,  ## probability of site one leaving
    migprob2 = 0 ## probability of site two leaving

  )

})


myTree <- ape::read.tree(text = res1$tree)
myTree$Nnode
plot(myTree, cex = .2)

length(res1$matrix_list[[1]])
length(res1$matrix_list[[2]])


plot(dropExtinct(myTree), cex= .2)

X1 <- c()
X2 <- c()
for (i in 1:length(res1$mig)){
  X1 <- append(X1,length(res1$mig[i][[1]][[1]]))
  X2 <- append(X2,length(res1$mig[i][[1]][[2]]))
}

plot(X2, typ = "l", col = "red", ylim =c(0,1000), main = "species richness")
lines(X1, typ = "l", col = "blue")


x1 <- c()
x2 <- c()
for (i in 1:length(res1$mig)){
  x1 <- append(x1,sum(res1$mig[i][[1]][[1]][,1]))
  x2 <- append(x2,sum(res1$mig[i][[1]][[2]][,1]))
}

dev.off()
plot(x2, typ = "l", main = "abundances", col = 'red')
lines(x1, typ = "l", col = "blue")


hist(res1$matrix_list[[1]], col = "red")
hist(res1$matrix_list[[2]], col = "blue")


####small helper funciton
'%!in%' <- function(x,y)!('%in%'(x,y))


###NODE LABELS
#### indexing the time steps they exist in

##will error if the first site speciates on the first step
sites <- c()
for (j in 2:length(res1$tree$node.label)) {
  t = 1
  while (res1$tree$node.label[j] %!in% c(unlist(row.names(res1$mig[[t]][[1]])), "t001", "t002",
                                         unlist(row.names(res1$mig[[t]][[2]]))))
                                       #  unlist(row.names(res1$mig[[t]][[3]]))))
    {
    t = t + 1
  }
  sites <- append(sites, t)
}



###NODE LABELS
#### indexing the site they exist in
R <- c(.5)
for (i in 2:length(res1$tree$node.label)) {
  tt <- sites[i-1]
  if (res1$tree$node.label[i] %in% unlist(row.names(res1$mig[[tt]][[1]]))) {
    R <- append(R, 0)
  }
  if (res1$tree$node.label[i] %in% unlist(row.names(res1$mig[[tt]][[2]]))) {
    R <- append(R, 1)
  }
 # if (res1$tree$node.label[i] %in% unlist(row.names(res1$mig[[tt]][[3]]))) {
   # R <- append(R, 2)
 # }

}



###NODE LABELS
#### indexing the time steps they exist in
Tsites <- c()
for (j in 1:length(res1$tree$tip.label)) {
  t = 1
  while (res1$tree$tip.label[j] %!in% c(unlist(row.names(res1$mig[[t]][[1]])), unlist(row.names(res1$mig[[t]][[2]])))) {
    t = t + 1
  }
  Tsites <- append(Tsites, t)
}


###NODE LABELS
#### indexing the site they exist in
TR <- c()
for (i in 1:length(res1$tree$tip.label)) {
  tt <- Tsites[i]
  if (res1$tree$tip.label[i] %in% unlist(row.names(res1$mig[[tt]][[1]]))) {
    TR <- append(TR, 0)
  }
  if (res1$tree$tip.label[i] %in% unlist(row.names(res1$mig[[tt]][[2]]))) {
    TR <- append(TR, 1)
  }

}


####################### try to track in tree
library(dplyr)

x <- res1$tree
d <- data.frame(label=x$tip.label, var1 = TR)
tree <- full_join(x, d, by='label')

var1 = R
var2 = +(!R)
var2[1] <- 0.5
ancstats <- as.data.frame(var1)
ancstats <- cbind(ancstats, var2)
ancstats$node <- 1:res1$tree$Nnode+Ntip(res1$tree)

nam <- res1$tree$tip.label
xx <- TR
names(xx) <- nam


p12 <- ggtree(res1$tree)

# add heatmap
p12 <- gheatmap(p12, data.frame(factor(xx)), offset=0, width=0.2, colnames = F, font.size=2) +
  scale_fill_manual(values=c("red","blue"))



#pies <- nodepie(ancstats, cols = 1:2)
#pies <- lapply(pies, function(g) g+scale_fill_manual(values=c("blue","red")))
#p2 <- p12 + geom_inset(pies, width = .05, height = .05)

#plot(p2, guides='collect', tag_levels='A')


library(epm)
LL<- DRstat(res1$tree)

trait.bi1 <- ifelse(LL > .2, 3, 2)

df <- data.frame(factor(xx),factor(trait.bi1))

p12 <- gheatmap(p12, df, offset=.00, width=0.2, colnames = F, font.size=2) +
  #  scale_fill_manual(values=c("red","blue"))
  scale_fill_manual(values=c("0" = "red","1" = "blue", "2" = "pink","3" = "lightblue"))

plot(p12)



q <-match(xtree$tip.label,row.names(df))
trimed <- df[q,]



p12 <- ggtree(xtree)

p12 <- gheatmap(p12, trimed, offset=.00, width=0.2, colnames = F, font.size=2) +
  #  scale_fill_manual(values=c("red","blue"))
  scale_fill_manual(values=c("0" = "red","1" = "blue", "2" = "blue","3" = "red"))
p12




tropic <- c(trimed$factor.xx.=="0")
tropic <- trimed[tropic,]
sum(as.numeric(tropic$factor.trait.bi1.)-1)
length(tropic[,1])
sum(as.numeric(tropic$factor.trait.bi1.)-1)/length(tropic[,1])


temp <- c(trimed$factor.xx.=="1")
temper <- trimed[temp,]
sum(as.numeric(temper$factor.trait.bi1.)-1)
sum(as.numeric(temper$factor.trait.bi1.)-1)/length(temper[,1])



plot(density(tropic$LL), col = "red")
lines(density(temper$LL), col = "blue")







