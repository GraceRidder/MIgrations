
#organizing files
Jmax = c(800, 1200, 1600, 2000)
sym = c(-2, -4, -6, -8 ,-10, -20)
model.directory = "~/Desktop/ETBDsim"
fileindex

fileindex <- c()
for (i in Jmax) {
  filename = file.path(
    model.directory,
    paste0(
      "flongLG_",
      as.character(Jmax),
      "data.",
      "RData"
    )
  )
  fileindex <- append(fileindex, filename)
}

fileindex<- fileindex[1:4]

fileindex <- "~/Desktop/SYM1200/LG_depCHANGESP_1600.conEX_RES.RData"
fileindex
name <- c()
finaldatalist <- list()
tiplist <- list()
afterlist <- list()
beforelist <- list()
fulllist <- list()
for (o in 1:length(fileindex)){
load(fileindex[o])
# name <-  append(name, paste("LG7", Jmax[o], "final", sep = ""))
 finaldatalist[[o]] <- finaldata
 tiplist[[o]] <- TIPS
 afterlist[[o]] <- AFTER
 beforelist[[o]] <- BEFORE
 fulllist[[o]] <- FULLTREE

}


plot(finaldatalist[[4]][,9], typ = "l", xlab = "speciation rate", ylab = "Time to equilibroum", col = "red", main = "time to equilibrium")
lines(finaldatalist[[3]][,9], typ = "l", col = "blue")
lines(finaldatalist[[2]][,9], typ = "l", col = "orange")
lines(finaldatalist[[1]][,9], typ = "l", col = "green")
legend("topright", title = "Jmax", c("1600", "1200", "1000", "800"),
       fill = c("red", "blue", "orange", "green"))



plot(finaldatalist[[4]][,1], typ = "l", xlab = "speciation rate", ylab = "Equilibroum", col = "red", main = "Equilibrium")
lines(finaldatalist[[3]][,1], typ = "l", col = "blue")
lines(finaldatalist[[2]][,1], typ = "l", col = "orange")
lines(finaldatalist[[1]][,1], typ = "l", col = "green")
legend("topleft", title = "Jmax", c("1600", "1200", "1000", "800"),
       fill = c("red", "blue", "orange", "green"))





eqs <- cbind(finaldatalist[[1]][1], finaldatalist[[2]][1], finaldatalist[[3]][1], finaldatalist[[4]][1])


plot(as.numeric(eqs[3,]), typ)



plot(as.numeric(eqs[6,]), typ = "l", xlab = "Jmax", ylab = "Equilibroum", col = "red", main = "Equilibrium", ylim = c(0,500))
lines(as.numeric(eqs[5,]), typ = "l", col = "blue")
lines(as.numeric(eqs[4,]), typ = "l", col = "orange")
lines(as.numeric(eqs[3,]), typ = "l", col = "green")
legend("topleft", title = "Speciation", c("1600", "1200", "1000", "800"),
       fill = c("red", "blue", "orange", "green"))


S <- c(0, 20, 300, 487.9687)
plot(S, typ = "l", xlab = "Jmax", ylab = "Equilibroum", col = "white", main = "Equilibrium", log = "xy")
lines(as.numeric(eqs[6,]), typ = "l", col = "red", log = "xy")
lines(as.numeric(eqs[3,]), typ = "l", col = "blue", log = "xy")


S <- c(0, 20, 300, 487.9687)




finaldata

# 1 is node numenr 2 is gamma 3 is beta
X <- 3
o  <- 1


{
e <- as.numeric(tiplist[[o]][X,])
w <- as.numeric(afterlist[[o]][X,])
q <- as.numeric(beforelist[[o]][X,])
s <- as.numeric(fulllist[[o]][X,])

ee <- as.numeric(tiplist[[o]][X+3,])
ww <- as.numeric(afterlist[[o]][X+3,])
qq <- as.numeric(beforelist[[o]][X+3,])
ss <- as.numeric(fulllist[[o]][X+5,])

eee <- as.numeric(tiplist[[o]][X+6,])
www<- as.numeric(afterlist[[o]][X+6,])
qqq <- as.numeric(beforelist[[o]][X+6,])
sss <- as.numeric(fulllist[[o]][X+10,])


eeee <- as.numeric(tiplist[[o]][X+9,])
wwww <- as.numeric(afterlist[[o]][X+9,])
qqqq <- as.numeric(beforelist[[o]][X+9,])
ssss <- as.numeric(fulllist[[o]][X+15,])

eeeee <- as.numeric(tiplist[[o]][X+12,])
wwwww <- as.numeric(afterlist[[o]][X+12,])
qqqqq <- as.numeric(beforelist[[o]][X+12,])
sssss <- as.numeric(fulllist[[o]][X+20,])

eeeeee <- as.numeric(tiplist[[o]][X+15,])
wwwwww <- as.numeric(afterlist[[o]][X+15,])
qqqqqq <- as.numeric(beforelist[[o]][X+15,])
ssssss <- as.numeric(fulllist[[o]][X+25,])

}


boxplot(e,ee, eee, eeee,eeeee, eeeeee)
boxplot(w, ww, www, wwww, wwwww, wwwwww)
boxplot(s, ss, sss, ssss, sssss, ssssss)

j2000N <- c(mean(na.omit(q)), mean(na.omit(qq)),mean(na.omit(qqq)),mean(na.omit(qqqq)),mean(na.omit(qqqq)),mean(na.omit(qqqqq)))
j2000Y <- c(mean(na.omit(q)), mean(na.omit(qq)),mean(na.omit(qqq)),mean(na.omit(qqqq)),mean(na.omit(qqqq)),mean(na.omit(qqqqq)))

j1600N <- c(mean(na.omit(w)), mean(na.omit(ww)),mean(na.omit(www)),mean(na.omit(wwww)),mean(na.omit(wwww)),mean(na.omit(wwwww)))
j1600Y <-c(mean(na.omit(w)), mean(na.omit(ww)),mean(na.omit(www)),mean(na.omit(wwww)),mean(na.omit(wwww)),mean(na.omit(wwwww)))

j1200N  <-c(mean(na.omit(e)), mean(na.omit(ee)),mean(na.omit(eee)),mean(na.omit(eeee)),mean(na.omit(eeeee)),mean(na.omit(eeeeee)))
j1200Y <- c(mean(na.omit(e)), mean(na.omit(ee)),mean(na.omit(eee)),mean(na.omit(eeee)),mean(na.omit(eeeee)),mean(na.omit(eeeeee)))

j800N<- c(mean(na.omit(s)), mean(na.omit(ss)),mean(na.omit(sss)),mean(na.omit(ssss)),mean(na.omit(sssss)),mean(na.omit(ssssss)))
j800Y <- c(mean(na.omit(s)), mean(na.omit(ss)),mean(na.omit(sss)),mean(na.omit(ssss)),mean(na.omit(sssss)),mean(na.omit(ssssss)))

finaldatalist[[1]]$exrate
plot(-(finaldatalist[[1]]$exrate), c(1, 2, 3, 4, 5,6), typ = "l")

lines(finaldatalist[[1]]$exrate)

plot(-(finaldatalist[[1]]$exrate), finaldatalist[[1]]$SPEQ, typ ="l"
     , main = "found extinction rate",
     xlab = "found extinction rate",
     ylab = "species richness at EQ")
finaldatalist[[1]]$SPEQ

plot(finaldatalist[[1]]$exrate, finaldatalist[[1]]$FTY, ylim = c(-10,6))
lines(finaldatalist[[1]]$exrate, finaldatalist[[1]]$tipsY)
lines(finaldatalist[[1]]$exrate, finaldatalist[[1]]$beforeEQY)
lines(finaldatalist[[1]]$exrate, finaldatalist[[1]]$afterEQY)

plot(finaldatalist[[1]]$FTY, typ = "l", ylim = c(-10,6))
lines(finaldatalist[[1]]$tipsY)
lines(finaldatalist[[1]]$beforeEQY)
lines(finaldatalist[[1]]$afterEQY)

plot(j2000N, j2000Y, typ = "l", col = "red",
     main = "gamma before EQ", xlab = "node #",
     ylab = "gamma", xlim = c(0,20), ylim = c(-1,50))
lines(j1600N, j1600Y, col = "blue")
lines(j1200N, j1200Y, col = "darkgreen")
lines(j800N, j800Y, col = "black")
legend("topleft", title = "jmax", c("2000", "1600", "1200", "800"),
       fill = c("red", "blue", "darkgreen", "black")
       , cex = .6)


plot(j2000Y, typ = "l", col = "red",
     main = "gamma 800", xlab = "high to low extionction", ylab = "gamma", xlim = c(1,7),
     ylim = c(-10,5.5))
lines(j1600Y, col = "blue")
lines(j1200Y, col = "darkgreen")
lines(j800Y, col = "black")
legend("bottomleft", title = "trims", c("before", "after", "tips", "full"),
       fill = c("red", "blue", "darkgreen", "black")
       , cex = .6)



plot(j2000N, typ = "l", col = "red",
     main = "node# 2000", xlab = "high to low extionction", ylab = "node#", xlim = c(1,7),
     ylim = c(0,150))
lines(j1600N, col = "blue")
lines(j1200N, col = "darkgreen")
lines(j800N, col = "black")
legend("topleft", title = "trims", c("before", "after", "tips", "full"),
       fill = c("red", "blue", "darkgreen", "black")
       , cex = .6)


#plot(finaldata$FTY, typ = "l", col = "blue")


#boxplot(t, a, b, f, tt, aa, bb, ff, ttt, aaa, bbb, fff, tttt, aaaa, bbbb, ffff,
#boxplot(l, p, u, j, ll, pp, uu, jj, lll, ppp, uuu, jjj, llll, pppp, uuuu, jjjj,
boxplot(
        e, w, q, s,
        ee, ww, qq, ss,
        eee, www, qqq, sss,
        eeee, wwww, qqqq, ssss,
        eeeee, wwwww, qqqqq, sssss,
        eeeeee, wwwwww, qqqqqq, ssssss,
        # eeeeee, wwwwww, qqqqqq, ssssss,
        # eeeee, wwwww, qqqqq, sssss,
        # eeee, wwww, qqqq, ssss,
        # eee, www, qqq, sss,
        # ee, ww, qq, ss,
        # e, w, q, s,

        main = "beta Jmax (1600)
        constant extinction
        dependent speciation",
      #  names = c(".3", ".3", ".3", ".3",".4", ".4", ".4", ".4", ".5", ".5", ".5", ".5", ".6", ".6", ".6", ".6"),
        #las = 2,
        col = c("orange","pink", "brown", "tan"),
        horizontal = F,
        notch = F,
        ylab = "beta",
        xlab = "high to low constant extinction"
       # ylim = c(-5,11)
)

legend("topleft", c("tip", "aft", "bf", "full"),
       fill = c("orange", "pink", "brown", "tan")
       , cex = .5)



install.packages("ETBDsim")
library(ETBDsim)

res1 = ETBDspaceSYM(
  t = 100,
  psymp = 0,
  SRS = F,
  GEO = F,
  LOGNRM = T,
  watchgrow = F,
  SADmarg = .1,
  siteN = 1,
  JmaxV = 1600,
  exparm = -8,
  NegExpEx = T,
  constantEX = .005
)

res <- res1

X <- c()
Xx <- c()

  tpac <- c()
  apac <- c()
  for (time in 1:length(res$mig)) {
    tpac <- append(tpac, length(unlist(res$mig[time][[1]][1])))
    apac <- append(apac,  sum(unlist(res$mig[time][[1]][1])))
  }



plot(tpac, typ = "l", main="species richness through time", ylab = " species richness", xlab = "time")

plot(apac, typ = "l")


plot(res$trees[[1000]], cex = .2)

library("treebalance")
library("geiger")
library("apTreeshape")
library("paleotree")

Xtree <- drop.extinct(res$tree)

plot(Xtree, cex = .7)
tim


gammaStat(Xtree)
maxlik.betasplit(Xtree,confidence.interval="none")

TrueAfter <- drop.extinct(res$trees[[950]])
plot(TrueAfter, cex = .7, main = "Tree saved within Sim T = 500")


tim = 1000 - 950

troll <- timeSliceTree(res$tree, tim, plot = F, drop.extinct = T)
plot(troll, cex = .7, main = "Tree sliced to time 500")

gammaStat(troll)
gammaStat(TrueAfter)

d <- 0
dd <- c()
for ( i in 1:100){
  d = d +100
  dd <- append(dd, d)
}
length(dd)

gam <- c()
for (i in 1:length(dd)){
TrueAfter <- drop.extinct(res$trees[[dd[i]]])
gam <- append(gam, gammaStat(TrueAfter))
}

plot(gam, type = "l", main = "Gamma through time", xlab = "time")



library(parallel)


Jmaxvalues <- c(1600)
speciations <- c(2)

for (v in speciations){
  for (i in Jmaxvalues){
    # #run the simulation multiple times in parallel
    trials = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

    repsim <- function(trials){
      reslist <- list()


      for (u in 1:length(trials)) {
        res1 = ETBDspaceSYM(
          t = 500,
          psymp = 0,
          SRS = F,
          GEO = F,
          LOGNRM = T,
          watchgrow = F,
          SADmarg = .1,
          siteN = 1,
          JmaxV = 1600,
          exparm = -6,
          NegExpEx = T,
          constantEX = 0,
          ExpSpParm = 2,
          ExpSp = T,
          SPgrow = 0,
          splitparm = .25,
        )

        reslist[[u]] <- res1

      }

      #result object
      results=list( reslist = reslist
      )

      #result
      return(results)

    }


    rs= mcmapply(repsim, trials, SIMPLIFY = F, mc.cores = 1)

    filename <- paste("LG_EX7_",i,".SP", v, ".RData", sep = "")

    save(rs, file = filename)

  }

}

wapper <- c()
for (i in 1:10){
apple <- maxlik.betasplit(rs[[i]]$reslist[[1]]$tree)
app <- apple$max_lik
wapper <- append(wapper, app)
}



Lg100 <- wapper



boxplot(Lg100, Lg50, Lg25, Lg05, Lg005)









