#### Display items 1 ####

## choose working directory:
path <- "~/virushostgc/"

setwd(path)

reset <- function() {
  par(mfrow=c(1,1), oma=rep(0,4), mar=rep(0,4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

data <- read.csv("s1.csv"); dim(data)

data$a <- 100*(data$a.gnm)/data$n.gnm; sum(is.na(data$a))
data$c <- 100*(data$c.gnm)/data$n.gnm; sum(is.na(data$c))
data$g <- 100*(data$g.gnm)/data$n.gnm; sum(is.na(data$g))
data$t <- 100*(data$t.gnm)/data$n.gnm; sum(is.na(data$t))

data$gc <- 100*(data$g.gnm+data$c.gnm)/data$n.gnm; sum(is.na(data$a))

data$baltimore <- factor(data$baltimore,
                         levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT"))
sum(is.na(data$baltimore))
# levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT",""))
# levels(data$baltimore)[nlevels(df$baltimore)] <- "unknown"

# display.brewer.all()
# display.brewer.all(colorblindFriendly = TRUE)
display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))

baltimore <- c(RColorBrewer::brewer.pal(12, "Paired")[10],
               RColorBrewer::brewer.pal(12, "Paired")[8],
               RColorBrewer::brewer.pal(12, "Paired")[2],
               RColorBrewer::brewer.pal(12, "Paired")[4],
               RColorBrewer::brewer.pal(12, "Paired")[6],
               RColorBrewer::brewer.pal(8, "Set2")[1],
               RColorBrewer::brewer.pal(12, "Paired")[9])

#####
## Figs.1 (& Figs.S1)

par <- par()
par(family = "serif")
par(oma = c(0,0,0,0))
par(mar = c(3,3,2,3))# + 0.1)
par(mfrow = c(1,1))

# display.brewer.all()
# display.brewer.all(colorblindFriendly = TRUE)
display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))

# png("fig1.png", width = 18, height = 25, units = "cm", res = 300)
# par(mfrow=c(4,2))
## (A)
"a" -> h
"all viruses" -> i
x <- data$gc; gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.02))
mtext("Density (GC)", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385 
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
# round(gc$x[which(s == -2) + 1],0)
round(gc$x[which(s == -2) + 1],1)
# round(gc$x[which(s == -2) + 1],2)
round(range(x),1)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))

## (B)
"b" -> h
"dsDNA" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.01))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385 
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))

## (C)
"c" -> h
"ssDNA" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper("a"),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.03))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))
a <- density(na.omit(data$a[data$baltimore==i])); round(a$x[which(diff(sign(diff(a$y))) == -2) + 1],2)
c <- density(na.omit(data$c[data$baltimore==i])); round(c$x[which(diff(sign(diff(c$y))) == -2) + 1],2)
g <- density(na.omit(data$g[data$baltimore==i])); round(g$x[which(diff(sign(diff(g$y))) == -2) + 1],2)
t <- density(na.omit(data$t[data$baltimore==i])); round(t$x[which(diff(sign(diff(t$y))) == -2) + 1],2)
plot(-x,ylim=c(0,max(c(g$y,c$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,50,by=10),labels = c(seq(0,40,by=10),""))
mtext("50%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*10,0,j*10,1, col = adjustcolor("black",0.5), lty=3)}
lines(g,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),lwd=5)
lines(c,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5),lwd=5)
legend(0,max(c(g$y,c$y)), c("G", "C"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(2, at=seq(0,ceiling(max(c(g$y,c$y))),by = 0.05))
mtext("G/C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
par(new = TRUE)
plot(-x,ylim=c(0,max(c(a$y,t$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
lines(a,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),lwd=5)
# lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5),lwd=5)
lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5),lwd=5)
# legend(40,max(c(a$y,t$y)), c("A", "U"), bg = adjustcolor("white", 0),
legend(40,max(c(a$y,t$y)), c("A", "T"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),
               # adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5)),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(4, at=seq(0,ceiling(max(c(a$y,t$y))),by = 0.04))
# mtext("A/U density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")
mtext("A/T density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")

## (D)
"d" -> h
"dsRNA" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.02))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))

## (E)
"e" -> h
"+ssRNA" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper("b"),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.02))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))
a <- density(na.omit(data$a[data$baltimore==i])); round(a$x[which(diff(sign(diff(a$y))) == -2) + 1],2)
c <- density(na.omit(data$c[data$baltimore==i])); round(c$x[which(diff(sign(diff(c$y))) == -2) + 1],2)
g <- density(na.omit(data$g[data$baltimore==i])); round(g$x[which(diff(sign(diff(g$y))) == -2) + 1],2)
t <- density(na.omit(data$t[data$baltimore==i])); round(t$x[which(diff(sign(diff(t$y))) == -2) + 1],2)
plot(-x,ylim=c(0,max(c(g$y,c$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,50,by=10),labels = c(seq(0,40,by=10),""))
mtext("50%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*10,0,j*10,1, col = adjustcolor("black",0.5), lty=3)}
lines(g,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),lwd=5)
lines(c,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5),lwd=5)
legend(0,max(c(g$y,c$y)), c("G", "C"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(2, at=seq(0,ceiling(max(c(g$y,c$y))),by = 0.04))
mtext("G/C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
par(new = TRUE)
plot(-x,ylim=c(0,max(c(a$y,t$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
lines(a,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),lwd=5)
lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5),lwd=5)
# lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5),lwd=5)
legend(40,max(c(a$y,t$y)), c("A", "U"), bg = adjustcolor("white", 0),
       # legend(40,max(c(a$y,t$y)), c("A", "T"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5)),
       # adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(4, at=seq(0,ceiling(max(c(a$y,t$y))),by = 0.03))
mtext("A/U density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")
# mtext("A/T density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")

## (F)
"f" -> h
"-ssRNA" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper("c"),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.02))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))
a <- density(na.omit(data$a[data$baltimore==i])); round(a$x[which(diff(sign(diff(a$y))) == -2) + 1],2)
c <- density(na.omit(data$c[data$baltimore==i])); round(c$x[which(diff(sign(diff(c$y))) == -2) + 1],2)
g <- density(na.omit(data$g[data$baltimore==i])); round(g$x[which(diff(sign(diff(g$y))) == -2) + 1],2)
t <- density(na.omit(data$t[data$baltimore==i])); round(t$x[which(diff(sign(diff(t$y))) == -2) + 1],2)
plot(-x,ylim=c(0,max(c(g$y,c$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,50,by=10),labels = c(seq(0,40,by=10),""))
mtext("50%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*10,0,j*10,1, col = adjustcolor("black",0.5), lty=3)}
lines(g,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),lwd=5)
lines(c,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5),lwd=5)
legend(0,max(c(g$y,c$y)), c("G", "C"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(2, at=seq(0,ceiling(max(c(g$y,c$y))),by = 0.04))
mtext("G/C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
par(new = TRUE)
plot(-x,ylim=c(0,max(c(a$y,t$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
lines(a,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),lwd=5)
lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5),lwd=5)
# lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5),lwd=5)
legend(40,max(c(a$y,t$y)), c("A", "U"), bg = adjustcolor("white", 0),
       # legend(40,max(c(a$y,t$y)), c("A", "T"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5)),
       # adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(4, at=seq(0,ceiling(max(c(a$y,t$y))),by = 0.04))
mtext("A/U density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")
# mtext("A/T density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")

## (G)
"g" -> h
"+ssRNA-RT" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper("d"),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.02))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))
a <- density(na.omit(data$a[data$baltimore==i])); round(a$x[which(diff(sign(diff(a$y))) == -2) + 1],2)
c <- density(na.omit(data$c[data$baltimore==i])); round(c$x[which(diff(sign(diff(c$y))) == -2) + 1],2)
g <- density(na.omit(data$g[data$baltimore==i])); round(g$x[which(diff(sign(diff(g$y))) == -2) + 1],2)
t <- density(na.omit(data$t[data$baltimore==i])); round(t$x[which(diff(sign(diff(t$y))) == -2) + 1],2)
plot(-x,ylim=c(0,max(c(g$y,c$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,50,by=10),labels = c(seq(0,40,by=10),""))
mtext("50%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*10,0,j*10,1, col = adjustcolor("black",0.5), lty=3)}
lines(g,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),lwd=5)
lines(c,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5),lwd=5)
legend(0,max(c(g$y,c$y)), c("G", "C"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[6], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[2], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(2, at=seq(0,ceiling(max(c(g$y,c$y))),by = 0.04))
mtext("G/C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
par(new = TRUE)
plot(-x,ylim=c(0,max(c(a$y,t$y))),col="white",lwd=3,xlim=c(0,50),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
lines(a,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),lwd=5)
lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5),lwd=5)
# lines(t,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5),lwd=5)
legend(40,max(c(a$y,t$y)), c("A", "U"), bg = adjustcolor("white", 0),
       # legend(40,max(c(a$y,t$y)), c("A", "T"), bg = adjustcolor("white", 0),
       col = c(adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[4], 0.5),
               adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[8], 0.5)),
       # adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[12], 0.5)),
       lwd = 5, lty = 1, box.lty = 0)
axis(4, at=seq(0,ceiling(max(c(a$y,t$y))),by = 0.04))
mtext("A/U density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")
# mtext("A/T density", side = 4, line = 2, cex.lab = 1, las = 0, col = "black")

## (H)
"h" -> h
"dsDNA-RT" -> i
x <- na.omit(data$gc[data$baltimore==i]); gc <- density(x)
plot(gc,col="white",lwd=3,xlim=c(0,100),xaxt="n",yaxt="n",main="",xlab="",ylab="",zero.line=FALSE)
axis(1, at=seq(0,100,by=20),labels = c(seq(0,80,by=20),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
title(paste0("(",toupper(h),")"),adj=0, cex.main=1.5)
title(i,font.main=1, cex.main=1)
title(paste0("n = ",length(x)),adj=1,font.main=1, cex.main=1)
for(j in 0:5){segments(j*20,0,j*20,1, col = adjustcolor("black",0.5), lty=3)}
lines(gc,col=adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),lwd=5)
legend(0,max(gc$y), "GC", bg = adjustcolor("white", 0),
       col = adjustcolor(RColorBrewer::brewer.pal(12, "Paired")[10], 0.5),
       lwd = 5, lty = 1, box.lty = 0)
# mtext("GC (%)", side=1, line=2, cex.lab=1,las=0, col="black")
# mtext("Density", side=2, line=2, cex.lab=1,las=0, col="black")
axis(2, at=seq(0,ceiling(max(gc$y)),by=0.04))
mtext("G+C density", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# d$x[pastecs::turnpoints(ts(d$y))$peaks]
## https://gist.github.com/ramhiser/5316385
s <- diff(sign(diff(gc$y)))
print(paste0(i, ":"))
round(gc$x[which(s == -2) + 1],0)
# round(gc$x[which(s == -2) + 1],1)
round(gc$x[which(s == -2) + 1],2)
range(x)
quantile(x,c(0.025,0.975))
quantile(x,c(0.25,0.75))

par(par)
palette(pal)
