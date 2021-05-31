#####
## Table 4 and Figures 3 (A-C + legend)

par <- par()
pal <- palette()

## you may need:
# install.packages("RColorBrewer")

## choose working directory:
path <- "~/virushostgc/"

setwd(path)

data <- read.csv("data4.csv"); dim(data)

data$prokaryotic <- as.factor(data$prokaryotic)
data$eukaryotic <- as.factor(data$eukaryotic)
# data$subviral <- as.factor(data$subviral.)
data$host <- as.factor(paste0(data$prokaryotic, data$eukaryotic))

# RColorBrewer::display.brewer.all()
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))
# levels = c("", "animals", "archaea", "bacteria", "fungi", "plants", "protists")
host <- c("gray",
          RColorBrewer::brewer.pal(12, "Paired")[6],
          RColorBrewer::brewer.pal(12, "Paired")[10],
          RColorBrewer::brewer.pal(12, "Paired")[2],
          RColorBrewer::brewer.pal(12, "Paired")[12],
          RColorBrewer::brewer.pal(12, "Paired")[4],
          RColorBrewer::brewer.pal(12, "Paired")[8])

#####
## Figs.3

par <- par()
par(family = "serif")
par(oma = c(0,0,0,0))
par(mar = c(3,3,2,3))# + 0.1)
par(mfrow = c(1,1))

palette(adjustcolor(host,0.5))

## (A)
"" -> H
"a" -> h
"gc" -> i
x <- data[,toupper(i)]
y <- data[,i]
z <- complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("host GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext("viral GC", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y ~ x)
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
points(x[z], y[z], pch = 20, col = adjustcolor("black",0.33), cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
title(H,font.main=1,cex.main=1)
title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
(rho <- sprintf(round(cor(x, y, method = "spearman", use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (B)
"" -> H
"b" -> h
"gc" -> i
x <- data[,toupper(i)]
y <- data[,i]
z <- complete.cases(x) & complete.cases(y)
Z <- data$prokaryotic!="" & complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext("host GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("viral G+C", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y[Z] ~ x[Z])
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
points(x[z&data$prokaryotic==""], y[z&data$prokaryotic==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(data$prokaryotic)[-1]){points(x[z&data$prokaryotic==i], y[z&data$prokaryotic==i], pch = 20, col = which(levels(data$prokaryotic)==i), cex=1.5)}
for(i in names(sort(table(data$prokaryotic)[-1],decreasing = TRUE))){points(x[z&data$prokaryotic==i], y[z&data$prokaryotic==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
points(x[z&data$prokaryotic=="archaea"], y[z&data$prokaryotic=="archaea"], pch = 20, col = adjustcolor(which(levels(data$host)=="archaea"),0.33), cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
title(H,font.main=1,cex.main=1)
title(paste0("n = ",format(sum(Z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x[Z], y[Z], use = "complete.obs"), 2), fmt = "%#.2f"))
(rho <- sprintf(round(cor(x[Z], y[Z], method = "spearman", use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (C)
"" -> H
"c" -> h
"gc" -> i
x <- data[,toupper(i)]
y <- data[,i]
z <- complete.cases(x) & complete.cases(y)
Z <- data$eukaryotic!="" & complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("host GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("viral G+C", side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y[Z] ~ x[Z])
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
points(x[z&data$eukaryotic==""], y[z&data$eukaryotic==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(data$eukaryotic)[-1]){points(x[z&data$eukaryotic==i], y[z&data$eukaryotic==i], pch = 20, col = which(levels(data$eukaryotic)==i), cex=1.5)}
for(i in names(sort(table(data$eukaryotic)[-1],decreasing = TRUE))){points(x[z&data$eukaryotic==i], y[z&data$eukaryotic==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
points(x[z&data$eukaryotic=="protists"], y[z&data$eukaryotic=="protists"], pch = 20, col = adjustcolor(which(levels(data$host)=="protists"),0.33), cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
title(H,font.main=1,cex.main=1)
title(paste0("n = ",format(sum(Z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x[Z], y[Z], use = "complete.obs"), 2), fmt = "%#.2f"))
(rho <- sprintf(round(cor(x[Z], y[Z], method = "spearman", use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (legend)
plot.new()
legend("left", legend=levels(data$host)[-1], col=host[-1],
       box.lty=1, pch = 20, pt.cex = 1.5,cex = 1)

par(par)
palette(pal)

#####

codons <- names(data)[87:(87+63)]

data[,toupper(codons)] <- data[,toupper(codons)]/rowSums(data[,toupper(codons)])
data[,codons] <- data[,codons]/rowSums(data[,codons])

data$prokaryotic <- factor(data$prokaryotic)
data$eukaryotic <- factor(data$eukaryotic)
# data$subviral. <- factor(data$subviral.)

dataPk <- NULL
for(j in c("gc","gc1","gc2","gc3")){
  i=toupper(j)
  for(k in levels(data$prokaryotic)){
    if(k == ""){
      dataPk[[paste(j,i,sep="~")]] <- cbind(k="prokaryotic",
                                            # r=cor(data[data$prokaryotic!=k,i],data[data$prokaryotic!=k,j], use="complete.obs"),
                                          rho=cor(data[data$prokaryotic!=k,i], method="spearman",data[data$prokaryotic!=k,j], use="complete.obs"),
                                          n=sum(complete.cases(data[data$prokaryotic!=k,c(i,j)])))
    }else{
      dataPk[[paste(j,i,sep="~")]] <- rbind(dataPk[[paste(j,i,sep="~")]],
                                          cbind(k=k,
                                                # r=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="complete.obs"),
                                                rho=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], method="spearman",use="complete.obs"),
                                                n=sum(complete.cases(data[data$prokaryotic==k,c(i,j)]))))
    }
  }
}
dataPk

dataEk <- NULL
for(j in c("gc","gc1","gc2","gc3")){
  i=toupper(j)
  for(k in levels(data$eukaryotic)){
    if(k == ""){
      dataEk[[paste(j,i,sep="~")]] <- cbind(k="eukaryotic",
                                            # r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="complete.obs"),r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="complete.obs"),
                                            rho=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], method="spearman",use="complete.obs"),
                                            n=sum(complete.cases(data[data$eukaryotic!=k,c(i,j)])))
    }else{
      dataEk[[paste(j,i,sep="~")]] <- rbind(dataEk[[paste(j,i,sep="~")]],
                                          cbind(k=k,
                                                # r=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="complete.obs"),
                                                rho=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], method="spearman",use="complete.obs"),
                                                n=sum(complete.cases(data[data$eukaryotic==k,c(i,j)]))))
    }
  }
}
dataEk

#####
## Table 4

## (left half)
codPk<- NULL
for(j in codons){
  i=toupper(j)
  for(k in levels(data$prokaryotic)){
    if(k == ""){
      codPk[[paste(j,i,sep="~")]] <- cbind(k="prokaryotic",
                                           # r=cor(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], use="complete.obs"),
                                           r=cor(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], use="na.or.complete"),
                                           p=cor.test(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], use="na.or.complete")$p.value,
                                           rho=cor(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], method="spearman", use="na.or.complete"),
                                           p=cor.test(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], method="spearman", use="na.or.complete")$p.value,
                                           R2=summary(lm(data[data$prokaryotic!=k,j]~data[data$prokaryotic!=k,i]))$adj.r.squared#,
                                           # n=sum(complete.cases(data[data$prokaryotic!=k,c(i,j)]))
      )
    }else{
      codPk[[paste(j,i,sep="~")]] <- rbind(codPk[[paste(j,i,sep="~")]],
                                           cbind(k=k,
                                                 # r=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="complete.obs"),
                                                 r=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="na.or.complete"),
                                                 p=cor.test(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="na.or.complete")$p.value,
                                                 rho=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], method="spearman", use="na.or.complete"),
                                                 p=cor.test(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], method="spearman", use="na.or.complete")$p.value,
                                                 R2=summary(lm(data[data$prokaryotic==k,j]~data[data$prokaryotic==k,i]))$adj.r.squared#,
                                                 # n=sum(complete.cases(data[data$prokaryotic==k,c(i,j)]))
                                           ))
    }
  }
}
codpk<-do.call(rbind.data.frame, codPk)
# write.table(codpk[codpk$k=="prokaryotic",],"codPk.csv",quote = F,row.names = T,col.names = T, sep=",")
# write.table(codpk[codpk$k!="prokaryotic",],"codPks.csv",quote = F,row.names = ,col.names = T, sep=",")
codpk1 <- codpk[codpk$k=="prokaryotic",c("rho","R2")]
codpk1 <- apply(apply(codpk1, 2, as.character), 2, as.numeric)
row.names(codpk1) <- gsub("T", "U", toupper(codons))
round(codpk1, 2)

## (right half)
codEk<- NULL
for(j in codons){
  i=toupper(j)
  for(k in levels(data$eukaryotic)){
    if(k == ""){
      codEk[[paste(j,i,sep="~")]] <- cbind(k="eukaryotic",
                                           # r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="complete.obs"),
                                           r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="na.or.complete"),
                                           p=cor.test(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="na.or.complete")$p.value,
                                           rho=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], method="spearman", use="complete.obs"),
                                           p=cor.test(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], method="spearman", use="na.or.complete")$p.value,
                                           R2=summary(lm(data[data$eukaryotic!=k,j]~data[data$eukaryotic!=k,i]))$adj.r.squared#,
                                           # n=sum(complete.cases(data[data$eukaryotic!=k,c(i,j)]))
      )
    }else{
      codEk[[paste(j,i,sep="~")]] <- rbind(codEk[[paste(j,i,sep="~")]],
                                           cbind(k=k,
                                                 # r=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="complete.obs"),
                                                 r=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="na.or.complete"),
                                                 p=cor.test(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="na.or.complete")$p.value,
                                                 rho=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], method="spearman", use="na.or.complete"),
                                                 p=cor.test(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], method="spearman", use="na.or.complete")$p.value,
                                                 R2=summary(lm(data[data$eukaryotic==k,j]~data[data$eukaryotic==k,i]))$adj.r.squared#,
                                                 # n=sum(complete.cases(data[data$eukaryotic==k,c(i,j)])
                                           ))
    }
  }
}
codek<-do.call(rbind.data.frame, codEk)
# write.table(codek[codek$k=="eukaryotic",],"codEk.csv",quote = F,row.names = T,col.names = T, sep=",")
# write.table(codek[codek$k!="eukaryotic",],"codEks.csv",quote = F,row.names = T,col.names = T, sep=",")
codek1 <- codek[codek$k=="eukaryotic",c("rho","R2")]
codek1 <- apply(apply(codek1, 2, as.character), 2, as.numeric)
row.names(codek1) <- gsub("T", "U", toupper(codons))
round(codek1, 2)
