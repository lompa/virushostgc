#### Display items 3 ####

## choose working directory:
path <- "~/virushostgc/"

setwd(path)

data <- read.csv("s3.csv"); dim(data)


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

data$host <- factor(rep("",nrow(data)),
                    levels = c("", "animals", "archaea", "bacteria", "fungi", "plants", "protists"))

sort(table(data$host_cell))
# sort(table(data$host_type))
data$host[data$host_cell=="animal; fungus; plant"] <- "fungi"
data$host[data$host_cell=="fungus; plant"] <- "fungi"
data$host[data$host_cell=="animal; protist"] <- "protists"
data$host[data$host_cell=="animal; plant"] <- "plants"
data$host[data$host_cell=="protist"] <- "protists"
data$host[data$host_cell=="fungus"] <- "fungi"
data$host[data$host_cell=="plant"] <- "plants"
data$host[data$host_cell=="animal"] <- "animals"
# data$host[data$host_cell=="prokaryote"] <- "prokaryote"
sort(table(data[data$host_cell=="prokaryote",]$host_type))
data$host[data$host_type=="archaea"] <- "archaea"
data$host[data$host_type=="bacteria"] <- "bacteria"
# data[data$host_cell=="prokaryote"&data$host_type=="",]
levels(data$host)

# display.brewer.all()
# display.brewer.all(colorblindFriendly = TRUE)
display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))

host <- c("gray",
          RColorBrewer::brewer.pal(12, "Paired")[6],
          RColorBrewer::brewer.pal(12, "Paired")[10],
          RColorBrewer::brewer.pal(12, "Paired")[2],
          RColorBrewer::brewer.pal(12, "Paired")[12],
          RColorBrewer::brewer.pal(12, "Paired")[4],
          RColorBrewer::brewer.pal(12, "Paired")[8])

# display.brewer.all()
# display.brewer.all(colorblindFriendly = TRUE)
display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))
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
x <- vh[,toupper(i)]
y <- vh[,i]
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
title(paste0("n = ",sum(z)),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adj.R = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (B)
"" -> H
"b" -> h
"gc" -> i
x <- vh[,toupper(i)]
y <- vh[,i]
z <- complete.cases(x) & complete.cases(y)
Z <- vh$Prokaryotes!="" & complete.cases(x) & complete.cases(y)
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
points(x[z&vh$Prokaryotes==""], y[z&vh$Prokaryotes==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(vh$Prokaryotes)[-1]){points(x[z&vh$Prokaryotes==i], y[z&vh$Prokaryotes==i], pch = 20, col = which(levels(vh$Prokaryotes)==i), cex=1.5)}
for(i in names(sort(table(vh$Prokaryotes)[-1],decreasing = TRUE))){points(x[z&vh$Prokaryotes==i], y[z&vh$Prokaryotes==i], pch = 20, col = which(levels(vh$Prokaryotes)==i), cex=1.5)}
points(x[z&vh$Prokaryotes=="archaea"], y[z&vh$Prokaryotes=="archaea"], pch = 20, col = adjustcolor(which(levels(vh$Prokaryotes)=="archaea"),0.33), cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
title(H,font.main=1,cex.main=1)
title(paste0("n = ",sum(Z)),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x[Z], y[Z], use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adj.R = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (C)
"" -> H
"c" -> h
"gc" -> i
x <- vh[,toupper(i)]
y <- vh[,i]
z <- complete.cases(x) & complete.cases(y)
Z <- vh$Eukaryotes!="" & complete.cases(x) & complete.cases(y)
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
points(x[z&vh$Eukaryotes==""], y[z&vh$Eukaryotes==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(vh$Eukaryotes)[-1]){points(x[z&vh$Eukaryotes==i], y[z&vh$Eukaryotes==i], pch = 20, col = which(levels(vh$Eukaryotes)==i), cex=1.5)}
for(i in names(sort(table(vh$Eukaryotes)[-1],decreasing = TRUE))){points(x[z&vh$Eukaryotes==i], y[z&vh$Eukaryotes==i], pch = 20, col = which(levels(vh$Eukaryotes)==i), cex=1.5)}
points(x[z&vh$Eukaryotes=="protists"], y[z&vh$Eukaryotes=="protists"], pch = 20, col = adjustcolor(which(levels(vh$Eukaryotes)=="protists"),0.33), cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
title(H,font.main=1,cex.main=1)
title(paste0("n = ",sum(Z)),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x[Z], y[Z], use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adj.R = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)

## (legend)
plot(1:10, col = "white", axes = FALSE)
legend("left", legend=levels(vh$Prokaryotes)[-1], col=host[-1],
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
                                          r=cor(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], use="complete.obs"),
                                          n=sum(complete.cases(data[data$prokaryotic!=k,c(i,j)])))
    }else{
      dataPk[[paste(j,i,sep="~")]] <- rbind(dataPk[[paste(j,i,sep="~")]],
                                          cbind(k=k,
                                                r=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="complete.obs"),
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
                                          r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="complete.obs"),
                                          n=sum(complete.cases(data[data$eukaryotic!=k,c(i,j)])))
    }else{
      dataEk[[paste(j,i,sep="~")]] <- rbind(dataEk[[paste(j,i,sep="~")]],
                                          cbind(k=k,
                                                r=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="complete.obs"),
                                                n=sum(complete.cases(data[data$eukaryotic==k,c(i,j)]))))
    }
  }
}
dataEk

codPk<- NULL
for(j in codons){
  i=toupper(j)
  for(k in levels(data$prokaryotic)){
    if(k == ""){
      codPk[[paste(j,i,sep="~")]] <- cbind(k="prokaryotic",
                                           r=cor(data[data$prokaryotic!=k,i], data[data$prokaryotic!=k,j], use="complete.obs"),
                                           R2=summary(lm(data[data$prokaryotic!=k,j]~data[data$prokaryotic!=k,i]))$adj.r.squared#,
                                           # n=sum(complete.cases(data[data$prokaryotic!=k,c(i,j)]))
      )
    }else{
      codPk[[paste(j,i,sep="~")]] <- rbind(codPk[[paste(j,i,sep="~")]],
                                           cbind(k=k,
                                                 r=cor(data[data$prokaryotic==k,i], data[data$prokaryotic==k,j], use="complete.obs"),
                                                 R2=summary(lm(data[data$prokaryotic==k,j]~data[data$prokaryotic==k,i]))$adj.r.squared#,
                                                 # n=sum(complete.cases(data[data$prokaryotic==k,c(i,j)]))
                                           ))
    }
  }
}
codPk
(codpk<-do.call(rbind.data.frame, codPk))
# write.table(codpk[codpk$k=="prokaryotic",],"codPk.csv",quote = F,row.names = T,col.names = T, sep=",")
# write.table(codpk[codpk$k!="prokaryotic",],"codPks.csv",quote = F,row.names = ,col.names = T, sep=",")

codEk<- NULL
for(j in codons){
  i=toupper(j)
  for(k in levels(data$eukaryotic)){
    if(k == ""){
      codEk[[paste(j,i,sep="~")]] <- cbind(k="eukaryotic",
                                           r=cor(data[data$eukaryotic!=k,i], data[data$eukaryotic!=k,j], use="complete.obs"),
                                           R2=summary(lm(data[data$eukaryotic!=k,j]~data[data$eukaryotic!=k,i]))$adj.r.squared#,
                                           # n=sum(complete.cases(data[data$eukaryotic!=k,c(i,j)]))
      )
    }else{
      codEk[[paste(j,i,sep="~")]] <- rbind(codEk[[paste(j,i,sep="~")]],
                                           cbind(k=k,
                                                 r=cor(data[data$eukaryotic==k,i], data[data$eukaryotic==k,j], use="complete.obs"),
                                                 R2=summary(lm(data[data$eukaryotic==k,j]~data[data$eukaryotic==k,i]))$adj.r.squared#,
                                                 # n=sum(complete.cases(data[data$eukaryotic==k,c(i,j)])
                                           ))
    }
  }
}
codEk
(codek<-do.call(rbind.data.frame, codEk))
# write.table(codek[codek$k=="eukaryotic",],"codEk.csv",quote = F,row.names = T,col.names = T, sep=",")
# write.table(codek[codek$k!="eukaryotic",],"codEks.csv",quote = F,row.names = T,col.names = T, sep=",")
