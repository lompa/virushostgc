#####
## Table 3 and Figure 2 (A-C + legend)

par <- par()
pal <- palette()

## you may need:
# install.packages("RColorBrewer")

## choose working directory:
path <- "~/virushostgc/"

setwd(path)

data <- read.csv("data3.csv"); dim(data)

data$baltimore <- factor(data$baltimore,
                         levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT"))
sum(is.na(data$baltimore))
# levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT",""))
# levels(data$baltimore)[nlevels(df$baltimore)] <- "unknown"

# RColorBrewer::display.brewer.all()
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))

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

# RColorBrewer::display.brewer.all()
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))

host <- c("gray",
          RColorBrewer::brewer.pal(12, "Paired")[6],
          RColorBrewer::brewer.pal(12, "Paired")[10],
          RColorBrewer::brewer.pal(12, "Paired")[2],
          RColorBrewer::brewer.pal(12, "Paired")[12],
          RColorBrewer::brewer.pal(12, "Paired")[4],
          RColorBrewer::brewer.pal(12, "Paired")[8])

#####
## Figures 2

par <- par()
par(family = "serif")
par(oma = c(0,0,0,0))
par(mar = c(3,3,2,3))# + 0.1)
par(mfrow = c(1,1))

baltimore <- c(RColorBrewer::brewer.pal(12, "Paired")[10],
               RColorBrewer::brewer.pal(12, "Paired")[8],
               RColorBrewer::brewer.pal(12, "Paired")[2],
               RColorBrewer::brewer.pal(12, "Paired")[4],
               RColorBrewer::brewer.pal(12, "Paired")[6],
               RColorBrewer::brewer.pal(8, "Set2")[1],
               RColorBrewer::brewer.pal(12, "Paired")[9])

palette(adjustcolor(baltimore,0.5))

## (A)
"a" -> h
"gc1" -> i
x <- data$gc0
y <- data[, i]
z <- complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y ~ x)
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
points(x[z], y[z], pch = 20, col = data$baltimore[z], cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# title(toupper(i),font.main=1,cex.main=1)
title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
(r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# library(grDevices)
# exp <- expression(R^2*" = "r2)
# https://rpubs.com/brouwern/superscript
# https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# text(2, 85, exp, adj = 0, cex = 1)
text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
text(72.5, 9, "2", adj = 0, cex = 0.7)
text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)
palette(adjustcolor(baltimore,0.5))

## (B)
"b" -> h
"gc2" -> i
x <- data$gc0
y <- data[, i]
z <- complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y ~ x)
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
points(x[z], y[z], pch = 20, col = data$baltimore[z], cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# title(toupper(i),font.main=1,cex.main=1)
title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
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
"c" -> h
"gc3" -> i
x <- data$gc0
y <- data[, i]
z <- complete.cases(x) & complete.cases(y)
plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
paste0(toupper(i), " (%)")
eq <- lm(y ~ x)
abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
points(x[z], y[z], pch = 20, col = data$baltimore[z], cex=1.5)
title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# title(toupper(i),font.main=1,cex.main=1)
title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# text(2, 85, paste0("r = ", , adj = 0, cex = 1)
(r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
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
legend("left", legend=levels(data$baltimore), col=baltimore,
       box.lty=1, pch = 20, pt.cex = 1.5,cex = 1)

# #####
# ## Figures S2
# 
# par <- par()
# par(family = "serif")
# par(oma = c(0,0,0,0))
# par(mar = c(3,3,2,3))# + 0.1)
# par(mfrow = c(1,1))
# 
# palette(adjustcolor(host,0.5))
# 
# ## (A)
# "a" -> h
# "gc1" -> i
# x <- data$gc0
# y <- data[, i]
# z <- complete.cases(x) & complete.cases(y)
# plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
# axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# # mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
# axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# paste0(toupper(i), " (%)")
# eq <- lm(y ~ x)
# abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# # points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
# points(x[z&data$host==""], y[z&data$host==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(data$host)[-1]){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# # for(i in names(sort(table(data$host)[-1],decreasing = TRUE))){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# # title(toupper(i),font.main=1,cex.main=1)
# title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# # text(2, 85, paste0("r = ", , adj = 0, cex = 1)
# (r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
# (r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# # library(grDevices)
# # exp <- expression(R^2*" = "r2)
# # https://rpubs.com/brouwern/superscript
# # https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# # text(2, 85, exp, adj = 0, cex = 1)
# text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
# text(72.5, 9, "2", adj = 0, cex = 0.7)
# text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)
# 
# ## (B)
# "b" -> h
# "gc2" -> i
# x <- data$gc0
# y <- data[, i]
# z <- complete.cases(x) & complete.cases(y)
# plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
# axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
# axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# paste0(toupper(i), " (%)")
# eq <- lm(y ~ x)
# abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# # points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
# points(x[z&data$host==""], y[z&data$host==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(data$host)[-1]){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# # for(i in names(sort(table(data$host)[-1],decreasing = TRUE))){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# # title(toupper(i),font.main=1,cex.main=1)
# title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# # text(2, 85, paste0("r = ", , adj = 0, cex = 1)
# (r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
# (r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# # library(grDevices)
# # exp <- expression(R^2*" = "r2)
# # https://rpubs.com/brouwern/superscript
# # https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# # text(2, 85, exp, adj = 0, cex = 1)
# text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
# text(72.5, 9, "2", adj = 0, cex = 0.7)
# text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)
# 
# ## (C)
# "c" -> h
# "gc3" -> i
# x <- data$gc0
# y <- data[, i]
# z <- complete.cases(x) & complete.cases(y)
# plot(x[z], y[z], xlab = "", xlim = c(0, 100), ylim = c(0, 100), xaxt = "n", yaxt = "n", col = "white")
# axis(1, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 1, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# # mtext("non-coding GC", side = 1, line = 2, cex.lab = 1, las = 0, col = "black")
# axis(2, at=seq(0,100,by=50),labels = c(seq(0,50,by=50),""))
# mtext("100%", side = 2, adj=1, line = 1, cex.lab = 1, las = 0, col = "black")
# mtext(toupper(i), side = 2, line = 2, cex.lab = 1, las = 0, col = "black")
# paste0(toupper(i), " (%)")
# eq <- lm(y ~ x)
# abline(eq, lwd = 2, col = adjustcolor("black", 0.5), lty = 3)
# # points(x[z], y[z], pch = 20, col = data$host[z], cex=1.5)
# points(x[z&data$host==""], y[z&data$host==""], pch = 20, col = 1, cex=1.5)
# for(i in levels(data$host)[-1]){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# # for(i in names(sort(table(data$host)[-1],decreasing = TRUE))){points(x[z&data$host==i], y[z&data$host==i], pch = 20, col = which(levels(data$host)==i), cex=1.5)}
# title(paste0("(",toupper(h),")"),adj=0,cex.main=1.5)
# # title(toupper(i),font.main=1,cex.main=1)
# title(paste0("n = ",format(sum(z),big.mark=",")),adj=1,font.main=1,cex.main=1)
# # text(2, 85, paste0("r = ", , adj = 0, cex = 1)
# (r <- sprintf(round(cor(x, y, use = "complete.obs"), 2), fmt = "%#.2f"))
# (r2 <- sprintf(round(summary(eq)$adj.r.squared, 2), fmt = "%#.2f"))
# # library(grDevices)
# # exp <- expression(R^2*" = "r2)
# # https://rpubs.com/brouwern/superscript
# # https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
# # text(2, 85, exp, adj = 0, cex = 1)
# text(55, 5, label=paste0("adjR  = ",r2), adj = 0, cex = 1)
# text(72.5, 9, "2", adj = 0, cex = 0.7)
# text(55, 15, paste0("slope = ", sprintf(round(eq$coefficients[2], 2), fmt = "%#.2f")),adj=0, cex = 1)
# 
# ## (legend)
# plot.new()
# legend("left", legend=levels(data$host)[-1], col=host[-1],
#        box.lty=1, pch = 20, pt.cex = 1.5,cex = 1)
# 
# par(par)
# palette(pal)

#####
## Correlations

with(data, cor(gc, gc0, use="complete.obs"))
with(data, cor(gc, gc123, use="complete.obs"))
with(data, cor(gc0, gc123, use="complete.obs"))

with(data,cor(gc, gc1, use="complete.obs"))
with(data,cor(gc, gc2, use="complete.obs"))
with(data,cor(gc, gc3, use="complete.obs"))

with(data,cor(gc0, gc1, use="complete.obs"))
with(data,cor(gc0, gc2, use="complete.obs"))
with(data,cor(gc0, gc3, use="complete.obs"))

with(data, cor(gc, gc0, method="spearman", use="complete.obs"))
with(data, cor(gc, gc123, method="spearman", use="complete.obs"))
with(data, cor(gc0, gc123, method="spearman", use="complete.obs"))

with(data,cor(gc, gc1, method="spearman", use="complete.obs"))
with(data,cor(gc, gc2, method="spearman", use="complete.obs"))
with(data,cor(gc, gc3, method="spearman", use="complete.obs"))

with(data,cor(gc0, gc1, method="spearman", use="complete.obs"))
with(data,cor(gc0, gc2, method="spearman", use="complete.obs"))
with(data,cor(gc0, gc3, method="spearman", use="complete.obs"))

rBaltimore <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rBaltimore[[paste(j,i,sep="~")]] <- do.call(rbind,
                                              lapply(
                                                split(data,data$baltimore),
                                                function(x) data.frame(x=i,
                                                                       y=j,
                                                                       group=x$baltimore[1],
                                                                       r=cor(x[,i], x[,j], use="complete.obs"),
                                                                       p=cor.test(x[,i], x[,j], use="complete.obs")$p.value,
                                                                       rho=cor(x[,i], x[,j], method="spearman", use="complete.obs"),
                                                                       p=cor.test(x[,i], x[,j], method="spearman", use="complete.obs")$p.value,
                                                                       n=sum(complete.cases(x[,c(i,j)])))))
}
rBaltimore
# write.table(do.call(rbind.data.frame,rBaltimore),"table2.csv",quote = F, row.names=T, col.names = T, sep=",")

#####
## Table 3

rBaltimore[[1]]

rHostCell <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostCell[[paste(j,i,sep="~")]] <- do.call(rbind,
                                             lapply(
                                               split(data,data$host_cell),
                                               function(x) data.frame(x=i,
                                                                      y=j,
                                                                      group=x$host_cell[1],
                                                                      r=cor(x[,i], x[,j], use="complete.obs"),
                                                                      rho=cor(x[,i], x[,j], method="spearman", use="complete.obs"),
                                                                      p=cor.test(x[,i], x[,j], method="spearman", use="complete.obs")$p.value,
                                                                      n=sum(complete.cases(x[,c(i,j)])))))
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][row.names(rHostCell[[paste(j,i,sep="~")]])!="1",]
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][order(rHostCell[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostCell

rHostType <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostType[[paste(j,i,sep="~")]] <- do.call(rbind,
                                             lapply(
                                               split(data,data$host_type),
                                               function(x) data.frame(x=i,
                                                                      y=j,
                                                                      group=x$host_type[1],
                                                                      r=cor(x[,i], x[,j], use="complete.obs"),
                                                                      rho=cor(x[,i], x[,j], method="spearman", use="complete.obs"),
                                                                      p=cor.test(x[,i], x[,j], method="spearman", use="complete.obs")$p.value,
                                                                      n=sum(complete.cases(x[,c(i,j)])))))
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][row.names(rHostType[[paste(j,i,sep="~")]])!="1",]
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][order(rHostType[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostType

rHost <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHost[[paste(j,i,sep="~")]] <- do.call(rbind,
                                         lapply(
                                           split(data,data$host),
                                           function(x) data.frame(x=i,
                                                                  y=j,
                                                                  group=x$host[1],
                                                                  r=cor(x[,i], x[,j], use="complete.obs"),
                                                                  rho=cor(x[,i], x[,j], method="spearman", use="complete.obs"),
                                                                  p=cor.test(x[,i], x[,j], method="spearman", use="complete.obs")$p.value,
                                                                  n=sum(complete.cases(x[,c(i,j)])))))
  rHost[[paste(j,i,sep="~")]] <- rHost[[paste(j,i,sep="~")]][row.names(rHost[[paste(j,i,sep="~")]])!="1",]
  rHost[[paste(j,i,sep="~")]] <- rHost[[paste(j,i,sep="~")]][order(rHost[[paste(j,i,sep="~")]]$n,
                                                                   decreasing = TRUE),]
}
rHost

dtPk <- data[data$host_cell == "prokaryote",]

with(dtPk,cor(gc, gc1, use = "complete.obs"))
with(dtPk,cor(gc, gc2, use = "complete.obs"))
with(dtPk,cor(gc, gc3, use = "complete.obs"))

with(dtPk,cor(gc, gc1, use = "complete.obs"))
with(dtPk,cor(gc, gc2, use = "complete.obs"))
with(dtPk,cor(gc, gc3, use = "complete.obs"))

rBaltimorePk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rBaltimorePk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                                lapply(
                                                  split(dtPk,dtPk$baltimore),
                                                  function(x) data.frame(x=i,
                                                                         y=j,
                                                                         group=x$baltimore[1],
                                                                         # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                         r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                         p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                         rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                         p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                         n=sum(complete.cases(x[,c(i,j)])))))
}
rBaltimorePk

rHostCellPk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostCellPk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                               lapply(
                                                 split(dtPk,dtPk$host_cell),
                                                 function(x) data.frame(x=i,
                                                                        y=j,
                                                                        group=x$host_cell[1],
                                                                        # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                        r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                        rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                        n=sum(complete.cases(x[,c(i,j)])))))
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][row.names(rHostCell[[paste(j,i,sep="~")]])!="1",]
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][order(rHostCell[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostCellPk

rHostTypePk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostTypePk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                               lapply(
                                                 split(dtPk,dtPk$host_type),
                                                 function(x) data.frame(x=i,
                                                                        y=j,
                                                                        group=x$host_type[1],
                                                                        # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                        r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                        rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                        n=sum(complete.cases(x[,c(i,j)])))))
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][row.names(rHostType[[paste(j,i,sep="~")]])!="1",]
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][order(rHostType[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostTypePk

dtEk <- data[data$host_cell != "prokaryote",]

with(dtEk,cor(gc, gc1, use = "complete.obs"))
with(dtEk,cor(gc, gc2, use = "complete.obs"))
with(dtEk,cor(gc, gc3, use = "complete.obs"))

with(dtEk,cor(gc0, gc1, use = "complete.obs"))
with(dtEk,cor(gc0, gc2, use = "complete.obs"))
with(dtEk,cor(gc0, gc3, use = "complete.obs"))

rBaltimoreEk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rBaltimoreEk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                                lapply(
                                                  split(dtEk,dtEk$baltimore),
                                                  function(x) data.frame(x=i,
                                                                         y=j,
                                                                         group=x$baltimore[1],
                                                                         # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                         r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                         p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                         rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                         p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                         n=sum(complete.cases(x[,c(i,j)])))))
}
rBaltimoreEk

rHostCellEk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostCellEk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                               lapply(
                                                 split(dtEk,dtEk$host_cell),
                                                 function(x) data.frame(x=i,
                                                                        y=j,
                                                                        group=x$host_cell[1],
                                                                        # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                        r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                        rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                        n=sum(complete.cases(x[,c(i,j)])))))
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][row.names(rHostCell[[paste(j,i,sep="~")]])!="1",]
  rHostCell[[paste(j,i,sep="~")]] <- rHostCell[[paste(j,i,sep="~")]][order(rHostCell[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostCellEk

rHostTypeEk <- NULL
"gc0" -> i
for(j in c("gc1","gc2","gc3")){
  rHostTypeEk[[paste(j,i,sep="~")]] <- do.call(rbind,
                                               lapply(
                                                 split(dtEk,dtEk$host_type),
                                                 function(x) data.frame(x=i,
                                                                        y=j,
                                                                        group=x$host_type[1],
                                                                        # r=cor(x[,i], x[,j], use="complete.obs"),
                                                                        r=cor(x[,i], x[,j], use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], use="na.or.complete")$p.value,
                                                                        rho=cor(x[,i], x[,j], method="spearman", use="na.or.complete"),
                                                                        p=cor.test(x[,i], x[,j], method="spearman", use="na.or.complete")$p.value,
                                                                        n=sum(complete.cases(x[,c(i,j)])))))
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][row.names(rHostType[[paste(j,i,sep="~")]])!="1",]
  rHostType[[paste(j,i,sep="~")]] <- rHostType[[paste(j,i,sep="~")]][order(rHostType[[paste(j,i,sep="~")]]$n,
                                                                           decreasing = TRUE),]
}
rHostTypeEk
