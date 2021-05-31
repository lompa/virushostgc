#####
## Table 2

par <- par()
pal <- palette()

## you may need:
# install.packages("RColorBrewer")

## choose working directory:
path <- "~/virushostgc/"

setwd(path)

data <- read.csv("data2.csv"); dim(data)

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
## Table 2


df <- t(data.frame(c(nrow(data), table(data$Host))))
rownames(df) <- ""
colnames(df)[1] <- "Total"
df

# pie(table(data$Host), cex = 1,
#     col = adjustcolor(host[-1], 1))
