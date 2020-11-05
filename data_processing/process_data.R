rm(list=ls())
library(here)
source(file.path(here(), 'Functions.R')) # source function scripts
setwd(file.path(here(), 'path/to/results'))

# load data
load(file.path(here(), 'path/to/data.RData'))

# catch output log
sink('log.txt')

cat("Initial number of traces: ", nrow(data), "\n")

# blankfiltering - removing contaminants (instrumental)
blanks=specgrep(data,"Blank")
data <- removegrep(data, "Blank")

to.remove=BlankFilter(blanks,data[,-c(1,2)],0.01)
data=data[-to.remove,]
cat("Blankfiltering", "\n", "Removed by BlankFilter: ",length(to.remove), "\n")

# dilution based filtering
d <- specgrep(data,"D")
data <- removegrep(data,"D")

conc <- c(0.5,1,2,4,8,16,32)

# calculate correlation between QC concentration and feature
c = list()
for (i in 1:dim(d)[1]) {
  x <- as.numeric(d[i,])
  y <- conc
  ind <- which(is.na(x))
  if (length(ind)>0) {
    x <- x[-ind]
    y <- y[-ind]
  }
  if (length(x)>2) {
    c[[i]] = cor.test(x,y,use="p")
  } else {
    c[[i]] = 0
  }
}

# Find features with a correlation significance < 0.05
significant = list()
to.keep = c()
n = 1
for (j in 1:length(c)) {
  if (length(c[[j]])>1){
    if (c[[j]]$p.value < 0.05 && !is.na(c[[j]]$p.value)) {
      significant[[n]] <- c[[j]]
      to.keep[n] <- j
      n = n + 1
    }
  }
}
cat("Dilution based filtering", "\n", "Features which correlate with the dilution serie: ", length(to.keep), "\n")
cat("Absolute correlation range: ", round(range(abs(unlist(lapply(significant, function(x) x$estimate)))), digits = 2), "\n")
data = data[to.keep,]

# 75% coverage cutoff
ind <- which(apply(data[,-c(1,2)], 1, function(x) sum(!is.na(x))/length(x))>=0.75)
data <- data[ind,]
cat("After 75% coverage threshold: ", nrow(data), "\n")

# log2 transformation
cat("log2 transformation", "\n")
data[,-c(1,2)] = log2(data[,-c(1,2)])
data.before <- data

# normalize by reg. LOESS
cat("Normalize using customized LOESS", "\n")
library(limma)
data.norm <- order.norm(data[,-c(1,2)], 0.2)
data <- data.frame(cbind(data[,c(1,2)], data.norm))

# QC outlier exclusion
cat("QC Outlier exclusion", "\n")
qcs <- specgrep(data,"QC")

TIC<-colSums(qcs, na.rm=T)
qcs.oe <- qcs[,!(colnames(qcs) %in% c(colnames(qcs)[TIC<mean(TIC)*0.75]))]

cat("Number of samples: ", ncol(qcs), "\n")
cat("Non-outlier samples: ", ncol(qcs.oe),"\n")
cat("Removed sample: ", outersect(names(qcs),names(qcs.oe)), "\n")
if (length(outersect(names(qcs),names(qcs.oe)))>0) {
  data <- data[,-which(names(data)==outersect(names(qcs),names(qcs.oe)))] }
data.before <- data.before[,names(data)]
qcs <- 2^qcs.oe

# calculate CVs in QCs
qc.cv <- data.frame(matrix(NA, nrow=nrow(qcs),ncol=3))
names(qc.cv) <- c("mz","rt", "CV")

qc.cv[,1] <- data$mz_cf
qc.cv[,2] <- data$rt_cf
qc.cv[,3] <- apply(qcs, 1, function(x) sd(as.numeric(x), na.rm=T)/mean(as.numeric(x), na.rm=T))

# filter on CV
cat("Filter on CV", "\n")
stable <- which(qc.cv$CV<0.20)
qc.cv <- qc.cv[stable,]
data <- data[stable,]
data.before <- data.before[stable,]
write.table(qc.cv, file="CoV.csv")
cat("After CV filter: ", nrow(data), "\n")

# plot the result of the normalization on 4 random metabolites
library(ggthemes)
library(ggplot2)
library(cowplot)
library(gridExtra)

show <- sample(1:nrow(data),4)
n <- 1
gg <- list()
for (i in 1:4) {
  # plot before normalization
  plotdata <- data.frame(names(data.before)[-c(1,2)], t(data.before[show[i],-c(1,2)]))
  names(plotdata) <- c("Sample", "Intensity")
  plotdata$Sample <- as.character(plotdata$Sample)
  plotdata$Sample <- factor(plotdata$Sample, levels = plotdata$Sample)
  print(range(plotdata$Intensity, na.rm=T))
  
  cols <- rep("#377EB8", nrow(plotdata))
  cols[grep("QC", plotdata$Sample)] <- "#C1002C"
  gg[[n]] <- ggplot(plotdata, aes(x=Sample, y=Intensity, col=Sample)) + geom_point() + 
    scale_color_manual(values=cols) + theme(legend.position = "none", axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank()) + ylim(10,29) +
    ylab("log2(Intensity)") + xlab("Sample runorder") + 
    ggtitle(paste("Raw mz", 
                  round(data.before[show[i],1], digits = 2), "rt", 
                  round(data.before[show[i],2], digits = 2)))
  n <- n + 1
  rm(plotdata)
  
  # plot after normalization
  plotdata <- data.frame(names(data)[-c(1,2)], t(data[show[i],-c(1,2)]))
  names(plotdata) <- c("Sample", "Intensity")
  plotdata$Sample <- as.character(plotdata$Sample)
  plotdata$Sample <- factor(plotdata$Sample, levels = plotdata$Sample)
  print(range(plotdata$Intensity, na.rm=T))
  
  gg[[n]] <- ggplot(plotdata, aes(x=Sample, y=Intensity, col=Sample)) + geom_point() + 
    scale_color_manual(values=cols) + theme(legend.position = "none", axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank()) + ylim(10,29) +
    ylab("log2(Intensity)") + xlab(paste("CV", round(qc.cv[show[i],2], digits = 3))) + 
    ggtitle(paste("Normalized mz", 
                  round(data[show[i],1], digits = 2), "rt", 
                  round(data[show[i],2], digits = 2)))
  n <- n + 1
  rm(plotdata)
}

plot <- list(gg[[1]], gg[[3]], gg[[5]], gg[[7]],
             gg[[2]], gg[[4]], gg[[6]], gg[[8]])

pdf(file = "randMetabos.pdf", width=19, height=9)
grid.arrange(grobs = plot, ncol = 4)
dev.off()

data <- removegrep(data,"QC")

# sample outlier exclusion
cat("Outlier exclusion", "\n")
data.oe <- data[,c(1,2)]
data <- data[,-c(1,2)]

TIC<-colSums(data, na.rm=T)
par(mfrow=c(1,2))
boxplot(TIC, ylab="Extracted peaks area")
hist(TIC, breaks=100, main="", xlab="Extracted peaks area")

cat("Number of samples: ", ncol(data), "\n")
cat("Non-outlier samples: ", ncol(data[,!(colnames(data) %in% c(colnames(data)[TIC<mean(TIC)*0.75]))]),"\n")

data.oe <- data.frame(cbind(data.oe, data[,!(colnames(data) %in% c(colnames(data)[TIC<mean(TIC)*0.75]))]))
names(data.oe) <- gsub("X", "", names(data.oe))
cat("Removed sample: ", outersect(names(data),names(data.oe)[-c(1,2)]), "\n")

TIC<-colSums(data.oe[,-c(1,2)], na.rm=T)
par(mfrow=c(1,2))
boxplot(TIC, ylab="Extracted peaks area")
hist(TIC, breaks=100, main="", xlab="Extracted peaks area")

data <- data.oe

# save processed data
save.image("workspace_processed.RData")
save(data, file="processed_data.RData")
cat("Final number of features: ", nrow(data), "\n")
sink()

