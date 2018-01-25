source("http://bioconductor.org/biocLite.R")
biocLite("cqn")
library(cqn)
library(gplots)

df <- read.table("./peaks_coverage/allPeaksCoverage.txt", sep="\t", header = T)

# Import GC table
GC <- df[,5]

# CQN normalisation
# Args 1 - matrix of counts, 2 - cofactor (GC), 3 - lengths (end - start)
tblCQN <- cqn(
  as.matrix(df[,-c(1:5)]),
  x=GC,
  lengths = (df[,3] - df[,2]), 
  lengthMethod="fixed"
)

# Display normalised Library Sizes before and after
libs <- data.frame(apply(df[,-c(1:5)],2,sum), apply(tblCQN$y, 2, sum), apply(tblCQN$y+tblCQN$offset, 2, sum))
colnames(libs) <- c("Un normalised", "Normalised", "Norm + offset")
print(libs)

#Correlation
x <- as.matrix(tblCQN$y)
rownames(x) <- df[,4]
y <- as.matrix(data.frame(df[,1:5], tblCQN$y))
write.table(y, "./peaks_coverage/allPeaksCoverage_CQN.txt", quote=F, row.names=F, sep="\t")

CorData <- cor(x, method="pearson") # t(x) for correlation across peaks

pdf("./peaks_coverage/PeakCorrelation.pdf",8,8)
	heatmap.2(
	  CorData,
	  trace= "none",
	  Rowv = T,
	  Colv = T,
	  dendrogram = "both",
	  scale = "none",
	  # col = hmcols,
	  main = "Heatmap - correlation"
	)
dev.off()

sampleDists <- as.matrix(dist(t(x)))

pdf("./peaks_coverage/PeakDistance.pdf",8,8)
	heatmap.2(
	  sampleDists,
	  trace= "none",
	  Rowv = T,
	  Colv = T,
	  dendrogram = "both",
	  scale = "none",
	  # col = hmcols,
	  main = "Heatmap - Distance Matrix"
	)
dev.off()
