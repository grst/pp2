library("ggplot2")

data <- read.csv("Übung03/DeepLoc.human_proteins_prediction.tsv", sep="\t", header=TRUE)
df <- as.data.frame(table(data[,2]))

png("deeploc.prediction_distribution.png", 1000,500, res=120)
p <- ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq, fill=Var1))
p <- p + geom_bar(stat="identity",position="dodge", colour="black")
p <- p + coord_flip()
p <- p + xlab("")
p <- p + ylab("Number of predictions")
p <- p + guides(fill=FALSE)
p <- p + theme(	axis.text.x = element_text(angle = 0, vjust = 0.5, size=12), 
                axis.text.y = element_text(size=12), 
                axis.title.y = element_text(size=12))
p + ggtitle("Deeploc data set", subtitle = "Location prediction distribution")
dev.off()

pred <- data[,-1:-2]
png("deeploc.predictionScore_distribution.png", 500,500, res=100)
hist(as.numeric(unlist(pred)), breaks=25, col=rainbow(100), xlab="Prediction score", main="Prediction score distribution")
dev.off()

pred_rowSum <- apply(pred, 1, sum)
png("deeploc.predictionScoreSum_distribution.png", 500,500, res=100)
hist(pred_rowSum, breaks=50, col=rainbow(100), xlab="Sum of all prediction scores", main="Pred. score sum for each protein")
dev.off()