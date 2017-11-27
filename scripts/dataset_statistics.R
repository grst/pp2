library("ggplot2")

data <- read.csv("../data/DeepLoc.human_proteins_prediction.txt", sep="\t", header=TRUE)
df <- as.data.frame(table(data[,2]))

#Plot class prediction distribution
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


#Plot all class prediction score distribution
pred <- data[,-1:-3]
png("deeploc.predictionScore_distribution.all_classes.png", 500,500, res=100)
hist(as.numeric(unlist(pred)), breaks=25, col=rainbow(100), xlab="Prediction score", main="All class prediction score distribution")
dev.off()

#Plot membrane bound prediction score distribution
png("deeploc.predictionScore_distribution.membrane_bound.png", 500,500, res=100)
hist(as.numeric(data$Membrane.bound), breaks=25, col=rainbow(100), xlab="Prediction score", main="Membrane bound prediction score distribution")
dev.off()

#Plot highest scoring class prediction score distribution - all
png("deeploc.HighestpredictionScore_distribution.all_classes.png", 500,500, res=100)
hist(apply(pred, 1, max), breaks=25, col=rainbow(100), xlab="Prediction score", main="Highest prediction score distribution")
dev.off()

#Plot highest scoring class prediction score distribution - density for each class
classes <- colnames(pred)[apply(pred, 1, which.max)]
scores <- apply(pred, 1, max)
df <- data.frame(class=classes, score=scores)
p <- ggplot(df, aes(score, fill=class))
p <- p + geom_density(alpha=.5)
p <- p + theme(	axis.text.x = element_text(size=12), 
                axis.text.y = element_text(size=12), 
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12))
p <- p + scale_fill_discrete(name = "Subcellular locations")
p <- p + xlab("Prediction score")
p <- p + ylab("Density")
png("deeploc.HighestpredictionScore_distribution.density.png", 750,500, res=100)
p + ggtitle("Highest prediction score for each protein", subtitle = NULL)
dev.off()

#Plot highest scoring class prediction score distribution - Boxplot for each class
p <- ggplot(data=df, aes(x=reorder(class, score), y=score, fill=class))
p <- p + geom_boxplot()
p <- p + scale_fill_discrete(name = "Subcellular locations")
p <- p + xlab("Subcellular locations")
p <- p + ylab("Prediction score")
p <- p + coord_flip()
p <- p + theme(	axis.text.x = element_text(angle = 0, vjust = 0.5, size=12), 
                axis.text.y = element_text(size=12), 
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12))
png("deeploc.HighestpredictionScore_distribution.boxplot.png", 750,500, res=100)
p + ggtitle("Highest prediction score for each protein", subtitle = NULL)
dev.off()


#Plot membrane-bound score associated with predicted location
df <- data.frame(class=classes, membrane_score=as.numeric(data$Membrane.bound))
p <- ggplot(data=df, aes(x=reorder(class, membrane_score), y=membrane_score, fill=class))
p <- p + geom_boxplot()
p <- p + scale_fill_discrete(name = "Subcellular locations")
p <- p + xlab("Subcellular locations")
p <- p + ylab("Membrane bound prediction score")
p <- p + coord_flip()
p <- p + theme(	axis.text.x = element_text(angle = 0, vjust = 0.5, size=12), 
                axis.text.y = element_text(size=12), 
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12))
png("deeploc.MembraneBoundScore.PredictedLocation_MembraneBoundScore.boxplot.png", 1000,500, res=100)
p + ggtitle("Membrane bound predicted score compared to the predicted location", subtitle = NULL)
dev.off()

#visualise the distribution of predicted compartments for different reliability indices
classes <- colnames(pred)[apply(pred, 1, which.max)]
scores <- apply(pred, 1, max)
df <- data.frame(class=classes, score=scores)
p <- ggplot(df, aes(x=class, y=score, fill=class))
p <- p + geom_violin(trim=TRUE)
p <- p + coord_flip()
p <- p + scale_fill_discrete(name = "Subcellular locations")
p <- p + xlab("Subcellular locations")
p <- p + ylab("Prediction score")
p <- p + ylim(0,1)
p <- p + theme(	axis.text.x = element_text(angle = 0, vjust = 0.5, size=12), 
                axis.text.y = element_text(size=12), 
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12))
png("deeploc.HighestpredictionScore_distribution.violin.png", 1000,500, res=100)
p + ggtitle("Predicted compartment distribution for different reliability indices", subtitle = NULL)
dev.off()
