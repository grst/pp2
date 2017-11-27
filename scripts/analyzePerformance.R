args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
	stop("<deeploc_dataset>")
}

#Read data
swissprot <- read.csv("../results/swissprot_filtered.tsv", sep="\t", header=TRUE)
hpa <- read.csv("../results/hpa_filtered.tsv", sep="\t", header=TRUE)
deeploc <- read.csv(args[1], sep="\t", header=TRUE)

#Read mapping + merge mapping with SwissProt & HPA
mapping <- read.csv("../results/location_mapping.tsv", sep="\t", header=TRUE)
swissprot_merge <- merge(swissprot, mapping[mapping$target_source=="swissprot",][,1:2], by.x="subcellular_location", by.y="target_location")
hpa_merge <- merge(hpa, mapping[mapping$target_source=="hpa",][,1:2], by.x="subcellular_location", by.y="target_location")

#Data set exploration
print("Number of subcellular locations annotated in SwissProt data set:")
length(unique(swissprot$subcellular_location))
print("Number of swissprot subcellular locations in the location_mapping.tsv:")
length(unique(mapping[mapping$target_source=="swissprot",][,1:2]$target_location))
print("Data set - Swissprot proteins:")
length(unique(swissprot$uniprot_id))
print("Data set - Swissprot_Mapping merge proteins:")
length(unique(swissprot_merge$uniprot_id))

#Deeploc filtering
deeploc <- deeploc[deeploc$Prediction!="no_accurate_estimation",]

#Performance calculation - SwissProt
##Nur Proteine mit einem 1:1 Annotation (1 Prediction + 1 Annotation durch Mapping) werden betrachtet
performance_swissprot <- merge(deeploc, swissprot_merge[c("hgnc","deeploc_location")], by.x="hgnc", by.y="hgnc")
set <- table(performance_swissprot$hgnc)==1
subset <- performance_swissprot[performance_swissprot$hgnc %in% names(set[set]),]
subset$Prediction <- as.character(subset$Prediction)
subset$deeploc_location <- as.character(subset$deeploc_location)
subset$deeploc_location <- gsub(" ", "_", subset$deeploc_location)
all <- length(subset[,1])
true <- table(subset$Prediction==subset$deeploc_location)[2]
false <- table(subset$Prediction==subset$deeploc_location)[1]
cat("Swissprot Performance (1:1 annotation):", true/all, "\n",  sep="\t")
cat("All:", all, "\t", "True:", true, "\t", "False:", false, "\n",  sep=" ")

##Alle Proteine mit einem 1:N Annotation (1 Prediction + N Annotationen durch Mapping) werden betrachtet
performance_swissprot$Prediction <- as.character(performance_swissprot$Prediction)
performance_swissprot$deeploc_location <- as.character(performance_swissprot$deeploc_location)
performance_swissprot$deeploc_location <- gsub(" ", "_", performance_swissprot$deeploc_location)
performance_swissprot <- performance_swissprot[!duplicated(performance_swissprot),]
all <- length(unique(performance_swissprot$hgnc))
true <- table(performance_swissprot$Prediction==performance_swissprot$deeploc_location)[2]
false <- all - true
cat("Swissprot Performance (1:N annotation):", true/all, "\n",  sep="\t")
cat("All:", all, "\t", "True:", true, "\t", "False:", false, "\n", sep=" ")
multiMapping <- table(performance_swissprot$hgnc)
cat("Average number of annotations per protein:", mean(multiMapping[multiMapping!=0]), "\n", sep=" ")

#Performance calculation - HPA
##Nur Proteine mit einem 1:1 Annotation (1 Prediction + 1 Annotation durch Mapping) werden betrachtet
performance_hpa <- merge(deeploc, hpa_merge[c("hgnc","deeploc_location")], by.x="hgnc", by.y="hgnc")
set <- table(performance_hpa$hgnc)==1
subset <- performance_hpa[performance_hpa$hgnc %in% names(set[set]),]
subset$Prediction <- as.character(subset$Prediction)
subset$deeploc_location <- as.character(subset$deeploc_location)
subset$deeploc_location <- gsub(" ", "_", subset$deeploc_location)
all <- length(subset[,1])
true <- table(subset$Prediction==subset$deeploc_location)[2]
false <- table(subset$Prediction==subset$deeploc_location)[1]
cat("HPA Performance (1:1 annotation):", true/all, "\n",  sep="\t")
cat("All:", all, "\t", "True:", true, "\t", "False:", false, "\n",  sep=" ")

##Alle Proteine mit einem 1:N Annotation (1 Prediction + N Annotationen durch Mapping) werden betrachtet
performance_hpa$Prediction <- as.character(performance_hpa$Prediction)
performance_hpa$deeploc_location <- as.character(performance_hpa$deeploc_location)
performance_hpa$deeploc_location <- gsub(" ", "_", performance_hpa$deeploc_location)
performance_hpa <- performance_hpa[!duplicated(performance_hpa),]
all <- length(unique(performance_hpa$hgnc))
true <- table(performance_hpa$Prediction==performance_hpa$deeploc_location)[2]
false <- all - true
cat("HPA Performance (1:N annotation):", true/all, "\n",  sep="\t")
cat("All:", all, "\t", "True:", true, "\t", "False:", false, "\n",  sep=" ")
multiMapping <- table(performance_hpa$hgnc)
cat("Average number of annotations per protein:", mean(multiMapping[multiMapping!=0]), "\n", sep=" ")
