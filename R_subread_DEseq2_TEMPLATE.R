#########################################################################################
#####Run Rsubread for counts
library(Rsubread);

files=c(bamlist)

gtf=("initial_annotation_file")

###Sense Unique
senseunique=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
nthreads = cpus,
strandSpecific = 1)

write.table(x=data.frame(
senseunique$annotation[,c("GeneID","Length")],
senseunique$counts,stringsAsFactors=FALSE),
file="outputprefix/Rsubread_sense.txt",
quote=FALSE,sep="\t",row.names=FALSE)

###Anti-sense unique
antisenseunique=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
nthreads = cpus,
strandSpecific = 2)

write.table(x=data.frame(
antisenseunique$annotation[,c("GeneID","Length")],
antisenseunique$counts,stringsAsFactors=FALSE),
file="outputprefix/Rsubread_antisense.txt",
quote=FALSE,sep="\t",row.names=FALSE)


#########################################################################################
#########################################################################################
#######DESeq2

library(DESeq2)
#Read in the count matrix
countData <- as.matrix(read.csv("outputprefix/Rsubread_sense.txt",sep="\t",row.names=1))
colData <- read.table("COL_DATA")
countData <- subset(countData, select=-c(Length)) #Get rid of length column
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition)

dds <- dds[ rowSums(counts(dds)) > 1,]   #Take everything over 1 because we added one in the normalization python steps
dds$condition <- factor(dds$condition, levels=c("C","T"))
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file="outputprefix/DESeq2_sense.csv")
#############################################


countData <- as.matrix(read.csv("outputprefix/Rsubread_antisense.txt",sep="\t",row.names=1))
countData <- subset(countData, select=-c(Length))

colData <- read.table("COL_DATA")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition)

dds <- dds[ rowSums(counts(dds)) > 1,]   #Take everything over 2 because we added one in the normalization python steps
dds$condition <- factor(dds$condition, levels=c("C","T"))
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file="outputprefix/DESeq2_antisense.csv")
#############################################

