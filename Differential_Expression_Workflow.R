#This code performs differential expression analysis for a Pre/Post Sleep Deprivation Microarray U133A gene chip dataset
#Written by Bryan Dighera


#Code Examples and Arguments: http://www.bioconductor.org/packages/devel/bioc//manuals/limma/man/limma.pdf

#Installs required libraries
#Bioconductor Installation
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

#U133a Genechip Annotation Library installation
if (!requireNamespace('BiocManager', quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install('hgu133a.db', version='3.8')

#GO - Gene Ontology Package installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GO.db")


#Loads necessary libraries
library(affy)
library(limma)
library(hgu133a.db)
library(annotate)
library(GO.db)

#*****Commented code block below moves into my .CEL file directory, loads the files, and normalizes using rma()*************

#Sets directory where .CEL files are located

datadir <- file.path('/Users/joewu/Desktop/Pycharm/Sleep_Dep_Controls')
setwd(datadir)

#Read raw .CEL file data contained in cd, normalize data using rma()
data <- ReadAffy()
eset <- rma(data)

check sample normalization
plotDensities(eset, legend=FALSE)

#***********************************************************************************************************************

#*****Commented code block below creates a complete annotation file included in the content directory I sent you along with expression set which you will use to construct the eset object again without .CEL files****

#Add Gene Annotation to Normalized Expression Set Output
my_frame <- data.frame(exprs(eset))
Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))

#Merge dataframes together
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)

#writes raw annotation data file including missing data (NA)
#write.table(all, file='annotated_data.txt', sep='/t')

#removes all rows that do not have any data associated with them (NA)
all <- all[complete.cases(all), ]

#writes cleaned annotation data file
#write.table(all, file='cleaned_annotation_file.tsv', sep='\t')

#writes expression set file
#exprs <- expres(eset)
#write.table(exprs, file='exprs_object.tsv, sep='\t')

#Corrects features in eset object so that it only has the features with their corresponding gene names
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "hgu133a.db")
fData(eset) <- data.frame(Symbol=Symbol)
HasSymbol <- !is.na(fData(eset)$Symbol)
eset <- eset[HasSymbol,]


#****************************************************************************************************************************************************************************************************************

#This code has been added for your convenience, typically I would use the previously created eset object in code above. However, since I have not included the raw .CEL files I am building new eset object from included files I sent
#sets directory
#cd <- getwd()
#setwd(cd)

#read in expression set
#exprs_input <- as.matrix(data.frame(read.delim(file='expression_set.tsv', sep='\t')))


#phenotype data
#I have included phenotypic data in file called targets.tsv, however I have included how I build the dataframe in R as well
#pheno_input <- data.frame(read.delim(file='pheno_data.tsv', sep='\t'))

patient <- factor(rep(c(1,2,3,4,5,6,7,8), each=2)) # patient ID
condition <- factor(rep(c('Post', 'Pre'), 8)) #Pre or Post Treatment
gender <- factor(c(rep('F', 8), rep('M', 8))) #gender
SSS <- c(6, 2, 3, 1, 3, 2, 5, 2, 5, 1, 3, 2, 2, 1, 5, 3) #response variable, NOT factored (continuous variable)
PVT <- c(339.67, 254.56, 423.33, 316.09, 640.13, 358.82, 321.15, 491.67, 338.99, 288.09, 261.96, 246.69, 276.48, 250.11, 267.14, 249.67) #Another response variable, not utilized in linear model (YET)

frame <- as.matrix(data.frame(patient, condition, gender, SSS, PVT))


#Build eset object from data files
#Initialize object with expression set
eset <- ExpressionSet(assayData = exprs_input)

#retrieve features (probe IDs) from initialized Expression Set object
ID <- featureNames(eset)

#retrieve gene symbol connecting to probe ID
Symbol <- getSYMBOL(ID, "hgu133a.db")

#apply gene symbol to feature data container in Expression Set Object
fData(eset) <- data.frame(Symbol=Symbol)

#Remove NA symbols
HasSymbol <- !is.na(fData(eset)$Symbol)

#Apply to Expression Set Object
eset <- eset[HasSymbol,]

#Can use this code below if the Entrez ID and Gene Name are also desired in output, feed INFO into pData(eset)
#then adds it to the feature container in the eset object so that it appears in the DE genes output
INFO <- select(hgu133a.db, ID, c('SYMBOL', 'ENTREZID', 'GENENAME')
fData(eset) <- INFO

#Apply phenotypic data to eset pData container
#pData(eset) <- pheno_input

#Build the design matrix
design <- model.matrix(~ patient + gender + condition * SSS)
design

#Fit eset to linear model
fit <- lmFit(eset, design)

#apply eBayes function
fit <- eBayes(fit)

#Dsiplay top differentially expressed genes, selecting for contrast conditionPre:SSS2, and applying Bonferroni correction
#Not currently sure if I am selecting the correct contrast
topTable(fit, coef='conditionPre:SSS', adjust='BH',)

#*****Everything below is a analysis tool for the DE gene data*******

#A volcano plot displays log fold changes on the x-axis versus a measure of statistical significance on the y-axis. Here the significance measure can be -log(p-value) or the B-statisticswhich give the posterior log-odds of differential expression.
volcanoplot(fit, coef = 'conditionPre:SSS', style = "p-value", highlight = 10, names = fit$genes$SYMBOL, hl.col="blue", xlab = "Log2 Fold Change", ylab = 'Log Odds', pch=16, cex=0.1)

#Subset the eset object, collect only the expression values of the DE genes

#Save the topTable Output into a variable
DE_Genes <- topTable(fit, coef='conditionPre:SSS', adjust='BH', 300)
#Take the gene IDS from the topTable ***(Requires run of code on line 119 to get the IDS in topTable output)
DE_IDS <- DE_Genes[1]

#Iterates through list of DE_IDS and saves all of the individual expressions mapped to that probe in a structure
for(id in DE_IDS){
	DE_subset <- exprs(eset[id])
	gene_sym <- getSYMBOL(id, 'hgu133a.db')
}

#New column names
patient_column <- c('P01_Post', 'P01_Pre', 'P02_Post', 'P02_Pre','P03_Post', 'P03_Pre', 'P04_Post', 'P04_Pre','P05_Post', 'P05_Pre', 'P06_Post', 'P06_Pre','P07_Post', 'P07_Pre', 'P08_Post', 'P08_Pre')

#Replaces the probe row names with gene row names 
row.names(DE_subset) <- gene_sym
#Renames the column names so it doesn't look as messy
colnames(DE_subset) <- patient_column


#Run heatmap analysis
coolmap(DE_subset, cluster.by='de pattern', linkage.row="complete", linkage.col="complete", show.dendrogram="both")

#Run Gene Ontology Analysis (GO Analysis)
#FDR has to be set at 0.06 because none of the DE genes made the threshold of 0.05
go <- goana(fit, coef='conditionPre:SSS', FDR=0.06, geneid='ENTREZID')
topGO(go)

#Building data file for GO Analysis
Symbol <- topTable['Symbol']
P.Value <- topTable['P.Value']

#Need to unlist the two variables and concatenate them back together
Symbol <- unlist(Symbol)
P.Value <- unlist(P.Value)

#Ready to be pumped into GO pipeline
genes <- setNames(P.Value, Symbol)

require(topGO)
require(org.Hs.eg.db)

#Code below retrieved from Kevin Blighe on BioStars, Link: https://www.biostars.org/p/350710/
selection <- function(allScore){ return(allScore < 0.055)} # function that returns TRUE/FALSE for p-values<0.05, changed to 0.055 b/c no genes were below 0.05 :(
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new('topGOdata', ontology='BP', allGenes=genes, annot=annFUN.GO2genes, GO2genes=allGO2genes, geneSel=selection, nodesize=10) #GO processing

#use the topGO vignette to run further analysis, plug and play
#in generative the p-value.elim vs. =-value.class comparison chart need the following function for var col map
# link to where I found function below: https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}



#Run KEGG analysis
kegg <- kegga(fit, coef='conditionPre:SSS', FDR=0.06, geneid='ENTREZID')
topKEGG(kegg)





