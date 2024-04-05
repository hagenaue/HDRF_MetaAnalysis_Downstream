#Gene Set Enrichment Analysis using the fGSEA package (fast GSEA) because the GSEA algorithm itself can't handle many permutations.
#Using effect size output from the HDRF hippocampal meta-analysis (male-only) - estimate ("average" hedge's g across depression models - the magnitude of effect divided by variability in the sample (sort-of) - so basically in units of standard deviation)
#This is a directional-based enrichment analysis - so we'll see if genes from a gene set tend to generally be upregulated or down-regulated using this analysis. If they show extreme effects in *both* directions we will not pick it up.

#Yusra Sannah & Megan Hagenauer
#2021-06-28

############################

#Required code packages:

#First requires installation - I tend to do this through the RStudio GUI

#Then load.

library(fgsea)
library(plyr)

#######################

#Set working directory - a folder where your files related to this project will be located:
#You can do this through the Rstudio GUI, but as it prints out the code for it in the console save it so you remember

#example code:
setwd("~/Documents/Microarray Gen/Angela_HRLR_DepressionModel_Stress/HC_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")

#######################

# HDRF hippocampal meta-analysis results: read them in

#Ideally, I would place them in your working directory
#It helps to save the results as a .csv or tab-delimited text (requires different code)

#Replace file name in this code with your actual file name.

HC_MetaResults<-read.csv("HC_MetaResults.csv", header=TRUE, stringsAsFactors = FALSE)
dim(HC_MetaResults)
#Tells you the dimensions of the dataset (nrow, ncol)

colnames(HC_MetaResults)
#Tells you the column names

str(HC_MetaResults)
#Tells you the structure of the dataset - double-check that the variables are defined properly (e.g., numeric variables haven't become character or factors)


#######################

#First step is to pull out rows with missing gene symbols:
#For this code to work, replace "gene_symbol" with the name of the column that has gene symbols OR rename the column to gene_symbol


#Dealing with date genes:

DateGeneDecoder<-read.delim("DateGeneDecoder.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
str(DateGeneDecoder)
# 'data.frame':	28 obs. of  3 variables:
# $ OldSymbol: chr  "Sept1" "Sept2" "Sept3" "Sept4" ...
# $ NewSymbol: chr  "Septin1" "Septin2" "Septin3" "Septin4" ...
# $ ExcelDate: chr  "1-Sep" "2-Sep" "3-Sep" "4-Sep" ...

#How many of the gene symbols in the dataset are outdated "date gene" names:
sum(HC_MetaResults$gene_symbol%in%DateGeneDecoder[,1])

HC_MetaResults$gene_symbol<-recode(HC_MetaResults$gene_symbol, Sept1='Septin1',
                                                         Sept2='Septin2',
                                                         Sept3='Septin3',
                                                         Sept4='Septin4',
                                                         Sept5='Septin5',
                                                         Sept6='Septin6',
                                                         Sept7='Septin7',
                                                         Sept8='Septin8',
                                                         Sept9='Septin9',
                                                         Sept10='Septin10',
                                                         Sept11='Septin11',
                                                         Sept12='Septin12',
                                                         Sept13='Septin13',
                                                         Sept14='Septin14',
                                                         Sept7p2='Septin7p2',
                                                         March1='Marchf1',
                                                         March2='Marchf2',
                                                         March3='Marchf3',
                                                         March4='Marchf4',
                                                         March5='Marchf5',
                                                         March6='Marchf6',
                                                         March7='Marchf7',
                                                         March8='Marchf8',
                                                         March9='Marchf9',
                                                         March10='Marchf10',
                                                         March11='Marchf11',
                                                         Dec1='Delec1',
                                                         Nov='Ccn3')

#Saving the results with the date genes fixed:
write.csv(SHC_MetaResults, "HC_MetaResults_DateGenesFixed.csv")

#How many rows are missing gene symbols:
sum(is.na(HC_MetaResults$gene_symbol))
#[1] 766

#This takes them out:
HC_MetaResults_noNA<-HC_MetaResults[is.na(HC_MetaResults$gene_symbol)==FALSE,]
dim(HC_MetaResults_noNA)
#[1] 16999    28

#This determines if there duplicated symbols:
sum(duplicated(HC_MetaResults_noNA$gene_symbol))
#[1] 143


#Dealing with multiple gene mappings:

#First: Which of the gene symbols are things we don't recognize?

#from a glance:
#we have multi-mappings written out as mouse gene symbol - rat gene symbol like this:  H2-D1 - RT1-CE4 
#...but unfortunately, many of the actual rat gene symbols also contain dashes, like this: RT1-CE4, so we can't just delimit by hyphen, grab the second entry, and be done with it

#Reading in a database of up-to-date gene symbols:
Orthologs_Mice_Rats<-read.delim("Orthologs_Mice_Rats_Symbols_NoNA_InformaticsJax_OrthologDatabase_20210228.txt", sep="\t", header=TRUE, stingsAsFactors=FALSE)

#How many symbols aren't standard rat symbols?
dim(HC_MetaResults_noNA$gene_symbol[which(HC_MetaResults_noNA$gene_symbol%in%Orthologs_Mice_Rats$Symbol_Rat)==FALSE])

#Let's take a look in excel at them:

write.csv(HC_MetaResults_noNA$gene_symbol[which(HC_MetaResults_noNA$gene_symbol%in%Orthologs_Mice_Rats$Symbol_Rat)==FALSE], "HC_MetaResults_noNA_UnmappedSymbols.csv")


###Stop and consider what to do next...


###################################


#This averages the results (effect estimate) for duplicated symbols:
#For this code to work, replace "Coef.DepressionModel" with whatever the name of the column is with the meta-analysis effect estimate

HC_DepressionModel_Betas_forGSEA<-tapply(X=HC_MetaResults_noNA$Coef.DepressionModel, INDEX=HC_MetaResults_noNA$gene_symbol, FUN=mean)
#This pulls out the gene symbols again and associates them with the averaged results (in the order of the results)
names(HC_DepressionModel_Betas_forGSEA)<-names(table(HC_MetaResults_noNA$gene_symbol))

#How long the results are now after averaging...
length(HC_DepressionModel_Betas_forGSEA)
#[1] 16856

HC_DepressionModel_Betas_forGSEARanked<-HC_DepressionModel_Betas_forGSEA[order(HC_DepressionModel_Betas_forGSEA)]
head(HC_DepressionModel_Betas_forGSEARanked)



#########################


#Updated GMT with Gemma and GeneWeaver derived gene sets:


GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 15914
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                                                                                                                                                                       : chr [1:2082] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                                                                                                                                                                           : chr [1:2082] "Ada" "Gnpda1" "Zglp1" "Syce1l" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_DepressionModel_Stress/HC_Remove2outliers/11_FGSEA/GO_C5_wBrainCellTypeFunctionGemmaGeneWeaver")

#Note - I also increased the # permutations for this, because the output from the others was looking like there might not be enough.

HC_DepressionModels_FGSEAResults<-fgsea(GMT_ForRats, HC_DepressionModel_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(HC_DepressionModels_FGSEAResults)

HC_DepressionModels_FGSEAResults$leadingEdge<-vapply(HC_DepressionModels_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))

write.csv(HC_DepressionModels_FGSEAResults, "HC_DepressionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunctionGemmaGeneWeaver.csv")

rm(HC_DepressionModels_FGSEAResults)



###################################################



###################################################

#All the code after this point is fancy code to pull out results of particular types:


#Gene sets to pull out:

HypothesisPathways_toExtract<-read.csv("HypothesisPathways_toExtract3.csv", header=TRUE, stringsAsFactors = FALSE)
str(HypothesisPathways_toExtract)

HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults<-join(HypothesisPathways_toExtract, HC_DepressionModels_FGSEAResults, by="pathway", type="left")
write.csv(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults, "HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults.csv")

CellTypePathways_HC_toExtract<-read.csv("CellTypePathways_HC_toExtract.csv", header=TRUE, stringsAsFactors = FALSE)
str(CellTypePathways_STR_toExtract)

CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults<-join(CellTypePathways_HC_toExtract, HC_DepressionModels_FGSEAResults, by="pathway", type="left")
write.csv(CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults, "CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults.csv")



