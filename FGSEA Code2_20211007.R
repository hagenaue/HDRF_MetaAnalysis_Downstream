#Gene Set Enrichment Analysis using the fGSEA package (fast GSEA) because the GSEA algorithm itself can't handle many permutations.
#Using effect size output from the HDRF hippocampal meta-analysis (male-only) - estimate ("average" hedge's g across depression models - the magnitude of effect divided by variability in the sample (sort-of) - so basically in units of standard deviation)
#This is a directional-based enrichment analysis - so we'll sDepressionModel if genes from a gene set tend to generally be upregulated or down-regulated using this analysis. If they show extreme effects in *both* directions we will not pick it up.

#Yusra Sannah & Megan Hagenauer
#2021-06-28

############################

#Required code packages:

#First requires installation - I tend to do this through the RStudio GUI

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

install.packages(c("fgsea", "plyr"))


#Then load.

library(fgsea)
library(plyr)

#######################

#Set working directory - a folder where your files related to this project will be located:
#You can do this through the Rstudio GUI, but as it prints out the code for it in the console save it so you remember

#example code:
setwd("~/Desktop/FGSEA")

#######################

# HDRF hippocampal meta-analysis results: read them in

#Ideally, I would place them in your working directory
#It helps to save the results as a .csv or tab-delimited text (requires different code)

#Replace file name in this code with your actual file name.

HC_MetaResults<-read.csv("HC_MetaResults_Fixed.csv", header=TRUE, stringsAsFactors = FALSE)
dim(HC_MetaResults)
#Tells you the dimensions of the dataset (nrow, ncol)
colnames(HC_MetaResults)
#Tells you the column names
# [1] "X"           "Symbol"      "Estimate"    "SE"          "Pvalue"      "FDR"         "MouseSymbol"
# [8] "RatSymbol"  

str(HC_MetaResults)
#Tells you the structure of the dataset - double-check that the variables are defined properly (e.g., numeric variables haven't become character or factors)
#data.frame':	14770 obs. of  12 variables:
#$ X          : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Symbol     : chr  "Slc27a1" "Pamr1" "Aloxe3" "Ephx2" ...
#$ Estimate   : num  -0.92 -0.784 -0.753 -1.455 -0.73 ...
#$ SE         : num  0.173 0.157 0.158 0.311 0.157 ...
# $ Pvalue     : num  1.03e-07 5.70e-07 1.79e-06 2.84e-06 3.21e-06 3.70e-06 5.33e-06 8.23e-06 1.17e-05 1.38e-05 ...
# $ FDR        : num  0.00151 0.00421 0.0088 0.0091 0.0091 ...
# $ MouseSymbol: chr  "Slc27a1" "Pamr1" "Aloxe3" "Ephx2" ...
# $ RatSymbol  : chr  "Slc27a1" "Pamr1" "Aloxe3" "Ephx2" ...

#######################

#First step is to pull out rows with missing gene symbols:
#For this code to work, replace "gene_symbol" with the name of the column that has gene symbols OR rename the column to gene_symbol

#Dealing with date genes: (DEALT BY HAND)

# 
# DateGeneDecoder<-read.delim("DateGeneDecoder.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
# str(DateGeneDecoder)
# # 'data.frame':	28 obs. of  3 variables:
# # $ OldSymbol: chr  "Sept1" "Sept2" "Sept3" "Sept4" ...
# # $ NewSymbol: chr  "Septin1" "Septin2" "Septin3" "Septin4" ...
# # $ ExcelDate: chr  "1-Sep" "2-Sep" "3-Sep" "4-Sep" ...
# 
# #How many of the gene symbols in the dataset are outdated "date gene" names:
# sum(HC_MetaResults$Symbol%in%DateGeneDecoder[,1])
# 
# #How many of the gene symbols in the dataset are NOT "date gene names:
# sum((HC_MetaResults$Symbol%in%DateGeneDecoder[,1])==FALSE)
# 
# #checking that the date genes are actually encoded as excel's version of dates:
# sum(HC_MetaResults$GeneSymbol%in%DateGeneDecoder[,3])
# 
# 
# HC_MetaResults$gene_symbol<-recode(HC_MetaResults$Symbol, Sept1='Septin1',
#                                    Sept2='Septin2',
#                                    Sept3='Septin3',
#                                    Sept4='Septin4',
#                                    Sept5='Septin5',
#                                    Sept6='Septin6',
#                                    Sept7='Septin7',
#                                    Sept8='Septin8',
#                                    Sept9='Septin9',
#                                    Sept10='Septin10',
#                                    Sept11='Septin11',
#                                    Sept12='Septin12',
#                                    Sept13='Septin13',
#                                    Sept14='Septin14',
#                                    Sept7p2='Septin7p2',
#                                    March1='Marchf1',
#                                    March2='Marchf2',
#                                    March3='Marchf3',
#                                    March4='Marchf4',
#                                    March5='Marchf5',
#                                    March6='Marchf6',
#                                    March7='Marchf7',
#                                    March8='Marchf8',
#                                    March9='Marchf9',
#                                    March10='Marchf10',
#                                    March11='Marchf11',
#                                    Dec1='Delec1',
#                                    Nov='Ccn3')

# #Saving the results with the date genes fixed:
# write.csv(HC_MetaResults, "HC_MetaResults_DateGenesFixed.csv")

#How many rows are missing gene symbols:
sum(is.na(HC_MetaResults$RatSymbol))
#0
sum(is.na(HC_MetaResults$MouseSymbol))
#0

#This takes them out:
HC_MetaResults_noNA<-HC_MetaResults[is.na(HC_MetaResults$RatSymbol)==FALSE,]
dim(HC_MetaResults_noNA)

#This determines if there duplicated symbols:
sum(duplicated(HC_MetaResults_noNA$RatSymbol))
#[1] 438



###################################


#This averages the results (effect estimate) for duplicated symbols:
#For this code to work, replace "Coef.DepressionModel" with whatever the name of the column is with the meta-analysis effect estimate

HC_DepressionModel_Betas_forGSEA<-tapply(X=HC_MetaResults_noNA$Estimate, INDEX=HC_MetaResults_noNA$RatSymbol, FUN=mean)
#This pulls out the gene symbols again and associates them with the averaged results (in the order of the results)
names(HC_DepressionModel_Betas_forGSEA)<-names(table(HC_MetaResults_noNA$RatSymbol))

head(HC_DepressionModel_Betas_forGSEA)

#How long the results are now after averaging...
length(HC_DepressionModel_Betas_forGSEA)
#[1] 14332

HC_DepressionModel_Betas_forGSEARanked<-HC_DepressionModel_Betas_forGSEA[order(HC_DepressionModel_Betas_forGSEA)]
head(HC_DepressionModel_Betas_forGSEARanked)

1925/14332

#########################


#Updated GMT with Gemma and GeneWeaver derived gene sets:

library(fgsea)
library(plyr)

GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrthologForYusra.gmt.txt")
str(GMT_ForRats)
# List of 15914 originally, then when we took some of them out so now it is 15892 
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                                                                                                                                                                   

#Note - I also increased the # permutations for this, because the output from the others was looking like there might not be enough.

HC_DepressionModels_FGSEAResults<-fgsea(GMT_ForRats, HC_DepressionModel_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(HC_DepressionModels_FGSEAResults)
str(HC_DepressionModel_Betas_forGSEARanked)
HC_DepressionModels_FGSEAResults$leadingEdge<-vapply(HC_DepressionModels_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))

write.csv(HC_DepressionModels_FGSEAResults, "HC_DepressionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunctionGemmaGeneWeaver.csv")

rm(HC_DepressionModels_FGSEAResults)
HC_DepressionModels_FGSEAResults<-read.csv("HC_DepressionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunctionGemmaGeneWeaver.csv", header=TRUE, stringsAsFactors = FALSE)
str(HC_DepressionModels_FGSEAResults)

###################################################



###################################################

#code to pull out results of particular types:


#Gene sets to pull out:

HypothesisPathways_toExtract<-read.csv("HypothesisPathways_toExtract3.csv", header=TRUE, stringsAsFactors = FALSE)
str(HypothesisPathways_toExtract)

HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults<-join(HypothesisPathways_toExtract, HC_DepressionModels_FGSEAResults , by="pathway", type="left")
write.csv(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults, "HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults.csv")

CellTypePathways_HC_toExtract<-read.csv("CellTypePathways_HC_toExtract.csv", header=TRUE, stringsAsFactors = FALSE)

CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults<-join(CellTypePathways_HC_toExtract, HC_DepressionModels_FGSEAResults, by="pathway", type="left")
write.csv(CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults, "CellTypePathways_HC_toExtract_vs_HC_DepressionModels_FGSEAResults.csv")












