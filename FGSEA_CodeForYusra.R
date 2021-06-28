#Gene Set Enrichment Analysis using the fGSEA package (fast GSEA) because the GSEA algorithm itself can't handle many permutations.
#Using effect size output from the HDRF hippocampal meta-analysis (male-only) - estimate ("average" hedge's g across depression models - the magnitude of effect divided by variability in the sample (sort-of) - so basically in units of standard deviation)
#This is a directional-based enrichment analysis - so we'll sDepressionModel if genes from a gene set tend to generally be upregulated or down-regulated using this analysis. If they show extreme effects in *both* directions we will not pick it up.

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

######################################

#Pulling out the common leading genes for a type of pathway from the hypothesis-driven pathways:

table(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway)

# Activity       Aggression       Enrichment      FearAnxiety FearConditioning   Glucocorticoid         Mood      OtherStress 
# 20                4                6               21               25               27               22               30 
# Social     SocialStress 
# 9               23

DepressionModel_LeadingEdge_FearConditioning<-HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$DepressionModelMain_leadingEdge[HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"]

DepressionModel_LeadingEdge_FearConditioning_Genes<-strsplit(stri_join_list(strsplit(DepressionModel_LeadingEdge_FearConditioning, ","), sep=" ", collapse=" "), " ")[[1]]

table(DepressionModel_LeadingEdge_FearConditioning_Genes)[order(table(DepressionModel_LeadingEdge_FearConditioning_Genes), decreasing=TRUE)]

TempTable<-table(DepressionModel_LeadingEdge_FearConditioning_Genes)[order(table(DepressionModel_LeadingEdge_FearConditioning_Genes), decreasing=TRUE)]
str(TempTable) 
names(TempTable)
TempTable[[1]]

sum(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning")
#[1] 25

TempTable_AsPercent<-(TempTable/25)*100

TempTable_DF<-data.frame("gene_symbol"=names(TempTable), "TotalGeneSets"=TempTable, "PercentGeneSets"=TempTable_AsPercent, stringsAsFactors = FALSE)

str(HC_MetaResults_noNA)

TempTable_DF_w_HC_MetaResults_noNA<-join(TempTable_DF, HC_MetaResults_noNA, by="gene_symbol", type="left")

write.csv(TempTable_DF_w_HC_MetaResults_noNA, "Table_DepressionModel_LeadingEdge_FearConditioning_Genes.csv")


#Or for a combination of types of pathways:



SD_LeadingEdge_AllStress<-HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$SDMain_leadingEdge[HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress"]

SD_LeadingEdge_AllStress_Genes<-strsplit(stri_join_list(strsplit(SD_LeadingEdge_AllStress, ","), sep=" ", collapse=" "), " ")[[1]]

table(SD_LeadingEdge_AllStress_Genes)[order(table(SD_LeadingEdge_AllStress_Genes), decreasing=TRUE)]

TempTable<-table(SD_LeadingEdge_AllStress_Genes)[order(table(SD_LeadingEdge_AllStress_Genes), decreasing=TRUE)]

sum(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress")

TempTable_AsPercent<-(TempTable/sum(HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress"))*100

TempTable_DF<-data.frame("gene_symbol"=names(TempTable), "TotalGeneSets"=TempTable, "PercentGeneSets"=TempTable_AsPercent, stringsAsFactors = FALSE)

str(HC_MetaResults_noNA)

TempTable_DF_w_HC_MetaResults_noNA<-join(TempTable_DF, HC_MetaResults_noNA, by="gene_symbol", type="left")

write.csv(TempTable_DF_w_HC_MetaResults_noNA, "Table_SD_LeadingEdge_AllStress_Genes.csv")

#Or using type of pathway and direction of effect:

#All stress, directional, from GMT

GMT_ForRats_AllStressUp<-GMT_ForRats[(names(GMT_ForRats)%in%HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$pathway[ (HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Up"])]

str(GMT_ForRats_AllStressUp)

GMT_ForRats_AllStressUp_Genes<-strsplit(stri_join_list(GMT_ForRats_AllStressUp, sep=" ", collapse=" "), " ")[[1]]

table(GMT_ForRats_AllStressUp_Genes)[order(table(GMT_ForRats_AllStressUp_Genes), decreasing=TRUE)]

TempTable<-table(GMT_ForRats_AllStressUp_Genes)[order(table(GMT_ForRats_AllStressUp_Genes), decreasing=TRUE)]

sum((HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Up")

TempTable_AsPercent<-(TempTable/sum((HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Up"))*100

TempTable_DF<-data.frame("gene_symbol"=names(TempTable), "TotalGeneSets"=TempTable, "PercentGeneSets"=TempTable_AsPercent, stringsAsFactors = FALSE)

str(HC_MetaResults_noNA)

TempTable_DF_w_HC_MetaResults_noNA<-join(TempTable_DF, HC_MetaResults_noNA, by="gene_symbol", type="left")

write.csv(TempTable_DF_w_HC_MetaResults_noNA, "Table_GMT_ForRats_AllStressUp_Genes.csv")


GMT_ForRats_AllStressDown<-GMT_ForRats[(names(GMT_ForRats)%in%HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$pathway[ (HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Down"])]

str(GMT_ForRats_AllStressDown)

GMT_ForRats_AllStressDown_Genes<-strsplit(stri_join_list(GMT_ForRats_AllStressDown, sep=" ", collapse=" "), " ")[[1]]

table(GMT_ForRats_AllStressDown_Genes)[order(table(GMT_ForRats_AllStressDown_Genes), decreasing=TRUE)]

TempTable<-table(GMT_ForRats_AllStressDown_Genes)[order(table(GMT_ForRats_AllStressDown_Genes), decreasing=TRUE)]

sum((HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Down")

TempTable_AsPercent<-(TempTable/sum((HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="SocialStress"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="FearConditioning"|HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$TypeOfPathway=="OtherStress") & HypothesisPathways_toExtract_vs_HC_DepressionModels_FGSEAResults$Direction=="Down"))*100

TempTable_DF<-data.frame("gene_symbol"=names(TempTable), "TotalGeneSets"=TempTable, "PercentGeneSets"=TempTable_AsPercent, stringsAsFactors = FALSE)

str(HC_MetaResults_noNA)

TempTable_DF_w_HC_MetaResults_noNA<-join(TempTable_DF, HC_MetaResults_noNA, by="gene_symbol", type="left")

write.csv(TempTable_DF_w_HC_MetaResults_noNA, "Table_GMT_ForRats_AllStressDown_Genes.csv")



#######################

#Pulling out which gene sets contain the top genes:
#Needs top genes entered instead of the ones that I have:

GeneSetsWithDepressionModelGenes<-lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5", "Pcdhga2", "Pcdhb6", "Pcdhb8", "Scn11a", "Dhrs2", "AABR07027137.1", "Megf8", "Mbd1", "Csf3r", "AABR07066700.1", "Ccl3")))
str(GeneSetsWithDepressionModelGenes)
write.csv(data.frame(names=names(GMT_ForRats), simplify2array(GeneSetsWithDepressionModelGenes)), "GeneSetsWithDepressionModelGenes.csv")

GeneSetsWithDepressionModelGenes_wInteractionModel<-lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5", "Pcdhga2", "Pcdhb6", "Pcdhb8", "Scn11a", "Dhrs2", "AABR07027137.1", "Megf8", "Mbd1", "Csf3r", "AABR07066700.1", "Ccl3", "Dele1", "Pcdhb7", "AC120486.2")))
str(GeneSetsWithDepressionModelGenes_wInteractionModel)
write.csv(data.frame(names=names(GMT_ForRats),simplify2array(GeneSetsWithDepressionModelGenes_wInteractionModel)), "GeneSetsWithDepressionModelGenes_wInteractionModel.csv")

DepressionModel_Pcdhb5<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5"))))
DepressionModel_Pcdhga2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhga2"))))
DepressionModel_Pcdhb6<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb6"))))
DepressionModel_Pcdhb8<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb8"))))
DepressionModel_Scn11a<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Scn11a"))))
DepressionModel_Dhrs2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Dhrs2"))))
DepressionModel_AABR07027137<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AABR07027137.1"))))
DepressionModel_Megf8<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Megf8"))))
DepressionModel_Mbd1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Mbd1"))))
DepressionModel_Csf3r<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Csf3r"))))
DepressionModel_AABR07066700<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AABR07066700.1"))))
DepressionModel_Ccl3<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Ccl3"))))
DepressionModel_Dele1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Dele1"))))
DepressionModel_Pcdhb7<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb7"))))
DepressionModel_AC120486<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AC120486.2"))))

head(data.frame(DepressionModel_Pcdhb5,DepressionModel_Pcdhga2,DepressionModel_Pcdhb6,DepressionModel_Pcdhb8,DepressionModel_Scn11a,DepressionModel_Dhrs2,DepressionModel_AABR07027137,DepressionModel_Megf8,DepressionModel_Mbd1,DepressionModel_Csf3r,DepressionModel_AABR07066700,DepressionModel_Ccl3, DepressionModel_Dele1,DepressionModel_Pcdhb7,DepressionModel_AC120486))

write.csv(data.frame(names(GMT_ForRats),DepressionModel_Pcdhb5,DepressionModel_Pcdhga2,DepressionModel_Pcdhb6,DepressionModel_Pcdhb8,DepressionModel_Scn11a,DepressionModel_Dhrs2,DepressionModel_AABR07027137,DepressionModel_Megf8,DepressionModel_Mbd1,DepressionModel_Csf3r,DepressionModel_AABR07066700,DepressionModel_Ccl3, DepressionModel_Dele1,DepressionModel_Pcdhb7,DepressionModel_AC120486), "GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene.csv")

GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene<-data.frame(pathway=names(GMT_ForRats),DepressionModel_Pcdhb5,DepressionModel_Pcdhga2,DepressionModel_Pcdhb6,DepressionModel_Pcdhb8,DepressionModel_Scn11a,DepressionModel_Dhrs2,DepressionModel_AABR07027137,DepressionModel_Megf8,DepressionModel_Mbd1,DepressionModel_Csf3r,DepressionModel_AABR07066700,DepressionModel_Ccl3, DepressionModel_Dele1,DepressionModel_Pcdhb7,DepressionModel_AC120486)

GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene_wDepressionModelSDResults<-join(HC_DepressionModels_FGSEAResults, GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene, by="pathway", type="left")

write.csv(GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene_wDepressionModelSDResults, "HC_GeneSetsWithDepressionModelGenes_wInteractionModel_ByGene_wDepressionModelSDResults.csv")



