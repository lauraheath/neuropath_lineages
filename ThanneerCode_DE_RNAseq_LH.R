# ---
#   title: "Covariate and differential expression analysis of ROSMAP reprocessed counts (with CQN normalisation)"
# author: "Thanneer Perumal"
# output: html_documents
# ---
#   Date of analysis update: "`r date()`"
# 
devtools::install_github('th1vairam/CovariateAnalysis@dev')
devtools::install_github('brian-bot/githubr')
### Load Libraries

## It is assumed your working directory is where this file
## Load required libraries
library(CovariateAnalysis) 
library(data.table)
library(plyr)
library(tidyverse)
library(psych)
library(limma)
library(edgeR)
library(biomaRt)
library(RColorBrewer)
library(cqn)
library(synapser) 
library(githubr) 
synLogin()
library(doParallel)
library(foreach)
library(knitr)





### Data download
#Obtain count matrix and metadata from synapse

# Download expression data
COUNT_ID <- 'syn8691134';
ALL_USED_IDs <- COUNT_ID
COUNT_OBJ <- synGet(COUNT_ID)
COUNT <- read.table(COUNT_OBJ$path, header=T, sep='\t', check.names = F, row.names = 1)
COUNT[,grep('150_120419', colnames(COUNT))[2]] = NULL
# Convert rownames of counts from tracking id to ensemble gene id
tmp = data.frame(Gene.ID = rownames(COUNT)) %>%
  dplyr::mutate(ID = Gene.ID) %>%
  tidyr::separate(ID, c('ensembl_gene_id', 'position'), sep = '\\.')
rownames(tmp) = tmp$Gene.ID
rownames(COUNT) = tmp[rownames(COUNT), 'ensembl_gene_id']
# Get clinical metadata
METADATA.CLINICAL_ID <- 'syn3191087'
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.CLINICAL_ID
METADATA.CLINICAL_OBJ <- synGet(METADATA.CLINICAL_ID)
METADATA.CLINICAL <- read.table(METADATA.CLINICAL_OBJ$path,sep=',',header=T)
# Get clinical metadata with uncensored ages
METADATA.CLINICAL_ID1 <- 'syn7116000'
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.CLINICAL_ID1
METADATA.CLINICAL_OBJ1 <- synGet(METADATA.CLINICAL_ID1)
METADATA.CLINICAL1 <- read.table(METADATA.CLINICAL_OBJ1$path,sep=',',header=T)
# Get technical covariates
METADATA.TECH_ID <- 'syn4300313'
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.TECH_ID
METADATA.TECH_OBJ <- synGet(METADATA.TECH_ID)
METADATA.TECH <- read.table(METADATA.TECH_OBJ$path,sep='\t',header=T)

# Get picard metrics from synapse
METADATA.PICARD_ID <- 'syn8698240';
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.PICARD_ID
METADATA.PICARD <- synGet(METADATA.PICARD_ID)$path %>%
  data.table::fread() %>%
  dplyr::rename(Sampleid = sample)
colnames(METADATA.PICARD) = gsub('AlignmentSummaryMetrics__','',colnames(METADATA.PICARD))
colnames(METADATA.PICARD) = gsub('RnaSeqMetrics__','',colnames(METADATA.PICARD))
#two columns in METADATA.PICARD have the same name. Rename one:
colnames(METADATA.PICARD)[10] <- "PF_ALIGNED_BASE2"

# Fix error in technical covariates data
KEY_ID <- 'syn3382527'
ALL_USED_IDs[length(ALL_USED_IDs)+1] = KEY_ID  
KEY <- synGet(KEY_ID)$path %>%
  read.csv %>% 
  dplyr::filter(!is.na(rnaseq_id)) %>%
  dplyr::select(projid, rnaseq_id) %>%
  tidyr::separate(rnaseq_id, c('a','b','batch'), sep = '_') %>% 
  unite(Sampleid, a, b) %>%
  dplyr::select(-batch) %>%
  unique
# Match technical and clinical covariates
METADATA <- METADATA.TECH %>%
  dplyr::left_join(METADATA.PICARD) %>%
  dplyr::select(-projid) %>%
  dplyr::left_join(KEY) %>%
  dplyr::left_join(METADATA.CLINICAL) %>%
  dplyr::select(-age_first_ad_dx, -age_death, -age_at_visit_max) %>%
  dplyr::left_join(METADATA.CLINICAL1)
# Pick higher quality RIN batch for sample 492_120515
METADATA <- METADATA %>%
  dplyr::group_by(Sampleid) %>%
  dplyr::top_n(1, RINcontinuous)
# Get gene specific parameters from synapse
GENE.PARAM = synGet('syn8449369')$path %>%
  data.table::fread(data.table = FALSE)
ALL_USED_IDs = c(ALL_USED_IDs, 'syn8449369')
GENE.LEN = dplyr::select(GENE.PARAM, ensembl_gene_id, gene.length) %>% 
  unique() 
rownames(GENE.LEN) = GENE.LEN$ensembl_gene_id
GENE.GC = dplyr::select(GENE.PARAM, ensembl_gene_id, percentage_gc_content) %>% 
  unique() 
rownames(GENE.GC) = GENE.GC$ensembl_gene_id 


### Data preprocessing
# Remove samples with no cogdx, RIN, PMI scores and age_death
METADATA <- METADATA %>%
  ungroup %>%
  dplyr::filter(Sampleid %in% colnames(COUNT)) %>%
  dplyr::filter(!is.na(cogdx), !is.na(braaksc), !is.na(ceradsc)) %>%
  dplyr::filter(!is.na(RINcontinuous)) %>%
  dplyr::filter(!is.na(pmi)) %>%
  dplyr::filter(!is.na(PCT_INTRONIC_BASES)) %>%
  dplyr::filter(!is.na(age_death)) %>%
  as.data.frame()
# Add harmonised case-control status
METADATA$Diagnosis = 'OTHER'
METADATA$Diagnosis[METADATA$cogdx == 1 & METADATA$braaksc <= 3 & METADATA$ceradsc >= 3] = 'CONTROL'
METADATA$Diagnosis[METADATA$cogdx == 4 & METADATA$braaksc >= 4 & METADATA$ceradsc <= 2] = 'AD'
#recode braaksc for three levels
METADATA$braak_category = 'MED'
METADATA$braak_category[METADATA$braaksc <= 2] = 'LOW'
METADATA$braak_category[METADATA$braaksc  >= 5] = 'HIGH'
#recode CERAD score 
METADATA$cerad_category = 'AD'
METADATA$cerad_category[METADATA$ceradsc >= 3] = 'NO_AD'
# Add sex variable 
METADATA$Sex = 'FEMALE'
METADATA$Sex[METADATA$msex == 1] = 'MALE'
# Add apoe4 genotype (0, 1, 2)
METADATA$APOE4 = 0
METADATA$APOE4[METADATA$apoe_genotype %in% c(24, 34)] = 1
METADATA$APOE4[METADATA$apoe_genotype %in% c(44)] = 2
# METADATA$APOE4[is.na(METADATA$apoe_genotype)] = NA
# Get square of RIN
METADATA$RINcontinuous2 = METADATA$RINcontinuous^2
# Match covariates to expression data
indToRetain = intersect(METADATA$Sampleid, colnames(COUNT))
removedIDs = setdiff(colnames(COUNT), METADATA$Sampleid)
COUNT = COUNT[,indToRetain]
rownames(METADATA) = METADATA$Sampleid
METADATA = METADATA[indToRetain,]

#Dorsolateral prefrontal cortex of `r dim(COUNT)[2]` subjects from the ROS and MAP cohorts are used for the analysis. Following sample are removed due to missing metadata `r paste(removedIDs, collapse = ',')`




### Covariate clustering
#Determine relationship between covariates
primaryVariable <- c("cogdx", "Diagnosis", "APOE4", "braak_category", "ceradsc")
FactorCovariates <- c("Batch", "Sex", "race", "spanish", "cogdx", "Diagnosis", "APOE4", "braak_category", "ceradsc")
ContCovariates <- c("RINcontinuous", "RINcontinuous2", "age_death", "pmi", "educ", 
                    "PCT_PF_READS_ALIGNED", "PCT_CODING_BASES",
                    "PCT_INTERGENIC_BASES", "PCT_INTRONIC_BASES", 
                    "PCT_RIBOSOMAL_BASES")
# Find inter relation between factor covariates
COVARIATES = METADATA[,c(FactorCovariates,ContCovariates),drop=F]
COVARIATES <- data.frame(lapply(COVARIATES,function(x){x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
rownames(COVARIATES) <- METADATA$Sampleid
# Convert factor covariates to factors
COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.character)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.numeric)


#Correlation/association between covariates at an FDR <= 0.1
COVARIATES.CORRELATION = CovariateAnalysis::getAssociationStatistics(COVARIATES, PVAL = 0.05)
draw(COVARIATES.CORRELATION$plot, heatmap_legend_side = 'left', padding  = unit(c(18,2,2,18), 'mm'))


### Explore metadata
# RIN
p = list()
p[[1]] = ggplot(COVARIATES, aes(x = Diagnosis, y = RINcontinuous)) + geom_boxplot()
p[[1]] = p[[1]] + ggtitle('RIN') + theme(legend.position = 'top')
# AgeAtDeath
p[[2]] = ggplot(COVARIATES, aes(x = Diagnosis, y = age_death)) + geom_boxplot()
p[[2]] = p[[2]] + ggtitle('AgeOfDeath') + theme(legend.position = 'top')
# PMI
p[[3]] = ggplot(COVARIATES, aes(x = Diagnosis, y = pmi)) + geom_boxplot()
p[[3]] = p[[3]] + ggtitle('PMI') + theme(legend.position = 'top')
# Education
p[[4]] = ggplot(COVARIATES, aes(x = Diagnosis, y = educ)) + geom_boxplot()
p[[4]] = p[[4]] + ggtitle('Education') + theme(legend.position = 'top')
# Intronic bases
p[[5]] = ggplot(COVARIATES, aes(x = Diagnosis, y = PCT_INTRONIC_BASES)) + geom_boxplot()
p[[5]] = p[[5]] + ggtitle('Fraction Intronic Bases') + theme(legend.position = 'top')
# Ribosomal bases
p[[6]] = ggplot(COVARIATES, aes(x = Diagnosis, y = PCT_RIBOSOMAL_BASES)) + geom_boxplot()
p[[6]] = p[[6]] + ggtitle('Fraction Ribosomal Bases') + theme(legend.position = 'top')
multiplot(plotlist = p, cols = 2)


### Filter genes
#* Remove genes that have less than 1 cpm counts in at least 50% of samples per Diagnosis
#* Remove genes with missing gene length and percentage GC content

genesToAnalyze = COVARIATES %>%
  rownameToFirstColumn('Sampleid') %>%
  dlply(.(Diagnosis), .fun = function(mtd, count){
    processed.counts = getGeneFilteredGeneExprMatrix(count[,mtd$Sampleid],
                                                     MIN_GENE_CPM=1, 
                                                     MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.5)
    processed.counts$filteredExprMatrix$genes
  }, COUNT)
genesToAnalyze = unlist(genesToAnalyze) %>% 
  unique() %>%
  intersect(GENE.GC$ensembl_gene_id[!is.na(GENE.GC$percentage_gc_content)]) %>%
  intersect(GENE.LEN$ensembl_gene_id[!is.na(GENE.LEN$gene.length)]) %>%
  setdiff(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))
PROCESSED_COUNTS = getGeneFilteredGeneExprMatrix(COUNT[genesToAnalyze, ], 
                                                 MIN_GENE_CPM=0, 
                                                 MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0)
# Check gene biotype
## Define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "dec2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
## Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                       filters = "ensembl_gene_id", 
                       values = PROCESSED_COUNTS$filteredExprMatrix$genes,
                       mart = mart)
summary(factor(Ensemble2HGNC$gene_biotype)) %>%
  rownameToFirstColumn('Biotype') %>%
  dplyr::rename(fraction = DF) %>%
  dplyr::mutate(fraction = fraction/dim(PROCESSED_COUNTS$filteredExprMatrix$genes)[1]) %>%
  dplyr::filter(fraction >= 0.01) %>%
  kable

  #Processing `r dim(PROCESSED_COUNTS$filteredExprMatrix)[1]` genes in `r dim(PROCESSED_COUNTS$filteredExprMatrix)[2]` samples

### Library Normalisation
#Library normalisation is performed using cqn (conditional quantile normalisation)

# Compute offset for gene length and gc content
CQN.GENE_EXPRESSION = cqn(PROCESSED_COUNTS$filteredExprMatrix$counts, 
                          x = GENE.GC[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'percentage_gc_content'],
                          lengths = GENE.LEN[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'gene.length'],
                          lengthMethod = "smooth", 
                          verbose = TRUE)
CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$y + CQN.GENE_EXPRESSION$offset


### Outlier Analysis
#### Sample outliers
#Outlier analysis is performed before library normalisation with raw cpm counts

indToRemove = c('380_120503', '500_120515')
# Find principal components of expression to plot
PC <- prcomp(CQN.GENE_EXPRESSION$E, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])
plotdata <- dplyr::left_join(plotdata, rownameToFirstColumn(COVARIATES, 'SampleID')) %>%
  dplyr::mutate(label = SampleID)
plotdata$label[!(plotdata$SampleID %in% indToRemove)] = ''
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Batch, shape=Diagnosis, size=RINcontinuous))
p <- p + theme_bw() + theme(legend.position="right")
p <- p + geom_text(aes(label= label), size=4, hjust=0)
p
# Plot abberent distribution of logcpm counts
tmp1 = CQN.GENE_EXPRESSION$E %>%
  rownameToFirstColumn('Gene.ID') %>%
  tidyr::gather(SampleID, logCPM, -Gene.ID) %>%
  dplyr::left_join(COVARIATES %>%
                     rownameToFirstColumn('SampleID'))
p = ggplot(tmp1 %>%
             dplyr::filter(SampleID %in% indToRemove),
           aes(x = logCPM, color = SampleID)) + geom_density() 
p = p + theme(legend.position = 'top')
p
indToRetain = setdiff(colnames(PROCESSED_COUNTS$filteredExprMatrix$counts), indToRemove)
PROCESSED_COUNTS$filteredExprMatrix$counts = PROCESSED_COUNTS$filteredExprMatrix$counts[,indToRetain]
CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$E[,indToRetain]
COVARIATES = COVARIATES[indToRetain,]
tmp = COVARIATES %>%
  dplyr::group_by(Diagnosis, cogdx) %>%
  dplyr::summarise(count = n()) %>%
  tidyr::spread(Diagnosis, count)

#Processing `r dim(PROCESSED_COUNTS$filteredExprMatrix)[1]` genes in `r dim(PROCESSED_COUNTS$filteredExprMatrix)[2]` samples

#Based on the expression pattern following samples were tagged as outliers: `r paste(indToRemove, collapse = ', ')`

#Distribution of samples are: `r kable(tmp)`

#### Gene outliers
#Assign NA values to genes that are above and below 3 std deviation of its distribution

# Set gene counts in specific samples that are deviating 3 sd from other samples to 3SD limit
LOG.CPM = apply(CQN.GENE_EXPRESSION$E, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  
  x[x < (mn-3*std.dev)] = NA
  x[x > (mn+3*std.dev)] = NA
  return(x)
}) %>% t
CQN.GENE_EXPRESSION$E = LOG.CPM
CQN.GENE_EXPRESSION$E.no.na = CQN.GENE_EXPRESSION$E
CQN.GENE_EXPRESSION$E.no.na[is.na(CQN.GENE_EXPRESSION$E.no.na)] = 0
LIB.SIZE = colSums(PROCESSED_COUNTS$filteredExprMatrix$counts)
NEW.COUNTS = (2^LOG.CPM) * t(replicate(dim(LOG.CPM)[1], LIB.SIZE))/1e6


### Sample clustering
#PCA based clustering of samples
# Find principal components of expression to plot
PC <- prcomp(CQN.GENE_EXPRESSION$E.no.na, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])
plotdata <- left_join(plotdata, rownameToFirstColumn(COVARIATES, 'SampleID'))
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Batch, shape=cogdx, size=RINcontinuous))
p <- p + theme_bw() + theme(legend.position="right")
p


#Tree based clustering of samples
# Eucledian tree based analysis
COVARIATES.tmp = data.matrix(COVARIATES[,c("Batch", "Sex", "Diagnosis")])
COVARIATES.tmp[is.na(COVARIATES.tmp)] = 0
tree = hclust(as.dist(t(CQN.GENE_EXPRESSION$E.no.na)))
cols = WGCNA::labels2colors(COVARIATES.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(COVARIATES.tmp))


### Distribution of samples (log cpm)
# Plot abberent distribution of logcpm counts
tmp1 = CQN.GENE_EXPRESSION$E %>%
  rownameToFirstColumn('Gene.ID') %>%
  tidyr::gather(SampleID, logCPM, -Gene.ID) %>%
  left_join(COVARIATES %>%
              rownameToFirstColumn('SampleID'))
p = ggplot(tmp1, aes(x = logCPM, color = SampleID)) + geom_density() 
p = p + theme(legend.position = 'NONE') + facet_grid(.~Diagnosis, scale = 'free')
p


#Coexpression of genes 
cr = cor(t(CQN.GENE_EXPRESSION$E.no.na))
hist(cr, main = 'Distribution of correlation between genes', xlab = 'Correlation')


### Significant Covariates
#Correlation between pca of unadjusted mRNA expression and covariates are used to find significant covariates

# Find correlation between PC's of gene expression with covariates
preAdjustedSigCovars = runPCAandPlotCorrelations(CQN.GENE_EXPRESSION$E.no.na, 
                                                 COVARIATES,
                                                 'NULL design(voom-normalized)', 
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)


#Significant covariates to adjust at FDR 0.1 are `r preAdjustedSigCovars$significantCovars`
preAdjustedSigCovars[["PC_res"]][[2]]$plotData


### Normalisation (iterative design)
# Since many covariates are correlated, re-normalising and re-adjusting COUNTS with an iterative design matrix
# 1. Adding Batch and Sex a priori to variable selection
# 2. Primary variable of interest Diagnosis is excluded from the pool of available covariates for selection

# Primary variable of interest
postAdjustCovars = c('Batch', 'Sex');
# Assign residual covariates
residualCovars = setdiff(preAdjustedSigCovars$significantCovars, c(postAdjustCovars, primaryVariable))
residualSigCovars = preAdjustedSigCovars
covariatesEffects = preAdjustedSigCovars$Effects.significantCovars[residualCovars]
postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects))) %>% unique()
loopCount = 0 
while(length(residualSigCovars$significantCovars)!=0 && loopCount <= 20){
  writeLines(paste('Using following covariates in the model:',
                   paste(postAdjustCovars, collapse=', '),
                   'as fixed effects'))
  
  # Post adjusted design matrix
  DM1 = getDesignMatrix(COVARIATES[,postAdjustCovars,drop=F],Intercept = F)
  DM1$design = DM1$design[,linColumnFinder(DM1$design)$indepCols]
  
  # Estimate voom weights for dispersion control
  cnts = NEW.COUNTS
  cnts[is.na(cnts)] = 0
  VOOM.GENE_EXPRESSION = voom(cnts, 
                              design=DM1$design, 
                              plot=F)#,
                              #na.rm = T)
  
  # Fit linear model using new weights and new design
  VOOM.ADJUSTED.FIT = lmFit(CQN.GENE_EXPRESSION$E,
                            design = DM1$design,
                            weights = VOOM.GENE_EXPRESSION$weights)
  
  # Residuals after normalisation
  RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(VOOM.ADJUSTED.FIT,
                                                CQN.GENE_EXPRESSION$E)
  
  # Residual covariates to choose from
  residCovars <- setdiff(c(FactorCovariates,ContCovariates), c(postAdjustCovars, primaryVariable))
  
  # Find PC of residual gene expression and significant covariates that are highly correlated with PCs
  expr = RESIDUAL.GENE_EXPRESSION
  expr[is.na(expr)] = 0
  residualSigCovars = runPCAandPlotCorrelations(expr, 
                                                COVARIATES[, residCovars, drop=F], 
                                                'adjusted design(voom-normalized)',
                                                isKeyPlot=TRUE)
  
  # Add postadjusted covariates (if any)
  residCovars = setdiff(residualSigCovars$significantCovars, c(postAdjustCovars, primaryVariable))
  covariatesEffects = residualSigCovars$Effects.significantCovars[residCovars]
  
  postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects)))
  loopCount = loopCount + 1
}
modelStr <- paste(paste(gsub('_','\\\\_',postAdjustCovars), collapse=', '),
                  'as fixed effects')
tmp <- paste('Using following covariates in the final model:', modelStr)


### Sanity check
# Find PC of residual gene expression and significant covariates that are highly correlated with PCs
residualSigCovars = runPCAandPlotCorrelations(expr, 
                                              COVARIATES,
                                              'adjusted design(voom-normalized)',
                                              isKeyPlot=TRUE)
residualSigCovars[["PC_res"]][[2]]$plotData

#Coexpression of genes 
cr = cor(t(expr))
hist(cr, main = 'Distribution of correlation between genes', xlab = 'Correlation')

#PCA of residual data
# Find principal components of expression to plot
PC <- prcomp(expr, scale.=T, center = T)
# Plot first 4 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])
plotdata <- left_join(plotdata, rownameToFirstColumn(COVARIATES, 'SampleID'))
p <- ggplot(plotdata, aes(x=PC1, y=PC2)) 
p <- p + geom_point(aes(color=Batch, shape=Diagnosis, size=RINcontinuous))
p <- p + theme_bw() + theme(legend.position="right")
p

#Tree based clustering of residual data
# Eucledian tree based analysis
COVARIATES.tmp = data.matrix(COVARIATES[,c('Batch', 'Sex', primaryVariable)])
COVARIATES.tmp[is.na(COVARIATES.tmp)] = 0
tree = hclust(as.dist(t(expr)))
cols = WGCNA::labels2colors(COVARIATES.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(COVARIATES.tmp))


### Adjust data with covariates for Network Analysis
#Identified covariates are regressed out from the expression matrix for network analysis
# Get design matrix
DESIGN.NET = getDesignMatrix(COVARIATES[, postAdjustCovars, drop = F], Intercept = F)
DESIGN.NET = DESIGN.NET$design[,linColumnFinder(DESIGN.NET$design)$indepCols]
# Estimate voom weights for dispersion control
cnts = NEW.COUNTS
cnts[is.na(cnts)] = 0
VOOM.NET.WEIGHTS = voom(cnts, design=DESIGN.NET, plot=F)
# Fit linear model using new weights and new design
VOOM.NET.FIT = lmFit(CQN.GENE_EXPRESSION$E,
                     design = DESIGN.NET,
                     weights = VOOM.NET.WEIGHTS$weights)
# Residuals after normalisation
RESIDUAL.NET.GENE_EXPRESSION = residuals.MArrayLM(VOOM.NET.FIT,
                                                  CQN.GENE_EXPRESSION$E)





### Differential expression analysis (with cogdx as primary variable)
#Differential expression is performed on the primary variable by controlling for covariates identified above

# Interpretation of cogdx scores
# 1. NCI, No cognitive impairment (No impaired domains)
# 2. MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
# 3. MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
# 4. AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
# 5. AD, Alzheimer's disease AND another cause of CI (NINCDS POSS AD)
# 6. Other dementia. Other primary cause of dementia
# 
# Genes that are differentially expressed at an FDR <= 0.05 are

# Get design matrix

# DESIGN = getDesignMatrix(COVARIATES[, c(primaryVariable[1], postAdjustCovars), drop = F], Intercept = F)
# DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# # Estimate voom weights for dispersion control
# VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# # Fit linear model using new weights and new design
# FIT = lmFit(CQN.GENE_EXPRESSION$E,
#             design = DESIGN$design,
#             weights = VOOM.WEIGHTS$weights)
# # Fit contrast
# contrast = makeContrasts(contrasts=c("cogdx2-cogdx1",
#                                      "cogdx4-cogdx1",
#                                      "cogdx4-cogdx2"),
#                          levels = colnames(FIT$coefficients))
# FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
# FIT.CONTR = eBayes(FIT.CONTR)
# # Get differnetial expression
# DE = lapply(1:3, function(i, FIT){
#   topTable(FIT, coef=i, number = 50000, confint = T) %>%
#     rownameToFirstColumn('ensembl_gene_id')
# }, FIT.CONTR)
# names(DE) = colnames(contrast)
# DE = DE %>% 
#   rbindlist(idcol = 'Comparison') %>%
#   dplyr::left_join(GENE.PARAM %>%
#                      dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gc_content, gene.length) %>%
#                      unique()) %>%
#   dplyr::mutate(Study = 'ROSMAP', Region = 'DLPFC',
#                 Direction = logFC/abs(logFC),
#                 Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
#                 Direction = as.character(Direction))
# DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
# writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
# tmp = DE %>%
#   dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
#   group_by(Comparison, Direction) %>%
#   dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
#   spread(Direction, FDR_0_05_FC_1.2) 
# kable(tmp)
# p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
# p = p + scale_color_manual(values = c('green','grey','red'))
# p = p + facet_grid(.~Comparison, scales = 'fixed')
# p
# all.diff.exp = list(cogdx = DE)
# all.fit = list(cogdx = FIT)
# 
# # Store differential expression results
# rbindlist(all.diff.exp, use.names = T, fill = T, idcol = 'Model') %>%
#   write.table(file = 'cogdx_results.tsv', sep = '\t', row.names=F, quote=F)
# 
# 
# cogdx_results <- read.delim(file="cogdx_results.tsv")


### Differential expression analysis (with braak category as primary variable)
#Differential expression is performed on the primary variable by controlling for covariates identified above

# Interpretation of braak category
# LOW: braak stage 0, 1, or 2
# MED: braak stage 3 or 4
# HIGH: braak stage 5 or 6


# Get design matrix
DESIGN = getDesignMatrix(COVARIATES[, c(primaryVariable[4], postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(CQN.GENE_EXPRESSION$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("braak_categoryLOW-braak_categoryMED",
                                     "braak_categoryLOW-braak_categoryHIGH",
                                     "braak_categoryMED-braak_categoryHIGH"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differential expression
DE = lapply(1:3, function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = colnames(contrast)
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  dplyr::left_join(GENE.PARAM %>%
                     dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gc_content, gene.length) %>%
                     unique()) %>%
  dplyr::mutate(Study = 'ROSMAP', Region = 'DLPFC',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = list(braak_category = DE)
all.fit = list(braak_category = FIT)

# # Store differential expression results
# rbindlist(all.diff.exp, use.names = T, fill = T, idcol = 'Model') %>%
#   write.table(file = 'braak_results.tsv', sep = '\t', row.names=F, quote=F)
# 
# 
# braak_results <- read.delim(file="braak_results.tsv")




### Differential expression analysis (with CERAD score as primary variable)
# Differential expression is performed on the primary variable by controlling for covariates identified above
# 
# Interpretation of CERAD score
#  Definite AD: ceradsc = 1
#  Probable AD: ceradsc = 2
#  Possible AD: ceradsc = 3
#  No AD: ceradsc = 4
# 

# Get design matrix
DESIGN = getDesignMatrix(COVARIATES[, c(primaryVariable[5], postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(CQN.GENE_EXPRESSION$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("ceradsc1-ceradsc2",
                                     "ceradsc1-ceradsc3",
                                     "ceradsc1-ceradsc4",
                                     "ceradsc2-ceradsc3",
                                     "ceradsc2-ceradsc4",
                                     "ceradsc3-ceradsc4"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:6, function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = colnames(contrast)
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  dplyr::left_join(GENE.PARAM %>%
                     dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gc_content, gene.length) %>%
                     unique()) %>%
  dplyr::mutate(Study = 'ROSMAP', Region = 'DLPFC',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = c(all.diff.exp, list(ceradsc = DE))
all.fit = list(ceradsc = FIT)

# Store differential expression results
rbindlist(all.diff.exp, use.names = T, fill = T, idcol = 'Model') %>%
  write.table(file = 'neuropath_results.tsv', sep = '\t', row.names=F, quote=F)


neuropath_results <- read.delim(file="neuropath_results.tsv")



#store differential expression results in synapse space for neuropath-based lineage
write.csv(neuropath_results, file='RNAseq_DLPFCneuropath_DE_genes.csv', row.names = FALSE)
file <- synapser::File(path='RNAseq_DLPFCneuropath_DE_genes.csv', parentId='syn27254960')
file <- synapser::synStore(file)













