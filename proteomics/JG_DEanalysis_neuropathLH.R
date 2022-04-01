devtools::install_github("kassambara/ggpubr")
library(ggpubr)


## Import Data
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
#Normalized <- TMT_Express_Load('syn21266453', 1)
Meta <- TMT_Express_Load('syn21323404', 1)
# - Public Facing BioSpecimin Data: syn21323366
# - Staged BioSpecimin Data: syn23583548
BioSpecimin <- TMT_Express_Load('syn21323366', 0)
if( "assay" %in% colnames(BioSpecimin) ){
}else{
  BioSpecimin <- TMT_Express_Load('syn23583548', 0)
}
BioSpecimin <- BioSpecimin[ BioSpecimin$assay == 'TMT quantitation', ]

Meta$specimenID <- row.names(Meta)
Meta <- dplyr::left_join(Meta, BioSpecimin, by = 'specimenID' )
Clinical <- TMT_Express_Load( 'syn3191087',0 )
Meta <- dplyr::left_join(Meta, Clinical, by = 'individualID' )
Meta <- Meta[ ,colSums(is.na(Meta))<nrow(Meta) ] 
row.names( Meta ) <- Meta$batchChannel
Meta <- Meta[ colnames(Log2_Normalized), ]
Meta <- Meta[ ,colnames(Meta)[ (colnames(Meta) %in%'controlType' )==F] ]
# Harmonize case-control status
Meta$braaksc <- as.numeric(Meta$braaksc)
Meta$ceradsc <- as.numeric(Meta$ceradsc)
Meta$cogdx <- as.numeric(Meta$cogdx)
# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"
kableExtra::kable( table(Meta$diagnosis) )
Meta$braak_category = 'MED'
Meta$braak_category[Meta$braaksc <= 2] = 'LOW'
Meta$braak_category[Meta$braaksc  >= 5] = 'HIGH'
## Add Ages over 90 for modeling
Mast <- TMT_Express_Load('syn23573928', 0)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')
#Convert APOE
Meta$apoe_genotype <- as.numeric( Meta$apoe_genotype )
APOS <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOS) <- APOS
APOS[names( table(Meta$apoe_genotype) )] <- 0
APOS[grepl("4",names(APOS))] <- 1
APOS[grepl("44",names(APOS))] <- 2
Meta$APOE <- as.numeric( APOS[ as.character(Meta$apoe_genotype) ] )
## Winzorize Expression Data
sink_preWinzor<- Log2_Normalized
for( i in 1:dim(Log2_Normalized)[1] ){
  Log2_Normalized[i,] <- DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) 
}
row.names(Meta) <- Meta$batchChannel
#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta$APOE <- as.factor(Meta$APOE)
Meta$ceradsc <- as.factor(Meta$ceradsc)
####Meta for Diagnosis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
Meta_D$diagnosis <- factor(Meta_D$diagnosis, levels = c("control", "AD"))

#update gene | uniprot identifiers (per Jake: some genes changed after checking ensembl gene ids)
p <- synapser::synGet('syn24216770')
correct_geneIDs <- read.csv(p$path)
Log2_Normalized$OldPeptideID <- rownames(Log2_Normalized)
Log2_Normalized <- dplyr::left_join(Log2_Normalized, correct_geneIDs, by="OldPeptideID")
rownames(Log2_Normalized) <- Log2_Normalized$NewPeptideID
Log2_Normalized$Old_Gene<-NULL
Log2_Normalized$Old_Pep<-NULL
Log2_Normalized$OldPeptideID<-NULL
Log2_Normalized$New_Gene<-NULL
Log2_Normalized$New_Pep<-NULL
Log2_Normalized$NewPeptideID<-NULL
Log2_Normalized$ENSG<-NULL

# save the log2-normalized matrix and metadata for downstream analyses:
saveRDS(Log2_Normalized, file="~/neuropath_lineages/proteomics/Dx_analysis/Log2_Normalized_winsorized.rds")
saveRDS(Meta, file="~/neuropath_lineages/proteomics/Dx_analysis/Meta.rds")


##
## Diagnosis Differential expression
Diagnosis_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_D) ],
                                                                                MET = Meta_D[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE' ) ],
                                                                                Vars = c( "pmi", "diagnosis", "APOE" ),
                                                                                Diag = "diagnosis",
                                                                                SampleID = 'batchChannel'
)
)
}) })


#sink<-Diagnosis_PVals

Diagnosis_Associations <- Cleaner(Diagnosis_PVals)



#for comparison to old analysis:
#upload JG's Non collapsed comparisons in synapse (from January 29, 2021)
p <- synapser::synGet('syn24216765')
JGdegs <- read.csv(p$path)
JGdegs_sig <- subset(JGdegs, JGdegs$PVal<0.05)
dim(JGdegs_sig)
dim(JGdegs)



# add ensemblIDs back to results
ensembls <- subset(correct_geneIDs, select=c(New_Pep, ENSG))
names(ensembls)[names(ensembls) == "New_Pep"] <- "ProtID"
Diagnosis_Associations <- dplyr::left_join(Diagnosis_Associations, ensembls)



corrs <- subset(Diagnosis_Associations, select=c(GeneID, Coefficient))

names(JGdegs)[names(JGdegs) == "Coefficient"] <- "Coefficient_orig"
corrs$Coefficient_orig <- JGdegs$Coefficient_orig
x <- corrs$Coefficient
y <- corrs$Coefficient_orig

cor(x, y, method="pearson")
cor.test(x,y,method="pearson")
cor(x, y, method="spearman")
cor.test(x,y,method="spearman")
ggscatter(corrs, x = "Coefficient", y = "Coefficient_orig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "new", ylab = "old")


x <- Diagnosis_Associations$PVal
y <- JGdegs$PVal
cor.test(x,y,method="pearson")
cor.test(x,y,method="spearman")


#overlap in genes selected:
oldgenes <- subset(JGdegs, JGdegs$PVal<0.05)
newgenes <- subset(Diagnosis_Associations, Diagnosis_Associations$PVal<0.05)
newgenes2 <- subset(Diagnosis_Associations, Diagnosis_Associations$FDR_PVal<0.1)

old1 <- oldgenes$Peptide
new1 <- newgenes$Peptide
new2 <- newgenes2$Peptide

#library(VennDiagram)
venn.diagram(
  x=list(old1,new2),
  category.names = c("Original", "Update"),
  filename = 'DEgene_overlap.png',
  output=TRUE
)

#save the DE results for trajectory analysis:
#store pseudotimes and states in synapse space for neuropath-based lineage
write.csv(Diagnosis_Associations, file="~/neuropath_lineages/proteomics/Dx_analysis/protDE_by_diagnosis.csv", row.names=FALSE)
file <- synapser::File(path='~/neuropath_lineages/proteomics/Dx_analysis/protDE_by_diagnosis.csv', parentId='syn28712340')
file <- synapser::synStore(file)



#### Determine differentially expressed proteins by neuropath: braak and cerad
#Braak comparisons: class into LOW, MED, HIGH (see below) and run all possible pairwise comparisons (LOW vs MED,
# MED vs HIGH, LOW vs HIGH)
#CERAD: run all comparisons for four categories (1 vs 2, 1 vs 3, 1 vs 4, 2 vs 3, 2 vs 4, 3 vs 4)

### BRAAK DE ####
#create Braak pairwise comparison sub-datasets
Meta_LOWMED <- Meta[ Meta$braak_category %in% c('LOW','MED'), ]
table(Meta_LOWMED$braak_category)
Meta_LOWMED$braak_category <- factor(Meta_LOWMED$braak_category, levels = c("LOW", "MED"))

Meta_LOWHIGH <- Meta[ Meta$braak_category %in% c('LOW','HIGH'), ]
table(Meta_LOWHIGH$braak_category)
Meta_LOWHIGH$braak_category <- factor(Meta_LOWHIGH$braak_category, levels = c("LOW", "HIGH"))

Meta_MEDHIGH <- Meta[ Meta$braak_category %in% c('MED','HIGH'), ]
table(Meta_MEDHIGH$braak_category)
Meta_MEDHIGH$braak_category <- factor(Meta_MEDHIGH$braak_category, levels = c("MED", "HIGH"))

## Diagnosis Differential expression
braak_LOWMED_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_LOWMED) ],
                                                                                MET = Meta_LOWMED[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                Vars = c( "pmi", "braak_category", "APOE" ),
                                                                                Diag = "braak_category",
                                                                                SampleID = 'batchChannel'
)
)
}) })


braak_LOWMED_Associations <- Cleaner(braak_LOWMED_PVals)
braak_LOWMED_Associations$category <- 'braak_category'
braak_LOWMED_Associations$comparison <- 'braak_LOWMED'
neuropath_DEgenes <- braak_LOWMED_Associations

braak_LOWHIGH_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_LOWHIGH) ],
                                                                                   MET = Meta_LOWHIGH[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                   Vars = c( "pmi", "braak_category", "APOE" ),
                                                                                   Diag = "braak_category",
                                                                                   SampleID = 'batchChannel'
)
)
}) })


braak_LOWHIGH_Associations <- Cleaner(braak_LOWHIGH_PVals)
braak_LOWHIGH_Associations$category <- 'braak_category'
braak_LOWHIGH_Associations$comparison <- 'braak_LOWHIGH'
neuropath_DEgenes <- rbind(neuropath_DEgenes, braak_LOWHIGH_Associations)

braak_MEDHIGH_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_MEDHIGH) ],
                                                                                   MET = Meta_MEDHIGH[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                   Vars = c( "pmi", "braak_category", "APOE" ),
                                                                                   Diag = "braak_category",
                                                                                   SampleID = 'batchChannel'
)
)
}) })


braak_MEDHIGH_Associations <- Cleaner(braak_MEDHIGH_PVals)
braak_MEDHIGH_Associations$category <- 'braak_category'
braak_MEDHIGH_Associations$comparison <- 'braak_MEDHIGH'
neuropath_DEgenes <- rbind(neuropath_DEgenes, braak_MEDHIGH_Associations)

#### CERAD score DE analysis
# Interpretation of CERAD score
#  Definite AD: ceradsc = 1
#  Probable AD: ceradsc = 2
#  Possible AD: ceradsc = 3
#  No AD: ceradsc = 4
#reset levels so that 4 is lowest, 1 is highest, to ensure the less-affected cerad score is the referent
Meta_2vs1 <- Meta[ Meta$ceradsc %in% c('1','2'), ]
table(Meta_2vs1$ceradsc)
Meta_2vs1$ceradsc <- factor(Meta_2vs1$ceradsc, levels = c("2", "1"))

Meta_3vs1 <- Meta[ Meta$ceradsc %in% c('1','3'), ]
table(Meta_3vs1$ceradsc)
Meta_3vs1$ceradsc <- factor(Meta_3vs1$ceradsc, levels = c("3", "1"))

Meta_4vs1 <- Meta[ Meta$ceradsc %in% c('1','4'), ]
table(Meta_4vs1$ceradsc)
Meta_4vs1$ceradsc <- factor(Meta_4vs1$ceradsc, levels = c("4", "1"))

Meta_3vs2 <- Meta[ Meta$ceradsc %in% c('3','2'), ]
table(Meta_3vs2$ceradsc)
Meta_3vs2$ceradsc <- factor(Meta_3vs2$ceradsc, levels = c("3", "2"))

Meta_4vs2 <- Meta[ Meta$ceradsc %in% c('4','2'), ]
table(Meta_4vs2$ceradsc)
Meta_4vs2$ceradsc <- factor(Meta_4vs2$ceradsc, levels = c("4", "2"))

Meta_4vs3 <- Meta[ Meta$ceradsc %in% c('4','3'), ]
table(Meta_4vs3$ceradsc)
Meta_4vs3$ceradsc <- factor(Meta_4vs3$ceradsc, levels = c("4", "3"))

cerad_2vs1_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_2vs1) ],
                                                                                    MET = Meta_2vs1[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                    Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                    Diag = "ceradsc",
                                                                                    SampleID = 'batchChannel'
)
)
}) })


cerad_2vs1_Associations <- Cleaner(cerad_2vs1_PVals)
cerad_2vs1_Associations$category <- 'ceradsc'
cerad_2vs1_Associations$comparison <- 'cerad2vs1'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_2vs1_Associations)

cerad_3vs1_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_3vs1) ],
                                                                                 MET = Meta_3vs1[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                 Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                 Diag = "ceradsc",
                                                                                 SampleID = 'batchChannel'
)
)
}) })


cerad_3vs1_Associations <- Cleaner(cerad_3vs1_PVals)
cerad_3vs1_Associations$category <- 'ceradsc'
cerad_3vs1_Associations$comparison <- 'cerad3vs1'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_3vs1_Associations)

cerad_4vs1_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_4vs1) ],
                                                                                 MET = Meta_4vs1[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                 Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                 Diag = "ceradsc",
                                                                                 SampleID = 'batchChannel'
)
)
}) })


cerad_4vs1_Associations <- Cleaner(cerad_4vs1_PVals)
cerad_4vs1_Associations$category <- 'ceradsc'
cerad_4vs1_Associations$comparison <- 'cerad4vs1'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_4vs1_Associations)

cerad_3vs2_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_3vs2) ],
                                                                                 MET = Meta_3vs2[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                 Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                 Diag = "ceradsc",
                                                                                 SampleID = 'batchChannel'
)
)
}) })


cerad_3vs2_Associations <- Cleaner(cerad_3vs2_PVals)
cerad_3vs2_Associations$category <- 'ceradsc'
cerad_3vs2_Associations$comparison <- 'cerad3vs2'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_3vs2_Associations)

cerad_4vs2_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_4vs2) ],
                                                                                 MET = Meta_4vs2[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                 Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                 Diag = "ceradsc",
                                                                                 SampleID = 'batchChannel'
)
)
}) })


cerad_4vs2_Associations <- Cleaner(cerad_4vs2_PVals)
cerad_4vs2_Associations$category <- 'ceradsc'
cerad_4vs2_Associations$comparison <- 'cerad4vs2'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_4vs2_Associations)

cerad_4vs3_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(1:dim(Log2_Normalized)[1], Logistic_Model, EXP=Log2_Normalized[ ,row.names(Meta_4vs3) ],
                                                                                 MET = Meta_4vs3[, c('SampleID', 'batchChannel', 'pmi', 'ceradsc', 'cogdx', 'dcfdx_lv', 'diagnosis', 'age_death', 'APOE', 'braak_category' ) ],
                                                                                 Vars = c( "pmi", "ceradsc", "APOE" ),
                                                                                 Diag = "ceradsc",
                                                                                 SampleID = 'batchChannel'
)
)
}) })


cerad_4vs3_Associations <- Cleaner(cerad_4vs3_PVals)
cerad_4vs3_Associations$category <- 'ceradsc'
cerad_4vs3_Associations$comparison <- 'cerad4vs3'
neuropath_DEgenes <- rbind(neuropath_DEgenes, cerad_4vs3_Associations)


#save the DE results for trajectory analysis:
#store pseudotimes and states in synapse space for neuropath-based lineage
write.csv(neuropath_DEgenes, file="~/neuropath_lineages/proteomics/prot_DLPFCneuropath_DE_genes.csv", row.names=FALSE)
file <- synapser::File(path='~/neuropath_lineages/proteomics/prot_DLPFCneuropath_DE_genes.csv', parentId='syn28712340')
file <- synapser::synStore(file)
