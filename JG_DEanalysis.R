devtools::install_github("kassambara/ggpubr")
library(ggpubr)


## Import Data
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
Normalized <- TMT_Express_Load('syn21266453', 1)
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
####Meta for Diagnosis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)

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

Diagnosis_Associations_new <- Cleaner(Diagnosis_PVals)

#upload JG's Non collapsed comparisons in synapse (from January 29, 2021)
p <- synapser::synGet('syn24216765')
JGdegs <- read.csv(p$path)
JGdegs_sig <- subset(JGdegs, JGdegs$PVal<0.05)
dim(JGdegs_sig)
dim(JGdegs)

corrs <- subset(Diagnosis_Associations_new, select=c(GeneID, Coefficient))

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


x <- Diagnosis_Associations_new$PVal
y <- JGdegs$PVal
cor.test(x,y,method="pearson")
cor.test(x,y,method="spearman")











# BRAAK - PMI, Diagnosis, APOE, 
## Braak Stage is a semiquantitative measure of severity of neurofibrillary tangle (NFT) pathology. 
BRAKK_PVals <-  suppressMessages({ suppressWarnings({ do.call( rbind, lapply(rownames(Log2_Normalized), NeuroPath_Calc, Exp=Log2_Normalized,
                                                                             MET=Meta,
                                                                             Path='braaksc'
)
)
}) })
BRAKK_Association <- Cleaner( BRAKK_PVals )



#upload JG's Non collapsed comparisons in synapse (from January 29, 2021)
p <- synapser::synGet('syn24216766')
JGdegs_br <- read.csv(p$path)
corrs <- subset(BRAKK_Association, select=c(GeneID, Coefficient))

names(JGdegs_br)[names(JGdegs_br) == "Coefficient"] <- "Coefficient_orig"
corrs$Coefficient_orig <- JGdegs_br$Coefficient_orig
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





# - pmi,diagnosis ( 13, 9)
# - pmi,diagnosis,APOE (12  9)
# - pmi,diagnosis,APOE,Age_death (   9 )
#tmp<-BRAKK_Association
# CERAD
## CERAD score is a semiquantitative measure of neuritic plaques.
CERAD_PVals <-  suppressMessages({ suppressWarnings({ do.call( rbind, lapply(rownames(Log2_Normalized), NeuroPath_Calc, Exp=Log2_Normalized,
                                                                             MET=Meta,
                                                                             Path='ceradsc'
)
)
}) })
CERAD_Association <- Cleaner( CERAD_PVals )
#168 without Age_of_death
# CogDx
## Final consensus cognitive diagnosis
# Value Coding
# 1     NCI: No cognitive impairment (No impaired domains)
# 2     MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
# 3     MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
# 4     AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)
# 5     AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
# 6     Other dementia: Other primary cause of dementia
## Model Uses 1, 2, and 4 only ( 374 out of 400 )
COGDX_PVals <- suppressMessages({ suppressWarnings({ 
  do.call( rbind, lapply(rownames(Log2_Normalized), NeuroPath_Calc, Exp=Log2_Normalized,
                         MET=Meta,
                         Path='cogdx'
  )
  )
}) })
COGDX_Association <- Cleaner( COGDX_PVals )
# 600+ without age_of_death
# 0   with age_of_death
# Dcfdx
## Clinical diagnosis of cognitive status [ lv = Last Valid Score ]
# 1 NCI: No cognitive impairment
# 2 MCI: Mild cognitive impairment, no other condition contributing to CI
# 3 MCI+: Mild cognitive impairment AND another condition contributing to CI
# 4 AD: Alzheimer’s dementia, no other condition contributing to CI (NINCDS/ADRDA Probable AD)
# 5 AD+: Alzheimer’s dementia AND other condition contributing to CI (NINCDS/ADRDA Possible AD)
# 6 Other dementia: Other primary cause of dementia, no clinical evidence of Alzheimer’s dementia
## - Model Uses 1,2,4
dcfdx_PVals <- suppressMessages({ suppressWarnings({ do.call( rbind, lapply(rownames(Log2_Normalized), NeuroPath_Calc, Exp=Log2_Normalized,
                                                                            MET=Meta,
                                                                            Path='dcfdx_lv'
)
)
}) })
dcfdx_Association <- Cleaner( dcfdx_PVals )



##### PULL ENSGs and ENSTs
name_trans <- synapser::synTableQuery( query='SELECT ENSG, GeneName FROM syn24168007' , resultsAs="csv" )$asDataFrame()[,c('ENSG', 'GeneName')]
head( dcfdx_Association[ !(dcfdx_Association$GeneID %in% name_trans$GeneName ),])
base <- dcfdx_Association[,  c("Peptide", "GeneID", "ProtID") ]
#### Pulls 7755 out of 8817
name_trans_slim <- name_trans[ as.character(name_trans$GeneName) %in% as.character(dcfdx_Association$GeneID), ]
name_trans_slim[ name_trans_slim$GeneName %in% 'DHRS11', ]$ENSG <- 'ENSG00000278535'
table(name_trans_slim$GeneName)[ table(name_trans_slim$GeneName) > 2 ]
name_trans_slim$Unis <- paste0( name_trans_slim$ENSG, '_', name_trans_slim$GeneName)
name_trans_slim <- name_trans_slim[ !duplicated( name_trans_slim$Unis ),]
name_trans_slim <- name_trans_slim[ !(name_trans_slim$ENSG == 'ENSG00000184779'), ]
table(table(name_trans_slim$GeneName))
table(table(name_trans_slim$ENSG))
Probs <- dcfdx_Association[ !( as.character(dcfdx_Association$GeneID) %in% as.character(name_trans$GeneName) ), 1:3 ]
#mart=biomaRt::useDataset("uniprot",mart=uniProt)
uniProt <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
attempt <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name","uniprotswissprot"),filter="uniprotswissprot",values=Probs$ProtID
                          ,mart=uniProt)
row.names(Probs) <- paste0( Probs$GeneID,'|', Probs$ProtID )
attempt$ID <- paste0( attempt$external_gene_name,'|', attempt$uniprotswissprot )
#attempt[ attempt$ID =='VARS1|P26640', ]$ensembl_gene_id <- 'ENSG00000204394'
#attempt[ attempt$ID =='HNRNPCL4|P0DMR1', ]$ensembl_gene_id <- 'ENSG00000179412'
#attempt[ attempt$ID =='ADGRL1|O94910', ]$ensembl_gene_id <- 'ENSG00000072071'
attempt[ attempt$ID =='AK6|Q9Y3D8', ]$ensembl_gene_id <- 'ENSG00000085231'
attempt[ attempt$ID =='ARHGAP23|Q9P227', ]$ensembl_gene_id <- 'ENSG00000275832'
#attempt[ attempt$ID =='ATP23|Q9Y6H3', ]$ensembl_gene_id <- 'ENSG00000166896'
#attempt[ attempt$ID =='ATP5IF1|Q9UII2', ]$ensembl_gene_id <- 'ENSG00000130770'
#attempt[ attempt$ID =='BMERB1|Q96MC5', ]$ensembl_gene_id <- 'ENSG00000166780'
#attempt[ attempt$ID =='BORCS5|Q969J3', ]$ensembl_gene_id <- 'ENSG00000165714'
#attempt[ attempt$ID =='CARMIL3|Q8ND23', ]$ensembl_gene_id <- 'ENSG00000186648'
#attempt[ attempt$ID =='DENND11|A4D1U4', ]$ensembl_gene_id <- 'ENSG00000257093'
attempt[ attempt$ID =='DUSP14|O95147', ]$ensembl_gene_id <- 'ENSG00000276023'
#attempt[ attempt$ID =='GON7|Q9BXV9', ]$ensembl_gene_id <- 'ENSG00000170270'
#attempt[ attempt$ID =='H2BU1|Q8N257', ]$ensembl_gene_id <- 'ENSG00000196890'
attempt[ attempt$ID =='IGHV1-69|P01742', ]$ensembl_gene_id <- 'ENSG00000211973'
attempt[ attempt$ID =='IGHV3-11|P01762', ]$ensembl_gene_id <- 'ENSG00000211941'
attempt[ attempt$ID =='IGHV3-15|A0A0B4J1V0', ]$ensembl_gene_id <- 'ENSG00000211943'
attempt[ attempt$ID =='IGHV3-74|A0A0B4J1X5', ]$ensembl_gene_id <- 'ENSG00000224650'
attempt[ attempt$ID =='IGHV5-51|A0A0C4DH38', ]$ensembl_gene_id <- 'ENSG00000211966'
attempt[ attempt$ID =='LENG9|Q96B70', ]$ensembl_gene_id <- 'ENSG00000275183'
attempt[ attempt$ID =='LGALS7|P47929', ]$ensembl_gene_id <- 'ENSG00000205076'
attempt[ attempt$ID =='LGALS7B|P47929', ]$ensembl_gene_id <- 'ENSG00000178934'
attempt[ attempt$ID =='MRM1|Q6IN84', ]$ensembl_gene_id <- 'ENSG00000278619'
#attempt[ attempt$ID =='MRNIP|Q6NTE8', ]$ensembl_gene_id <- 'ENSG00000161010'
#attempt[ attempt$ID =='NSMCE3|Q96MG7', ]$ensembl_gene_id <- 'ENSG00000185115'
attempt[ attempt$ID =='PIP4K2B|P78356', ]$ensembl_gene_id <- 'ENSG00000276293'
attempt[ attempt$ID =='PLSCR3|Q9NRY6', ]$ensembl_gene_id <- 'ENSG00000187838'
attempt[ attempt$ID =='PSMB3|P49720', ]$ensembl_gene_id <- 'ENSG00000277791'
#attempt[ attempt$ID =='PYCR3|Q53H96', ]$ensembl_gene_id <- 'ENSG00000104524'
#attempt[ attempt$ID =='STING1|Q86WV6', ]$ensembl_gene_id <- 'ENSG00000184584'
attempt[ attempt$ID =='SYNRG|Q9UMZ2', ]$ensembl_gene_id <- 'ENSG00000275066'
attempt[ attempt$ID =='TAF9|Q16594', ]$ensembl_gene_id <- 'ENSG00000273841'
#attempt[ attempt$ID =='UTP4|Q969X6', ]$ensembl_gene_id <- 'ENSG00000141076'
attempt <- attempt[!duplicated(attempt),]
# table(table(attempt$external_gene_name))
row.names( attempt  )  <- attempt$external_gene_name
row.names( name_trans_slim  )  <- name_trans_slim$GeneName  
dcfdx_Association$ENSG <- NA
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$GeneID %in% row.names(name_trans_slim) ){
    dcfdx_Association[ i, ]$ENSG <- name_trans_slim[dcfdx_Association[ i, ]$GeneID,]$ENSG
  }else{}
  if( dcfdx_Association[ i, ]$GeneID %in% row.names(attempt) ){
    dcfdx_Association[ i, ]$ENSG <- attempt[dcfdx_Association[ i, ]$GeneID,]$ensembl_gene_id
  }else{}
  
}
another <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name","uniprotswissprot"),filter="uniprotswissprot",
                          values=dcfdx_Association[ is.na(dcfdx_Association$ENSG), ]$ProtID,
                          mart=uniProt
)
#another[ another$uniprotswissprot =='P62805',]$ensembl_gene_id <- 'ENSG00000197238'
#another[ another$uniprotswissprot =='P62805',]$external_gene_name <- 'H4C11'
another[ another$uniprotswissprot =='P26640',]$ensembl_gene_id <- 'ENSG00000204394'
another[ another$uniprotswissprot =='P0DMR1',]$ensembl_gene_id <- 'ENSG00000179412'
another[ another$uniprotswissprot =='A4D1U4',]$ensembl_gene_id <- 'ENSG00000257093'
another[ another$uniprotswissprot =='O94910',]$ensembl_gene_id <- 'ENSG00000072071'
another[ another$uniprotswissprot =='P01742',]$ensembl_gene_id <- 'ENSG00000211973'
another[ another$uniprotswissprot =='P01762',]$ensembl_gene_id <- 'ENSG00000211941'
another[ another$uniprotswissprot =='P84243',]$ensembl_gene_id <- 'ENSG00000163041'
another[ another$uniprotswissprot =='P84243',]$external_gene_name <- 'H3-3A'
another[ another$uniprotswissprot =='Q53H96',]$ensembl_gene_id <- 'ENSG00000104524'
another[ another$uniprotswissprot =='Q6FI13',]$ensembl_gene_id <- 'ENSG00000203812'
another[ another$uniprotswissprot =='Q6FI13',]$external_gene_name <- 'H2AC18'
another[ another$uniprotswissprot =='Q6NTE8',]$ensembl_gene_id <- 'ENSG00000161010'
another[ another$uniprotswissprot =='Q86WV6',]$ensembl_gene_id <- 'ENSG00000184584'
another[ another$uniprotswissprot =='Q8N257',]$ensembl_gene_id <- 'ENSG00000196890'
another[ another$uniprotswissprot =='Q8ND23',]$ensembl_gene_id <- 'ENSG00000186648'
another[ another$uniprotswissprot =='Q8NI60',]$ensembl_gene_id <- 'ENSG00000163050'
another[ another$uniprotswissprot =='Q8NI60',]$external_gene_name <-'COQ8A'
another[ another$uniprotswissprot =='Q969J3',]$ensembl_gene_id <- 'ENSG00000165714'
another[ another$uniprotswissprot =='Q969X6',]$ensembl_gene_id <- 'ENSG00000141076'
another[ another$uniprotswissprot =='Q96MC5',]$ensembl_gene_id <- 'ENSG00000166780'
another[ another$uniprotswissprot =='Q96MG7',]$ensembl_gene_id <- 'ENSG00000185115'
another[ another$uniprotswissprot =='Q9BXV9',]$ensembl_gene_id <- 'ENSG00000170270'
#nother[ another$uniprotswissprot =='Q9NRY6',]$ensembl_gene_id <- 'ENSG00000187838'
another[ another$uniprotswissprot =='Q9UII2',]$ensembl_gene_id <- 'ENSG00000130770'
another[ another$uniprotswissprot =='Q9Y6H3',]$ensembl_gene_id <- 'ENSG00000166896'
another <- another[ !duplicated(another),]
row.names(another) <- another$uniprotswissprot
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(another) ){
    dcfdx_Association[ i, ]$ENSG <- another[dcfdx_Association[ i, ]$ProtID,]$ensembl_gene_id
    dcfdx_Association[ i, ]$GeneID <- another[dcfdx_Association[ i, ]$ProtID,]$external_gene_name
  }else{}
}
biomaRt::biomartCacheClear()
another <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", 'uniprot_isoform'),filter="uniprot_isoform",
                          values=dcfdx_Association[ is.na(dcfdx_Association$ENSG), ]$ProtID,
                          mart=uniProt
)
another[ another$uniprot_isoform =='P49589-3',]$ensembl_gene_id <- 'ENSG00000110619'
another <- another[!duplicated(another),]
row.names(another) <- another$uniprot_isoform
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(another) ){
    dcfdx_Association[ i, ]$ENSG <- another[dcfdx_Association[ i, ]$ProtID,]$ensembl_gene_id
    dcfdx_Association[ i, ]$GeneID <- another[dcfdx_Association[ i, ]$ProtID,]$external_gene_name
  }else{}
}
#another <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", "uniprot_gn_id"),filter="uniprot_gn_id",
#                          values=dcfdx_Association[ is.na(dcfdx_Association$ENSG), ]$ProtID,
#                          mart=uniProt
#                        )
another <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", "uniprot_gn_id"),filter="uniprotsptrembl",
                          values=dcfdx_Association[ is.na(dcfdx_Association$ENSG), ]$ProtID,
                          mart=uniProt
)
another <- another[ !(another$ensembl_gene_id %in% 'ENSG00000275125'),]
another <- another[ !(another$ensembl_gene_id %in% 'ENSG00000284874'),]
##_-##
another <- another[complete.cases(another),]
another <- another[ !(another$uniprot_gn_id == ""),]
row.names(another) <- another$uniprot_gn_id
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(another) ){
    dcfdx_Association[ i, ]$ENSG <- another[dcfdx_Association[ i, ]$ProtID,]$ensembl_gene_id
    dcfdx_Association[ i, ]$GeneID <- another[dcfdx_Association[ i, ]$ProtID,]$external_gene_name
  }else{}
}
another <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", "uniprotsptrembl"),filter="uniprotsptrembl",
                          values=dcfdx_Association[ is.na(dcfdx_Association$ENSG), ]$ProtID,
                          mart=uniProt
)
another <- another[ !(another$ensembl_gene_id %in% 'ENSG00000278823'),]
row.names(another) <- another$uniprotsptrembl
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(another) ){
    dcfdx_Association[ i, ]$ENSG <- another[dcfdx_Association[ i, ]$ProtID,]$ensembl_gene_id
    dcfdx_Association[ i, ]$GeneID <- another[dcfdx_Association[ i, ]$ProtID,]$external_gene_name
  }else{}
}
Rep1 <- as.data.frame( rbind(
  c("ENSG00000178605", "H0Y2S1", "GTPBP6"),
  c("ENSG00000211949", "P01779", "IGHV3-23"),
  c("ENSG00000239975", "P01613", "IGKV1D-33"),
  c("ENSG00000211598", "P01625", "IGKV4-1"),
  c("ENSG00000211956", "A0A0A0MS12", "IGHV4-34"),
  c("ENSG00000270550", "A0A0B4J2B7", "IGHV3-30"),
  c("ENSG00000225698", "A0A087WW89", "IGHV3-72"),
  c("ENSG00000211658", "A0A075B6J3", "IGLV3-27"),
  c("ENSG00000099625", "K7EJP2", "CBARP"),
  c("ENSG00000211679", "A0A075B6L0", "IGLC3"),
  c("ENSG00000270550", "P01769", "IGHV3-30"),
  c("ENSG00000183291", "A0A0B4J1S4", "SELENOF"),
  c("ENSG00000211644", "A0A075B6I5", "IGLV1-51"),
  c("ENSG00000174917", "A0A0B4J2A5", "MICOS13"),
  c("ENSG00000211598", "P06313", "IGKV4-1"),
  c("ENSG00000243238", "A0A075B6S3", "IGKV2-30"),
  c("ENSG00000243466", "P01598", "IGKV1-5"),
  c("ENSG00000239951", "P01623", "IGKV3-20"),
  c("ENSG00000239264", "Q86UY0", "TXNDC5"),
  c("ENSG00000182484", "Q9NQA3", "WASH6P"),
  c("ENSG00000261716", "Q6DN03", "H2BC20P"),
  c("ENSG00000123009", "O60361", "NME2P1"),
  c("ENSG00000242802", "A4D1Z4", "AP5Z1"),
  c("ENSG00000212643", "Q15695", "ZRSR2P1"),
  c("ENSG00000204434", "Q9BYX7", "POTEKP"),
  c("ENSG00000182487", "A6NI72", "NCF1B")
), stringsAsFactors = F)  
row.names(Rep1) <- Rep1$V2
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(Rep1) ){
    dcfdx_Association[ i, ]$ENSG <- Rep1[dcfdx_Association[ i, ]$ProtID,]$V1
    dcfdx_Association[ i, ]$GeneID <- Rep1[dcfdx_Association[ i, ]$ProtID,]$V3
  }else{}
}

Rep2 <- as.data.frame( rbind(
  c("ENSG00000070010", "Q92890-1", "UFD1"),
  c("ENSG00000221823", "Q86YQ0", "PPP3R1"),
  c("ENSG00000146955", "C9JJQ5", "RAB19"),
  c("ENSG00000198843", "A0A087WTN3", "SELENOT"),
  c("ENSG00000146556", "A0A096LP75", "WASHP2"),
  c("ENSG00000117600", "A0A0A6YYH7", "PLPPR4"),
  c("ENSG00000198356", "A0A087WXS7", "GET3"),
  c("ENSG00000078618", "B1AKJ5", "NRDC"),
  c("ENSG00000130396", "P55196-2", "AFDN"),
  c("ENSG00000167302", "Q96N21-2", "TEPSIN"),
  c("ENSG00000263563", "J3QRK5", "UBBP4"),
  c("ENSG00000058453", "H7C117", "CROCC"),
  c("ENSG00000183199", "Q58FF7", "HSP90AB3P"),
  c("ENSG00000205100", "Q58FG1", "HSP90AA4P"),
  c("ENSG00000157654", "Q9Y2D5-4", "PALM2AKAP2"),
  c("ENSG00000099290", "Q5SNT6", "WASHC2A"),
  c("ENSG00000122545", "Q16181-2", "SEPTIN7"),
  c("ENSG00000157654", "Q8IXS6-2", "PALM2"),
  c("ENSG00000182774", "P0CW22", "RPS17"),
  c("ENSG00000239975", "P01608", "IGKV1D-33"),
  c("ENSG00000127314", "A6NIZ1", "RAP1B")), 
  stringsAsFactors = F) 
row.names(Rep2) <- Rep2$V2
for( i in 1:dim(dcfdx_Association)[1] ){
  if( dcfdx_Association[ i, ]$ProtID %in% row.names(Rep2) ){
    dcfdx_Association[ i, ]$ENSG <- Rep2[dcfdx_Association[ i, ]$ProtID,]$V1
    dcfdx_Association[ i, ]$GeneID <- Rep2[dcfdx_Association[ i, ]$ProtID,]$V3
  }else{}
}
Peptide_Translation <- as.data.frame( cbind( OldPeptideID=dcfdx_Association$Peptide, NewPeptideID= paste0(dcfdx_Association$GeneID, '|', dcfdx_Association$ProtID)),
                                      stringsAsFactors = F )
row.names( dcfdx_Association ) <- dcfdx_Association$Peptide
row.names(BRAKK_Association ) <- BRAKK_Association$Peptide
CERAD_Association <- CERAD_Association[ !(is.na(CERAD_Association$ENSG)) ,]
row.names(CERAD_Association ) <- CERAD_Association$Peptide
row.names(COGDX_Association ) <- COGDX_Association$Peptide
row.names(Diagnosis_Associations ) <- Diagnosis_Associations$Peptide
dcfdx_Association$Peptide <- paste0(dcfdx_Association$GeneID, '|', dcfdx_Association$ProtID )
BRAKK_Association[ row.names(dcfdx_Association),]$Peptide <- dcfdx_Association$Peptide
BRAKK_Association[ row.names(dcfdx_Association),]$GeneID <- dcfdx_Association$GeneID
BRAKK_Association[ row.names(dcfdx_Association),]$ProtID <- dcfdx_Association$ProtID
BRAKK_Association$ENSG <- NA
BRAKK_Association[ row.names(dcfdx_Association),]$ENSG <- dcfdx_Association$ENSG
CERAD_Association[ row.names(dcfdx_Association),]$Peptide <- dcfdx_Association$Peptide
CERAD_Association[ row.names(dcfdx_Association),]$GeneID <- dcfdx_Association$GeneID
CERAD_Association[ row.names(dcfdx_Association),]$ProtID <- dcfdx_Association$ProtID
CERAD_Association$ENSG <- NA
CERAD_Association[ row.names(dcfdx_Association),]$ENSG <- dcfdx_Association$ENSG
COGDX_Association[ row.names(dcfdx_Association),]$Peptide <- dcfdx_Association$Peptide
COGDX_Association[ row.names(dcfdx_Association),]$GeneID <- dcfdx_Association$GeneID
COGDX_Association[ row.names(dcfdx_Association),]$ProtID <- dcfdx_Association$ProtID
COGDX_Association$ENSG <- NA
COGDX_Association[ row.names(dcfdx_Association),]$ENSG <- dcfdx_Association$ENSG
Diagnosis_Associations[ row.names(dcfdx_Association),]$Peptide <- dcfdx_Association$Peptide
Diagnosis_Associations[ row.names(dcfdx_Association),]$GeneID <- dcfdx_Association$GeneID
Diagnosis_Associations[ row.names(dcfdx_Association),]$ProtID <- dcfdx_Association$ProtID
Diagnosis_Associations$ENSG <- NA
Diagnosis_Associations[ row.names(dcfdx_Association),]$ENSG <- dcfdx_Association$ENSG
row.names( dcfdx_Association ) <- dcfdx_Association$Peptide
row.names(BRAKK_Association ) <- BRAKK_Association$Peptide
row.names(CERAD_Association ) <- CERAD_Association$Peptide
row.names(COGDX_Association ) <- COGDX_Association$Peptide
row.names(Diagnosis_Associations ) <- Diagnosis_Associations$Peptide
#Fix the Palm2Gene Name:
Diagnosis_Associations[ Diagnosis_Associations$GeneID == 'PALM2', ]$GeneID <- 'PALM2AKAP2'
BRAKK_Association[ BRAKK_Association$GeneID == 'PALM2', ]$GeneID <- 'PALM2AKAP2'
CERAD_Association[ CERAD_Association$GeneID == 'PALM2', ]$GeneID <- 'PALM2AKAP2'
COGDX_Association[ COGDX_Association$GeneID == 'PALM2', ]$GeneID <- 'PALM2AKAP2'
dcfdx_Association[ dcfdx_Association$GeneID == 'PALM2', ]$GeneID <- 'PALM2AKAP2'
Non_Collapsed_Comparisons <- list(
  Dignosis = Diagnosis_Associations[ ,c("Peptide",  "GeneID", "ProtID", "ENSG", "Coefficient", "OR", "CI_L", "CI_H", "PVal", "FDR_PVal") ],
  BRAAK = BRAKK_Association[ ,c("Peptide",  "GeneID", "ProtID", "ENSG", "Coefficient", "OR", "CI_L", "CI_H", "PVal", "FDR_PVal") ],
  CERAD = CERAD_Association[ ,c("Peptide",  "GeneID", "ProtID", "ENSG", "Coefficient", "OR", "CI_L", "CI_H", "PVal", "FDR_PVal") ],
  COGDX = COGDX_Association[ ,c("Peptide",  "GeneID", "ProtID", "ENSG", "Coefficient", "OR", "CI_L", "CI_H", "PVal", "FDR_PVal") ],
  dcfdx = dcfdx_Association[ ,c("Peptide",  "GeneID", "ProtID", "ENSG", "Coefficient", "OR", "CI_L", "CI_H", "PVal", "FDR_PVal") ]
)  
Collapsed_Comparisons <- list()
for( Title in names(Non_Collapsed_Comparisons) ){
  eval(parse( text = paste0( " foo <- Non_Collapsed_Comparisons$", Title )  ))
  Scroll <- names(table(foo$GeneID)[table(foo$GeneID)>1])
  for( Gene in Scroll ){
    KeepPoint <- as.numeric( min( foo[ foo$GeneID %in% Gene, ]$FDR_PVal ) )
    if( var( foo[ foo$GeneID %in% Gene, ]$FDR_PVal ) >0 ){
      foo[ foo$GeneID %in% Gene, ][ foo[ foo$GeneID %in% Gene, ]$FDR_PVal > KeepPoint, ]$ENSG <- "TOSS"
    }
  }
  deDUP <- foo[ !(foo$ENSG %in% 'TOSS'), ]
  #Were there ties?
  if( length(names(table(deDUP$GeneID)[table(deDUP$GeneID)>1])) > 0 ){
    Scroll2 <- names(table(deDUP$GeneID)[table(deDUP$GeneID)>1])
    for( Gene in Scroll2 ){
      KeepPoint <- as.numeric( min( deDUP[ deDUP$GeneID %in% Gene, ]$OR ) )
      deDUP[ deDUP$GeneID %in% Gene, ][ deDUP[ deDUP$GeneID %in% Gene, ]$OR > KeepPoint , ]$ENSG <- "TOSS"
    }
    deDUP <- deDUP[ !(deDUP$ENSG %in% 'TOSS'), ]
  }
  eval(parse( text = paste0( "Collapsed_Comparisons$", Title, " <- deDUP" )  ))
}
```

```{r TMT_Coexpression}
install.packages('vbsr', repos='http://cran.us.r-project.org')
#Log2_Normalized
#att <- pvbsrBootstrap( y=y, x=X, nsamp=10, cores=1 )
Cor_Input <- data.table::transpose(Log2_Normalized)
row.names(Cor_Input) <- colnames(Log2_Normalized)
colnames(Cor_Input) <- row.names(Log2_Normalized)
Cor_Input <- as.matrix( Cor_Input )
Cor_Object <- Hmisc::rcorr( Cor_Input, type=c("spearman") )
Cors <- Cor_Object$r
Cor_PVal <- Cor_Object$P
#38% of intereactions are si co-exp with out corrections for multiple comps
# 7.23% of intereactions are sig co-exp with Bonferroni
# 27.86 of intereactions are sig co-exp with fdr
Cor_PVal_fdr <- as.data.frame(Cor_Object$P )
Cor_PVal_fdr <- apply( Cor_PVal_fdr, 1,  p.adjust, method = 'fdr', n=8817 )
#### Clean CoExpression to only FDR significant interactions;
cor( as.vector(as.numeric(Log2_Normalized[1,])), as.vector(as.numeric(Log2_Normalized[3,])), method = "spearman")
#mark<-Sys.time()
#FOO <- psych::corr.test(x, method = "spearman")
#Sys.time()-mark

#Prep for partial correlation detection
Final <- bcv::impute.svd( Cor_Input, k = 100 )
Fin <- Final$x
colnames( Fin ) <- colnames( Cor_Input )
row.names( Fin ) <- row.names( Cor_Input )
#Run Partial correlations for each gene
RUNNe <- function( i=i, p=p){
  source('TMT_Proteomics/code/01_Analysis_Functions.R')
  OBS <- i
  y <- as.matrix(p[,OBS]) 
  colnames(y) <- i
  X <- p[,(colnames(p) %in% OBS) == F ]
  att <- pvbsrBootstrap( y=y, x=X, nsamp=10, cores=1 )
  eval( parse(text = paste0( 'names( att )[ names( att ) == \'intercept\' ] <- \'', OBS, '\'') ))
  return(att[colnames(p)])
}
core <- parallel::detectCores()-2 
cl <- parallel::makePSOCKcluster(core)
doParallel::registerDoParallel(cl)
LIST <- colnames( Fin )
mark <- Sys.time()
FOO <- t(  parallel::parApply(cl, as.matrix( LIST[1:length(LIST)] ), 1, RUNNe, p=Fin) )
message( paste0( Sys.time()-mark )) 
write.table( FOO, file='TMT_Proteomics/vsbr_PartialCor_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
row.names( FOO ) <- colnames( FOO )
FOO_SIF <- reshape2::melt( FOO )
FOO_SIF$Var1 <- as.character(FOO_SIF$Var1)
FOO_SIF$Var2 <- as.character(FOO_SIF$Var2)
FOO_SIF$value <- as.numeric( as.character(FOO_SIF$value) )
write.table( FOO_SIF,
             file='TMT_Proteomics/vsbr_PartialCor_SIF.tsv',
             row.names = F, col.names = T,
             quote=F, sep='\t' 
)
```


```{R Clean_RedundantENSGs }
#FOO <- read.table( synapser::synGet( 'syn24201483' )$path, sep = '\t', header = T) 
## FOO <- read.table( file = '~/Downloads/vsbr_PartialCor_Matrix.tsv', sep = '\t', header = T) 
#FOO <- data.table::fread( file = '~/Downloads/vsbr_PartialCor_Matrix.tsv', sep = '\t', header = T) 
#colnames(FOO) <- gsub( '[.]', '|', colnames(FOO) )
#FOO <- as.data.frame( FOO )
#row.names( FOO ) <- colnames( FOO )
#row.names(FOO) <- colnames(FOO) 
SINK_FOO <- FOO
# - Expand the Peptide 
Peptide_Translation$Old_Gene <- do.call( rbind, strsplit( Peptide_Translation$OldPeptideID, '[|]' ) )[ ,1 ]
Peptide_Translation$Old_Pep <- do.call( rbind, strsplit( Peptide_Translation$OldPeptideID, '[|]' ) )[ ,2 ]
Peptide_Translation$New_Gene <- do.call( rbind, strsplit( Peptide_Translation$NewPeptideID, '[|]' ) )[ ,1 ]
Peptide_Translation$New_Pep <- do.call( rbind, strsplit( Peptide_Translation$NewPeptideID, '[|]' ) )[ ,2 ]
row.names(Peptide_Translation) <- Peptide_Translation$OldPeptideID 
# - Rename co-expression matrix
colnames(FOO) <- Peptide_Translation[ colnames(FOO), ]$NewPeptideID
row.names(FOO) <- Peptide_Translation[ row.names(FOO), ]$NewPeptideID
# - Cleans
head( row.names( FOO)[ !( row.names(FOO) %in% Collapsed_Comparisons$Dignosis$Peptide ) ] )
RepdGenes <- names( table(Peptide_Translation$New_Gene )[ table(Peptide_Translation$New_Gene )>1 ] )
#Nam <- 'AAK1'
#Average partial Corelations across Gene-Annotation types***
Remove <- NULL
for( Nam in RepdGenes ){
  RWS <- Peptide_Translation[ Peptide_Translation$New_Gene == Nam, ]$NewPeptideID
  Remove <- c( Remove, RWS[ !(RWS %in% Collapsed_Comparisons$Dignosis$Peptide ) ] )
  temp <- apply( FOO[ RWS, ], 2, mean)
  for( i in RWS ){
    FOO[ i, ] <- temp
  }
}
Filt_Foo_Path <- FOO[ !( row.names(FOO) %in% Remove), !( colnames(FOO) %in% Remove) ]
FOO <- SINK_FOO
#Convert the square Matrix to SIF format
####First zero out everything below the diagonal
#trial<- Filt_Foo_Path
#for( i in 1:(dim(trial)[2]-1) ){
#  trial[ (i+1):dim(trial)[1] ,i] <- NA
#}
Filt_Foo_SIF_Path <- reshape2::melt( as.matrix(Filt_Foo_Path) )
```


```{r FormatTheSpearman }
Path_fdr <- Cor_PVal_fdr
Path_cor <- Cors
for( i in 1:dim(Path_cor)[1] ){
  Path_cor[ i, which( !(Cor_PVal_fdr[i,] < 0.05) ) ] <- NA
}
for( i in 1:dim(Path_cor)[1] ){
  Path_cor[ i, i ] <- NA
}
colnames(Path_cor) <- Peptide_Translation[ colnames(Path_cor), ]$NewPeptideID
row.names(Path_cor) <- Peptide_Translation[ row.names(Path_cor), ]$NewPeptideID
Path_cor_SIF <- reshape2::melt( Path_cor )
Path_cor_SIF <- Path_cor_SIF[ !(is.na(Path_cor_SIF$value)) ,]
Path_cor_SIF$absvalue <- abs( Path_cor_SIF$value)
Path_cor_SIF$GeneA <- do.call(rbind,strsplit(as.character(Path_cor_SIF$Var1), '[|]'))[,1]
Path_cor_SIF$GeneB <- do.call(rbind,strsplit(as.character(Path_cor_SIF$Var2), '[|]'))[,1] 
Path_cor_SIF <- Path_cor_SIF[, c("Var1", "Var2",
                                 "GeneA", "GeneB",
                                 "GeneB", "absvalue" ) ]
```

```{r PushToSynapse}
#Push to Synapse
parentId ="syn21917175"
folderName = 'TMT_Proteomics'
CODE <- synapser::Folder(name = folderName, parentId = parentId)
CODE <- synapser::synStore(CODE)
activityName = 'TMT VSBR Partial Correlation Stats'
activityDescription = 'TMT Proteomics Partial Correlation Statistics of Variable Bayes Spike Regression'
thisFileName <- 'code/02_DE_Analysis.Rmd'
# Github link
thisRepo <- githubr::getRepo(repository = "jgockley62/TMT_Proteomics", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0(thisFileName))
#Set Used SynIDs For Provenance
#Need to push covariates to synapse to automate above metadata parts and then add to provenance
used <- c( 'syn21266454', 'syn21266453', 'syn21323404', 'syn21323366', 'syn23583548', 'syn3191087', 'syn23573928' )
# Set annotations
all.annotations = list(
  dataType = 'TMT',
  dataSubType = 'geneExp',
  summaryLevel = 'Peptide',
  assay	 = 'TMT',
  tissueTypeAbrv	= 'DLPFC', 
  study = 'ROSMAP', 
  organism = 'HomoSapiens',
  consortium	= 'AMP-AD',
  normalizationStatus	= TRUE,
  normalizationType	= 'TAMPOR'
)
################################################################################################
## Push Non-Collapsed vbsr files to synapse
ENRICH_OBJ <- synapser::File( path='TMT_Proteomics/vsbr_PartialCor_Matrix.tsv', name = 'vsbr PartialCor Statistics Matrix', parentId=CODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)

ENRICH_OBJ <- synapser::File( path='TMT_Proteomics/vsbr_PartialCor_SIF.tsv', name = 'vsbr PartialCor Statistics Table SIF Format', parentId=CODE$properties$id )
all.annotations$dataSubType = 'corelation table'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
################################################################################################
## Push Peptide Translator to Synapse
write.csv( Peptide_Translation, file="Peptide_Translation.csv", row.names=F, quote=F)
ENRICH_OBJ <- synapser::File( path='Peptide_Translation.csv', name = 'Translation File for Peptides Poset Ensemble Pull', parentId=CODE$properties$id )
all.annotations$dataSubType = 'ID Translation'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
################################################################################################
## Push Non-Collapsed Comparisons to Synapase
subFolder <- synapser::synFindEntityId( name=folderName , parent=parentId )
subCODE <- synapser::Folder(name = "NonCollapsed Comparison Stats", parentId = subFolder)
subCODE <- synapser::synStore(subCODE)
#Push Association Comparisons - Non Collapsed
for( Data in names( Non_Collapsed_Comparisons ) ){
  PATH <- paste0( 'Non_Collapsed_', Data, '_Comps.csv' )
  NAME <- paste0( 'Non Collapsed ', Data, ' Comparisons' )
  OBJ <- paste0( 'Non_Collapsed_Comparisons$', Data)
  
  eval(parse( text = paste0( 'write.csv( ', OBJ, ', file= \'', PATH,'\', row.names=F, quote=F)' ) ))
  ENRICH_OBJ <- synapser::File( path=PATH, name = NAME, parentId=subCODE$properties$id )
  all.annotations$dataSubType = 'AssociationStats'
  ENRICH_OBJ$annotations = all.annotations
  synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
}
################################################################################################
## Push Collapsed Comparisons to Synapase
subFolder <- synapser::synFindEntityId( name=folderName , parent=parentId )
subCODE <- synapser::Folder(name = "Collapsed Comparison Stats", parentId = subFolder)
subCODE <- synapser::synStore(subCODE)
#Push Association Comparisons - Non Collapsed
for( Data in names( Collapsed_Comparisons ) ){
  PATH <- paste0( 'Collapsed_', Data, '_Comps.csv' )
  NAME <- paste0( 'Collapsed ', Data, ' Comparisons' )
  OBJ <- paste0( 'Collapsed_Comparisons$', Data)
  
  eval(parse( text = paste0( 'write.csv( ', OBJ, ', file= \'', PATH,'\', row.names=F, quote=F)' ) ))
  ENRICH_OBJ <- synapser::File( path=PATH, name = NAME, parentId=subCODE$properties$id )
  all.annotations$dataSubType = 'AssociationStats'
  ENRICH_OBJ$annotations = all.annotations
  synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
}
################################################################################################
## Push Pathway App Filtered VBSR comps to synapse
subFolder <- synapser::synFindEntityId( name=folderName , parent=parentId )
subCODE <- synapser::Folder(name = "Pathway Tracing CoExp Data", parentId = subFolder)
subCODE <- synapser::synStore(subCODE)
write.table( Filt_Foo_Path, file='PathwayTrace_Input_vsbr_PartialCor_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
write.table( Filt_Foo_SIF_Path,
             file='PathwayTrace_Input_vsbr_PartialCor_SIF.tsv',
             row.names = F, col.names = T,
             quote=F, sep='\t' 
)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_vsbr_PartialCor_Matrix.tsv', name = 'vsbr PartialCor Statistics Matrix', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_vsbr_PartialCor_SIF.tsv', name = 'vsbr PartialCor Statistics Table SIF Format', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation table'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
################################################################################################
## Push Pathway App Filtered Spearman comps to synapse
subFolder <- synapser::synFindEntityId( name=folderName , parent=parentId )
subCODE <- synapser::Folder(name = "Pathway Tracing CoExp Data", parentId = subFolder)
subCODE <- synapser::synStore(subCODE)
write.table( Cors, file='PathwayTrace_Input_Spearman_Cor_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
write.table( Cor_PVal, file='PathwayTrace_Input_Spearman_pval_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
write.table( Cor_PVal_fdr, file='PathwayTrace_Input_Spearman_pvalFDR_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
write.table( Path_cor, file='PathwayTrace_Input_Spearman_FiltCor_Matrix.tsv', row.names = F, col.names = T, quote=F, sep='\t' )
write.table( Path_cor_SIF,
             file='PathwayTrace_Input_Spearman_CorFilt_SIF.tsv',
             row.names = F, col.names = T,
             quote=F, sep='\t' 
)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_Spearman_Cor_Matrix.tsv', name = 'Spearman Correlation Matrix Matrix', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_Spearman_pval_Matrix.tsv', name = 'Spearman P-Values Matrix Matrix', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_Spearman_pvalFDR_Matrix.tsv', name = 'Spearman FDR Corected P-Value  Matrix', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_Spearman_FiltCor_Matrix.tsv', name = 'Spearman Filtered Correlation Matrix Matrix', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation matrix'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
ENRICH_OBJ <- synapser::File( path='PathwayTrace_Input_Spearman_CorFilt_SIF.tsv', name = 'Filtered Spearman Correlation Table', parentId=subCODE$properties$id )
all.annotations$dataSubType = 'corelation table'
ENRICH_OBJ$annotations = all.annotations
synapser::synStore( ENRICH_OBJ, used = used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
