synapser::synLogin()

TMT_Express_Load <- function( SynID, names ){
  #'@SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  grepl('You do not have READ permission for the requested entity', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
    return("User does not have READ access")
  }
  if( !grepl('.csv$', synapser::synGet(SynID)$path ) ){
    return("File to import must be a CSV File")
  }
  if( !is.numeric(names) ){
    return("names is not a number")
  }
  
  if( names == 0 ){
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ))
  }else{
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ), row.names = 1)
  }
  return( import )
}

Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
#Normalized <- TMT_Express_Load('syn21266453', 1)
Meta <- TMT_Express_Load('syn21323404', 1)

p2 <- synapser::synGet('syn23583548')
BioSpecimin <- read.csv(p2$path)

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
Meta$diagnosis <- "Other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "Control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"

## Add Ages over 90 for modeling
p <- synapser::synGet('syn23573928')
Mast <- read.csv(p$path)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')

#Convert APOE
APOE <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOE) <- APOE
APOE[names( table(Meta$apoe_genotype) )] <- 0
APOE[grepl("4",names(APOE))] <- 1
APOE[grepl("44",names(APOE))] <- 2
Meta$APOE <- as.numeric( APOE[ as.character(Meta$apoe_genotype) ] )

#impute missing pmi (2 missing: fill in with median pmi)
paste('Imputing PMI to:',median(Meta$pmi[!is.na(Meta$pmi)]))
#add this back into the metadata file
Meta$pmi[is.na(Meta$pmi)] <- 6.5

row.names(Meta) <- Meta$batchChannel

#need to reorder Meta rows to match the protein data column order
meta_ordered <- match(colnames(Log2_Normalized), rownames(Meta))
Meta <- Meta[meta_ordered,]


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
saveRDS(Log2_Normalized, file="~/neuropath_lineages/proteomics/Dx_analysis/Log2_Normalized.rds")
saveRDS(Meta, file="~/neuropath_lineages/proteomics/Dx_analysis/Meta.rds")
