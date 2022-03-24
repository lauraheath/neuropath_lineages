### extract the gene list from the rna-bulk analysis and proteomics analysis and see how many genes overlap

#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
#load ROSMAP filtered counts logCPM:
dlpfcCPMObj <- synapser::synGet('syn8456638')

Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
dlpfcCovObj <- synapser::synGet('syn11024258')
Dat2 <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)

foobar <- synapser::synGet('syn11695124')
foobar2 <- data.table::fread(foobar$path,data.table=F)

#load ROSMAP_DLPFC_DiffExpression.tsv
de_file <- synapser::synGet('syn8456721')
de1 <- data.table::fread(de_file$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL')
AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]

#Normalize all columns 
GeneNames <- Dat$ensembl_gene_id
GeneNamesAD <- AMP_mods$GeneID

Names <- colnames(Dat)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat) <- Names
cNames <- Dat2$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat <- Dat[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
Dat2 <- Dat2[In,]

ColNorm <- function(Dat){
  
  M = max(colSums(Dat))
  l <- length(colnames(Dat))
  
  for( i in 1:l){
    
    Dat[,i] = Dat[,i]*(M/sum(Dat[,i]))
    
  }
  
  return(Dat)
}

DatNorm <- ColNorm(Dat)
In_genes <- which(GeneNames %in% GeneNamesAD)
DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]

#removing bad batches
DatNorm2 <- DatNorm2[,Dat2$Batch<7]
Dat2 <- Dat2[Dat2$Batch<7,] 



DatNorm3 <- DatNorm2
Dat3 <- Dat2

#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(Dat3$msex == 0)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]

# In_cov <- which(Cov$ID %in% Dat4$Donor_ID)
# Cov <- Cov[In_cov,]
# In_cov <- c()
# for(i in 1:length(Dat4$Donor_ID)){
#   temp <- which(Cov$ID == Dat4$Donor_ID[i])
#   In_cov <- c(In_cov,temp[1])
# }
# Cov <- Cov[In_cov,]
# 
# for (i in 23:26){
#   Cov[,i] <- (Cov[,i] - min(Cov[,i]))/(max(Cov[,i])-min(Cov[,i]))
# }

source('LineageFunctions.R')
temp <- DatNorm4
temp2 <- Dat4
#temp2$APOE4 <- as.character(temp2$APOE4)
#temp2$braaksc <- as.character(temp2$braaksc)
#temp2$ceradsc <- as.character(temp2$ceradsc)
#temp2$cogdx.1 <- as.character(temp2$cogdx.1)

#converting ENSG to gene symbols
#convert to gene symbol
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='useast.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
gene_short_name <- Make.Gene.Symb(GeneNamesAD)
rownames(temp)<-NULL
rownames(temp2)<-NULL

#save gene list for later use:
rnaseq_genes <- as.data.frame(gene_short_name)
write.csv(rnaseq_genes, file="~/prot-lineage/data/rnaseq_genes.csv", row.names=FALSE)


#Run Monocle2: (ignore warning messages that occur)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

table(MonRun$State)

#MonRun2@assayData$exprs <- ScaledDat
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6

MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
tiff(file='~/prot-lineage/figures/FEMALE_RNAseqDLPFC_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.0)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()






############ Recreate monocle object for mayo (tcx) data #################

RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

### extract the gene list from the rna-bulk analysis and proteomics analysis and see how many genes overlap
#load rosmap filtered counts logCPM:
tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file: syn8466814
tcxCovObj <- synapser::synGet('syn8466814')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)

#subsetting genes based on differential expression
#synapse id of DE file: syn8468023?
de_file <- synapser::synGet('syn8468023')
de_file2 <- synapser::synGet('syn18475579')

de1 <- data.table::fread(de_file$path,data.table=F)
de2 <- data.table::fread(de_file2$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL',Tissue.ref=='TCX')
#AMP_mods <-  read.csv('Data/TCX_DE.csv')
#AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- de2[,-1]
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]




#Normalize all columns 
GeneNames <- Dat$ensembl_gene_id
GeneNamesAD <- AMP_mods$GeneID

Names <- colnames(Dat)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat) <- Names
cNames <- Dat2$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
Dat <- Dat[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
Dat2 <- Dat2[In,]

ColNorm <- function(Dat){
  
  M = max(colSums(Dat))
  l <- length(colnames(Dat))
  
  for( i in 1:l){
    
    Dat[,i] = Dat[,i]*(M/sum(Dat[,i]))
    
  }
  
  return(Dat)
}

DatNorm <- ColNorm(Dat)
DatNorm <- ColNorm(Dat)
In_genes <- which(GeneNames %in% GeneNamesAD)
DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]

#Subsetting based on brain region
In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
DatNorm3 <- DatNorm2[,In_BR]
Dat3 <- Dat2[In_BR,]

#subsetting based on gender
Sex <- 'FEMALE'
In_S <- which(Dat3$Sex == Sex)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]

temp <- DatNorm4
#temp2 <- cbind(Dat4,Cov)
temp2 <- Dat4

temp2$Diagnosis <- temp2$Tissue.SourceDiagnosis

temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.AD'] <- 'AD'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.CONTROL'] <- 'Control'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PATH_AGE'] <- 'PA'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PSP'] <- 'PSP'


#converting ENSG to gene symbols
#convert to gene symbol
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='www.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
gene_short_name <- Make.Gene.Symb(GeneNamesAD)

rownames(temp) <- NULL
colnames(temp) <- NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 2)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

head(MonRun$Pseudotime)


table(MonRun$State)

MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 1
MonRun$State2[MonRun$State == 1] <- 6
MonRun$State2[MonRun$State == 6] <- 4


MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)

g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 2)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g


tiff(file='~/prot-lineage/figures/FEMALE_RNAseqTCX_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.0)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()