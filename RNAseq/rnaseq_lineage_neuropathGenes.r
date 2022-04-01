#reinstall monocle if uninstalled for DE analysis
BiocManager::install("monocle")
install.packages("VennDiagram")
library(monocle)
library(VennDiagram)

#generate trajectories based on genes DE by braak and cerad (two trajectories)

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
  #reverse for CERAD lineage
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

#load neuropath_results.tsv (DE genes by braak and cerad from 'ThanneerCode_DE_RNAseq_LH.R')
#de1 <- read.delim(file="neuropath_results.tsv")
neuroDEobj <- synapser::synGet('syn27255271')
de1 <- read.csv(neuroDEobj$path)
#subset for all braak comparisons (LOW vs MED, LOW vs HIGH, MED vs HIGH)
de3 <- dplyr::filter(de1,Model=='braak_category')

#subset for all ceradsc comparisons (1-2, 1-3, 1-4, 2-3, 2-4, 3-4)
#de3 <- dplyr::filter(de1,Model=='ceradsc')


AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]

#Normalize all columns 
GeneNames <- Dat$ensembl_gene_id
GeneNamesAD <- unique(AMP_mods$GeneID)

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

rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)
rosmap <- dplyr::select(rosmap,braaksc,ceradsc,projid)
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)

rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
rosmapRNAid <- rosmapRNAid[!duplicated(rosmapRNAid),]

Dat4<-dplyr::left_join(Dat4,rosmapRNAid,by=c('SampleID'='rnaseq_id'))

#recode braaksc for three levels
Dat4$braak_category = 'MED'
Dat4$braak_category[Dat4$braaksc <= 2] = 'LOW'
Dat4$braak_category[Dat4$braaksc  >= 5] = 'HIGH'

Dat4$braaksc <- factor(Dat4$braaksc,levels = c(0:6))
Dat4$ceradsc <- factor(Dat4$ceradsc,levels = c(1:4))
Dat4$cogdxNew <- rep(NA,nrow(Dat4))
Dat4$cogdxNew[Dat4$cogdx==1] <- 'NCI'
Dat4$cogdxNew[Dat4$cogdx==2] <- 'MCI'
Dat4$cogdxNew[Dat4$cogdx==4] <- 'LOAD'
Dat4$cogdxNew <- factor(Dat4$cogdxNew,levels = c('NCI','MCI','LOAD'))

#source('LineageFunctions.R')
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
write.csv(rnaseq_genes, file="~/prot-lineage/rnaseq_genes_braak.csv", row.names=FALSE)

rnaseq_genes2 <- as.data.frame(gene_short_name)
write.csv(rnaseq_genes2, file="~/prot-lineage/rnaseq_genes_cerad.csv", row.names=FALSE)

#Run Monocle2: (ignore warning messages that occur)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#check: is tree ordered in the right direction?
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
table(MonRun$State)
#collapse some states due to small numbers and rearrange

#BRAAK:female samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 11] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 10] <- 3
MonRun$State2[MonRun$State == 8] <- 3
MonRun$State2[MonRun$State == 9] <- 4
MonRun$State2[MonRun$State == 7] <- 5
MonRun$State2[MonRun$State == 4] <- 6
MonRun$State2[MonRun$State == 5] <- 6
MonRun$State2[MonRun$State == 6] <- 7
plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1)
MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)


#CERAD:female samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 9] <- 6
MonRun$State2[MonRun$State == 8] <- 7
MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)

###### FIGURES #########
tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()

tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
dev.off()



tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "ceradsc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CERAD Score")
g
dev.off()


tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

tiff(file='~/prot-lineage/figures/BRAAK_FEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_FEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()

MonRun$cogdx <- as.factor(MonRun$cogdx)
tiff(file='~/prot-lineage/figures/BRAAK_ FEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/CERAD_ FEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()


#compare the genes DE by case/control and DE by neuropathology

foobar <- synapser::synGet('syn11695124')
foobar2 <- data.table::fread(foobar$path,data.table=F)

#load ROSMAP_DLPFC_DiffExpression.tsv (AD case-control DE analysis)
#save braak genes:
braakGenes <- GeneNamesAD
#ceradGenes <- GeneNamesAD

de_file <- synapser::synGet('syn8456721')
de1 <- data.table::fread(de_file$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL')
AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]
GeneNamesAD <- unique(AMP_mods$GeneID)

venn.diagram(
  x=list(braakGenes,GeneNamesAD),
  category.names = c("Braak Genes", "AD Genes"),
  filename = 'braak_AD_venn.png',
  output=TRUE
)

venn.diagram(
  x=list(ceradGenes,GeneNamesAD),
  category.names = c("Cerad Genes", "AD Genes"),
  filename = 'cerad_AD_venn.png',
  output=TRUE
)




######## for stats & other figures, create a dataframe with all relevant covariates & pseudotime & state

x <- list()
x$SampleID <- MonRun$SampleID
x$State <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$Diagnosis
x$braaksc <- MonRun$braaksc
x$braak_category <- MonRun$braak_category
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$apoe_genotype
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$batch <- MonRun$Batch

#rename and create a scaled pseudotime variable
covars <- as.data.frame(x)
covars$pseudotime_sc <- scale(covars$Pseudotime, center=F)

#save variables file for later
write.csv(covars, file="~/prot-lineage/BraakLineage_pseudotimes_states.csv", row.names=FALSE)
#store pseudotimes and states in synapse space for neuropath-based lineage
write.csv(covars, file='RNAseq_DLPFCbraak_lineage_pseudotimesF.csv', row.names = FALSE)
file <- synapser::File(path='RNAseq_DLPFCbraak_lineage_pseudotimesF.csv', parentId='syn27254960')
file <- synapser::synStore(file)


#save variables file for later
write.csv(covars, file="~/prot-lineage/CeradLineage_pseudotimes_states.csv", row.names=FALSE)
#store pseudotimes and states in synapse space for neuropath-based lineage
write.csv(covars, file='RNAseq_DLPFCcerad_lineage_pseudotimesF.csv', row.names = FALSE)
file <- synapser::File(path='RNAseq_DLPFCcerad_lineage_pseudotimesF.csv', parentId='syn27254960')
file <- synapser::synStore(file)




#run logistic regression comparing pseudotiem between cases and controls only
casecontrolF <- subset(covars, covars$diagnosis=='AD'|covars$diagnosis=='CONTROL')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)

#summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))
#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(casecontrolF,aes(x=diagnosis,
                             y=pseudotime_sc,
                             color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
#dev.off()


#store Monocle object for DE analysis
saveRDS(MonRun, file='~/prot-lineage/BRAAK_MonocleObject.R')
saveRDS(MonRun, file='~/prot-lineage/CERAD_MonocleObject.R')












