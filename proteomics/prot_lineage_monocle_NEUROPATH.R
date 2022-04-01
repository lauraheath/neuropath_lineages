
#first run "get_proteins_metadata.r" script and save the protein matrix and metadata file in your directory


#reinstall monocle if uninstalled for DE analysis
#BiocManager::install("monocle")
#install.packages("VennDiagram")
library(monocle)
library(VennDiagram)


Log2_Normalized <- readRDS(file="~/neuropath_lineages/proteomics/Dx_analysis/Log2_Normalized.rds")
Meta <- readRDS(file="~/neuropath_lineages/proteomics/Dx_analysis/Meta.rds")

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
  #will need to reverse for braak & cerad female samples
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

Dat <- Log2_Normalized
Dat[is.na(Dat)] <- 0


#get list of proteins that are differentially expressed between neuropath states 
p <- synapser::synGet('syn28712770')
Neuro_genes <- read.csv(p$path)

#subset for all braak comparisons (LOW vs MED, LOW vs HIGH, MED vs HIGH)
#DEgenes <- dplyr::filter(Neuro_genes,category=='braak_category')

#subset for all ceradsc comparisons (1-2, 1-3, 1-4, 2-3, 2-4, 3-4)
DEgenes <- dplyr::filter(Neuro_genes,category=='ceradsc')

DEgenes1 <- subset(DEgenes, DEgenes$FDR_PVal<0.1)
DEgenes1 <- unique(DEgenes1$Peptide)

#subset the protein matrix
genes2<-c()
#for (gene in unique(c(as.vector(DEgenes1$Peptide)))){
for (gene in unique(c(as.vector(DEgenes1)))){
  if (gene %in% rownames(Dat)){
    genes2 <- c(genes2,which(rownames(Dat)==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)



#msex==0 is female, msex==1 is male; run separately for sex-specific analysis
#In_S <- which(Meta$msex == 0)
In_S <- which(Meta$msex == 1)
Dat2 <- Dat2[,In_S]
Meta2 <- Meta[In_S,]

gene_short_name <- rownames(Dat2)

#save the protein matrix with rownames for later
#saveRDS(Dat2, file="~/neuropath_lineages/proteomics/BRAAKFemale_matrix_rownames.rds")
#saveRDS(Dat2, file="~/neuropath_lineages/proteomics/BRAAKMale_matrix_rownames.rds")
#saveRDS(Dat2, file="~/neuropath_lineages/proteomics/CERADFemale_matrix_rownames.rds")
saveRDS(Dat2, file="~/neuropath_lineages/proteomics/CERADMale_matrix_rownames.rds")
temp <- Dat2
temp2 <- Meta2


temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))

rownames(temp)<-NULL
rownames(temp2)<-NULL


#Run Monocle2: (ignore warning messages that occur)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

#check to see if State designations make sense
g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

#reorder with state 3 as root, and rerun Monocle:
#for cerad samples, root at state 3
#MonRun <- orderCells(MonRun, root_state = 3)


#check to see if pseudotime is ordered correctly
plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)

table(MonRun$State)
#collapse some states due to small numbers and rearrange
plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)

#braak female samples
# MonRun$State2 <- MonRun$State
# MonRun$State2[MonRun$State == 2] <- 3
# MonRun$State2[MonRun$State == 5] <- 2
# MonRun$State2[MonRun$State == 7] <- 2
# MonRun$State2[MonRun$State == 6] <- 2
# MonRun$State2[MonRun$State == 4] <- 3
# MonRun$State2[MonRun$State == 3] <- 4
# 
# #braak male samples
# MonRun$State2 <- MonRun$State
# MonRun$State2[MonRun$State == 9] <- 2
# MonRun$State2[MonRun$State == 2] <- 3
# MonRun$State2[MonRun$State == 8] <- 3
# MonRun$State2[MonRun$State == 3] <- 4
# MonRun$State2[MonRun$State == 4] <- 4
# MonRun$State2[MonRun$State == 5] <- 4
# MonRun$State2[MonRun$State == 7] <- 5
# MonRun$State2[MonRun$State == 6] <- 6

#cerad female samples
# MonRun$State2 <- MonRun$State
# MonRun$State2[MonRun$State == 3] <- 1
# MonRun$State2[MonRun$State == 11] <- 2
# MonRun$State2[MonRun$State == 1] <- 2
# MonRun$State2[MonRun$State == 4] <- 3
# MonRun$State2[MonRun$State == 10] <- 3
# MonRun$State2[MonRun$State == 5] <- 4
# MonRun$State2[MonRun$State == 6] <- 4
# MonRun$State2[MonRun$State == 7] <- 5
# MonRun$State2[MonRun$State == 9] <- 5
# MonRun$State2[MonRun$State == 8] <- 6

#cerad male samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 5] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 3] <- 4
MonRun$State2[MonRun$State == 4] <- 5


plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)



MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

#save the monocle object for later:
#saveRDS(MonRun, file="~/neuropath_lineages/proteomics/BRAAKFemale_monocleObject.rds")
#saveRDS(MonRun, file="~/neuropath_lineages/proteomics/BRAAKMale_monocleObject.rds")
#saveRDS(MonRun, file="~/neuropath_lineages/proteomics/CERADFemale_monocleObject.rds")
saveRDS(MonRun, file="~/neuropath_lineages/proteomics/CERADMale_monocleObject.rds")

###### FIGURES #########
#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)

g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()

#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='neuropath_lineages/proteomics/BRAAKMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='neuropath_lineages/proteomics/CERADMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
dev.off()

MonRun$APOE <- factor(MonRun$APOE,levels=c(0,1,2))
#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "APOE",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="APOE e4 Dosage")
g
dev.off()


#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "ceradsc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CERAD Score")
g
dev.off()

#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "cogdx",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Cognitive Diagnosis")
g
dev.off()


#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()


#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()

#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()





######## for stats & other figures, create a dataframe with all relevant covariates & pseudotime & state

x <- list()
x$SampleID <- MonRun$batchChannel
x$individualID <- MonRun$individualID
x$rnaseqID <- MonRun$rnaseq_id
x$State2 <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$APO
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$batch <- MonRun$batch
#x$mmse <- MonRun$cts_mmse30_lv
#x$age_death <- MonRun$age_death
x$rna_seq_sample <- MonRun$rnaseq
x$SampleID <- as.character(x$SampleID)

#rename and create a scaled pseudotime variable
prot_pstime_covars_F <- as.data.frame(x)
prot_pstime_covars_F$pseudotime_sc <- scale(prot_pstime_covars_F$Pseudotime, center=F)

#save variables file for later
# write.csv(prot_pstime_covars_F, file="~/neuropath_lineages/proteomics/Prot_DLPFCbraak_lineage_pseudotimesF.csv", row.names=FALSE)
# file <- synapser::File(path='~/neuropath_lineages/proteomics/Prot_DLPFCbraak_lineage_pseudotimesF.csv', parentId='syn28712340')
# file <- synapser::synStore(file)

# write.csv(prot_pstime_covars_F, file="~/neuropath_lineages/proteomics/Prot_DLPFCbraak_lineage_pseudotimesM.csv", row.names=FALSE)
# file <- synapser::File(path='~/neuropath_lineages/proteomics/Prot_DLPFCbraak_lineage_pseudotimesM.csv', parentId='syn28712340')
# file <- synapser::synStore(file)

# write.csv(prot_pstime_covars_F, file="~/neuropath_lineages/proteomics/Prot_DLPFCcerad_lineage_pseudotimesF.csv", row.names=FALSE)
# file <- synapser::File(path='~/neuropath_lineages/proteomics/Prot_DLPFCcerad_lineage_pseudotimesF.csv', parentId='syn28712340')
# file <- synapser::synStore(file)

write.csv(prot_pstime_covars_F, file="~/neuropath_lineages/proteomics/Prot_DLPFCcerad_lineage_pseudotimesM.csv", row.names=FALSE)
file <- synapser::File(path='~/neuropath_lineages/proteomics/Prot_DLPFCcerad_lineage_pseudotimesM.csv', parentId='syn28712340')
file <- synapser::synStore(file)


#run female and male analyses separately
Fvariables <- prot_pstime_covars_F
#Fvariables <- prot_pstime_covars_M

#run logistic regression comparing pseudotiem between cases and controls only
casecontrolF <- subset(Fvariables, Fvariables$diagnosis=='AD'|Fvariables$diagnosis=='Control')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))
#tiff(file='~/neuropath_lineages/proteomics/BRAAKFEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/BRAAKMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/neuropath_lineages/proteomics/CERADFEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/CERADMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(casecontrolF,aes(x=diagnosis,
                        y=pseudotime_sc,
                        color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()

#run proportional odds logistic regression for neuropath/cognitive endpoints:
braakfit <- MASS::polr(braaksc ~ pseudotime_sc,Fvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,Fvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,Fvariables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')



#look for correlations with GWAS LOAD genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

Dat3 <- Dat2
Dat3$gene_names <- rownames(Dat3)
Dat3$gene_short_name <- gsub("\\|.*", "", Dat3$gene_names)


# dlpfcCPMObj <- synapser::synGet('syn8456638')
# Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat2)#[-120]
#sampleIds <- sampleIds[-120]
sampleIds
geneIds <- Dat3$gene_short_name
Dat3$gene_short_name<-NULL
Dat3$gene_names<-NULL
Dat3 <- t(Dat3)
colnames(Dat3) <- geneIds
Dat3 <- data.frame(Dat3,stringsAsFactors=F)
Dat3$sampleId <- sampleIds
dlpfc <- dplyr::left_join(Fvariables,Dat3,by=c('SampleID'='sampleId'))
dlpfc2 <- dlpfc[,14:1850]

corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor))

mean(corDfdlpfc$cor)
mean(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor)


corDfdlpfc$adGwas <- corDfdlpfc$geneid %in% ad_gwas
colnames(corDfdlpfc)[3] <- 'LOADGWASGene'
corDfdlpfc$LOADGWASGene2 <- ifelse(corDfdlpfc$LOADGWASGene==FALSE, "NOT GWAS GENE", "GWAS GENE")

#tiff(file='~/neuropath_lineages/proteomics/Dx_analysis/FEMALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/neuropath_lineages/proteomics/Dx_analysis/MALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfdlpfc,ggplot2::aes(x=LOADGWASGene2,y=cor,fill=LOADGWASGene2))
g <- g + ggplot2::geom_boxplot() + theme(legend.position="none")
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'LOAD GWAS STATUS',y='Correlation with pseudotime')
g
dev.off()




#how much overlap between AD genes and neuropath genes?
#get list of proteins that are differentially expressed between AD case & control (rerun of jake's analysis March 2022)
p <- synapser::synGet('syn28559312')
ADgenes_prot <- read.csv(p$path)
ADgenes_prot1 <- subset(ADgenes_prot, ADgenes_prot$FDR_PVal<0.1)
dim(ADgenes_prot1)
ADgenes <- unique(ADgenes_prot1$Peptide)

braakGenes <- DEgenes1
#ceradGenes <- DEgenes1




venn.diagram(
  x=list(braakGenes,ADgenes),
  category.names = c("Braak Genes", "AD Genes"),
  filename = 'braak_AD_venn.png',
  output=TRUE
)

venn.diagram(
  x=list(ceradGenes,ADgenes),
  category.names = c("Cerad Genes", "AD Genes"),
  filename = 'cerad_AD_venn.png',
  output=TRUE
)
