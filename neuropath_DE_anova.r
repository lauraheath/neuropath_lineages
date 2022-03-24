#Loading pre-computed lineage from rnasea_lineage_neuropathGenes.r
MonRun <- readRDS(file='~/prot-lineage/CERAD_MonocleObject.R')
#get gene list from same script
GeneNamesCERAD <- read.csv(file='~/prot-lineage/rnaseq_genes_cerad.csv')


#getting gene symbol
gene_short_name <- GeneNamesCERAD$gene_short_name

#check number of states
table(MonRun$State2)

#Create differential expression based on state 
MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)


#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))
l2$p_7 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  l2$p_7[i] <- tk$s[6,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  l2$d_7[i] <- tk$s[6,1]
  
}

#Save data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

write.csv(df3,file='~/prot-lineage/CERAD_DE_by_state_genes.csv',quote=F,row.names=F)
#store differentially expressed genes in synapse space
file <- synapser::File(path='~/prot-lineage/CERAD_DE_by_state_genes.csv', parentId='syn27254960')
file <- synapser::synStore(file)



######## repeat for BRAAK lineage states ######

MonRun <- readRDS(file='~/prot-lineage/BRAAK_MonocleObject.RDS')
#get gene list from same script
GeneNamesCERAD <- read.csv(file='~/prot-lineage/rnaseq_genes_braak.csv')


#getting gene symbol
gene_short_name <- GeneNamesCERAD$gene_short_name

#check number of states
table(MonRun$State2)

#Create differential expression based on state 
MonRun2 = MonRun


#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))
l2$p_7 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  l2$p_7[i] <- tk$s[6,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  l2$d_7[i] <- tk$s[6,1]
  
}

#Save data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

write.csv(df3,file='~/prot-lineage/BRAAK_DE_by_state_genes.csv',quote=F,row.names=F)
#store differentially expressed genes in synapse space
file <- synapser::File(path='~/prot-lineage/BRAAK_DE_by_state_genes.csv', parentId='syn27254960')
file <- synapser::synStore(file)
