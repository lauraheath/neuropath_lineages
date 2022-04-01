######### branch-specific differential expression analysis

#need Monrun object from protein_lineage_monocle2 and gene_short_name vector
#scripts in order by neuropath & sex (4 states in braak females, 6 in braak males, 6 in cerad females, 5 in cerad males)

MonRun <- readRDS(file="~/neuropath_lineages/proteomics/BRAAKFemale_monocleObject.rds")
temp <- readRDS(file="~/neuropath_lineages/proteomics/BRAAKFemale_matrix_rownames.rds")
table(MonRun$State2)



#there are 4 states in the braak female tree
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- rownames(temp)
gene_short_name <- rownames(temp)
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

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
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
names(df4)[names(df4) == "gene_names"] <- "peptide"

write.csv(df4,file='~/neuropath_lineages/proteomics/BRAAK_DE_by_stateF.csv',quote=F,row.names=F)
file <- synapser::File(path='~/neuropath_lineages/proteomics/BRAAK_DE_by_stateF.csv', parentId='syn28712340')
file <- synapser::synStore(file)





#run Male monocle object with script below:
MonRun <- readRDS(file="~/neuropath_lineages/proteomics/BRAAKMale_monocleObject.rds")
temp <- readRDS(file="~/neuropath_lineages/proteomics/BRAAKMale_matrix_rownames.rds")
table(MonRun$State2)
#for braak males there are 6 states
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- rownames(temp)
gene_short_name <- rownames(temp)
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

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
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
names(df4)[names(df4) == "gene_names"] <- "peptide"

write.csv(df4,file='~/neuropath_lineages/proteomics/BRAAK_DE_by_stateM.csv',quote=F,row.names=F)
file <- synapser::File(path='~/neuropath_lineages/proteomics/BRAAK_DE_by_stateM.csv', parentId='syn28712340')
file <- synapser::synStore(file)






MonRun <- readRDS(file="~/neuropath_lineages/proteomics/CERADFemale_monocleObject.rds")
temp <- readRDS(file="~/neuropath_lineages/proteomics/CERADFemale_matrix_rownames.rds")
table(MonRun$State2)
#there are 6 states in the cerad female tree
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- rownames(temp)
gene_short_name <- rownames(temp)
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

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
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
names(df4)[names(df4) == "gene_names"] <- "peptide"

write.csv(df4,file='~/neuropath_lineages/proteomics/CERAD_DE_by_stateF.csv',quote=F,row.names=F)
file <- synapser::File(path='~/neuropath_lineages/proteomics/CERAD_DE_by_stateF.csv', parentId='syn28712340')
file <- synapser::synStore(file)






MonRun <- readRDS(file="~/neuropath_lineages/proteomics/CERADMale_monocleObject.rds")
temp <- readRDS(file="~/neuropath_lineages/proteomics/CERADMale_matrix_rownames.rds")
table(MonRun$State2)
#there are 6 states in the cerad male tree
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- rownames(temp)
gene_short_name <- rownames(temp)
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

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
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
names(df4)[names(df4) == "gene_names"] <- "peptide"

write.csv(df4,file='~/neuropath_lineages/proteomics/CERAD_DE_by_stateM.csv',quote=F,row.names=F)
file <- synapser::File(path='~/neuropath_lineages/proteomics/CERAD_DE_by_stateM.csv', parentId='syn28712340')
file <- synapser::synStore(file)

