#### Differential Proteomics for ROSMAP TMT Samples

#setwd('~/Desktop/Projects/TMT_Proteomics/')
#usethis::use_testthat()


#' SynID <- 'syn23583548'
#' A function
#' 
#' @param SynID A vector
#' @param names A numeric
#' 
#' @return The dataframe pulled from synapse
#' @importFrom synapser synGet
#' @importFrom data.table fread
#' @export

TMT_Express_Load <- function( SynID, names ){
  #'@param  SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@param  names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  'You do not have READ permission for the requested entity' %in% as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) ) ){
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

#' A function
#' 
#' @param Input A matrix
#' 
#' @return A vector
#' @importFrom stats p.adjust
#' @export
Cleaner <- function( Input ){
  #Test if its a matrix
  if( !is.matrix(Input) ){
    return("Input is not a Matrix")
  }
  #test if it has 6 columns
  if( !(dim(Input)[2] == 6) ){
    return("Input not six columns")
  }
  #test if column 6 are p-values
  if( !( as.numeric(summary(as.numeric(Input[,6]))['Max.']) <= 1 & as.numeric(summary(as.numeric(Input[,6]))['Max.']) <= 1 ) ){
    return("Column Six Doesn't Appear to be P-Values")
  }
  
  colnames(Input) <- c( "Peptide", 'Coefficient', 'OR', 'CI_L', 'CI_H', "PVal" )
  Input <- as.data.frame(Input)
  Input$Peptide <- as.character(Input$Peptide)
  Input$Coefficient <- as.numeric(as.character(Input$Coefficient))
  Input$OR <- as.numeric(as.character(Input$OR))
  Input$CI_L <- as.numeric(as.character(Input$CI_L))
  Input$CI_H <- as.numeric(as.character(Input$CI_H))
  Input$PVal <- as.numeric(as.character(Input$PVal))
  
  
  Input$GeneID <- do.call( rbind,strsplit( Input$Peptide,'\\|' ) )[,1]
  Input$ProtID <- do.call( rbind,strsplit( Input$Peptide,'\\|' ) )[,2]
  Input <- Input[,c( 'Peptide', 'GeneID', 'ProtID', 'Coefficient', 'OR', 'CI_L', 'CI_H', 'PVal' ) ]
  Input$FDR_PVal <- p.adjust(Input$PVal, method = 'fdr', n = dim(Input)[1] )  
  return( Input )
}


#' A function
#' 
#'@param i Row Number to use in expression data
#'@param EXP - Expression data frame object colnames must equal MET rownames
#'@param MET - MetaData data frame object rownames must equal EXP colnames
#'@param Vars - Character vector of metadata columns to use in the model eg c("diagnosis", "msex", "APOE", "age_death", "pmi")
#'@param Diag - Character vector of metadata column to use as response variable must be binary factor eg. "diagnosis"
#'@param SampleID - Character vector of Meta Column name to use as sample ID eg. 'batchChannel'
#' 
#' @return A vector
#' @importFrom dplyr left_join
#' @importFrom stats confint
#' @importFrom stats coef
#' @export
#' 
Logistic_Model <- function( i, EXP, MET, Vars, Diag, SampleID ){
  #'@param i Row Number to use in expression data
  #'@param EXP - Expression data frame object colnames must equal MET rownames
  #'@param MET - MetaData data frame object rownames must equal EXP colnames
  #'@param Vars - Character vector of metadata columns to use in the model eg c("diagnosis", "msex", "APOE", "age_death", "pmi")
  #'@param Diag - Character vector of metadata column to use as response variable must be binary factor eg. "diagnosis"
  #'@param SampleID - Character vector of Meta Column name to use as sample ID eg. 'batchChannel'
  
  #Are the tables filered correctly
  if( !(as.numeric(table(row.names(MET)==colnames(EXP))[TRUE]) == length(row.names(MET))) ){
    return( "Error Meta Data and Expression Sample Row and Column Names Don't Align")
  }
  #Is entry a row
  if( !(as.numeric( i ) )){
    return( "i must be a row number")
  }
  #test if row exists
  if( i > dim(EXP)[1] ){
    return( "i must be a valid row number")
  }
  #Test Vars - Character Vector
  if( !(is.character(Vars) ) ){
    #Test Vars - Null
    if( (is.null(Vars) ) ){
    }else{
      return( "Vars Must be Character Vector or NULL")
    }
  }else{
    #Test Vars - Present in Meta Data
    if( !(is.null(Vars)) & !( as.numeric(table(Vars %in% colnames(MET))['TRUE'])== length(Vars) )  ){
      return( "Error Meta Data Vars Not Column Names in Meta Data Dataframe")
    }
  }
  #Test dataframes 
  if( !(is.data.frame(EXP) & is.data.frame(MET) ) ){
    return( "Error Meta Data and Expression Objects Must be Data Frames")
  }
  #Test Vars - Character Vector
  if( !(is.character(Diag) ) ){
    return( "Diag Must be Character")
  }
  #Test Vars - Character Vector
  if( !(length(Diag) == 1 ) ){
    return( "Diag Must Only 1 Value")
  }
  if( !(is.character(SampleID) ) ){
    return( "SampleID's Must be Character")
  }
  
  
  #Vars<-NULL
  #Form the Model data frame object
  EXP$batchChannel <- row.names(EXP)
  if( is.null(Vars)){
    Model_df <- dplyr::left_join(MET[ ,c(SampleID, Diag)], as.data.frame( t(rbind( batchChannel=colnames(EXP), GeneExp=EXP[i,] ))), by = SampleID )
  }else{
    Model_df <- dplyr::left_join(MET[ ,c(SampleID,Vars)], as.data.frame( t(rbind( batchChannel=as.character(colnames(EXP)), GeneExp=EXP[i,] ))), by = SampleID )
  }
  Model_df$GeneExp <- as.numeric(as.character(Model_df$GeneExp)) 
  
  if( is.null(Vars)){
    mylogit <- eval( parse( text=paste0( 'glm( ',Diag, ' ~ GeneExp, data = Model_df, family = \'binomial\')' )))
  }else{
    mylogit <- eval( parse( text=paste0( 'glm( ',Diag, ' ~ GeneExp ', paste0( c( " ", Vars[Vars != Diag]), collapse = " + "), ', data = Model_df, family = \'binomial\')' )))
  }
  
  Logit_Sum <- summary(mylogit)
  
  Gene <- row.names(EXP)[i]
  PVal <- Logit_Sum$coefficients['GeneExp','Pr(>|z|)']
  OR <- as.character( exp(coef(mylogit))['GeneExp'] )
  Coef <- Logit_Sum$coefficients['GeneExp','Estimate']
  CI_L <- confint(mylogit)[ 'GeneExp','2.5 %' ]
  CI_H <- confint(mylogit)[ 'GeneExp','97.5 %' ]
  
  return( c( Gene, Coef, OR, CI_L, CI_H, PVal) )
}

#' A function
#' 
#'@param GN a character ENSG Gene name eg 'VAMP1|P23763'
#'@param Path the character string of Neuropath column to use for model eg. EITHER: 'ceradsc', 'braaksc', 'cogdx', 'dcfdx_lv'
#'@param EXP Expression data frame
#'@param MET Metadata data frame
#' 
#' @return A vector
#' @importFrom dplyr left_join
#' @importFrom stats confint
#' @importFrom stats coef
#' @importFrom stats pnorm
#' @export
#'
NeuroPath_Calc <- function( GN,Path,Exp, MET ){
  #'@param GN a character ENSG Gene name eg 'VAMP1|P23763'
  #'@param Path the character string of Neuropath column to use for model eg. EITHER: 'ceradsc', 'braaksc', 'cogdx', 'dcfdx_lv'
  #'@param EXP Expression data frame
  #'@param MET Metadata data frame
  
  #Are the tables filered correctly
  if( !(as.numeric(table(row.names(MET)==colnames(Exp))[TRUE]) == length(colnames(Exp)) ) ){
    return( "Error Meta Data and Expression Sample Row and Column Names Don't Align")
  }
  #test if row exists
  if( !(GN %in% row.names(Exp)) ){
    return( "GN must be a valid row name in the expression dataset")
  }
  
  #Test Path - Character Vector
  if( !(is.character(Path) ) ){
    return( "Path Must be Character Vector")
  }
  #Test Path - Present in Meta Data
  if( !( Path %in% colnames(MET)) ){
    return( "Path Not a Valid Meta Data Column name")
  }
  #Test Path - length
  if( !( length(Path) == 1) ){
    return( "Path Can only be one variable at the moment")
  }
  
  #Test dataframes 
  if( !(is.data.frame(Exp) & is.data.frame(MET) ) ){
    return( "Error Meta Data and Expression Objects Must be Data Frames")
  }
  
  #if( !(is.character(SampleID) ) ){
  #  return( "SampleID's Must be Character")
  #}
  
  message(paste0("Processing: ", GN))
  
  Dat <- as.data.frame( cbind( MET, Gene = scale( as.numeric(Exp[GN,]) ) ), stringsAsFactors = F )
  
  if( Path %in% c( 'cogdx', 'dcfdx_lv') ){
    if( Path == 'cogdx'){
      Dat <- Dat[ Dat$cogdx %in% c(1,2,4), ]
    }else{
      Dat <- Dat[ Dat$dcfdx_lv %in% c(1,2,4), ]
    }
  }
  
  Dat <- Dat[ !is.na(Dat$Gene),]
  Dat$diagnosis <- as.factor(Dat$diagnosis)
  Dat$ceradsc <- as.factor(Dat$ceradsc)
  Dat$braaksc <- as.factor(Dat$braaksc)
  Dat$cogdx <- as.factor(Dat$cogdx)
  Dat$dcfdx_lv <- as.factor(Dat$dcfdx_lv)
  Dat$APOE <- as.factor(as.character(Dat$APOE))
  
  #Dat$pmi <- scale(log2(Dat$pmi))
  if( Path %in% c( 'cogdx', 'dcfdx_lv') ){
    Dat <- Dat[ !is.na(Dat$age_death), ] 
    Dat$APOE <- as.factor(as.character(Dat$APOE))
    Dat <- Dat[ , c('Gene','pmi', 'APOE', 'age_death', 'diagnosis','ceradsc','braaksc','cogdx','dcfdx_lv') ]
    #m <- eval(parse(text=paste0( 'MASS::polr(', Path, '~ Gene + pmi + APOE + diagnosis , data = Dat, na.action = na.omit, Hess=TRUE)' )))
    m <- eval(parse(text=paste0( 'MASS::polr(', Path, '~ Gene + pmi + diagnosis , data = Dat, na.action = na.omit, Hess=TRUE)' )))
  }else{
    Dat <- Dat[ , c('Gene','pmi', 'APOE', 'diagnosis','ceradsc','braaksc','cogdx','dcfdx_lv') ]
    m <- eval(parse(text=paste0( 'MASS::polr(', Path, '~ Gene + pmi + APOE + diagnosis, data = Dat, na.action = na.omit, Hess=TRUE)' )))
  }
  ctable <- coef(summary(m))
  ctable <- coef(summary(m))
  
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  ctable <- cbind(ctable, "p value" = p)
  
  ci <- confint( m, parm=as.vector("Gene") )
  
  Coeff <- m$coefficients['Gene']
  names( GN ) <- 'Gene'
  PVal <- p['Gene']
  names( PVal ) <- 'PVal'
  OR <- exp(coef(m))['Gene']
  names( OR ) <- 'OR'
  CI_L <- as.vector(ci[1])
  CI_H <- as.vector(ci)[2]
  return( as.vector(c( GN, Coeff, OR, CI_L, CI_H, PVal ) ))
}



#' A function
#' 
#'@param y 
#'@param x 
#'@param nsamp a numeric
#'@param cores a numeric
#' 
#' @return A matrix
#' @importFrom dplyr left_join
#' @importFrom spike fastlmbeta2
#' @importFrom parallel makeCluster
#' @importFrom parallel parApply
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @importFrom vbsr vbsr
#' @importFrom spike fastlmbeta2
#' @export
#'
spvbsrBootstrap = function(y,x,nsamp=100,cores=8){
  library(dplyr)
  library(parallel)
  library(doParallel)
  n <- length(y)
  replicateMatrix = sample(1:(n*nsamp),replace=TRUE) %>%
    matrix(n,nsamp)
  replicateMatrix = replicateMatrix%%n +1
  #print(replicateMatrix[1:5,])
  fxn1 <- function(shuf,y,x){
    #library(utilityFunctions)
    #return(utilityFunctions::fastlmbeta(y[shuf],x[shuf,]))
    #library(vbsr)
    #y=y[shuf]
    #x=x[shuf,]
    res <- vbsr::vbsr(y[shuf],x[shuf,])
    
    ###identify significant features
    ###return unpenalized betas
    baz = rep(0,ncol(x)+1)
    names(baz) = c('intercept',colnames(x))
    whichSig = which(res$pval < 0.05/ncol(x))
    if(length(whichSig)>0){
      #write a fastlmbeta function that also returns correlation
      baz2 = c()
      try(baz2 <- spike::fastlmbeta2(y[shuf],x[shuf,whichSig],colnames(x)[whichSig]),silent=T)
      if(length(baz2)>0){
        baz[names(baz2)] <- baz2
      }
    }
    
    #return(c(res$alpha,res$beta))
    return(baz)
  }
  
  cl <- parallel::makeCluster(cores)
  registerDoParallel(cl)
  #betaMatrix <- foreach(i=1:nsamp,.combine='rbind') %dopar% fxn1(replicateMatrix[,i],y,data.matrix(x))
  betaMatrix <- t(parallel::parApply(cl,replicateMatrix,2,fxn1,y,data.matrix(x)))
  #betaMatrix <- t(apply(replicateMatrix,2,fxn1,y,data.matrix(x)))
  parallel::stopCluster(cl)
  #colnames(betaMatrix) <- c(colnames(x))
  return(apply(betaMatrix!=0,2,mean))
  #return(betaMatrix)
  
}

