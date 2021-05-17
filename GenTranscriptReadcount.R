# Ted's version of ExCluster R script for randomizing transcript levels
library(stringi)
setwd('/Users/tperkins/Mine/MyProjects/Matt_Excluster_Project/CodeTJP1/')

# install.packages('DirichletReg')
library(DirichletReg)

#################################
# COMMAND LINE ARGUMENT PARSING #
# This function reads the command line arguments and returns them in a list
# where each is accessible by name. It fills in default values where 
# possible. 
GetArguments <- function() {
  # Default values
  GGTFile <- NA
  OutStem <- "./RT"
  RSeed <- 1
  ReadsPath <- NA
  Reads1 <- NA
  Reads2 <- NA
  NumSpliced <- 2 # Number of genes to be differentially spliced
  MinTransSpliced <- 2 # Minimum number of transcripts for a differentially spliced gene
  MinExprSpliced <- 200 # Minimum expression in rpm for a differentially spliced gene
  NumTransExprLo <- 2 # Minimum number of transcripts expressed per gene
  NumTransExprHi <- 4 # Maximum number of transcripts expressed per gene
  Dispersion <- 1000 # Dispersion parameter for biological variability
  
  
  # Get the arguments
  #args = commandArgs(trailingOnly=TRUE)
  args = c('GGTFile=TestGGT2.txt',
           'ReadsPath=/Users/tperkins/Mine/MyProjects/Matt_Excluster_Project/CodeTJP1',
           'Reads1=TestGE2_C1_Rep',
           'Reads2=TestGE2_C2_Rep');
  
  # Parse the arguments
  for(i in 1:length(args)) {
    arg = args[i]
    ArgRecognized <- FALSE
    
    # Testing for GGTFile
    if (nchar(arg)>=7) {
      if (stri_cmp_eq(substr(arg,1,7),'GGTFile')) {
        GGTFile <- substr(arg,9,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for OutStem
    if (nchar(arg)>=7) {
      if (stri_cmp_eq(substr(arg,1,7),'OutStem')) {
        OutStem <- substr(arg,9,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for RSeed
    if (nchar(arg)>=5) {
      if (stri_cmp_eq(substr(arg,1,5),'RSeed')) {
        RSeed <- as.numeric(substr(arg,7,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for ReadsPath
    if (nchar(arg)>=9) {
      if (stri_cmp_eq(substr(arg,1,9),'ReadsPath')) {
        ReadsPath <- substr(arg,11,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for Reads1
    if (nchar(arg)>=6) {
      if (stri_cmp_eq(substr(arg,1,6),'Reads1')) {
        Reads1 <- substr(arg,8,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for Reads2
    if (nchar(arg)>=6) {
      if (stri_cmp_eq(substr(arg,1,6),'Reads2')) {
        Reads2 <- substr(arg,8,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for NumSpliced
    if (nchar(arg)>=10) {
      if (stri_cmp_eq(substr(arg,1,10),'NumSpliced')) {
        NumSpliced <- as.numeric(substr(arg,12,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for MinTransSpliced
    if (nchar(arg)>=15) {
      if (stri_cmp_eq(substr(arg,1,15),'MinTransSpliced')) {
        MinTransSpliced <- as.numeric(substr(arg,17,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for MinExprSpliced
    if (nchar(arg)>=14) {
      if (stri_cmp_eq(substr(arg,1,14),'MinExprSpliced')) {
        MinExprSpliced <- as.numeric(substr(arg,16,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for NumTransExprLo
    if (nchar(arg)>=14) {
      if (stri_cmp_eq(substr(arg,1,14),'NumTransExprLo')) {
        NumTransExprLo <- as.numeric(substr(arg,16,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for NumTransExprHi
    if (nchar(arg)>=14) {
      if (stri_cmp_eq(substr(arg,1,14),'NumTransExprHi')) {
        NumTransExprHi <- as.numeric(substr(arg,16,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for Dispersion
    if (nchar(arg)>=10) {
      if (stri_cmp_eq(substr(arg,1,10),'Dispersion')) {
        Dispersion <- as.numeric(substr(arg,12,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    if (!(ArgRecognized)) {
      print("Warning: Command line argument not recognized: ")
      print(arg)
    }
  }
  
  Args <- list("GGTFile" = GGTFile,
               "OutStem" = OutStem,
               "RSeed" = RSeed,
               "ReadsPath" = ReadsPath,
               "Reads1" = Reads1,
               "Reads2" = Reads2,
               "NumSpliced" = NumSpliced,
               "MinTransSpliced" = MinTransSpliced,
               "MinExprSpliced" = MinExprSpliced,
               "NumTransExprLo" = NumTransExprLo,
               "NumTransExprHi" = NumTransExprHi,
               "Dispersion" = Dispersion
  )
  return(Args)
}

########################################################
# READING READ COUNTS FROM FILES FOR A SINGLE CONDITION
# This function reads in the all the files in ReadsPath where the file
# begins with the string FNameStart. Such files, one or more, are 
# assumed to be text files with first column containing gene names, and
# second column containing non-negative integer counts. If there are
# multiple files (replicates/patients/etc.) they are assumed to list the
# same genes and in the same order.
# NOTE: In the data frame returned, the gene names are the row names and
# 1, 2, 3, ... are the column names.
# NOTE: The first line of each input file is assumed to be a header line,
# and is ignored.
# Returns NULL if no such files can be found
ConditionReadCounts <- function(ReadsPath,FNameStart) {
  
  # Figure out all the file names in the right directory, and starting
  # in the right way.
  FNames <- list.files(path=ReadsPath,pattern=FNameStart,full.names = TRUE)
  
  # If no files found, return NULL
  if (length(FNames)==0) {
    return(NULL)
  }
  
  # Read in the first file
  RC <- read.table(file=FNames[1],header=TRUE,stringsAsFactors=FALSE)
  RC <- data.frame(RC[,2],row.names=RC[,1],stringsAsFactors = FALSE)
  
  # Read in the remaining files, appending their gene counts in additional columns
  # Record also the column names
  for (i in 2:length(FNames)) {
    TempRC <- read.table(file=FNames[i],header=TRUE,stringsAsFactors=FALSE)
    RC <- cbind(RC,TempRC[,2])
  }
  
  # Fill in the column names
  colnames(RC) <- 1:length(FNames)
  
  return(RC)
}

#############################################################
# GETTING READ COUNTS FOR BOTH CONDITIONS AND ALL REPLICATES
# We also compute, for later convenience, reads per million (RPM)
# mean RPM in each condition, and the mean of that across conditions.
GetRCandRPM <- function(ReadsPath,Reads1,Reads2) {
  
  # Read counts
  RC1 <- ConditionReadCounts(ReadsPath,Reads1)
  RC2 <- ConditionReadCounts(ReadsPath,Reads2)
  RC <- cbind(RC1,RC2)
  
  # Dimensions (for NGenes and NConditions)
  D1 <- dim(RC1)
  D2 <- dim(RC2)
  
  # Reads per million
  RPM1 <- RC1
  for (j in 1:D1[2]) {
    RPM1[,j] <- RPM1[,j]*1000000/sum(RPM1[,j])
  }
  RPM2 <- RC2
  for (j in 1:D2[2]) {
    RPM2[,j] <- RPM2[,j]*1000000/sum(RPM2[,j])
  }
  RPM <- cbind(RPM1,RPM2)
  
  # Mean RPM for each condition
  MRPM1 <- RPM1[,1]
  for (j in 2:D1[2]) {
    MRPM1 <- MRPM1 + RPM1[,2]
  }
  MRPM1 <- MRPM1 / D1[2]
  MRPM2 <- RPM2[,1]
  for (j in 2:D2[2]) {
    MRPM2 <- MRPM2 + RPM2[,2]
  }
  MRPM2 <- MRPM2 / D2[2]
  
  # Overall mean RPM
  MRPM <- (MRPM1+MRPM2)/2
  
  # Return value
  L <- list("NGenes"=D1[1],
            "NC1"=D1[2],
            "NC2"=D2[2],
            "RC1"=RC1, 
            "RC2"=RC2,
            "RC"=RC,
            "RPM1"=RPM1, 
            "RPM2"=RPM2,
            "RPM"=RPM,
            "MRPM1"=MRPM1,
            "MRPM2"=MRPM2,
            "MRPM"=MRPM)
}

###########################################
# GETTING GENE AND TRANSCRIPT DEFINITIONS #
# This function reads in a "GGT" table with three columns: gene ID,
# gene name, transcript ID. The GGT table may have additional columns, but they
# are ignored. The function then processes these data into two separate tables
# for convenience in future calculations. Specifically, it provides
# Genes, which has a column of unique gene IDs, gene name, number of
#  transcripts, and start and end index of transcripts in the Transcripts table.
# Transcripts, which has a column of gene IDs, gene names, and transcript names --
#  basically just a copy of the first three columns of the input file, but 
#  potentially re-ordered for consistency with gene table sorting.
ReadGGT <- function(GGTFName) {
  
  # Read in the file
  F <- read.table(file=GGTFName,header=TRUE,stringsAsFactors=FALSE)
  
  # Find out the unique gene IDs
  UniqueGeneIDs <- sort(as.character(unique(F[,1])))
  NGenes <- length(UniqueGeneIDs)
  NTrans <- dim(F)[1]
  
  # Create the gene table
  GeneTable <- data.frame(matrix(NA,length(UniqueGeneIDs),5))
  colnames(GeneTable) <- c("GeneID","GeneName","NTrans","TransStart","TransEnd")
  
  # Create the transcript table
  TransTable <- data.frame(matrix(NA,length(nrow(F)),3))
  colnames(TransTable) <- c("GeneID","GeneName","TransID")
  TransIndex <- 1
  
  # Loop through
  for (gi in 1:NGenes) {
    
    # What gene are we talking about?
    GeneID <- UniqueGeneIDs[gi]
    GeneTable[gi,1] <- GeneID
    GeneTable[gi,3] <- 0
    
    # Loop through the table looking for that GeneID
    FoundFirstTrans = FALSE
    for (ti in 1:NTrans) {
      # Checking for same GeneID
      if (stri_cmp_eq(GeneID,F[ti,1])) {
        # Is this the first time we're finding it? If so, set gene name
        # and transcript start 
        if (!FoundFirstTrans) {
          FoundFirstTrans <- TRUE
          GeneTable[gi,2] <- F[ti,2]
          GeneTable[gi,4] <- TransIndex
        }
        # Account for the newly found transcript
        GeneTable[gi,3] <- GeneTable[gi,3]+1
        GeneTable[gi,5] <- TransIndex
        #print(TransIndex)
        #print(GeneID)
        TransTable[TransIndex,1] <- GeneID
        TransTable[TransIndex,2] <- GeneTable[gi,2]
        TransTable[TransIndex,3] <- F[ti,3]
        TransIndex <- TransIndex+1
      }
    }
  }
  return(list("GeneTable"=GeneTable,"TransTable"=TransTable))
}

#######################################
# PICKING DIFFERENTIALLY SPLICED GENES
# Here we randomly pick a specified number of genes to be differentially
# spliced between conditions 1 and 2. Those genes must satisfy certain
# properties, particularly, having a minimum expression level, and
# having a minimum number of distinct transcripts. These parameters
# are all specified as arguments. The function returns a logical (TRUE/FALSE)
# vector of length equal to the number of genes in the GGT table.
# WARNING: If the user asks for a number of differentially spliced
# genes that is impossible, (say, because minimum expression or transcript 
# number thresholds are too high) then we print a warning and return NULL.
ChooseDiffSplicedGenes <- function(GGT,MeanExpression,NumSpliced,MinNTranscripts,MinExpression) {
  # How many genes are there?
  NGenes <- nrow(GGT$GeneTable)
  
  # Initializing output to a vecto of FALSE
  SplicedOrNo <- rep(FALSE,NGenes)
  
  # Create a random permutation of the genes
  RP <- sample(NGenes)
  
  # March through the permutation until we find enough genes satisfying the conditions
  NFound <- 0
  NextTryIndex <- 1
  while ((NFound<NumSpliced) & (NextTryIndex<=NGenes)) {
    
    # Which gene are we testing?
    NextTry <- RP[NextTryIndex]
    
    # Check if this gene passes criteria
    if ((MeanExpression[NextTry]>=MinExpression) & (GGT$GeneTable$NTrans[NextTry]>=MinNTranscripts)) {
      # Great, we found one!
      SplicedOrNo[NextTry] <- TRUE
      NFound <- NFound+1
    }
    
    # Increment to next gene
    NextTryIndex <- NextTryIndex+1
  }
  
  # Check if we found enough genes
  if (NFound<NumSpliced) {
    print('WARNING: Not enough genes could be selected for differential splicing!')
    SplicedOrNo <- NULL
  }
  
  return(SplicedOrNo)
}



#######################################
# ASSIGNING PERCENTAGES TO TRANSCRIPTS
# Given the GGT table, we do the following to randomly assign percentages
# to each transcript in each gene, producing a "transcript fraction
# table" (TFT). The TFT has rows corresponding to each transcript in the GGT,
# two columns stating the transcript percentages for condition 1 and 2,
# and columns corresponding to each condition. In each of those is listed
# a fraction between 0 and 1, representing the fraction of the gene's
# expression accounted for by each transcript.
# For either expressed or non-expressed genes, we first randomly choose
# NumTransExprLo/Hi transcripts to be expressed (i.e. >0 fraction of expression).
# The number of transcripts to be expressed is chosen randomly in the range
# min(NumTransExprLo,Z):min(NumTransExprHi,Z), where Z is the total number of
# transcripts for the gene.
# If a gene is not-spliced, fractions will be assigned to those transcripts,
# the same for conditions 1 and 2, else fractions will be chosen randomly for
# both conditions.
# Per-condition fractions are then slightly perturbed versions of those.
MakeTFT <- function(GGT,SplicedOrNo,NumTransExprLo,NumTransExprHi,NCond1,NCond2,Dispersion) {
  # Set up the column names for the data frame
  ColNames <- c('Cond1','Cond2')
  for (i in 1:NCond1) {
    ColNames <- c(ColNames,paste('Cond1_Rep',as.character(i),sep=''))
  }
  for (i in 1:NCond2) {
    ColNames <- c(ColNames,paste('Cond2_Rep',as.character(i),sep=''))
  }
  
  # Make datafame
  NTrans <- nrow(GGT$TransTable)
  Fracs <- matrix(0,nrow=NTrans,ncol=2+NCond1+NCond2)
  TFT <- data.frame(Fracs)
  colnames(TFT) <- ColNames
  
  # Fill it in, gene by gene!
  NGenes <- nrow(GGT$GeneTable)
  for (g in 1:NGenes) {
    NTranscripts <- GGT$GeneTable$NTrans[g]
    GeneStartIndex <- GGT$GeneTable$TransStart[g]
    
    # Choose how many transcripts to express
    NTElo <- min(NumTransExprLo,NTranscripts)
    NTEhi <- min(NumTransExprHi,NTranscripts)
    if (NTElo==NTEhi) {
      NTE = NTElo
    } else {
      NTE = as.integer(round(NTElo+(NTEhi-NTElo)*runif(1)))
    }
    
    # Choose the particular transcripts to express
    RandPerm <- sample(NTranscripts)
    ToExpress <- RandPerm[1:NTE]
    
    # Choose condition 1 expression levels
    Cond1E <- runif(NTE)
    Cond1E <- Cond1E / sum(Cond1E)
    for (j in 1:NTE) {
      TFT[GeneStartIndex+ToExpress[j]-1,1] <- Cond1E[j]
    }
    
    # Choose condition 2 expression levels -- if gene is spiced, we choose new ones, 
    # else we keep the condition 1 ones.
    if (SplicedOrNo[g]) {
      Cond2E <- runif(NTE)
      Cond2E <- Cond2E / sum(Cond2E)
    } else {
      Cond2E <- Cond1E
    }
    for (j in 1:NTE) {
      TFT[GeneStartIndex+ToExpress[j]-1,2] <- Cond2E[j]
    }
    
    # Choose the condition 1 replicate expression levels
    for (c in 1:NCond1) {
      RepE <- rdirichlet(1,Cond1E*Dispersion)
      for (j in 1:NTE) {
        TFT[GeneStartIndex+ToExpress[j]-1,2+c] <- RepE[j]
      }
    }
    
    # Choose the condition 2 replicate expression levels
    for (c in 1:NCond2) {
      RepE <- rdirichlet(1,Cond2E*Dispersion)
      for (j in 1:NTE) {
        TFT[GeneStartIndex+ToExpress[j]-1,2+NCond1+c] <- RepE[j]
      }
    }
  }
  
  # Add gene and transcript info to the front
  TFT <- cbind(GGT$TransTable,TFT)
  
  return(TFT)
}

###########################################
# ASSIGNING READ COUNTS TO EACH TRANSCRIPT
# Taking as input tables of reads per gene and of expression fractions per
# transcript, we use multinomial sampling do determine the number of reads
# per transcript. Rows are in the same order as the transcript table.
# Read counts per gene are taken to be a random multinomial sample with
# the total number of reads being equal to the specified read count in

MakeTRT <- function(GGT,TFT,Expression) {
  NConds <- Expression$NC1 + Expression$NC2
  
  # Data frame of the right size, of all zeros
  TRT <- data.frame(matrix(0,nrow(GGT$TransTable),NConds))
  colnames(TRT) <- colnames(TFT)[6:(5+NConds)]
  
  # Loop over the gene and conditions, sampling reads
  NGenes <- nrow(GGT$GeneTable)
  for (g in 1:NGenes) {
    
    # How many transcripts in this gene, and start and stop indeces
    GeneStartIndex <- GGT$GeneTable$TransStart[g]
    GeneEndIndex <- GGT$GeneTable$TransEnd[g]
    
    # Loop over the conditions
    for (c in 1:NConds) {
      
      # If the gene has reads at all
      if (Expression$RC[g,c]>0) {
        
        # Multinomially sample from the transcript fractions, and slot them
        # into the TRT 
        Multi <- rmultinom(1,Expression$RC[g,c],TFT[GeneStartIndex:GeneEndIndex,5+c])
        TRT[GeneStartIndex:GeneEndIndex,c] <- Multi
        #for (j in 1:NTrans) {
        #  TRT[GeneStartIndex+j-1,c] <- Multi[j]
        #}
      }
    }
  }
  
  # Put gene and transcript info on the front
  TRT <- cbind(GGT$TransTable,TRT)
  
  return(TRT)
}

###############
# MAIN SCRIPT #
###############

# Read the command line arguments
print("Getting command line arguments...")
Args = GetArguments()
#print(Args)

# Set the random number generator seed
print("Setting random number seed...")
set.seed(as.numeric(Args$RSeed))

# Read in the GGT file
print("Reading GGT file...")
GGT <- ReadGGT(Args$GGTFile)
#print(GGT)

# Read in the read counts data
print("Readig read counts data...")
Expr <- GetRCandRPM(Args$ReadsPath,Args$Reads1,Args$Reads2)
#print(Expr)

# Choose differentially spliced genes
print("Choosing differentially spliced genes...")
S <- ChooseDiffSplicedGenes(GGT,Expr$MRPM,Args$NumSpliced,Args$MinTransSpliced,Args$MinExprSpliced)
print(S)

TFT <- MakeTFT(GGT,S,Args$NumTransExprLo,Args$NumTransExprHi,Expr$NC1,Expr$NC2,Args$Dispersion)
print(TFT)

TRT <- MakeTRT(GGT,TFT,Expr)
print(TRT)

# Saving the output to file
write.csv(TFT,paste(Args$OutStem,"_TFT.csv"),row.names = FALSE)
write.csv(TRT,paste(Args$OutStem,"_TRT.csv"),row.names = FALSE)

