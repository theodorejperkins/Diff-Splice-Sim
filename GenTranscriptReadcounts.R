# Ted's version of ExCluster R script for randomizing transcript levels
library(stringi)
setwd('/Users/tperkins/Mine/MyProjects/Matt_Excluster_Project/CodeTJP1/')

#################################
# COMMAND LINE ARGUMENT PARSING #
# This function reads the command line arguments and returns them in a list
# where each is accessible by name. It fills in default values where 
# possible. If a needed argument (i.e. one without a default value)
# is not specified, it prints and error message, and returns FALSE.
GetArguments <- function() {
  # Default values
  exonCov <- 0.5
  SD <- 2
  CDS_File <- NA
  #GFF_File <- NA
  #GTF_File <- NA
  GTE_File <- NA
  OutFolder <- "./"
  RSeed <- 1
  RepNum <- 1
  ReadsPath <- NA
  Reads1 <- NA
  Reads2 <- NA
  OutDir <- NA # How is this different from OutFolder?
  NumSpliced <- 1000
  MC <- 10
  SeqLength <- 100
  ReadType <- "pe"
  
  # Get the arguments
  #args = commandArgs(trailingOnly=TRUE)
  args = c('GTE_File=TestGTE.txt')
  
  # Parse the arguments
  for(i in 1:length(args)) {
    arg = args[i]
    ArgRecognized <- FALSE

    # Testing for exonCov    
    if (nchar(arg)>=7) {
      if (stri_cmp_eq(substr(arg,1,7),'exonCov')) {
        exonCov <- as.numeric(substr(arg,9,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for SD
    if (nchar(arg)>=2) {
      if (stri_cmp_eq(substr(arg,1,2),'SD')) {
        SD <- as.numeric(substr(arg,4,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for CDS_File
    if (nchar(arg)>=8) {
      if (stri_cmp_eq(substr(arg,1,8),'CDS_File')) {
        CDS_File <- substr(arg,10,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for GFF_File
    #if (nchar(arg)>=8) {
    #  if (stri_cmp_eq(substr(arg,1,8),'GFF_File')) {
    #    GFF_File <- substr(arg,10,nchar(arg))
    #    ArgRecognized <- TRUE
    #  }
    #}

    # Testing for GTF_File
    #if (nchar(arg)>=8) {
    #  if (stri_cmp_eq(substr(arg,1,8),'GTF_File')) {
    #    GTF_File <- substr(arg,10,nchar(arg))
    #    ArgRecognized <- TRUE
    #  }
    #}
    
    # Testing for GTE
    if (nchar(arg)>=8) {
      if (stri_cmp_eq(substr(arg,1,8),'GTE_File')) {
        GTE_File <- substr(arg,10,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for OutFolder
    if (nchar(arg)>=9) {
      if (stri_cmp_eq(substr(arg,1,9),'OutFolder')) {
        OutFolder <- substr(arg,11,nchar(arg))
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
    
    # Testing for RepNum
    if (nchar(arg)>=6) {
      if (stri_cmp_eq(substr(arg,1,6),'RepNum')) {
        RepNum <- as.numeric(substr(arg,8,nchar(arg)))
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
    
    # Testing for OutDir
    if (nchar(arg)>=6) {
      if (stri_cmp_eq(substr(arg,1,6),'OutDir')) {
        OutDir <- substr(arg,8,nchar(arg))
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
    
    # Testing for MC
    if (nchar(arg)>=2) {
      if (stri_cmp_eq(substr(arg,1,2),'MC')) {
        MC <- as.numeric(substr(arg,4,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for SeqLength
    if (nchar(arg)>=9) {
      if (stri_cmp_eq(substr(arg,1,9),'SeqLength')) {
        SeqLength <- as.numeric(substr(arg,11,nchar(arg)))
        ArgRecognized <- TRUE
      }
    }
    
    # Testing for ReadType
    if (nchar(arg)>=8) {
      if (stri_cmp_eq(substr(arg,1,8),'ReadType')) {
        ReadType <- substr(arg,10,nchar(arg))
        ArgRecognized <- TRUE
      }
    }
    
    if (!(ArgRecognized)) {
      print("Warning: Command line argument not recognized: ")
      print(arg)
    }
  }

  Args <- list("exonCov" = exonCov, 
               "SD" = SD,
               "CDS_File" = CDS_File,
               #"GFF_File" = GFF_File,
               #"GTF_File" = GTF_File,
               "GTE_File" = GTE_File,
               "OutFolder" = OutFolder,
               "RSeed" = RSeed,
               "RepNum" = RepNum,
               "ReadsPath" = ReadsPath,
               "Reads1" = Reads1,
               "Reads2" = Reads2,
               "OutDir" = OutDir,
               "NumSpliced" = NumSpliced,
               "MC" = MC,
               "SeqLength" = SeqLength,
               "ReadType" = ReadType
               )
  return(Args)
}

##################################
# READING READ COUNTS FROM FILES #
# This function reads in the all the files in ReadsPath where the file
# begins with the string FNameStart. Such files, one or more, are 
# assumed to be text files with first column containing gene names, and
# second column containing non-negative integer counts. If there are
# multiple files (replicates/patients/etc.) they are assumed to list the
# same genes and in the same order.
# NOTE: In the data frame return, the gene names are in the first column, 
# rather than being row names.
# Returns FALSE if no such files can be found
ReadReadCounts <- function(ReadsPath,FNameStart) {

  # Figure out all the file names in the right directory, and starting
  # in the right way.
  FNames <- list.files(path=ReadsPath,pattern=FNameStart,full.names = TRUE)
  
  # If no files found, return FALSE
  if (length(FNames)==0) {
    return(FALSE)
  }
  
  # Read in the first file
  RC <- read.table(file=FNames[1],header=TRUE,stringsAsFactors=FALSE)
  RC <- data.frame(RC,stringsAsFactors = FALSE)
  ColNames <- colnames(RC)
  
  # Read in the remaining files, appending their gene counts in additional columns
  # Record also the column names
  for (i in 2:length(FNames)) {
    TempRC <- read.table(file=FNames[i],header=TRUE,stringsAsFactors=FALSE)
    RC[,i+1] <- TempRC[,2]
    ColNames[i+1] <- colnames(TempRC)[2]
  }
  
  # Fill in the column names
  colnames(RC) <- ColNames
  return(RC)
}

##################################################
# GETTING GENE, TRANSCRIPT, AND EXON DEFINITIONS #
# This function reads in a "GTE" table with three columns: gene,
# transcript, and exon. The function then processes it in several 
# ways convenient for future calculations
ReadGTE <- function(GTEFName) {
  F <- read.table(file=GTEFName,header=TRUE,stringsAsFactors=FALSE)
  print(F)

  GTE <- NULL
  # Copy over gene names, uniquifying them
  GTE$Genes <- sort(as.character(unique(F[,1])))
  GTE$NGenes <- length(GTE$Genes)

  # Create a list of the transcripts for each gene. This could, and should!,
  # be done more efficiently.
  GTE$TransPerGene <- vector(mode='list',GTE$NGenes)
  GTE$NTransPerGene <- vector(mode='list',GTE$NGenes)
  for (gi in 1:GTE$NGenes) {
    # What gene are we talking about?
    GeneID <- GTE$Genes[gi] 
    # Create a new, empty list of its transcripts
    TransList <- c()
    for (li in 1:dim(F)[1]) {
      if (stri_cmp_eq(GeneID,F[li,1])) {
        TransList <- c(TransList,F[li,2])
      }
    }
    TransList <- unique(TransList)
    print(GeneID)
    print(TransList)
    GTE$TransPerGene[[gi]] <- TransList
    GTE$NTransPerGene[[gi]] <- length(TransList)
  }
  return(GTE)
}


# Matt does a ton of stuff manipulating different annotation files. Really
# that should be done once and for all, and we should have a clear definition
# for the input files that we actually need.

# Differential splicing simulation

# Randomly chooses some genes to be differentially spliced. They must have
# >= 200 average reads per condition, and at least three transcripts

# Then he has a really weird way of generating log2 fold changes from random numbers
# Why not just generate random normal numbers with absolute value >= 1?



# Tested somewhat as Dec 22, 2020
Args = GetArguments()
print(Args)

# Testing and working as of 10am, Dec 23, 2020
#ReadCounts1 <- ReadReadCounts(Args$ReadsPath,Args$Reads1)
#print(ReadCounts1)

# Testing GTE input. Working as of 1:40pm, Dec 31, 2020
GTE <- ReadGTE(Args$GTE_File)
print(GTE)