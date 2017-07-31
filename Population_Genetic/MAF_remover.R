library('seqinr')

FastatoMatrix <- function(infile) {
  indata <- read.fasta(infile)
  ###getting gene name
  ###fasta file should end up with .fas
  gene <- sub(".fas","",infile)
  print(gene)
  
  ###get no of individuals
  nidv=length(indata)
  seqlength=getLength(indata[1])
  print(nidv)
  print(seqlength)
  
  ###convert fasta to matrix for easier use
  seq_matrix<-matrix(data = NA, nrow = (nidv), ncol = seqlength, byrow = FALSE,dimnames = NULL)
  print(dim(seq_matrix))
  
  ###Filling matrix with sequence
  for (i in 1:nidv) {seq_matrix[i,]<-getSequence(indata[[i]],as.string=FALSE)}
  print(dim(seq_matrix))
  
  ###name for indata
  seq_id <- c()
  for (i in 1:nidv) {seq_id <- c(seq_id,getName(indata[[i]]))}
  
  return(list(gene=gene,id=seq_id,matrix=seq_matrix))
}


minoremover <- function(dnacol){
  tabseq <- as.data.frame(table(dnacol),stringsAsFactors = FALSE)
  majorallele <- tabseq[which.max(tabseq$Freq),1]
  minorallele <- tabseq[which(tabseq$Freq<=3 & as.vector(tabseq[,1]) != 'n'),1]
  if ((length(minorallele)) != 0){
    for(i in 1:length(minorallele)){
      dnacol[dnacol==minorallele[i]] <- majorallele 
    }
  }
  return(dnacol)
}


###Make parser for all fasta file in the directory 
#Identify all fasta
FastatoID <- function(fasta) {
  indata <- read.fasta(fasta)
  nidv=length(indata)
  seqlength=getLength(indata[1])
  seq_id <- c()
  for (i in 1:nidv) {seq_id <- c(seq_id,getName(indata[[i]]))}
  return(seq_id)
}


###Function to put back matrix to fasta format and save it 
MatrixtoFasta <- function(ID,seq_mat,gene) {
##Writing fasta
filename=paste(gene,'_minrem',sep='')
cat("",file=filename, sep="",append=FALSE)
   for(i in 1:nrow(seq_mat)){
   cat('>',ID[i],'\n',sep='',file=filename,append=TRUE)
   cat(seq_mat[i,],'\n', sep='',file=filename,append=TRUE)
   }
}



#Opening all files with fasta identifier
fasfile <- dir(pattern='*fas')



###Looping for each column and adding into rem_matrix
for (i in 1:length(fasfile)){
  ###Information for fasta ID
  fastaID <- FastatoID(fasta=fasfile[i])
  dpseloc <- grep('Dpse',fastaID)
  dmirloc <- grep('Dmir',fastaID)
  ###Getting all fasta sequence
  all_seq <-FastatoMatrix(infile=fasfile[i])
  working_sequence <- all_seq$matrix
  lenseq <- length(all_seq$matrix[1,])
  noseq <- length(all_seq$matrix[,1])
  ###Setting matrix for minor allele 
  rem_matrix <- matrix(data=NA,ncol=lenseq,nrow=noseq)
  ###Iterate over column
  for (i in 1:lenseq) {
  ###work on dpse column
  dpse_rem <- minoremover(dnacol=working_sequence[dpseloc,i])
  ###then on dmir column
  dmir_rem <-  minoremover(dnacol=working_sequence[dmirloc,i])
  ###combine together into single column
  col_rem <- c(dpse_rem,dmir_rem)
  ###put back into matrix
  rem_matrix[,i] <- col_rem
  MatrixtoFasta(fastaID,rem_matrix,all_seq$gene)
   }  
}

