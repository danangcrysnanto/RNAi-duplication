library('seqinr')
FastatoPHASE <- function(infile) {
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
indata_matrix<-matrix(data = NA, nrow = (nidv), ncol = seqlength, byrow = FALSE,dimnames = NULL)
print(dim(indata_matrix))

###Filling matrix with sequence
for (i in 1:nidv) {indata_matrix[i,]<-getSequence(indata[[i]],as.string=FALSE)}
print(dim(indata_matrix))

###name for indata
indata_id <- c()
for (i in 1:nidv) {indata_id <- c(indata_id,getName(indata[[i]]))}
###Writetofile
filename=paste(gene,'_id',sep="")
write.table(indata_id,file=filename)

####Create random pseudohaplotype
####haplotype split
hapsplit<-matrix(data = NA, nrow = (2*nidv), ncol = seqlength, byrow = FALSE,dimnames = NULL)

###Filling all with the same nuc but twice 
for (i in 1:nidv) {
hapsplit[2*i-1,] <- indata_matrix[i,]
hapsplit[2*i,] <- indata_matrix[i,]
}

####Saving information for site with n (to be recovered after phase)
###This will perform on operation on cell-to-cell basis
hapsplitno <- hapsplit == 'n'
###saved in a table
filename=paste(gene,'_hapsplitno',sep='')
write.table(hapsplitno,file=filename)

###Setting random filling for heterozygous site
for (i in 1:nidv){
for (j in 1:seqlength){
if (indata_matrix[i,j]=="y"){
hapsplit[2*i-1,j] <- 'c'
hapsplit[2*i,j] <- 't'
}
if (indata_matrix[i,j]=="r"){
hapsplit[2*i-1,j] <- 'a'
hapsplit[2*i,j] <- 'g'
}
if (indata_matrix[i,j]=="m"){
hapsplit[2*i-1,j] <- 'c'
hapsplit[2*i,j] <- 'a'
}
if (indata_matrix[i,j]=="s"){
hapsplit[2*i-1,j] <- 'c'
hapsplit[2*i,j] <- 'g'
}
if (indata_matrix[i,j]=="k"){
hapsplit[2*i-1,j] <- 'g'
hapsplit[2*i,j] <- 't'
}
if (indata_matrix[i,j]=="w"){
hapsplit[2*i-1,j] <- 'a'
hapsplit[2*i,j] <- 't'
}
if (indata_matrix[i,j]=="n"){
hapsplit[2*i-1,j] <- '?'
hapsplit[2*i,j] <- '?'
}
}
}


###Looking for polymorphic sites
#c() to create an empty vector
poly <- c()
#Create vector for count polymorphic sites but exclude ?
for (i in 1:seqlength){
if (is.element('?',hapsplit[,i])){
nounique<-length(unique(hapsplit[,i]))-1
poly <- c(poly,nounique)
}
else {
nounique<-length(unique(hapsplit[,i]))
poly <- c(poly,nounique)
}
}


###storing multisite haplo information
###saved in haplomulti
multisite <- which(poly>2)
haplomulti <- matrix(data=NA,nrow=length(hapsplit[,1]),ncol=length(multisite))
for (i in 1:length(multisite)){
haplomulti[,i] <- hapsplit[,multisite[i]] 
}

###Writetofile
filename=paste(gene,'_multisite',sep='')
write.table(multisite,file=filename)

###Writetofile
filename=paste(gene,'_haplomulti',sep='')
write.table(haplomulti,file=filename)



####masking multiallelic site (>2 alleles) with 1
for (i in 1:length(multisite)){
hapsplit[,multisite[i]] <- rep('1',length(hapsplit[,1]))
}


filename=paste(gene,'_unphased.inp',sep='')
###write single letter, cat for contatenate to vector need to use \n
write(nidv,file=filename,append=FALSE)
###write seqlength
write(seqlength,file=filename,append=TRUE)

###Loop for genotype
for (i in 1:nidv) {
#write the haplotypes 
cat(hapsplit[(2*i-1),],"\n",file=filename, sep="",append=TRUE)
cat(hapsplit[(2*i),],"\n",file=filename, sep="",append=TRUE)

	}

}


PHASEtoFasta <- function(infile) {
###Opening generated fasta
gene <- sub(".phased","",infile)
datout <- read.table(infile,sep='')  
dat2 <- as.matrix(datout)
multisite<-read.table(paste(gene,"_multisite",sep=""),stringsAsFactors = FALSE)
hapsplitno<-read.table(paste(gene,"_hapsplitno",sep=""),stringsAsFactors = FALSE)
haplomulti<-read.table(paste(gene,"_haplomulti",sep=""),stringsAsFactors = FALSE)
geneid <-read.table(paste(gene,"_id",sep=""),stringsAsFactors = FALSE)


###Reinsert multi-site
for (i in 1:length(multisite[,1])){
idx <- multisite[i,1]
dat2[,idx] <- haplomulti[,i]  
}
multisite
dim(dat2)


###Reinsert n
for (i in 1:nrow(hapsplitno)) {
idmis <- which(hapsplitno[i,] == TRUE)
dat2[i,idmis] <- 'n'
}

##Writing fasta
filename=paste(gene,'_phase.fas',sep='')
cat("",file=filename, sep="",append=FALSE)

for (i in 1:(nrow(dat2)/2)) {
cat('>',geneid[i,1],'_1','\n',sep='',file=filename,append=TRUE)
cat(dat2[(2*i-1),],'\n', sep='',file=filename,append=TRUE)
cat('>',geneid[i,1],'_2','\n',sep='',file=filename,append=TRUE)
cat(dat2[(2*i),],'\n',sep='',file=filename,append=TRUE)
}
}


##Trying mael1
####Run terminal inside r
genelist=c('arx.fas')
for (i in 1:length(genelist)) {
###Getting gene identifier
geneid<-sub(".fas","",genelist[i])
###Convert fasta to Phase format
FastatoPHASE(genelist[i])
###Running Fast-Phase
system(paste("./fastPHASE -o",geneid," -n ",geneid,"_unphased.inp",sep=""))
###Tidy-up the phase output
system(paste("tail -n +22 ",geneid,"_hapguess_switch.out |grep -v '#' | grep -v 'END'* > ",geneid,".phased",sep=""))
###Convert back from Phase to Fasta
PHASEtoFasta(paste(geneid,".phased",sep=""))
}







