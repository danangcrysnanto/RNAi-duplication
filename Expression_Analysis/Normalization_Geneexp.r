relex <- dat[,7]/dat[,6]
dat2 <- cbind(dat,relex)

#Getting only data from Dmir
dat3 <- dat2[which(dat2[,2]=='Dmir'),]
head(dat3)
reldmir <- dat3[,7]/dat3[,6]
#Getting position of rpl32
posrpl <- seq(16,nrow(dat3),by=16)
#Getting value of rpl32
valrpl <- relex[posrpl]
stdmir <- rep(valrpl,each=16)
#Getting normalized by rpl32
normdmir <- reldmir/stdmir
head(normdmir,n=20)

#Getting only data from Dobs
datobs <- dat2[which(dat2[,2]=='Dobs'),]
head(datobs)
reldobs <- datobs[,7]/datobs[,6]
#Getting position of rpl32
posrplobs <- seq(10,nrow(datobs),by=10)
#Getting value of rpl32
valrplobs <- reldobs[posrplobs]
stdobs <- rep(valrplobs,each=10)
#Getting normalized by rpl32
normdobs <- reldobs/stdobs
head(normdobs,n=20)


#Getting only data from Dpse
datdpse <- dat2[which(dat2[,2]=='Dpse'),]
head(datdpse)
relddpse <- datdpse[,7]/datdpse[,6]
#Getting position of rpl32
posrpldpse <- seq(18,nrow(datdpse),by=18)
#Getting value of rpl32
valrpldpse <- relddpse[posrpldpse]
stddpse <- rep(valrpldpse,each=18)
#Getting normalized by rpl32
normddpse <- relddpse/stddpse
head(normddpse,n=20)


#Joined all together
normalized <- c(normdmir,normdobs,normddpse)

#Joined with original dataset
nigmaptable<-cbind(dat2,normalized)
head(nigmaptable,n=30)
colnames(nigmaptable) <- c('accession','species','organ','sex','gene','genelength','relex','normalized')
write.table(nigmaptable,sep='\t',file='nigmaptable.tsv')
getwd()

#check infinite
is.infinite
is.na





