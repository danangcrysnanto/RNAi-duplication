###I used this script to do automatic mcmcGLMM analysis on Dpse expression dataset 


###Getting data from tabulation table
dat <- read.table('bigtable.tsv',header=TRUE)

###Getting-dpse data
expdpse <- dat %>% filter( species=='Dpse'& gene != 'Dpse_rpl32') %>% droplevels()
head(expdpse,n=100)
length(levels(expdpse$gene))


###Taking log of normalized dataset
lognorm <- log(expdpse$normalized)
lognorm[which(is.infinite(lognorm))]<--15
exp1 <- cbind(expdpse,lognorm)
hist(lognorm)
hist(expdpse$normalized)

###Filtering accesory gland
exp3 <- filter(exp1, organ !='accessory_glands') %>% droplevels()
levels(exp3$organ)

##Change ovaries into reproductive tissue
levels(exp3$organ)[9] <- 'reproductive_tissue'
levels(exp3$organ)[4]

###Change abdomen without gonad to abdomen
levels(exp3$organ)[4] <- 'abdomens'
levels(exp3$organ)

###Change testes to reproductive tissue
levels(exp3$organ)[10] <- 'reproductive_tissue'
levels(exp3$organ)

###Getting the gene name
gn1 <- sub('Dpse','',levels(exp3$gene))
gn2 <- sub('_','',gn1)
gn2
lev <- levels(exp3$gene)

###Doing mcmcglmm
cat('','dpse_MCMCglmm.txt',append=FALSE)
for (i in 1:length(lev)){
name<-lev[i]
cat(name,'\n',file='dpse_MCMCglmm.txt',append=TRUE)
selgen <- exp3 %>% filter(gene==lev[i]) %>% droplevels()
MCMCglmm(lognorm~sex, random=~organ,data=selgen)->fit
write.table(summary(fit)$solutions,append=TRUE,file='dpse_MCMCglmm.txt')
}
