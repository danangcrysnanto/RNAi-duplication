
##I used this script to create boxplot of expression divided between sex across genes in 3 different Drosophila species.
##I used ggplot2 for plotting

dat<-read.csv('bigtable.tsv',sep='\t')
head(dat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)
library(plotrix)
library(svglite)
###Getting Dpse and count accession number####
dpsedat <- unique(dat %>% filter(species=='Dpse') %>% select(accession,species,organ,sex))
nodataset <- nrow(unique(dpsedat))
###Dpse has 163 transcriptome datasets
write.table(dpsedat,file = 'dpsedat.tsv',sep='\t')

###Getting Dmir and count accession number####
dmirdat <- unique(dat %>% filter(species=='Dmir') %>% select(accession,species,organ,sex))
nodataset <- nrow(unique(dmirdat))
###Dpse has 163 transcriptome datasets
write.table(dmirdat,file = 'dmirdat.tsv',sep='\t')


###Getting Dmir and count accession number####
dobsdat <- unique(dat %>% filter(species=='Dobs') %>% select(accession,species,organ,sex))
nodataset <- nrow(unique(dobsdat))
###Dpse has 163 transcriptome datasets
write.table(dobsdat,file = 'dobsdat.tsv',sep='\t')


###Getting Dmir and count accession number####
dathdat <- unique(dat %>% filter(species=='Dath') %>% select(accession,species,organ,sex))
nodataset <- nrow(unique(dathdat))
###Dpse has 163 transcriptome datasets
write.table(dathdat,file = 'dathdat.tsv',sep='\t')


###Getting-dpse data
expdpse <- dat %>% filter( species=='Dpse'& gene != 'Dpse_rpl32') %>% droplevels()
head(expdpse)
length(levels(expdpse$gene))

###Renaming leves and reordering
levels(expdpse$gene)<-c('Armi','Arx-duplicate','Arx-ancestral','Cuff-ancestral','Cuff-duplicate','Mael-ancestral','Mael-duplicate','Tejas-ancestral','Tejas-duplicate','Vret-duplicate','Vret-ancestral','Ago2a1','Ago2a3','Ago2b','Ago2c','Ago2d','Ago2e')
expdpse$gene <- factor(expdpse$gene,levels=c('Ago2a1','Ago2a3','Ago2b','Ago2c','Ago2d','Ago2e','Arx-ancestral','Arx-duplicate','Cuff-ancestral','Cuff-duplicate','Armi','Mael-ancestral','Mael-duplicate','Tejas-ancestral','Tejas-duplicate','Vret-ancestral','Vret-duplicate'))
levels(expdpse$species)<-c('pseudoobscura')
length(levels(expdpse$gene))




####Ago Dpse boxplot####
agolevel <- levels(expdpse$gene)[1:6]
agodpse <- filter(expdpse, gene %in% agolevel) %>% droplevels()
plotab <- ggplot(agodpse,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() +
  scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position='none',axis.title.y=element_blank())
ggsave(plot=plotab,file='dpsegene1.svg',device="svg",height=2,width=7.5)
plotab
dat2 <- read.table(file='datjoin.tsv',sep='\t')

###Other genes, Armi,Arx, Cuff####
duplevel <- levels(expdpse$gene)[7:11]
dupdpse <- filter(expdpse, gene %in% duplevel) %>% droplevels()
plotab2 <- ggplot(dupdpse,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none")+labs(y=NULL,x=NULL)
ggsave(plot=plotab2,file='dupdpse1.svg',device="svg",height=2,width=7.5)
plotab2

###Other genes, vret,Mael, Tejas####
duplevel <- levels(expdpse$gene)[12:17]
dupdpse <- filter(expdpse, gene %in% duplevel) %>% droplevels()
plotab2 <- ggplot(dupdpse,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none")+labs(y=NULL,x=NULL)
ggsave(plot=plotab2,file='dupdpse2.svg',device="svg",height=2,width=7.5)
plotab2

head(expdpse)



###Getting-dmir data
expdmir <- dat2 %>% filter( species=='Dmir'& gene != 'rpl32Dmir') %>% droplevels()
head(expdmir)
levels(expdmir$gene)

###Renaming leves and reordering
levels(expdmir$gene)<-c('Ago2ab','Ago2c','Ago2d','Ago2e','Armi','Arx-ancestral','Arx-duplicate','Cuff-duplicate','Cuff-ancestral','Mael-duplicate','Mael-ancestral','Tejas-duplicate','Tejas-ancestral','Vret-duplicate','Vret-ancestral')
expdmir$gene <- factor(expdmir$gene,levels=c('Ago2ab','Ago2c','Ago2d','Ago2e','Arx-ancestral','Arx-duplicate','Cuff-ancestral','Cuff-duplicate','Armi','Mael-ancestral','Mael-duplicate','Tejas-ancestral','Tejas-duplicate','Vret-ancestral','Vret-duplicate'))
summary(expdmir$gene)
levels(expdmir$species)<-c('miranda')
length(levels(expdmir$gene))




####Ago Dmir boxplot####
agolevel <- levels(expdmir$gene)[1:4]
summary(expdmir$gene)
agodmir <- filter(expdmir, gene %in% agolevel) %>% droplevels()
plotab <- ggplot(agodmir,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position='none')+labs(x=NULL,y=NULL)
ggsave(plot=plotab,file='dmirargo.svg',device="svg",height=2,width=5)
plotab

#Dmirallplot########
#Renaminglevels
summary(expdmir$gene)
levels(expdmir$gene)<-c('Ago2ab','Ago2c','Ago2d','Ago2e','Armi','Arx1','Arx2','Cuff1','Cuff2','Mael1','Mael2','Tejas1','Tejas2','Vret1','Vret2')
ggplot(expdmir,aes(x=gene,y=log(normalized), fill = sex)) +
  geom_boxplot() +theme_bw() + scale_fill_npg()

plotexp <- ggplot(dat2,aes(x=sex,y=log(normalized)),fill=sex) +
  geom_boxplot(aes(fill=sex)) + facet_grid(.~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank())

###Other genes, Armi,Arx, Cuff####
duplevel <- levels(expdmir$gene)[5:9]
dupdpse <- filter(expdmir, gene %in% duplevel) %>% droplevels()
plotab2 <- ggplot(dupdpse,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none")+labs(y=NULL,x=NULL)
ggsave(plot=plotab2,file='dupdmir1.svg',device="svg",height=2,width=7.5)
plotab2



###Obscura


###Getting-obscura data
expdobs <- dat %>% filter( species=='Dobs'& gene != 'Dobs_rpl32') %>% droplevels()
head(expdobs)
levels(expdpse$gene)

###Renaming leves and reordering
levels(expdobs$gene)<-c('Ago2ab','Ago2e','Ago2f','Armi','Arx','Cuff','Mael','Tejas','Vret')
expdobs$gene <- factor(expdobs$gene,levels=c('Ago2ab','Ago2e','Ago2f','Arx','Cuff','Armi','Mael','Tejas','Vret'))
levels(expdobs$species)<-c('obscura')
length(levels(expdobs$gene))


####Ago Dpse boxplot####
agolevel <- levels(expdobs$gene)[1:3]
agodobs <- filter(expdobs, gene %in% agolevel) %>% droplevels()
plotab <- ggplot(agodobs,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() +
  scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position='none')+labs(x=NULL,y=NULL)
ggsave(plot=plotab,file='dobsago.svg',device="svg",height=2,width=3.75)
plotab


###Other genes, Armi,Arx, Cuff####
duplevel <- levels(expdobs$gene)[4:6]
dupdobs <- filter(expdobs, gene %in% duplevel) %>% droplevels()
plotab2 <- ggplot(dupdobs,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none")+labs(y=NULL,x=NULL)
ggsave(plot=plotab2,file='dupdobs1.svg',device="svg",height=2,width=3.75)
plotab2

###Other genes, vret,Mael, Tejas####
duplevel <- levels(expdobs$gene)[7:9]
dupdobs <- filter(expdobs, gene %in% duplevel) %>% droplevels()
plotab2 <- ggplot(dupdobs,aes(x=sex,y=log(normalized)),fill=sex) + geom_boxplot(aes(fill=sex)) + facet_grid(species~gene)+theme_bw() + scale_fill_npg()+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none")+labs(y=NULL,x=NULL)
ggsave(plot=plotab2,file='dupdobs2.svg',device="svg",height=2,width=3.75)
plotab2


View(dat2)

