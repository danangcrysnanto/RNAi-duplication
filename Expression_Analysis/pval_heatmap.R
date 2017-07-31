###I used this script to create heat map of the p-value from difference between two posterior distribution

###Open posterior distribution table 
totagd <- read.table(totagd,file='totagd.tsv')

###Getting levels and how many levels
lev <- levels(totagd$type)
lenlev <- length(levels(totagd$type))


###Set blank data frame to be filled with 
tabpval <- data.frame(matrix(NA,nrow=110,ncol=3))


###Taking pairwise posterior distribution and calculate p-value
itera <- 0
for (i in 1:lenlev){
  for (j in 1:lenlev){
    itera <- itera+1
    gene1 <- lev[i]
    gene2 <- lev[j]
    ###Getting posterior distribution for each gene and getting only the age column
    dist1 <- filter(totagd,totagd$type==gene1)$age
    dist2 <- filter(totagd,totagd$type==gene2)$age
    ###True or False for difference in distribution for more than 0
    diffdist <- (dist1-dist2) >0
    ###Taking sum and divided by number of observation
    x<-(sum(diffdist)/length(diffdist))
    pval<-2*min(x,(1-x))
    ###Filling the p-value table
    tabpval[itera,1]<-gene1
    tabpval[itera,2]<-gene2
    tabpval[itera,3]<-pval
  }
}
colnames(tabpval)<-c('gen1','gen2','pval')

write.table(tabpval,file='pvaldup.txt',sep='\t')


###Plot heatmap using geom_tile
plotval<-ggplot(tabpval, aes(x=gen1,y=gen2,fill=pval)) +
geom_tile(color='black') +
geom_text(aes(gen1, gen2, label = round(pval,3)), color = "white", size = 4)+theme_minimal() +
scale_fill_distiller(palette = "Spectral", trans = "reverse")

###View plot
plotval

###Save plot
ggsave(plot=plotval,file='pval.svg',device="svg")


