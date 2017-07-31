###getting all log files in directory
logfiles <- list.files(pattern = "*.log.txt")
head(logfiles)

##Loop to normalize the age based on the age of root
for (i in 1:length(logfiles)) {
  gene <- gsub('_root.log.txt','',logfiles[i])
  initdata <- read.table(logfiles[i],skip = 2,header = TRUE)
  ###Burnin(1000 initial, discard initial 1000 rows) and only taking the age column 
  burndata <- initdata[-c(1:1000),grep('tmrca',colnames(initdata))]
  normdata <- t(apply(burndata,1,
                      function(x) {x/x[grep('obsr',names(x))]}))
  assign(paste("norm",gene,sep = '_'),normdata)
}


###Converting into long format suitable for ggplot
agetab <- data.frame(type='type',age=0)
standard <- ls(pattern='norm')
for (i in 1:length(standard)) {
  longdata <- gather(as.data.frame(get(standard[i])),type,age,grep('paralog',colnames(get(standard[i]))))
  agedata <- longdata[,grep('type|age',colnames(longdata))]
  agetab <- rbind(agetab,agedata)
  
}

###Tidying up the result

agetab <- agetab[-1,]
head(agetab)
gene <- levels(agetab$type)[c(2,3,5,6,8,9)]
agetab2 <- filter(agetab,type %in% gene) %>% droplevels()
levels(agetab2$type) <- c('Armi','Arx','Cuff','Mael','Tejas','Vret')
levels(agetab2$type)
write.table(agetab2,file='agedup_stdroot.txt',sep='\t')

###Plotting ggplot
plotage <- ggplot(agetab2, aes(x=age,fill=type)) + geom_density(alpha=0.5)+
  theme_bw()+labs(x='normalized age of duplication',y='posterior density')+scale_fill_aaas()

##Saving the plot
ggsave(plot=plotage,file='dmel_dobs_root.svg',device="svg")


