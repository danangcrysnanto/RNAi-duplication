dpsedat <-expdmir
lgene <- as.character(levels(expdmir$gene))
dat7 <- dpsedat %>% group_by(sex,gene,organ)  %>% filter(gene %in% lgene[c(8)]) %>% droplevels()
dat8 <- summarise(dat7,av=mean(normalized),serror=std.error(normalized))
dat9 <- as.data.frame(dat8)
dat10 <- rbind(dat9,c('female',lgene[i],'testes',0,0))
dat11 <- rbind(dat10,c('male',lgene[i],'ovaries',0,0))
dat12 <- rbind(dat11,c('female',lgene[i],'accessory_glands',0,0))
dat12$av <- as.numeric(dat12$av)
dat12$serror<-as.numeric(dat12$serror)


explot <- ggplot(data=dat9, aes(x=organ, y=av,fill=sex)) +
  geom_bar(stat="identity",position='dodge',color='black')+
  theme_bw()+scale_fill_npg()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position='none')+
  facet_wrap(~gene,ncol=2)+
  theme(text= element_text(size=12), axis.text.x=element_text())+
  labs(y='normalized expression')+
  geom_errorbar(aes(ymin=av-serror,ymax=av+serror),width=0.25,position=position_dodge(0.9))


ggsave(explot,file='obs_1.svg',device='svg')
