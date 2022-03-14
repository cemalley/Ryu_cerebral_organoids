library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)

load("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData")





exp <- FetchData(is019, c('SLC17A7', 'GRIN2B','GABBR1', 'GABBR2', 'GAD2', 'GAD1', 'KCNJ6', 'NR4A2', 'TH'), slot = 'data')
exp$barcode <- row.names(exp)
exp <- as.data.table(exp)
all(exp$barcode == row.names(is019@meta.data))

meta.mini <- is019@meta.data
meta.mini <- subset(meta.mini, select=c('clusters_refactor'))
meta.mini$barcode <- row.names(meta.mini)
meta.mini <- as.data.table(meta.mini)
meta.mini

# plots----
exp <- merge(exp, meta.mini, by='barcode')

exp.melt <- melt(exp, id.vars=c('barcode','clusters_refactor'))

exp.melt <- exp.melt[order(clusters_refactor)]
perc.expr <- exp.melt %>% as_tibble() %>% group_by(clusters_refactor, variable) %>% dplyr::summarize(n_cells = n(),
                                                                                                n_expressed = sum(value > 0),
                                                                                                p_expressed = n_expressed / n_cells) %>% as.data.table()

perc.expr

exp.final <- merge(exp.melt, perc.expr, by=c('clusters_refactor','variable'))
exp.final

range(exp.final$p_expressed)

exp.final[,new_clusters := factor(clusters_refactor, levels=c('CEPT1','CEPT2', 'CEPT3','Y4','Y5','Y6','Y7'))]
length(unique(exp.final$barcode))


genelist <- data.table(GeneId=c('SLC17A7', 'GRIN2B','GABBR1', 'GABBR2', 'GAD2', 'GAD1', 'KCNJ6', 'NR4A2', 'TH'), Category=c(rep('Glutamatergic',2),rep('Gabaergic',4),rep('Dopaminergic',3)))
genelist

exp.final <- merge(exp.final, genelist, by.x='variable', by.y='GeneId')
exp.final

# fancy dot plot 02/18/20-----
install.packages('ggnewscale')
library(ggnewscale)

is.data.table(exp.final)
exp.final[,variable_reorder:= factor(variable, levels=c('SLC17A7','GRIN2B','GABBR1','GABBR2','GAD1','GAD2','KCNJ6','NR4A2','TH'))]
exp.final[,GeneID := variable_reorder]
exp.final$GeneID <- gsub('SLC17A7','C2_SLC17A7', exp.final$GeneID)
exp.final$GeneID <- gsub('GRIN2B','C1_GRIN2B', exp.final$GeneID)
exp.final$GeneID <- gsub('GABBR1','B4_GABBR1', exp.final$GeneID)
exp.final$GeneID <- gsub('GABBR2','B3_GABBR2', exp.final$GeneID)
exp.final$GeneID <- gsub('GAD1','B2_GAD1', exp.final$GeneID)
exp.final$GeneID <- gsub('GAD2','B1_GAD2', exp.final$GeneID)
exp.final$GeneID <- gsub('KCNJ6','A3_KCNJ6', exp.final$GeneID)
exp.final$GeneID <- gsub('NR4A2','A2_NR4A2', exp.final$GeneID)
exp.final$GeneID <- gsub('TH','A1_TH', exp.final$GeneID)

ggplot(exp.final[Category=='Glutamatergic',],
       aes(y=GeneID, x=new_clusters, size=p_expressed)) + 
  scale_size_continuous(range = c(1,5))+
  geom_point(aes(color=value))+
  scale_color_gradient2("Score1", low='black',mid='lightgrey', high='blue')+
  
  new_scale('color')+
  geom_point(data=exp.final[Category=='Gabaergic',], aes(color=value))+
  scale_color_gradient2('Score2', low='black', mid='lightgrey', high='green')+
  
  new_scale('color')+
  geom_point(data=exp.final[Category=='Dopaminergic',], aes(color=value))+
  scale_color_gradient2('Score3', low='black', mid='lightgrey', high='red')+
  
  theme_bw()+theme(axis.text = element_text(size=10), axis.title = element_text(size=14), panel.grid = element_blank())+
  labs(x='Cluster', y='Gene', size='Percent expressed')

#glutamatergic
DotPlot(is019, features=c(unique(exp.final$variable))) +
  coord_cartesian(xlim=c(2,5))+ scale_color_gradient2(low='grey', mid='red', high='darkred')


DotPlot(is019, features=c(unique(exp.final$variable)))

#gabaergic
DotPlot(is019, features=c(unique(exp.final$variable))) +
  coord_cartesian(xlim=c(6:9))+ scale_color_gradient2(low='grey', mid='lightgreen', high='darkgreen')

#dopaminergic
DotPlot(is019, features=c(unique(exp.final$variable))) +
  coord_cartesian(xlim=c(1,3:4))+ scale_color_gradient2(low='grey', mid='lightblue', high='darkblue')


DotPlot(is019, features=c('ERO1A','DDIT3','BCAP31','CASP3'
,'BAX','ATF4')) +coord_flip()+labs(x='Cluster',y='Gene')

DotPlot(is019, features=c('ERO1A','DDIT3','BCAP31','CASP3'
                          ,'BAX','ATF4')) +coord_flip()+labs(x='Cluster',y='Gene')
