load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')

cluster.markers.0.2_merged <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Res0.2_merged_cluster_DE.csv')
cluster.markers.0.2_merged <- cluster.markers.0.2_merged[p_val_adj <= 0.001,]
# # top 100 genes up in each cluster:
clusters_up <- setorder(setDT(cluster.markers.0.2_merged), -avg_logFC)[, head(.SD, 100), keyby = cluster]
# 


#recheck the findallmarkers


cluster.markers.0.2 <- FindAllMarkers(is019, logfc.threshold = 0.25, return.thresh=0.001)
cluster.markers.0.2
fwrite(cluster.markers.0.2, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Res0.2_merged_cluster_DE_new_names.csv')

clusters_up <- setorder(setDT(cluster.markers.0.2), -avg_logFC)[, head(.SD, 100), keyby = cluster]

clusters_up[1:100 & cluster=='CEPT1' & avg_logFC >0,c(gene)]

##
glycolysis <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=2)
erstress <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=1)

DotPlot(is019, features=glycolysis$Gene)+ coord_flip()
DotPlot(is019, features=erstress$Gene)+ coord_flip()


avgs.glyc <- AverageExpression(is019, features=glycolysis$Gene)$RNA
avgs.glyc
avgs.glyc$GeneId <- row.names(avgs.glyc)
avgs.glyc <- as.data.table(avgs.glyc)

avgs.erstress <- AverageExpression(is019, features=erstress$Gene[erstress$Gene %in% row.names(is019@assays$RNA@counts)])$RNA
avgs.erstress
avgs.erstress$GeneId <- row.names(avgs.erstress)
avgs.erstress <- as.data.table(avgs.erstress)

avgs.glyc <- melt(avgs.glyc, id.vars='GeneId')
avgs.erstress <- melt(avgs.erstress, id.vars='GeneId')
names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.erstress)[2:3] <- c('Cluster','Expression')

avgs.glyc <- avgs.glyc[Cluster != 'CEPT3' & Cluster != 'Y7',]
avgs.erstress <- avgs.erstress[Cluster != 'CEPT3' & Cluster != 'Y7',]


glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.erstress[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

library(cowplot)

plot_grid(glycolysis.plot, erstress.plot, ncol=2)

#

is019.neuronal <- subset(is019, idents=c('CEPT1',
                                         'CEPT2',
                                         'Y4',
                                         'Y5',
                                         'Y6'))

glyc.DE <- as.data.table(FindAllMarkers(is019.neuronal, only.pos=T, features=glycolysis$Gene, logfc.threshold = 0, min.pct = 0))
glyc.DE<- glyc.DE[p_val_adj < 0.001 & avg_logFC >= 0.25,]

erstress.DE <- as.data.table(FindAllMarkers(is019.neuronal, only.pos=T, features=avgs.erstress$GeneId, logfc.threshold = 0.25, min.pct = 0))
erstress.DE<- erstress.DE[p_val_adj < 0.001 & avg_logFC >= 0.25,]

fwrite(erstress.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/ER_stress_DE.csv')
fwrite(glyc.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Glycolysis_DE.csv')
erstress.DE <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/ER_stress_DE_clusters.csv')
glyc.DE <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Glycolysis_DE_clusters.csv')

# glycolysis wilcox test cluster by cluster----

avgs.glyc <- AverageExpression(is019.neuronal, features=unique(glyc.DE$gene) )$RNA
avgs.glyc
avgs.glyc$GeneId <- row.names(avgs.glyc)
avgs.glyc <- as.data.table(avgs.glyc)

avgs.erstress <- AverageExpression(is019.neuronal, features=unique(erstress.DE$gene))$RNA
avgs.erstress
avgs.erstress$GeneId <- row.names(avgs.erstress)
avgs.erstress <- as.data.table(avgs.erstress)

avgs.glyc <- melt(avgs.glyc, id.vars='GeneId')
avgs.erstress <- melt(avgs.erstress, id.vars='GeneId')
names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.erstress)[2:3] <- c('Cluster','Expression')

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y4'),Expression]),
            alternative='less')
# p-value = 0.1164

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y5'),Expression]),
            alternative='less')
# p-value = 0.2432

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y6'),Expression]),
            alternative='less')
# p-value = 0.1946

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y4'),Expression]),
            alternative='less')
# p-value = 0.6435

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y5'),Expression]),
            alternative='less')
# p-value = 0.6587

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y6'),Expression]),
            alternative='less')
# p-value = 0.6882

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y4'),Expression]),
            alternative='two.sided')
# p-value = 0.2328

#***
wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y5'),Expression]),
            alternative='two.sided')
# p-value = 0.4864

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y6'),Expression]),
            alternative='two.sided')
# p-value = 0.3892

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y4'),Expression]),
            alternative='two.sided')
# p-value = 0.7437

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y5'),Expression]),
            alternative='two.sided')
# p-value = 0.713

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y6'),Expression]),
            alternative='two.sided')
# p-value = 0.6529

# er stress wilcoxon test by cluster----
wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y4'),Expression]),
            alternative='less')
# p-value = 0.005293

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y5'),Expression]),
            alternative='less')
# p-value =0.01638

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y6'),Expression]),
            alternative='less')
# p-value = 0.0005847

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y4'),Expression]),
            alternative='less')
# p-value = 0.3715

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y5'),Expression]),
            alternative='less')
# p-value = 0.3565

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y6'),Expression]),
            alternative='less')
# p-value = 0.07834

#***
wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y4'),Expression]),
            alternative='two.sided')
# p-value = 0.01059


wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y5'),Expression]),
            alternative='two.sided')
# p-value = 0.03275

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT1'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y6'),Expression]),
            alternative='two.sided')
# p-value = 0.001169

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y4'),Expression]),
            alternative='two.sided')
# p-value = 0.743

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y5'),Expression]),
            alternative='two.sided')
# p-value = 0.7131

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT2'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y6'),Expression]),
            alternative='two.sided')
# p-value = 0.1567


# now plot only the genes in these two lists which are DE and at least 0.25 LFC------

glycolysis <- glyc.DE$gene
erstress <- erstress.DE$gene

avgs.glyc <- AverageExpression(is019.neuronal, features=unique(glycolysis) )$RNA
avgs.glyc
avgs.glyc$GeneId <- row.names(avgs.glyc)
avgs.glyc <- as.data.table(avgs.glyc)

avgs.erstress <- AverageExpression(is019.neuronal, features=unique(erstress) )$RNA
avgs.erstress
avgs.erstress$GeneId <- row.names(avgs.erstress)
avgs.erstress <- as.data.table(avgs.erstress)

avgs.glyc <- melt(avgs.glyc, id.vars='GeneId')
avgs.erstress <- melt(avgs.erstress, id.vars='GeneId')
names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.erstress)[2:3] <- c('Cluster','Expression')

fwrite(avgs.erstress, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/ER_stress_avgs.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Glycolysis_avgs.csv')

glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.erstress[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

library(cowplot)

plot_grid(glycolysis.plot, erstress.plot, ncol=2)
#

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT1','CEPT2'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y4','Y5', 'Y6'),Expression]))

# p = 0.679


# try test overall of CEPT vs Y----
Idents(is019.neuronal) <- 'sample'
glyc.DE.treatment <- as.data.table(FindAllMarkers(is019.neuronal, only.pos=T, features=glycolysis, logfc.threshold = 0, min.pct = 0))
glyc.DE.treatment<- glyc.DE.treatment[p_val_adj < 0.001 & avg_logFC >= 0.25,]

erstress.DE.treatment <- as.data.table(FindAllMarkers(is019.neuronal, only.pos=T, features=erstress, logfc.threshold = 0.25, min.pct = 0))
erstress.DE.treatment <- erstress.DE.treatment[p_val_adj < 0.001 & avg_logFC >= 0.25,]

fwrite(glyc.DE.treatment, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Glycolysis_DE_treatment.csv')
fwrite(erstress.DE.treatment, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/ER_stress_DE_treatment.csv')

# plot avgs by treatment for DE-----
glycolysis <- glyc.DE.treatment$gene
erstress <- erstress.DE.treatment$gene

avgs.glyc <- AverageExpression(is019.neuronal, features=unique(glycolysis) )$RNA
avgs.glyc
avgs.glyc$GeneId <- row.names(avgs.glyc)
avgs.glyc <- as.data.table(avgs.glyc)

avgs.erstress <- AverageExpression(is019.neuronal, features=unique(erstress) )$RNA
avgs.erstress
avgs.erstress$GeneId <- row.names(avgs.erstress)
avgs.erstress <- as.data.table(avgs.erstress)

avgs.glyc <- melt(avgs.glyc, id.vars='GeneId')
avgs.erstress <- melt(avgs.erstress, id.vars='GeneId')
names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.erstress)[2:3] <- c('Cluster','Expression')

fwrite(avgs.erstress, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/ER_stress_avgs_treatment.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Glycolysis_avgs_treatment.csv')


glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.erstress[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

library(cowplot)

plot_grid(glycolysis.plot, erstress.plot, ncol=2)

# wilcox test by treatment-----

wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y'),Expression]),
            alternative="less")

# p-value = 0.02857

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y'),Expression]),
            alternative="less")
# p-value = 0.01876


# two sided:
wilcox.test(x= c(avgs.glyc[Cluster %in% c('CEPT'),Expression]),
            y= c(avgs.glyc[Cluster %in% c('Y'),Expression]),
            alternative="two.sided")

# p-value = 0.05714

wilcox.test(x= c(avgs.erstress[Cluster %in% c('CEPT'),Expression]),
            y= c(avgs.erstress[Cluster %in% c('Y'),Expression]),
            alternative="two.sided")
# p-value = 0.03752


# ISB017 tests--------------------
load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017.DDS.CM.RData")
unique(dds@colData$condition)

res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/3DCEPT_vs_Y/DE.3D_CEPT.vs.3D_Y.csv')
res
glycolysis <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=2)
erstress <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=1)

res.glyc <- res[GeneId %in% glycolysis$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #3
res.er <- res[GeneId %in% erstress$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #13

fwrite(res.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/ISB017_Glycolysis_DE_treatment.csv')
fwrite(res.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/ISB017_ER_stress_DE_treatment.csv')


res.glyc <- res.glyc[,c('GeneId','3D_CEPT','3D_Y')]
avgs.glyc <- melt(res.glyc, id.vars = 'GeneId')

res.er <- res.er[,c('GeneId','3D_CEPT','3D_Y')]
avgs.er <- melt(res.er, id.vars = 'GeneId')

names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.er)[2:3] <- c('Cluster','Expression')

fwrite(avgs.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/ISB017_ER_stress_avgs_treatment.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/ISB017_Glycolysis_avgs_treatment.csv')



glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.er[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

plot_grid(glycolysis.plot, erstress.plot, ncol=2)


wilcox.test(x=avgs.er[Cluster=='3D_CEPT',Expression],
            y=avgs.er[Cluster=='3D_Y',Expression],
            alternative='less')
#p-value = 0.5
wilcox.test(x=avgs.er[Cluster=='3D_CEPT',Expression],
            y=avgs.er[Cluster=='3D_Y',Expression],
            alternative='two.sided')
#p-value = 0.1


wilcox.test(x=avgs.glyc[Cluster=='3D_CEPT',Expression],
            y=avgs.glyc[Cluster=='3D_Y',Expression],
            alternative='less')
#p-value = 0.35
wilcox.test(x=avgs.glyc[Cluster=='3D_CEPT',Expression],
            y=avgs.glyc[Cluster=='3D_Y',Expression],
            alternative='two.sided')
#p-value = 0.7

## using all counts without averaging------

norm.counts  <- as.data.frame(counts(dds, normalized=T))
norm.counts$GeneId <- row.names(norm.counts)
norm.counts <- as.data.table(norm.counts)
norm.counts <- norm.counts[,c(45,8:44)]
names(norm.counts)

all.glyc <- norm.counts[norm.counts$GeneId %in% res.glyc$GeneId,]
wilcox.test(x=unlist(all.glyc[,2:20], use.names=F),
            y=unlist(all.glyc[,21:38], use.names=F),
            alternative='less')
# p-value = 0.01287

all.glyc.for.plot <- melt(all.glyc, id.vars = 'GeneId')
all.glyc.for.plot[,Treatment:=tstrsplit(variable, '_')[1]]

ggplot(data=all.glyc.for.plot[value >0,], aes(x=Treatment, y=log(value, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))


# er
all.erstress <- norm.counts[norm.counts$GeneId %in% res.er$GeneId,]
wilcox.test(x=unlist(all.erstress[,2:20], use.names=F),
            y=unlist(all.erstress[,21:38], use.names=F),
            alternative='less')
#p-value = 0.5508

all.erstress.for.plot <- melt(all.erstress, id.vars = 'GeneId')
all.erstress.for.plot[,Treatment:=tstrsplit(variable, '_')[1]]

ggplot(data=all.erstress.for.plot[value >0,], aes(x=Treatment, y=log(value, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))


# GM23279-----
glycolysis <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=2)
erstress <- readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/ERstress and glycolysis gene sets.xlsx', sheet=1)
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/DE.CEPT_GM23279.vs.Y_GM23279.csv')
res.glyc <- res[GeneId %in% glycolysis$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #7
res.er <- res[GeneId %in% erstress$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #57

fwrite(res.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM23279/ISB017_GM23279_Glycolysis_DE_treatment.csv')
fwrite(res.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM23279/ISB017_GM23279_ER_stress_DE_treatment.csv')


res.glyc <- res.glyc[,c('GeneId','CEPT_GM23279','Y_GM23279')]
avgs.glyc <- melt(res.glyc, id.vars = 'GeneId')

res.er <- res.er[,c('GeneId','CEPT_GM23279','Y_GM23279')]
avgs.er <- melt(res.er, id.vars = 'GeneId')

names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.er)[2:3] <- c('Cluster','Expression')

fwrite(avgs.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM23279/ISB017_GM23279_ER_stress_avgs_treatment.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM23279/ISB017_GM23279_Glycolysis_avgs_treatment.csv')

wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM23279',Expression],
            y=avgs.glyc[Cluster=='Y_GM23279',Expression],
            alternative='less')
#p-value = 0.6448
wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM23279',Expression],
            y=avgs.glyc[Cluster=='Y_GM23279',Expression],
            alternative='two.sided')
#p-value =  0.8048


wilcox.test(x=avgs.er[Cluster=='CEPT_GM23279',Expression],
            y=avgs.er[Cluster=='Y_GM23279',Expression],
            alternative='less')
#p-value = 0.5451
wilcox.test(x=avgs.er[Cluster=='CEPT_GM23279',Expression],
            y=avgs.er[Cluster=='Y_GM23279',Expression],
            alternative='two.sided')
#p-value =  0.9143


# GM25256-----
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/DE.CEPT_GMG25256.vs.Y_GM25256.csv')
res.glyc <- res[GeneId %in% glycolysis$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #11
res.er <- res[GeneId %in% erstress$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #55

fwrite(res.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM25256/ISB017_GM25256_Glycolysis_DE_treatment.csv')
fwrite(res.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM25256/ISB017_GM25256_ER_stress_DE_treatment.csv')


res.glyc <- res.glyc[,c('GeneId','CEPT_GMG25256','Y_GM25256')]
names(res.glyc)[2] <- 'CEPT_GM25256'
avgs.glyc <- melt(res.glyc, id.vars = 'GeneId')

res.er <- res.er[,c('GeneId','CEPT_GMG25256','Y_GM25256')]
names(res.er)[2] <- 'CEPT_GM25256'
avgs.er <- melt(res.er, id.vars = 'GeneId')

names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.er)[2:3] <- c('Cluster','Expression')

fwrite(avgs.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM25256/ISB017_GM25256_ER_stress_avgs_treatment.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/GM25256/ISB017_GM25256_Glycolysis_avgs_treatment.csv')



glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2',
       x='Treatment and cell line')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.er[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2',
       x='Treatment and cell line')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

plot_grid(glycolysis.plot, erstress.plot, ncol=2)

wilcox.test(x=avgs.er[Cluster=='CEPT_GM25256',Expression],
            y=avgs.er[Cluster=='Y_GM25256',Expression],
            alternative='less')
#p-value = 0.7936
wilcox.test(x=avgs.er[Cluster=='CEPT_GM25256',Expression],
            y=avgs.er[Cluster=='Y_GM25256',Expression],
            alternative='two.sided')
#p-value =  0.4162


wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM25256',Expression],
            y=avgs.glyc[Cluster=='Y_GM25256',Expression],
            alternative='less')
#p-value = 0.4744
wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM25256',Expression],
            y=avgs.glyc[Cluster=='Y_GM25256',Expression],
            alternative='two.sided')
#p-value = 0.9487

# Lonza-----
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/DE.CEPT_Lonza.vs.Y_Lonza.csv')
res.glyc <- res[GeneId %in% glycolysis$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #2 #11
res.er <- res[GeneId %in% erstress$Gene & padj <= 0.001 & abs(log2FoldChange) >= 0.25,] #0

fwrite(res.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/Lonza/ISB017_Lonza_Glycolysis_DE_treatment.csv')
fwrite(res.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/Lonza/ISB017_Lonza_ER_stress_DE_treatment.csv')

res.glyc <- res.glyc[,c('GeneId','CEPT_Lonza','Y_Lonza')]
avgs.glyc <- melt(res.glyc, id.vars = 'GeneId')

res.er <- res.er[,c('GeneId','CEPT_Lonza','Y_Lonza')]
avgs.er <- melt(res.er, id.vars = 'GeneId')

names(avgs.glyc)[2:3] <- c('Cluster','Expression')
names(avgs.er)[2:3] <- c('Cluster','Expression')

fwrite(avgs.er, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/Lonza/ISB017_Lonza_ER_stress_avgs_treatment.csv')
fwrite(avgs.glyc, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Glycolysis_ER_stress/Lonza/ISB017_Lonza_Glycolysis_avgs_treatment.csv')



glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2',
       x='Treatment and cell line')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

glycolysis.plot


wilcox.test(x=avgs.glyc[Cluster=='CEPT_Lonza',Expression],
            y=avgs.glyc[Cluster=='Y_Lonza',Expression],
            alternative='less')
#p-value =  0.3333
wilcox.test(x=avgs.glyc[Cluster=='CEPT_Lonza',Expression],
            y=avgs.glyc[Cluster=='Y_Lonza',Expression],
            alternative='two.sided')
#p-value = 0.6667

glycolysis.plot <- ggplot(data=avgs.glyc[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='Glycolysis',
       y='Expression, log2',
       x='Treatment and cell line')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

erstress.plot <- ggplot(data=avgs.er[Expression >0,], aes(x=Cluster, y=log(Expression, base=2) ))+
  geom_boxplot()+
  labs(title='ER stress',
       y='Expression, log2',
       x='Treatment and cell line')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),
        title=element_text(size=15))

plot_grid(glycolysis.plot, erstress.plot, ncol=2)

wilcox.test(x=avgs.er[Cluster=='CEPT_GM23279',Expression],
            y=avgs.er[Cluster=='Y_GM23279',Expression],
            alternative='less')
#p-value = 0.5451
wilcox.test(x=avgs.er[Cluster=='CEPT_GM23279',Expression],
            y=avgs.er[Cluster=='Y_GM23279',Expression],
            alternative='two.sided')
#p-value = 0.9143


wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM23279',Expression],
            y=avgs.glyc[Cluster=='Y_GM23279',Expression],
            alternative='less')
#p-value = 0.6448
wilcox.test(x=avgs.glyc[Cluster=='CEPT_GM23279',Expression],
            y=avgs.glyc[Cluster=='Y_GM23279',Expression],
            alternative='two.sided')
#p-value = 0.8048