
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
