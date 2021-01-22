
d <- fread('data/Fulldata_communityanalysis.txt')


d.lm <- lm(Yield_Average~ sp1 + sp2st1 + sp2st2 + sp2st3 + sp3 + sp4 + sp5 + sp6 + sp7 + sp8, data=d)

d.lm <- lm(Yield_Average~ sp1, data=d)

nms <- c('sp1', 'sp2st1', 'sp2st2', 'sp2st3', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8')
w_yield <- c(with(d, wilcox.test(Yield_Average~sp1, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp2st1, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp2st2, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp2st3, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp3, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp4, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp5, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp6, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp7, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yield_Average~sp8, alternative = "two.sided"))$p.value)
names(w_yield) <- nms

w_yield <- setDT(as.list(w_yield))
fwrite(w_yield, 'w_yield.csv')

w_yeast <- c(with(d, wilcox.test(Yeast_Average~sp1, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp2st1, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp2st2, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp2st3, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp3, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp4, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp5, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp6, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp7, alternative = "two.sided"))$p.value,
             with(d, wilcox.test(Yeast_Average~sp8, alternative = "two.sided"))$p.value)
names(w_yeast) <- nms

w_yeast <- setDT(as.list(w_yeast))
fwrite(w_yeast, 'w_yeast.csv')



with(d, wilcox.test(Yeast_Average~sp1, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp2st1, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp2st2, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp2st3, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp3, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp4, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp5, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp6, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp7, alternative = "two.sided"))$p.value,
with(d, wilcox.test(Yeast_Average~sp8, alternative = "two.sided"))$p.value,

##
pdf('results/Amylovorans_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp1, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp1))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/Amylovorans_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp1, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp1))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()


##
pdf('results/sp2st1_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st1, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp2st1))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp2st1_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st1, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp2st1))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp2st2_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st2, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp2st2))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp2st2_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st2, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp2st2))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp2st3_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st3, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp2st3))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp2st3_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp2st3, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp2st3))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp3_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp3, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp3))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp3_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp3, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp3))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp4_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp4, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp4))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp4_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp4, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp4))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp5_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp5, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp5))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp5_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp5, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp5))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp6_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp6, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp6))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp6_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp6, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp6))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()



##
pdf('results/sp7_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp7, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp7))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp7_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp7, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp7))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()


##
pdf('results/sp8_yield.pdf', 3, 3)
ggplot(d, aes(x = factor(sp8, labels = c('Absent', 'Present')),
              y = Yield_Average,
              fill = factor(sp8))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Ethanol Yield (%)')
dev.off()

pdf('results/sp8_yeastCount.pdf', 3, 3)
ggplot(d, aes(x = factor(sp8, labels = c('Absent', 'Present')),
              y = Yeast_Average,
              fill = factor(sp8))) +
  geom_dotplot(binaxis = "y", stackdir = "center", method='histodot') +
  theme(panel.grid = element_blank(),
        legend.position='none') +
  labs(x = 'L. amylovorus', y='Yeast cell counts / ml')
dev.off()


