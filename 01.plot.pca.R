setwd("/Users/wuz6/Documents/Project/18.FGFR3-TACC3/01.cluster.revision1/")
rm(list=ls())
library(ggplot2)
library(umap)
library(ggrepel)

pair.col = RColorBrewer::brewer.pal(12, "Paired")
#
meta = data.frame(readxl::read_xlsx("../MS/v2/TableS1-2.xlsx", skip = 1), stringsAsFactors = F)
#meta = read.csv("../00.data/FGFR3-TACC3_cases.csv", head = T, stringsAsFactors = F)
head(meta)
meta$subtype = sapply(meta$Subclass..calibrated.score., function(x) unlist(strsplit(x, "[;:]"))[1])
meta$subtype.score = as.numeric(sapply(meta$Subclass..calibrated.score., function(x) unlist(strsplit(x, split = "[;:]"))[2]))
meta$subtype[which(meta$subtype.score <0.5)] = "No match"
table(meta$Tumor)
table(meta$subtype)

meta$subtype[meta$Tumor == "ANA PA"] = "ANA-PA"
meta$subtype[meta$Tumor == "DMG, K27"] = "DMG-K27"
meta$subtype[meta$Tumor == "GBM, G34"] = "GBM-G34"
meta$subtype[meta$Tumor == "LGG, GG"] = "Ganglioglioma"
meta$subtype[grep("CONTR", meta$subtype)] = "No match"
table(meta$subtype)

set.seed(123)

x <- meta$SampleID[meta$FGFR3.TACC3.fusion == "Y" & meta$Final_class == "GBM"]
#write.csv(x, file = "tmp.csv")

table(meta$FGFR3.TACC3.fusion, meta$Source)
table(meta$subtype, meta$FGFR3.TACC3.fusion)
table(meta$subtype, meta$Source)

nrow(meta)
table(meta$Final_class)

#RColorBrewer::display.brewer.all()
subtype.col =  RColorBrewer::brewer.pal(8, "Dark2")[c(2,3,4,5,6,7,8)]
subtype.col = c(pair.col[c(10, 1,2,9)], subtype.col, 'red')
names(subtype.col) = c("Ganglioglioma", "ANA-PA", "DMG-K27", "GBM-G34",
                       "GBM, MES", "GBM, MID", "GBM, MYCN", "GBM, RTK I", "GBM, RTK II", 
                       "GBM, RTK III", "No match", "GBM-F3T3-O")

meta$Mean.global.methylation = as.numeric(meta$Mean.global.methylation)
ggplot(meta, aes(x = reorder(Final_class,Mean.global.methylation), y = Mean.global.methylation)) +
  geom_boxplot()+ theme_bw()

aggregate(meta$Mean.global.methylation, list(meta$Final_class), mean)

wilcox.test(meta$Mean.global.methylation[which(meta$Final_class == "GBM-F3T3-O")], 
            meta$Mean.global.methylation[which(meta$Final_class == "GBM")], alternative = "greater")

#CDKN2A/B homo loss
fisher.test(matrix(c(84, 229-84, 0, 13), nrow = 2))

#oligo-area
fisher.test(matrix(c(8, 12-8, 3, 32-3), nrow = 2), alternative = "greater")
#calcifications
fisher.test(matrix(c(4, 12-4,6, 32-6), nrow = 2))
#endocrinoid vascular
fisher.test(matrix(c(5, 12-5,7, 32-7), , nrow = 2))
#Microvascular Proliferation (Y/N)
fisher.test(matrix(c(5, 12-5,20, 32-20), nrow = 2))
#necrosis
fisher.test(matrix(c(1, 12-1,23, 32-20), nrow = 2), alternative = "less")

#0.002744
#load("01.F3T3.preprocessQuantile.32k.pc.tsne.rda") #similar to preprocessFunnorm
load("01.F3T3.preprocessFunnorm.32k.pc.tsne.rda") #better than 20k, best
colnames(tsne) = c("tsne1", "tsne2")


N_non_trivial_pc #24 functionnorm, 29 processQuantile
custom.config = umap.defaults
custom.config$n_neighbors =10
#custom.config$random_state = 123
ref.umap = umap(betas.pca$x[,1:N_non_trivial_pc], config = custom.config)

um = ref.umap$layout; colnames(um) = c("umap1", "umap2")

rownames(betas.pca$x)[!rownames(betas.pca$x) %in% meta$SampleID]
meta$SampleID[!meta$SampleID %in% rownames(betas.pca$x)]
rownames(betas.pca$x)[!rownames(betas.pca$x) %in% meta$SampleID]
#meta = meta[meta$SampleID != "V197",]
meta = meta[match(rownames(tsne), meta$SampleID),]

identical(meta$SampleID, rownames(betas.pca$x))
identical(meta$SampleID, rownames(tsne))
identical(meta$SampleID, rownames(um))

res = cbind(meta, betas.pca$x[,1:4], tsne, um)
pc.var.explained = round(100*betas.pca$sdev^2/sum(betas.pca$sdev^2), 1)
#write.csv(res, file = "01.umap.csv", row.names = F)
res$Final_class[is.na(res$Final_class)] <- "new"
res$FGFR3.TACC3.fusion[is.na(res$FGFR3.TACC3.fusion)] <- "NA"

cancer.col <- c("#6A3D9A", "#CAB2D6", "#66A61E", "red", "gray")
names(cancer.col) <- c("Ganglioglioma", "GBM-G34", "GBM", "GBM-F3T3-O", "reference")

p4 = ggplot(res, aes(x = umap1, y = umap2, shape = FGFR3.TACC3.fusion, color = Final_class,label = SampleID)) + 
  scale_color_manual(values = cancer.col)+
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #geom_text_repel(data = res[res$Final_class == "GBM-F3T3-O",], max.overlaps = 30, show.legend = F)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())
p4


p3 = ggplot(res, aes(x = umap1, y = umap2, shape = FGFR3.TACC3.fusion, color = subtype)) + 
  scale_color_manual(values = subtype.col[names(subtype.col) %in% res$subtype])+
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #scale_color_brewer(palette = "Dark2")+
  #geom_text_repel(inherit.aes = F, data = gg.res, aes(x = umap1, y = umap2, label = SampleID),max.overlaps = 30)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1,panel.grid = element_blank(), axis.ticks = element_blank())
p3

pdf("01.F3T3.preprocessFunnorm.32k.umap.pdf", width = 10, height = 5, useDingbats = F)
#pdf("01.F3T3.processIllumina.32k.umap.pdf", width = 10, height = 5, useDingbats = F)
#pdf("01.F3T3.processQuantile.32k.umap.pdf", width = 10, height = 5, useDingbats = F)
cowplot::plot_grid(p3, p4)#
dev.off()

 #save(ref.umap, file = "F3T3.refUmap.rda")
##

p1 = ggplot(res, aes(x = PC1, y = PC2, shape = FGFR3.TACC3.fusion, color = subtype)) + 
  scale_color_manual(values = subtype.col) +
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #scale_color_brewer(palette = "Dark2")+
  labs(x = paste0("PC1 (", round(pc.var.explained[1], 1), " %)"), y = paste0("PC2 (", pc.var.explained[2], "%)"))+
  #geom_text_repel(inherit.aes = F, data = gg.res, aes(x = PC1, y = PC2, label = SampleID), max.overlaps = 20)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())
p1


p2 = ggplot(res, aes(x = tsne1, y = tsne2, shape = FGFR3.TACC3.fusion, color = subtype)) + 
  scale_color_manual(values = subtype.col)+
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #scale_color_brewer(palette = "Dark2")+
  #geom_text_repel(inherit.aes = F, data = gg.res, aes(x = tsne1, y = tsne2, label = SampleID),max.overlaps = 30)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())
p2




p6 = ggplot(res, aes(x = PC1, y = PC2, shape = FGFR3.TACC3.fusion, color = Final_class)) + 
  #scale_color_manual(values = list("Ganglioglioma" = "#6A3D9A", "Others" = "#66A61E", "NA" = "#666666"))+
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #scale_color_brewer(palette = "Dark2")+
  labs(x = paste0("PC1 (", round(pc.var.explained[1], 1), " %)"), y = paste0("PC2 (", pc.var.explained[2], "%)"))+
  #geom_text_repel(inherit.aes = F, data = gg.res, aes(x = PC1, y = PC2, label = SampleID), max.overlaps = 20)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())
p6


p7 = ggplot(res, aes(x = tsne1, y = tsne2, shape = FGFR3.TACC3.fusion, color = Final_class,label = SampleID)) + 
  #scale_color_manual(values = list("Ganglioglioma" = "#6A3D9A", "Others" = "#66A61E", "NA" = "#666666"))+
  scale_shape_manual(values = c("Y" = 17, "N" = 1, "NA" = 0))+
  #geom_text_repel(show.legend = F)+
  #scale_color_brewer(palette = "Dark2")+
  #geom_text_repel(data = res[res$Final_class == "GBM-F3T3-O" & res$umap2 >8,], max.overlaps = 30)+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())
p7

#pdf("01.F3T3.pca_tsne.gg.pdf", width = 15, height = 8, useDingbats = F)
cowplot::plot_grid(p1, p2, p3,p6, p5, p4, nrow =2)
dev.off()


###

load("01.F3T3.classifierv11b6.tsne.rda")
load("01.F3T3.classifierv11b6.umap.rda")
#test.umap <- data.frame(test.umap)
head(test.tsne)        
test.tsne$class <- meta$subtype[match(rownames(test.tsne), meta$SampleID)]
test.tsne$class1 <- meta$Final_class[match(rownames(test.tsne), meta$SampleID)]
test.tsne$class2 <- test.tsne$class
test.tsne$class2 <- ifelse(test.tsne$class1 == 'GBM-F3T3-O', test.tsne$class1, test.tsne$class2)
#test.umap$class <- meta$Final_class[match(rownames(test.umap), meta$SampleID)]

cancer.col <- c("#6A3D9A", "#CAB2D6", "#66A61E", "red", "gray")
names(cancer.col) <- c("Ganglioglioma", "GBM-G34", "GBM", "GBM-F3T3-O", "reference")

#pdf("03.F3T3.DKFZ.tsne.pdf", width = 10, height = 10, useDingbats = F)
ggplot(ref.tsne, aes(x = X1, y = X2)) +
  geom_point(color = "grey")+ theme_bw()+
  geom_point(inherit.aes = F, data = test.tsne, aes(x =X1, y = X2, color = class2), shape = 17)+
  scale_color_manual(values = subtype.col[names(subtype.col) %in% test.tsne$class2])+
  #geom_point(inherit.aes = F, data = test.tsne, aes(x =X1, y = X2, color = class1), shape = 17)+
  #scale_color_manual(values = cancer.col)+
  labs(x = "tsne 1", y = "tsne 2")+
  theme(panel.grid = element_blank(), aspect.ratio = 1)

dev.off()  


#pdf("03.F3T3.DKFZ.umap.pdf", width = 10, height = 10, useDingbats = F)
ggplot(ref.umap, aes(x = X1, y = X2)) +
  geom_point(color = "grey")+ theme_bw()+
  theme(panel.grid = element_blank())+
  geom_point(inherit.aes = F, data = test.umap, aes(x =X1, y = X2, color = class))+
  scale_color_manual(values = cancer.col)

dev.off()

subtype.col['reference'] <- "gray"
subtype.col

library(plotly)
p1 = plot_ly(x = ref.tsne$X1, y = ref.tsne$X2, text = y.ref,
             #shape = res$, 
             #color = y.ref, 
             color = rep("reference", nrow(ref.tsne)),
             type="scatter", 
             colors = subtype.col,
             mode = "markers", 
             marker = list(size = 3, line = list(width=1)) #symbols = c('dot')
) %>% 
  add_trace(x = test.tsne$X1, y = test.tsne$X2, text = paste0(rownames(test.tsne), test.tsne$class), color = test.tsne$class2, type="scatter") %>% 
  layout(showlegend = TRUE, 
         xaxis = list(title = "tsne 1", zeroline = FALSE),
         yaxis = list(title = "tsne 2", zeroline = FALSE))

p1


p2 = plot_ly(x = ref.umap$X1, y = ref.umap$X2, text = y.ref,
             #shape = res$, 
             #color = y.ref, 
             color = rep("reference", nrow(ref.umap)),
             type="scatter", 
             colors = cancer.col,
             mode = "markers", 
             marker = list(size = 3, line = list(width=1)) #symbols = c('dot')
) %>% 
  add_trace(x = test.umap$X1, y = test.umap$X2, text = paste0(rownames(test.umap), test.umap$class), color = test.umap$class, type="scatter") %>% 
  layout(showlegend = TRUE, 
         xaxis = list(title = "umap 1", zeroline = FALSE),
         yaxis = list(title = "umap 2", zeroline = FALSE))

p2

p1

##
load("03.DKFZ.NCI.preprocessFunnorm.betas.combat.32k.pc.tsne.rda")
meta1<- read.csv("DKFZ.NCI.samples.csv", head =T)
head(tsne)
N_non_trivial_pc #40, combat
custom.config = umap.defaults
custom.config$n_neighbors = 8
umap <- data.frame(umap::umap(betas.pca$x[,1:N_non_trivial_pc], config = custom.config)$layout)
umap$class <- meta1$methylation_class[match(rownames(umap), meta1$ID)]
umap$class[umap$class == "GBM, G34"] <- "GBM-G34"
umap$class[umap$class == "LGG, GG"] <- "Ganglioglioma"
umap$class2 <- meta$subtype[match(rownames(umap), meta$SampleID)]
umap$class <- ifelse(!is.na(umap$class2), umap$class2, umap$class)
umap$class2 <- ifelse(is.na(umap$class2), umap$class, umap$class2)
umap$class3 <- meta$Final_class[match(rownames(umap), meta$SampleID)]
umap$class2 <- ifelse(umap$class3== "GBM-F3T3-O" & !is.na(umap$class3), umap$class3, umap$class2)
umap$class <- ifelse(umap$class2 == "GBM-F3T3-O" & !is.na(umap$class3), umap$class3, umap$class)
umap$source <- meta1$Source[match(rownames(umap), meta1$ID)]


head(umap)
tail(umap)

unique(umap$class)[!unique(umap$class) %in% names(subtype.col)]
umap$class2 <- umap$class
umap$class2[umap$source == "DKFZ.reference"] <- "reference"

unique(umap$class)
unique(umap$class2)

subtype.col
cancer.col <- c("#6A3D9A", "#CAB2D6", "#66A61E", "red", "#888888", "#888888", "#000000")
names(cancer.col) <- c("Ganglioglioma", "GBM-G34", "GBM", "GBM-F3T3-O", "CONTR, REACT", "CONTR, INFLAM", "CONTR, WM")
cancer.col[5:7]
subtype.col <-c(subtype.col, cancer.col[5:7])


#pdf("03.F3T3.DKFZ.reNorm.umap.pdf")

pdf("03.DKFZ.NCI.preprocessFunnorm.betas.combat.32k.pdf")
ggplot(umap, aes(x = X1, y = X2, col = class, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggplot(umap, aes(x = X1, y = X2, col = class2, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

dev.off()


plot_ly(x = umap$X1, y = umap$X2, text = umap$class,
        #shape = umap$source, 
        #color = y.ref, 
        color = umap$class2,
        type="scatter", 
        colors = c(cancer.col, subtype.col),
        mode = "markers", 
        marker = list(size = 3, line = list(width=1)) #symbols = c('dot')
) %>% 
  #add_trace(x = test.umap$X1, y = test.umap$X2, text = paste0(rownames(test.umap), test.umap$class), color = test.umap$class, type="scatter") %>% 
  layout(showlegend = TRUE, 
         xaxis = list(title = "umap 1", zeroline = FALSE),
         yaxis = list(title = "umap 2", zeroline = FALSE))

###
#load("03.DKFZ.NCI.meffil.pc.tsne.rda")
rownames(betas.pca$x)

meta$idat2 <- meta$idat
meta$idat2[grep("^GSM", meta$idat, value = F)] <- sapply(grep("GSM", meta$idat, value = T), function(x) unlist(strsplit(x, "\\_"))[1])

rownames(betas.pca$x) <- meta$ID[match(rownames(betas.pca$x), meta$idat2)]

N_non_trivial_pc
custom.config = umap.defaults
#custom.config$n_neighbors = 10
umap <- data.frame(umap::umap(betas.pca$x[,1:N_non_trivial_pc], config = custom.config)$layout)
umap$class <- meta$methylation_class[match(rownames(umap), meta$ID)]

umap$class[umap$class == "GBM, G34"] <- "GBM-G34"
umap$class[umap$class == "LGG, GG"] <- "Ganglioglioma"

umap$source <- meta$Source[match(rownames(umap), meta$ID)]
umap$class2 <- umap$class
umap$class2[umap$source == "DKFZ.reference"] <- "reference"

head(umap)
tail(umap)

unique(umap$class)[!unique(umap$class) %in% names(subtype.col)]

#pdf("03.F3T3.DKFZ.reNorm.Meffil.umap.pdf")
ggplot(umap, aes(x = X1, y = X2, col = class, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggplot(umap, aes(x = X1, y = X2, col = class2, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

dev.off()


###preprocessQuantile
#load("03.DKFZ.NCI.preprocessQuantile.pc.tsne.rda")
load("03.DKFZ.NCI.ppreprocessQuantile.betas.combat.32k.pc.tsne.rda")

meta1 <- read.csv("DKFZ.NCI.samples.csv", head =T)
head(tsne)
N_non_trivial_pc #46, combat
custom.config = umap.defaults
custom.config$n_neighbors = 8
umap <- data.frame(umap::umap(betas.pca$x[,1:N_non_trivial_pc], config = custom.config)$layout)
umap$class <- meta1$methylation_class[match(rownames(umap), meta1$ID)]
umap$class[umap$class == "GBM, G34"] <- "GBM-G34"
umap$class[umap$class == "LGG, GG"] <- "Ganglioglioma"
umap$class2 <- meta$subtype[match(rownames(umap), meta$SampleID)]
umap$class <- ifelse(!is.na(umap$class2), umap$class2, umap$class)
umap$class2 <- ifelse(is.na(umap$class2), umap$class, umap$class2)
umap$class3 <- meta$Final_class[match(rownames(umap), meta$SampleID)]
umap$class2 <- ifelse(umap$class3== "GBM-F3T3-O" & !is.na(umap$class3), umap$class3, umap$class2)
umap$class <- ifelse(umap$class2 == "GBM-F3T3-O" & !is.na(umap$class3), umap$class3, umap$class)
umap$source <- meta1$Source[match(rownames(umap), meta1$ID)]
head(umap)
tail(umap)

unique(umap$class)[!unique(umap$class) %in% names(subtype.col)]


dev.off()

#pdf("03.F3T3.DKFZ.preprocessQuantile.umap.pdf")
#pdf("03.DKFZ.NCI.ppreprocessQuantile.betas.combat.32k.pdf")
ggplot(umap, aes(x = X1, y = X2, col = class, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggplot(umap, aes(x = X1, y = X2, col = class2, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

dev.off()

library(dplyr)
library(plotly)


plot_ly(x = umap$X1, y = umap$X2, text = paste0(umap$class),
        #shape = umap$source, 
        #color = y.ref, 
        color = umap$class2,
        type="scatter", 
        colors = subtype.col,
        mode = "markers", 
        marker = list(size = 3, line = list(width=1)) #symbols = c('dot')
) %>% 
  #add_trace(x = test.umap$X1, y = test.umap$X2, text = paste0(rownames(test.umap), test.umap$class), color = test.umap$class, type="scatter") %>% 
  layout(showlegend = TRUE, 
         xaxis = list(title = "umap 1", zeroline = FALSE),
         yaxis = list(title = "umap 2", zeroline = FALSE))






###preprocessIllumina
load("03.DKFZ.NCI.preprocessIllumina.pc.tsne.rda")
load("03.DKFZ.NCI.preprocessIllumina.combat.pc.tsne.rda")

meta <- read.csv("DKFZ.NCI.samples.csv", head =T)
head(tsne)
N_non_trivial_pc
custom.config = umap.defaults
custom.config$n_neighbors = 10
umap <- data.frame(umap::umap(betas.pca$x[,1:N_non_trivial_pc], config = custom.config)$layout)
umap$class <- meta$methylation_class[match(rownames(umap), meta$ID)]

umap$class[umap$class == "GBM, G34"] <- "GBM-G34"
umap$class[umap$class == "LGG, GG"] <- "Ganglioglioma"

umap$source <- meta$Source[match(rownames(umap), meta$ID)]


head(umap)
tail(umap)

unique(umap$class)[!unique(umap$class) %in% names(subtype.col)]
umap$class2 <- umap$class
umap$class2[umap$source == "DKFZ.reference"] <- "reference"

subtype.col
cancer.col <- c("#6A3D9A", "#CAB2D6", "#66A61E", "red", "#888888", "#888888", "#000000")
names(cancer.col) <- c("Ganglioglioma", "GBM-G34", "GBM", "GBM-F3T3-O", "CONTR, REACT", "CONTR, INFLAM", "CONTR, WM")



#pdf("03.F3T3.DKFZ.preprocessIllumina.umap.pdf")
#pdf("03.F3T3.DKFZ.preprocessIllumina.combat.umap.pdf")
ggplot(umap, aes(x = X1, y = X2, col = class, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggplot(umap, aes(x = X1, y = X2, col = class2, shape = source))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(subtype.col, cancer.col))+
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()


