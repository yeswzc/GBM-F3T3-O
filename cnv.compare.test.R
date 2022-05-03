setwd("/Users/wuz6/Documents/Project/18.FGFR3-TACC3/03.CNV")
rm(list=ls())
library(maftools)

meta <- data.frame(readxl::read_xlsx("../00.data/FGFR3-TACC3.cases.xlsx"), check.names = F)
#meta <- meta[meta$Tumor == "GBM, IDHwt",]
head(meta)


if(0){
postive.gistic <- readGistic(gisticAllLesionsFile = "F3T3.postive.gistic/all_lesions.conf_90.txt", 
                             gisticAmpGenesFile = "F3T3.postive.gistic/amp_genes.conf_90.txt",
                             gisticDelGenesFile = "F3T3.postive.gistic/del_genes.conf_90.txt",
                             gisticScoresFile = "F3T3.postive.gistic/scores.gistic",
                             isTCGA = F)
postive.gistic
pdf("F3T3.postive.gistic-1.pdf", width = 10, height = 5)
gisticChromPlot(gistic = postive.gistic, markBands = "all")
dev.off()
pdf("F3T3.postive.gistic-2.pdf", width = 10, height = 10)
gisticOncoPlot(gistic = postive.gistic, top = 20)
dev.off()

negative.gistic <- readGistic(gisticAllLesionsFile = "F3T3.negative.gistic/all_lesions.conf_90.txt", 
                             gisticAmpGenesFile = "F3T3.negative.gistic/amp_genes.conf_90.txt",
                             gisticDelGenesFile = "F3T3.negative.gistic/del_genes.conf_90.txt",
                             gisticScoresFile = "F3T3.negative.gistic/scores.gistic",
                             isTCGA = F)
negative.gistic
pdf("F3T3.negative.gistic-1.pdf", width = 10, height = 5)
gisticChromPlot(gistic = negative.gistic, markBands = "all")
dev.off()

pdf("F3T3.negative.gistic-2.pdf", width = 10, height = 10)
gisticOncoPlot(gistic = negative.gistic, top = 20)
dev.off()


rm(list=ls())
gistic <- readGistic(gisticAllLesionsFile = "GBM.gistic/all_lesions.conf_90.txt", 
                     gisticAmpGenesFile = "GBM.gistic/amp_genes.conf_90.txt",
                     gisticDelGenesFile = "GBM.gistic/del_genes.conf_90.txt",
                     gisticScoresFile = "GBM.gistic/scores.gistic",
                     isTCGA = F)

fb <- meta[,c("SampleID", "FGFR3-TACC3")]
colnames(fb)[1] = "Tumor_Sample_Barcode"
gistic@data[1:4,1:4]

gisticOncoPlot(gistic = gistic, clinicalData = fb, clinicalFeatures = "FGFR3-TACC3", 
               sortByAnnotation = T, top = 10)

}
##
#rm(gistic)

d <- read.table("GBM.n.GG.gistic/all_lesions.conf_90.txt", head =T, sep ="\t")

dim(d)

idx1 <- which(colnames(d) %in% meta$SampleID[meta$`FGFR3-TACC3`=="Y"])
idx2 <- which(colnames(d) %in% meta$SampleID[meta$`FGFR3-TACC3` == "N" & meta$Tumor == "GBM, IDHwt"])
length(idx1)
length(idx2)

f.test <- lapply(1:round(nrow(d)/2), function(k){

     n11 <- sum(d[k,idx1] >0)
     n12 <- length(idx1) - n11
     n21 <- sum(d[k,idx2] >0)
     n22 <- length(idx2) - n21
     f <- fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
     return(c(f$p.value, f$estimate, n11/length(idx1), n21/length(idx2)))
})

f.test <- data.frame(do.call(rbind, f.test))
f.test$name <- d[1:round(nrow(d)/2),2]
d[,2] <- gsub("\\s","", d[,2])
colnames(f.test) <- c("p.value", "odd.ratio", "freq.Y", "freq.N")
f.test$cnv <- sapply(d[1:round(nrow(d)/2),1], function(x) unlist(strsplit(x, "\\s+"))[1])
f.test$name <- paste0(f.test$cnv, ":", d[1:round(nrow(d)/2),2])
f.test$FDR <- p.adjust(f.test$p.value, method = "fdr")
sum(f.test$FDR <0.05)
f.test[f.test$p.value <0.05, ]

f.test[f.test$FDR <0.05, ]

x<- which(f.test$FDR <0.05)

paste0(d[x,2], collapse = "|")

gistic <- readGistic(gisticAllLesionsFile = "GBM.n.GG.gistic/all_lesions.conf_90.txt", 
                     gisticAmpGenesFile = "GBM.n.GG.gistic/amp_genes.conf_90.txt",
                     gisticDelGenesFile = "GBM.n.GG.gistic/del_genes.conf_90.txt",
                     gisticScoresFile = "GBM.n.GG.gistic/scores.gistic",
                     isTCGA = F)

fb0 <- meta[meta$SampleID %in% colnames(d),]
fb0$`FGFR3-TACC3`[fb0$Tumor == "LGG, GG"] = "LGG-GG"
fb <- fb0[,c("SampleID", "FGFR3-TACC3")]
colnames(fb)[1] = "Tumor_Sample_Barcode"

gistic@data[1:4,1:4]
keep.bands <- grep(paste0(d[x,2], collapse = "|"), rownames(gistic@cnMatrix), value = T)

pdf("01.F3T3.postive.vs.neg.cnv.pdf", width = 10, height = 8, useDingbats = F)
gisticOncoPlot(gistic = gistic, clinicalData = fb, clinicalFeatures = "FGFR3-TACC3", 
               bands = keep.bands, sortByAnnotation = T )
dev.off()

library(ggplot2)
library(ggrepel)

p <- ggplot(f.test, aes(x= odd.ratio, y = -log10(FDR), label = name)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  #geom_vline(xintercept = 1, linetype = "dashed")+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
  geom_text_repel(data = subset(f.test, FDR <0.05), size =3 )
p

pdf("01.F3T3.postive.vs.neg.cnv.seg.fisherTest.pdf", width = 7, height = 7, useDingbats = F)  
p
dev.off()
write.csv(f.test, file = "01.F3T3.postive.vs.neg.cnv.seg.fisherTest.csv", quote = F, row.names = F)




### GBM-GG v.s. other GBM
idx1 <- which(colnames(d) %in% meta$SampleID[meta$`FGFR3-TACC3`=="Y" & meta$Classifier.v12 == "Ganglioglioma"])
idx2 <- which(colnames(d) %in% meta$SampleID[meta$Tumor == "GBM, IDHwt" & 
                                               (meta$Classifier.v12 != "Ganglioglioma" | meta$`FGFR3-TACC3`!="Y")])
length(idx1)
length(idx2)

f.test <- lapply(1:round(nrow(d)/2), function(k){
  
  n11 <- sum(d[k,idx1] >0)
  n12 <- length(idx1) - n11
  n21 <- sum(d[k,idx2] >0)
  n22 <- length(idx2) - n21
  f <- fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
  return(c(f$p.value, f$estimate, n11/length(idx1), n21/length(idx2)))
})

f.test <- data.frame(do.call(rbind, f.test))
f.test$name <- d[1:round(nrow(d)/2),2]
d[,2] <- gsub("\\s","", d[,2])
colnames(f.test) <- c("p.value", "odd.ratio", "freq.Y", "freq.N")
f.test$cnv <- sapply(d[1:round(nrow(d)/2),1], function(x) unlist(strsplit(x, "\\s+"))[1])
f.test$name <- paste0(f.test$cnv, ":", d[1:round(nrow(d)/2),2])
f.test$FDR <- p.adjust(f.test$p.value, method = "fdr")
sum(f.test$FDR <0.05)
f.test[f.test$p.value <0.05, ]

f.test[f.test$FDR <0.05, ]

x<- which(f.test$FDR <0.05)

paste0(d[x,2], collapse = "|")


fb0$`FGFR3-TACC3`[fb0$Tumor == "LGG, GG"] = "LGG-GG"
fb <- fb0[,c("SampleID", "FGFR3-TACC3")]
colnames(fb)[1] = "Tumor_Sample_Barcode"
fb$`FGFR3-TACC3`[fb$Tumor_Sample_Barcode %in% fb0$SampleID[fb0$Classifier.v12 == "Ganglioglioma" &
                                                             fb0$`FGFR3-TACC3`== "Y"]] = "GBM-GG"

table(fb$`FGFR3-TACC3`)

keep.bands <- grep(paste0(d[x,2], collapse = "|"), rownames(gistic@cnMatrix), value = T)

pdf("03.GBM-GG.vs.otherGBM.cnv.pdf", width = 10, height = 5, useDingbats = F)
gisticOncoPlot(gistic = gistic, clinicalData = fb, clinicalFeatures = "FGFR3-TACC3", 
               bands = keep.bands, sortByAnnotation = T )
dev.off()

p <- ggplot(f.test, aes(x= odd.ratio, y = -log10(FDR), label = name)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  #geom_vline(xintercept = 1, linetype = "dashed")+
  geom_point() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
  geom_text_repel(data = subset(f.test, FDR <0.05), size =3 )
p

pdf("03.GBM-GG.vs.otherGBM.cnv.seg.fisherTest.pdf", width = 7, height = 7, useDingbats = F)  
p
dev.off()

write.csv(f.test, file = "03.GBM-GG.vs.otherGBM.cnv.seg.fisherTest.csv", quote = F, row.names = F)

