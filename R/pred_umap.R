library(umap)
library(ggplot2)
library(ggrepel)

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

subtype.col = RColorBrewer::brewer.pal(12, "Paired")[c(6,10,8, 2,4)]
names(subtype.col) = c("NewCase", "GBM-F3T3-O", "Ganglioglioma", "GBM-G34", "GBM")

###
load("F3T3.umap.model.RDA")
     
###

RGset <- read.metharray("BL97/205841320139_R02C01")

###  
if(ncol(RGset) > 1){
  GRset <- preprocessFunnorm(RGset)
  betas <- getBeta(GRset); rm(GRset, RGset)
}else{
  Mset <- preprocessIllumina(RGset,  bg.correct = TRUE, normalize = "controls") 
  #Mset <- preprocessRaw(RGset)
  betas <- getBeta(Mset); rm(Mset, RGset)
}


if( sum(!rownames(betas.pca$rotation) %in% rownames(betas)) >0 ){
  stop("Error: some probes required not found: ", rownames(pc.rotation)[!rownames(pc.rotation) %in% rownames(betas)])
}
betas.32k <- betas[rownames(betas.pca$rotation),,drop=F]
betas.32k[is.na(betas.32k)] = 0.5

#pca transform
pc.pred <- predict(betas.pca, newdata = t(betas.32k))
#UMAP
umap.pred = predict(ref.umap, data = pc.pred[,1:N_non_trivial_pc, drop = F]);

res <- data.frame(rbind(ref.umap$layout, umap.pred))
colnames(res) <- c("umap1", "umap2")
res$class <- c(ref.tumor.class, rep("NewCase", ncol(betas.32k)))
res$source <- c(rep("ref", nrow(ref.umap$layout)), rep("NewCase", ncol(betas.32k)))
res$ID = rownames(res)



p <- ggplot(res, aes(x = umap1, y = umap2, col = class, label = ID))+
  geom_point()+
  geom_text_repel(data = res[res$source == "NewCase",], 
                  aes(x = umap1, y = umap2, label = ID),
                  show.legend = F, nudge_y = -0.8)+
  scale_color_manual(values = subtype.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())

pdf("result.pdf")
p
dev.off()

