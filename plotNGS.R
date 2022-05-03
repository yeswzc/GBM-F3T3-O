setwd("/Users/wuz6/Documents/Project/18.FGFR3-TACC3/04.NGS/")
rm(list=ls())



data <- data.frame(readxl::read_xlsx("../MS/v2/TableS1-2.xlsx", skip = 1))
data$tumor <- data$Final_class

table(data$Tumor)
x <- lapply(data$Subclass..calibrated.score., function(k){
   xx <- unlist(strsplit(k, "[\\(:;\\)]"))
})
x <- do.call(rbind, x)
data$subtype <- x[,1]
data$subtype.score <- as.numeric(x[,2])
head(data)
min(data$subtype.score, na.rm = T)
data$subtype.score <- ifelse(data$subtype.score >0.3 & !is.na(data$subtype.score), data$subtype, "No match")

data$Molecule.result[1:4]
xx <- lapply(1:nrow(data), function(k){
   x <- data$Molecule.result[k]
   if(x == "NA" | is.na(x)) return(NULL)
   res <- unlist(strsplit(x, "[;,]"))
   res <- gsub("^\\s", "", res)
   res <- gsub("\\s$", "", res)
   if(length(res) >0){
      return(cbind(data$SampleID[k], data$tumor[k], data$FGFR3.TACC3.fusion[k], res, data$subtype[k]))
   }else{
      return(NULL)
   }
})

ngs.res <- data.frame(do.call(rbind, xx))
tail(ngs.res)

colnames(ngs.res) <- c("ID", "Tumor", "F3T3", "X", "subtype")
ngs.res = ngs.res[ngs.res$X != "",]
table(ngs.res$Tumor)
ngs.res <- ngs.res[ngs.res$Tumor %in% c("GBM", "GBM-F3T3-O"),]
ngs.res$tumor <- ngs.res$Tumor
ngs.res$tumor[ngs.res$Tumor == "GBM" & 
                ngs.res$F3T3 == "Y" & 
                ngs.res$Tumor != "GBM-F3T3-O"] = "GBM-F3T3"
ngs.res$tumor[ngs.res$tumor == "GBM"] = "GBM"
ngs.res$status = "MUT"
unique(ngs.res$tumor)

table(ngs.res$X)
#top.genes = names(sort(rowSums(table(sub.ngs$GENE, sub.ngs$ID)), decreasing = T)[1:4])
#ngs.res <- ngs.res[ngs.res$X %in% c("TERT", "TP53", "PTEN loss", "EGFR ampl", "FGFR3-TACC3", "CDKN2A loss", "CDKN2B loss"),]
select.genes <- c( "FGFR3-TACC3", "FGFR3-RABGAP1L", "CDKN2A/B loss","TERT", "TP53", '+7', '-10', "EGFR ampl", "MGMT methy")




ngs.res$to.order <- paste0(ngs.res$tumor, ":", ngs.res$subtype)
ngs.res <- ngs.res[order(ngs.res$to.order),]
ngs.res$ID <- factor(ngs.res$ID, levels = unique(ngs.res$ID))

f3t3 = data.frame(ID = unique(ngs.res$ID))
f3t3$tumor = ngs.res$tumor[match(f3t3$ID, ngs.res$ID)]
f3t3$subtype = ngs.res$subtype[match(f3t3$ID, ngs.res$ID)]
head(f3t3)

head(ngs.res)

fill.col =  RColorBrewer::brewer.pal(8, "Dark2")[c(1,2,3)]
names(fill.col) = c("GBM", "GBM-F3T3-O", "GBM-F3T3")

subtype.col =  RColorBrewer::brewer.pal(8, "Dark2")[c(2,3,4,5,6,7,8)]
subtype.col = c(subtype.col)
names(subtype.col) = c("GBM, MES", "GBM, MID", "GBM, MYCN", "GBM, RTK I", "GBM, RTK II", 
                       "GBM, RTK III", "No match")



p1 = ggplot(ngs.res, aes(x = ID, y = X)) + #reorder did not work here
  geom_bar(inherit.aes = F, data = f3t3, 
           aes(x= ID, y = 0.5, fill = subtype),stat = "identity")+
  geom_bar(inherit.aes = F, data = f3t3, 
           aes(x= ID, y = -1, fill = tumor),stat = "identity")+
  geom_tile(color = "white", na.rm = T) +
  geom_hline(yintercept = 0)+
  theme_bw() + 
  scale_x_discrete(expand=c(0,0))+ #limits=rev
  scale_y_discrete(limits = rev(select.genes), expand=c(0,0))+
  scale_fill_manual(values = c(fill.col, subtype.col))+
  theme(#axis.text.x = element_text(angle = 90), 
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  labs(x= "", y = "")

p1

pdf("molecular.pdf", width = 8, height = 3, useDingbats = F)
p1
dev.off()

