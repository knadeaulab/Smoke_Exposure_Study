library(EnhancedVolcano)
library(limma)
library(circlize)
library(annotatr)
library(GenomicRanges)
library(regioneR)
library(reshape)
library(randomcoloR)
library(RColorBrewer)
library(dplyr)
options(bedtools.path = "/opt/homebrew/bin/")
library(bedtoolsr)
library(bios2mds)
library(magrittr)
library(ggpubr)
library(gridBase)
library(grid)
library(ComplexHeatmap)
library("viridis")
library(minfi)
library(dplyr)
library(PCAtools)
library(ENmix)
library(mixOmics)
library(mixOmics)
library(BayesFactor)
library(bayestestR)
load("MethylationDataset.RData") ## Loading pre-processed dataset ##
load("Corrmat.Rdata")

#### Section-1: PCA analysis ####
colnames(Mdata) == rownames(metadata) ## Mdata variable contrins M-values for each CpG site
x_vars = Mdata
y_var = metadata[colnames(x_vars),c("age","Race","Group","Batch_number","BMI")]
rownames(y_var) == colnames(x_vars) ## should be all TRUE
p <- PCAtools::pca(x_vars, metadata = as.data.frame(y_var), removeVar = 0.1)
out = eigencorplot(p, components = getComponents(p, seq_len(5)),
                   metavars = c('age', 'Race', 'Group', 'Batch_number','BMI'), ##
                   main  = "PCA analysis", scale = F, corMultipleTestCorrection = "BH", 
                   signifCutpoints = c(0, 0.00001, 0.001, 0.005, 1))
pdf("SuppFig2.pdf", width = 9, height = 6)
print(out)
dev.off()

#### Section-2: Differential methylation analysis ####
phenoData <- new("AnnotatedDataFrame", data=metadata)
eSet <- Biobase::ExpressionSet(assayData=as.matrix(Mdata), phenoData=phenoData)
eSet$Group <- factor(eSet$Group)
tmp = ifelse(eSet$RaceInfo=="White", "White", "other")
tmp[is.na(tmp)]<- "other"
eSet$Race2 <- factor(tmp)
design <- model.matrix(~ eSet$Group + eSet$Race2 + eSet$Batch_number)
fit <- lmFit(eSet, design)
fit <- eBayes(fit)
res = topTable(fit, n = Inf)
res$logFC = res$eSet.Group1
## volcano plot 
vp = EnhancedVolcano::EnhancedVolcano(toptable = res, x = 'logFC', y = 'adj.P.Val', lab = rownames(res), 
                                      labSize=0,  FCcutoff = 0.5, pCutoff = 0.2, xlim = c(-2, 2.3), ylim = c(0,7.5),
                                      pointSize = 4)


pdf("Fig2A.pdf", width = 5, height = 6)
print(vp)
dev.off()

#### Section-3: Manhattan plot ####
mm = match(rownames(res), locinfo$probes)
res$chr = locinfo$chr[mm]
res$start = locinfo$start[mm]
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
         'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
         'chr18','chr19','chr20','chr21','chr22')

colrs = randomcoloR::randomColor(22, luminosity = "dark")
res$chr = factor(res$chr, levels = chrs )
manhat = ggplot(res, aes(x = chr, y = -log10(adj.P.Val), color = as.factor(chr))) +
  geom_jitter(alpha = 0.75, width = 0.4, size = 1.2) +
  scale_color_manual(values = colrs) +
  geom_hline(yintercept = -log10(0.001), color = "red", linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.35))

pdf("SuppFig1B.pdf", width = 11, height = 4)
print(manhat)
dev.off()



#### Section-4: MDS plot ####
res_sig = dplyr::filter(res, abs(logFC) > 0.5 & adj.P.Val < 0.2 )
Mdata_t <- t(Mdata[rownames(res_sig),]) 
Mdata_df <- as.data.frame(Mdata_t)
pls_lda_model <- plsda(X = Mdata_t, Y = metadata$Group, ncomp = 2)  # 2 components for 2D plot
lda_scores <- pls_lda_model$variates$X
colnames(lda_scores) <- c("Dim.1", "Dim.2")
lda_scores = as.data.frame(lda_scores)
mm = match(rownames(lda_scores), metadata$ID)
lda_scores$Group = metadata$Group[mm]
lda_scores$female = metadata$female[mm]

PLSLDAplot = ggplot(lda_scores, aes(x = Dim.1, y = Dim.2, color = as.factor(Group))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('#4363d8','#f58231')) +
  theme_classic(base_size = 16)

pdf("Fig2B.pdf", width = 5.5, height = 3)
print(PLSLDAplot)
dev.off()


### Section-5: GREAT tool for CpG site-> gene-ID association analysis ####
## Note: GREAT (http://great.stanford.edu/public/html/index.php) 
## can be executed using web-browser. The following section of code 
## will generate two bed files that can be used as input for GREAT 
## web-browser. The output will include distance from TSS also.
## If you want, you can skip this section, please directly move to next section

### Enrichment Analysis: abs(logFC) > 0.5 & adj.P.Val < 0.2 ###
res_sig = dplyr::filter(res, abs(logFC) > 0.5 & adj.P.Val < 0.2 )

set.seed(123)

dd <- toGRanges(locinfo)
annots = c('hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries')
annotations = build_annotations(genome = 'hg19', annotations = annots)
dm_annotated = annotate_regions(
                  regions = dd,
                  annotations = annotations,
                  ignore.strand = TRUE,
                  quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
dm_annsum = summarize_annotations(annotated_regions = dm_annotated,quiet = TRUE)
print(dm_annsum)
bedout = data.frame()
for (feature in unique(df_dm_annotated$annot.type)){
  message(feature)
  tmp = dplyr::filter(df_dm_annotated, annot.type == feature)
  tmp = tmp[,c("seqnames","start", "end")] %>% dplyr::distinct()
  ## merge regions that fall under 1000 bp window and make it a single entry
  tet= bt.cluster(i = tmp, d = 1000)
  a = tet[,c("V1","V2","V4")] %>% group_by(V1,V4) %>% dplyr::summarise_all(min)
  b = tet[,c("V1","V2","V4")] %>% group_by(V1,V4) %>% dplyr::summarise_all(max)
  tmp = cbind(a[,c("V1","V2")],a[,"V2"] + 1)
  colnames(tmp) = c("chr","start","end")
  tmp$group = gsub("hg19_", "", feature)
  bedout = rbind(bedout, tmp)
}
mm = match(rownames(res_sig), df_dm_annotated$probes)
res_sig$chr = df_dm_annotated$seqnames[mm]
res_sig$start = df_dm_annotated$start[mm]
res_sig$end = df_dm_annotated$end[mm]
a = res_sig[,c("chr", "start", "end")] ## foreground for pathway analysis 
b = locinfo[,c("chr", "start", "end")] ## background 
colnames(a) = NULL
colnames(b) = NULL
write.table(a, file = "significantSite.bed", sep="\t", quote = F, row.names = F)
write.table(b, file = "background.bed", sep="\t", quote = F, row.names = F)
## the above two files were used as input for http://great.stanford.edu/public/html/index.php
## from the result page of GREAT, I extracted the region-geneID association excel file- GREAT_output_association_table.xlsx
GL = read.table("IntermediateFiles/DataF2C.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "TableS1.tsv", sep="\t", quote = F, row.names = F)

### Enrichment Analysis: adj.P.Val < 0.001 ###
res_sig = dplyr::filter(res, adj.P.Val < 0.005 ) # using only q < 0.01 as threshold #
mm = match(rownames(res_sig), df_dm_annotated$probes)
res_sig$chr = df_dm_annotated$seqnames[mm]
res_sig$start = df_dm_annotated$start[mm]
res_sig$end = df_dm_annotated$end[mm]
a = res_sig[,c("chr", "start", "end")] ## foreground for new pathway analysis 
colnames(a) = NULL
write.table(a, file = "significantSite_2.bed", sep="\t", quote = F, row.names = F)
## the above file was also used as input for http://great.stanford.edu/public/html/index.php
## from the result page of GREAT, I extracted the region-geneID association excel file- GREAT_output_association_table.xlsx
GL1 = read.table("IntermediateFiles/DataF2D.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL1$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "DataFile-S02.tsv", sep="\t", quote = F, row.names = F)

#### Section-6: upload gene list from association file to metascape (https://metascape.org/gp/index.html)  ####
write.table(file="CpGassociatedGenes4Metascape.txt", GL$V1, quote = F, row.names = F)
write.table(file="CpGassociatedGenes4Metascape_withQ_less_than_0.05.txt", GL1$V1, quote = F, row.names = F)

## Upload the gene list (abs(log FC) > 0.5 & q < 0.2) to Metascape pathway analysis software and download the Summary of enrichment analysis 
ggdf = read.table("IntermediateFiles/Metascape_1.tsv", sep="\t", header = T)
ggdf$logQ = -1 * (ggdf$Log10.P.)

pdf("Fig2C.pdf", height = 4, width = 8)
ggplot(ggdf, aes(x = logQ, y = reorder(Description, logQ, decreasing = FALSE) , fill = logQ)) + 
  geom_bar(stat = "identity", width = 0.3) +
  geom_point(alpha=0.99, size = 3.0) + 
  theme_classic(base_size = 12) +
  scale_fill_continuous(type = "viridis", trans = 'reverse')
dev.off()

## Upload the gene list (q < 0.005) to Metascape pathway analysis software and download the Summary of enrichment analysis in DisGeNET as TSV file 
ggdf = read.table("IntermediateFiles/Metascape_2.tsv", sep="\t", header = T)
ggdf$logQ = -1 * (ggdf$Log10.P.)
insteringPathaways = c("Allergic Reaction","Adult onset asthma","Squamous cell carcinoma of lung","Juvenile arthritis",
                       "Graves Disease","Primary biliary cirrhosis","Allergic rhinitis (disorder)",
                       "Hypothyroidism","Respiratory Tract Diseases","Autoimmune thyroid disease (AITD)","Cholangitis, Sclerosing")
ggdf = ggdf[ggdf$Description %in% insteringPathaways, ]

pdf("Fig2D.pdf", height = 2.5, width = 5.8)
ggplot(ggdf, aes(x = logQ, y = reorder(Description, logQ, decreasing = FALSE) , fill = logQ)) + 
  geom_bar(stat = "identity", width = 0.3) +
  geom_point(alpha=0.99, size = 3.0) + 
  theme_classic(base_size = 12) +
  scale_fill_continuous(type = "viridis", trans = 'reverse')
dev.off()


#### circos plot ####
Mdata = as.matrix(Mdata)
FFid = rownames(metadata)[metadata$Group==1]
ctrlid = rownames(metadata)[metadata$Group!=1]
a = rowMedians(Mdata[,FFid])
b = rowMedians(Mdata[,ctrlid])
a1 = data.frame(Methsite = rownames(Mdata), median = a)
b1 = data.frame(Methsite = rownames(Mdata), median = b)
mm = match(a1$Methsite, locinfo$probes)
a1 = cbind(locinfo[mm,c("chr","start","end")], a1)
mm = match(b1$Methsite, locinfo$probes)
b1 = cbind(locinfo[mm,c("chr","start","end")], b1)
a1$Methsite = NULL
colnames(a1)  = c("chr","start","end","value1")
b1$Methsite = NULL
colnames(b1)  = c("chr","start","end","value1")
target = c("genes_promoters","genes_intergenic","genes_introns","genes_exons","genes_3UTRs","genes_5UTRs" )
colrs=randomcoloR::randomColor(count = length(target), luminosity="dark" )
colrs = c('#e6194b','#3cb44b','#f032e6','#4363d8','#f58231','#911eb4')
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))

pdf("SuppFig1A.pdf")
circos.clear()
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
#circos.genomicHeatmap(a1, col = col_fun, side = "outside",connection_height = NULL, heatmap_height = 0.07) ## FF 
circos.genomicHeatmap(b1, col = col_fun, side = "outside",connection_height = NULL, heatmap_height = 0.07) ## controls
for ( i in 1:length(colrs)){
  gene_type = dplyr::filter(bedout,group %in% c(target[i]))
  cnt = length(unique(gene_type$group))
  gene_type$group = as.numeric(factor(gene_type$group, labels = cnt))
  circos.genomicDensity(gene_type, col = colrs[i], track.height = 0.05)
}
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
upViewport()
h = dev.size()[2]
lgd_meth = Legend(title = "Methylation", col_fun = col_fun)
lgd_list = packLegend(lgd_meth, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")
circos.clear()
dev.off()


#### PFAS analysis ####
metadata00 = metadata[!is.na(metadata$TotalPFOS),]
X = Mdata[,metadata00$ID]
colnames(X) = metadata00$PPID

getcorrPFOA <- function(x2){
  
  corr = cor.test(TPFOA,as.numeric(x2), method = "s")
  rho = corr$estimate
  pval = corr$p.value
  result <- correlationBF(TPFOA,as.numeric(x2))
  res = describe_posterior(result)
  out = 0
  if(res$ROPE_Percentage < 0.10 & res$BF > 1.0){
    if(pval < 0.05){
      out = res$Median
      out = rho
    }
  }
  return(out)
}
getcorrPFOS <- function(x2){
  
  corr = cor.test(TPFOS,as.numeric(x2), method = "s")
  rho = corr$estimate
  pval = corr$p.value
  out = 0
  result <- correlationBF(TPFOS,as.numeric(x2))
  res = describe_posterior(result)
  if(res$ROPE_Percentage < 0.10 & res$BF > 1.0){
    if(pval < 0.05){
      out = res$Median
      out = rho
    }
  }
  return(out)
}

TPFOA = metadata00$TotalPFOA
TPFOS = metadata00$TotalPFOS

## following bayesian analysis is extremely time-consuming. however Outcome variables (corrmat1,corrmat2) are already loaded 
#corrmat1 = apply(X, MARGIN = 1, FUN=function(x2) getcorrPFOA(x2))
#corrmat2 = apply(X, MARGIN = 1, FUN=function(x2) getcorrPFOS(x2))
#save(corrmat1, corrmat2, file = "Corrmat.Rdata")


PFOA_epigen = corrmat1[abs(corrmat1) >= 0.60]
PFOS_epigen = corrmat2[abs(corrmat2) >= 0.60] 
common = intersect(names(PFOA_epigen), names(PFOS_epigen))
edgelist1 = as.data.frame(PFOA_epigen)
edgelist1$PFAS = "PFOA"
edgelist1$from = rownames(edgelist1)
colnames(edgelist1) = c("weight","from","to")
edgelist2 = as.data.frame(PFOS_epigen)
edgelist2$PFAS = "PFOS"
edgelist2$from = rownames(edgelist2)
colnames(edgelist2) = c("weight","from","to")
edgelist = rbind(edgelist1,edgelist2)
rownames(edgelist) =NULL
edgelist$weight = NULL

mm = match(common, locinfo$probes)
bedc = locinfo[mm,c("chr", "start", "end")]
colnames(bedc) = NULL
write.table(bedc, file = "significantSite_common.bed", sep="\t", quote = F, row.names = F)


## Circos plot pfos in SE or firefighters ##
colr = randomcoloR::randomColor(count = 2, luminosity = c("random"))
names(colr) = c("PFOS","PFOS")

colr = c(PFOS="red",PFOA="darkblue")

pdf("Fig2F.pdf")
circos.clear()
par(cex = 0.9, mar = c(0, 0, 0, 0))
#circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
chordDiagram(edgelist, transparency = 0, grid.col = colr, 
             annotationTrack = c("grid"), preAllocateTracks = list(
               track.height = 0.040,
               col = colr,
               track.margin = c(0.00, 0.00)))
dev.off()


## bed file for GREAT analysis ##
## for those cpG sites that are significantly 
## associated pfoa pfos 
mm = match(edgelist1$to, locinfo$probes)
bed1 = locinfo[mm,c("chr", "start", "end")]
colnames(bed1) = NULL
mm = match(edgelist2$to, locinfo$probes)
bed2 = locinfo[mm,c("chr", "start", "end")]
colnames(bed2) = NULL
write.table(bed1, file = "significantSite_PFOA.bed", sep="\t", quote = F, row.names = F)
write.table(bed2, file = "significantSite_PFOS.bed", sep="\t", quote = F, row.names = F)
## use these bed files represent CpG sites significantly associated with pfas (total pfoa & total pfos)
## along with previously created background.bed file to upload in great webserver ###
## also get thier distance from TSS and gene-ID association ##
## upload the gene list from gene-ID association file to metascape webserver to get affected biological processes/pathways 
ggdf = read.table("IntermediateFiles/Metascape_PFOS.tsv", sep="\t", header = T)
ggdf$logQ = -1 * (ggdf$Log10.P.)
insteringPathaways = c("Allergic Reaction","Juvenile arthritis", "Cholangitis, Sclerosing",
                       "Allergic rhinitis (disorder)","Hypothyroidism","Ovarian Carcinoma",)
ggdf = ggdf[ggdf$Description %in% insteringPathaways, ]

pdf("SuppFig6D.pdf", height = 2.0, width = 4.8)
ggplot(ggdf, aes(x = logQ, y = reorder(Description, logQ, decreasing = FALSE) , fill = logQ)) + 
  geom_bar(stat = "identity", width = 0.3) +
  geom_point(alpha=0.99, size = 3.0) + 
  theme_classic(base_size = 12) +
  scale_fill_continuous(type = "viridis", trans = 'reverse')
dev.off()

GL00 = read.table("IntermediateFiles/F7C.txt", sep="\t")
# Here is the PFOS output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = GL00$V1,
                                  columns = c("GENENAME"),
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "DataFile-S05.tsv", sep="\t", quote = F, row.names = F)
GL01 = read.table("IntermediateFiles/F8C.txt", sep="\t")
# Here is the PFOA output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = GL01$V1,
                                  columns = c("GENENAME"),
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "DataFile-S06.tsv", sep="\t", quote = F, row.names = F)





## Lets create a scatter plot between most  correlated sites and Total PFOA ##
BX2= ENmix::M2B(X) ## CONVERT M values to B values for plotting 
ggdf = as.data.frame(PFOA_epigen)
ggdf$MethSite = rownames(ggdf)
ggdf$cor = abs(ggdf$PFOA_epigen)
ord = order(ggdf$cor, decreasing = T)
ggdf = ggdf[ord,]
#methlation = as.data.frame(t(BX2[ggdf$MethSite[1], ]))
#pfas2$X == rownames(methlation) ## if all TRUE we are good to go 
tmpsource = data.frame()
plots = list()

main_Cpgs = c("cg07709773", "cg12533148", "cg14190488")

## following two for loops generate data for Figure 2E ##
fig2e = data.frame()
for (i in 1:length(ggdf$MethSite)){
  message(i)
  #si = "cg04432046_BC11"
  si = ggdf$MethSite[i]
  methlation = as.data.frame(t(BX2[si, ]))
  mydf = data.frame(X = methlation[,1], Y =  metadata00$TotalPFOA, G = "Firefighter")
  mod = cor.test(x = mydf$X, y = mydf$Y, method = "p")
  pval = mod$p.value
  mydf$CpGSite = si
  tmpsource = rbind(tmpsource, mydf)
  if(sum(grepl(sub("_.*","",si), main_Cpgs)) > 0){
    idx= grepl(sub("_.*","",si), main_Cpgs)
    fig2e = rbind(fig2e, mydf)
  }
  g = ggplot(mydf, aes(x = X, y = Y, fill = G)) + 
    geom_point(size = 2.0, aes(color= G)) + 
    geom_smooth(method = "lm", aes(fill=G), color = "#4D4D4D") + 
    xlab(label = colnames(methlation)[1] ) + ylab(label = "Total PFOA") +
    scale_fill_manual(values = c('#f58231')) +
    scale_color_manual(values = c('#f58231')) +
    theme(line= element_blank()) + theme_bw(base_size = 11)
  plots[[i]] = g
}
fig2e$PFAS = "PFOA"

pdf("Fig2E_1.pdf", width = 2.8, height = 2.0) ## scatterplot DNA methylation total PFOA
for(pl in plots){
  print(pl)
}
dev.off()

### Lets create a scatter plot between most correlated sites and Total PFOS ###
ggdf = as.data.frame(PFOS_epigen)
ggdf$MethSite = rownames(ggdf)
ggdf$cor = abs(ggdf$PFOS_epigen)
ord = order(ggdf$cor, decreasing = T)
ggdf = ggdf[ord,]
plots = list()
main_Cpgs = c("cg00191575", "cg20793599", "cg03162381")

for (i in 1:length(ggdf$MethSite)){
  message(i)
  #si = "cg04432046_BC11"
  si = ggdf$MethSite[i]
  methlation = as.data.frame(t(BX2[si, ]))
  mydf = data.frame(X = methlation[,1], Y =  metadata00$TotalPFOS, G = "Firefighter")
  mod = cor.test(x = mydf$X, y = mydf$Y, method = "p")
  pval = mod$p.value
  mydf$CpGSite = si
  mydf$PFAS = "PFOS"
  #tmpsource = rbind(tmpsource, mydf)
  if(sum(grepl(sub("_.*","",si), main_Cpgs)) > 0){
    idx= grepl(sub("_.*","",si), main_Cpgs)
    fig2e = rbind(fig2e, mydf[,c("X","Y","G","CpGSite","PFAS")])
  }
  g = ggplot(mydf, aes(x = X, y = Y, fill = G)) + 
    geom_point(size = 2.0, aes(color= G)) + 
    geom_smooth(method = "lm", aes(fill=G), color = "#4D4D4D") + 
    xlab(label = colnames(methlation)[1] ) + ylab(label = "Total PFOS") +
    scale_fill_manual(values = c('#f58231')) +
    scale_color_manual(values = c('#f58231')) +
    theme(line= element_blank()) + theme_bw(base_size = 11)
  plots[[i]] = g
  
}

pdf("Fig2E_2.pdf", width = 2.8, height = 2.0) ## scatterplot DNA methylation total PFOS
for(pl in plots){
  print(pl)
}
dev.off()


#### Section-9: FF experience vs DNAmethylation analysis ####
metadataFF = metadata[!is.na(metadata$exposure_yrs),]
FFset = Mdata[,metadataFF$ID]

# Check if IDs match
if (!all(colnames(FFset) == metadataFF$ID)) {
  stop("Column names of FFset do not match IDs in metadataFF")
}

row_sd <- function(x) {apply(x, 1, sd)}
row_variation <- row_sd(FFset)
FFset <- FFset[row_variation > 0.2, ] ## analyzing only variable CpG sites in SE group 

results <- data.frame(CpG_site = rownames(FFset), 
                      Estimate = numeric(nrow(FFset)),
                      Std_Error = numeric(nrow(FFset)),
                      t_value = numeric(nrow(FFset)),
                      p_value = numeric(nrow(FFset)))
plots = list()
# Run GLM for each CpG site
for (i in 1:nrow(FFset)) {
  #message(i)
  
  cpg_data <- data.frame(M_value = as.numeric(FFset[i, ]), 
                         Experience = metadataFF$exposure_yrs,
                         G = "SE")
  
  # Fit GLM
  glm_model <- glm(M_value ~ Experience, data = cpg_data, family = gaussian())
  
  # Extract results
  model_summary <- summary(glm_model)
  coef_summary <- coef(model_summary)
  
  results$Estimate[i] <- coef_summary['Experience', 'Estimate']
  results$Std_Error[i] <- coef_summary['Experience', 'Std. Error']
  results$t_value[i] <- coef_summary['Experience', 't value']
  results$p_value[i] <- coef_summary['Experience', 'Pr(>|t|)']
  pval = coef_summary['Experience', 'Pr(>|t|)']
  if (pval < 0.05){
    #message(i)
  g = ggplot(cpg_data, aes(x = M_value, y = Experience, fill=G)) + 
      geom_point(size = 2.0, aes(color= G)) + 
      geom_smooth(method = "lm", aes(fill=G), color = "#4D4D4D") + 
      xlab(label = rownames(FFset)[i] ) + ylab(label = "Experience") +
      scale_fill_manual(values = c('#f58231')) +
      scale_color_manual(values = c('#f58231')) +
      theme(line= element_blank()) + theme_bw(base_size = 11)
  plots[[paste(i)]] = g
  }

}


results$fdr = p.adjust(results$p_value, method = "fdr")
results2 = dplyr::filter(results, fdr < 0.06 )
write.csv(results, "cpg_glm_results.csv", row.names = FALSE)


## check common sites ##
results2 = dplyr::filter(results, p_value < 0.05 )
res_sig = dplyr::filter(res, adj.P.Val < 0.005)
com = intersect(rownames(res_sig), results2$CpG_site)
per = length(com)/nrow(res_sig) #174 of the sites assocaites with experience are also differentially methylated
a = res_sig[com,]

sum(a$eSet.Group1 < 0) ## 31 higher methylation (hyper-methylation) in firefighter
sum(a$eSet.Group1 > 0) ## 46 lower methylation (hypo-methylation) in firefighter

mm = match(com, df_dm_annotated$probes)
a$chr = df_dm_annotated$seqnames[mm]
a$start = df_dm_annotated$start[mm]
a$end = df_dm_annotated$end[mm]
a = a[,c("chr","start","end")]
colnames(a) = NULL
write.table(a, file = "significantSite_experience.bed", sep="\t", quote = F, row.names = F)
## use this bed file to represent CpG sites significantly Differentially methylated sites 
## associated with years of exposure, along with previously created background.bed file 
## to upload in great webserver to get their distance from TSS and gene-ID association ##
## upload the gene list from gene-ID association file to metascape webserver 
## to get affected biological processes/pathways 

GL3 = read.table("DataEXP.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL3$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "Gene_Description_Experience.tsv", sep="\t", quote = F, row.names = F)

ggdf = read.table("IntermediateFiles/Metascape_Exp.tsv", sep="\t", header = T)
ggdf$logQ = -1 * (ggdf$Log10.P.)
notinsteringPathaways = c("Eosinophil count procedure","White Blood Cell Count procedure",
                          "Blood basophil count (lab test)", "Ankylosing spondylitis", "Corpuscular Hemoglobin Concentration Mean",
                          "Platelet Count measurement","Hyperplastic Polyp","Single umbilical artery")
ggdf = ggdf[!ggdf$Description %in% notinsteringPathaways, ]
pdf("SuppFig5C.pdf", height = 3.0, width = 8.8)
ggplot(ggdf, aes(x = logQ, y = reorder(Description, logQ, decreasing = FALSE) , fill = logQ)) + 
  geom_bar(stat = "identity", width = 0.3) +
  geom_point(alpha=0.99, size = 3.0) + 
  theme_classic(base_size = 12) +
  scale_fill_continuous(type = "viridis", trans = 'reverse')
dev.off()


