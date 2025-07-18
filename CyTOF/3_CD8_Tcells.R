library(ggplot2)
library(dplyr)
library(reshape2)
library(FlowSOM)
library(flowCore)
library(flowVS)
library(flowStats)
library(pheatmap)
library(patchwork)
library(mclust)
library(Seurat)
library(MAST)
library(umap)
library(Biobase)
library(preprocessCore)

dir.create("Data")

md = read.csv("Files/FileInfoCD8.csv") ## Sample annotation 
cd3 = read.csv("Files/CD3Pos_cell_counts.csv") ## CD3+ cell count per sample 
finalSamples = c("1294","1804","2765","2766","5023","5024","5865","5866","FF-024","WF046","WF135","FF-013",
                 "WF034","FF-035","FF-022","WF008","WF055","WF045","WF128","WF110","WF081","FF-052","FF-023",
                 "WF015","WF057","FF-045","WF023","WF002","WF005","WF042","WF051","WF094","WF102","WF120",
                 "FF-042","FF-026","FF-027","FF-028","FF-029","FF-031","FF-033","FF-041","FF-043","FF-016",
                 "FF-018","FF-020","FF-021","WF116","FF-008","FF-014","FF-019","FF-038","W0043B","WF066","W0050A",
                 "WF006","W0053A","WF069","WF025","WF138")

toxicMetals= c("79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W", 
               "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")

fm = c("CD27","CD28","CD154", "CD294","CD185", "CD62L", "CD183", "CD194", "CD196", "CD197", "CD195","CD279",
        "CXCR5","CXCR3","CCR4","CCR6","CCR7","CCR5","PD.1", "CRTH2")
functional_markers = unique(c(fm,toxicMetals))

load("Data/CD8.Rdata")
load("Data/CD8_T_HandGated_subsets.Rdata")

colrs = c( '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
           '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', 
           '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

######### Functions ###########
extract_key_cells <- function(data, keymarkers, clus, md2, toxicMetals){
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  hcless = dplyr::filter(data, cluster %in% clus)
  mplots = list()
  msioutput = data.frame()
  tempout = data.frame()
  for (mark in keymarkers){
    message(mark)
    X = as.matrix(hcless[mark])
    if (mark %in% toxicMetals){
      model <- Mclust(X, G = 2)
      clusters <- model$classification  # Get the cluster assignments
      a = median(X[clusters == 1])
      b = median(X[clusters != 1])
      if(a > b){target = 1} else{target = 2}
      tmpdf2 = hcless[,c(mark,"sample_id")]
      tmpdf2 = tmpdf2[clusters==target,]
    } else {
      thrsh = quantile(X, probs = 0.90)
      target = ifelse(X > thrsh, T, F)
      tmpdf2  = hcless[target, c(mark,"sample_id")]
    }
    colnames(tmpdf2) = c("Value","sample_id")
    
    prop = as.data.frame(table(tmpdf2$sample_id))
    Total =  as.data.frame(table(hcless$sample_id))
    mm = match(prop$Var1, Total$Var1)
    prop$Total = Total$Freq[mm]
    prop$Prop = (prop$Freq/prop$Total)*100
    colnames(prop)[1] = c("sample_id")
    mm = match(prop$sample_id, md2$SampleID)
    prop$Group = md2$Group[mm]
    prop$Group = factor(prop$Group, levels = c("nonSE", "SE"))
    a = table(prop$Group)[1]
    b = table(prop$Group)[2]
    if(b > 8){ ## if there are sufficient values in SE group 
      if (a == 0){ ## if there is no finding in nonSE add dummy data 
        tmp = data.frame(sample_id = "XXX", Freq = 0, Total = 0, Prop = 0, Group = "nonSE", Batch = "R1")
        prop = rbind(prop, tmp)
      }
      ggdf = prop
      ggdf$Group = ifelse(ggdf$Group == "nonSE", "nonSE", "SE")
      ggdf$Marker = mark
      ggdf$Group = factor(ggdf$Group, levels = c("nonSE", "SE"))
      gg = ggplot(ggdf, aes(x = Group, y = Prop, color = Group)) +
        geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
        geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
        scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
        scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
        xlab(label = "" ) + ylab(label = "% cluster cells") +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 14),
              legend.position = "none")  +  facet_wrap(~Marker) 
      mplots[[paste(mark)]]= gg
    }
  }
  return(mplots)
}
extract_key_cells_seuret <- function(data, keymarkers=fm, clus, md2, toxicMetals){
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  hcless = dplyr::filter(data, cluster %in% clus)
  keymarkers = intersect(colnames(hcless), keymarkers)
  
  Mdata = as.data.frame(t(hcless[,keymarkers]))
  colnames(Mdata) = rownames(hcless)
  seurat_object <- CreateSeuratObject(counts = Mdata)
  
  # Adding a new assay with additional data forms
  seurat_object[["Abseq"]] <- CreateAssayObject(
    data = as.matrix(Mdata),  # If normalized data available
    scale.data = as.matrix(Mdata)  # If scaled data available
  )
  DefaultAssay(seurat_object) <- "Abseq" 
  seurat_cell_names <- colnames(seurat_object)
  
  mm = match(hcless$sample_id, md2$SampleID)
  hcless$Group = md2$Group[mm]
  hcless$Batch = md2$Batch[mm]
  
  metadata <- data.frame(
    CellID = rownames(hcless),
    Group =  hcless$Group,
    Batch = hcless$Batch)
  seurat_object <- AddMetaData(seurat_object, metadata)
  # Scale the data if needed
  ## seurat_object <- ScaleData(seurat_object) ## not needed 
  Idents(seurat_object) <- seurat_object@meta.data$Group
  seurat_object <- SetIdent(seurat_object, value = "Group") 
  
  # Example of differential expression analysis using MAST
  differential_markers <- FindMarkers(
    seurat_object,
    ident.1 = "nonSE",
    ident.2 = "SE",
    test.use = "MAST",
    latent.vars = "Batch"
  )
  
  return(differential_markers)
}
MSI_analysis <- function(data, clus, functional_markers, md2, toxicMetals){ 
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  ggdf = dplyr::filter(data, cluster == clus)
  mm = match(ggdf$sample_id, md2$SampleID)
  ggdf$Group = md2$Group[mm]
  ggdf$group2 = ifelse(ggdf$Group == "nonSE", 0 , 1)
  ggdf$Group = factor(ifelse(ggdf$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
  fm = intersect(functional_markers, colnames(ggdf))
  BTplot = list()
  for (mark in fm){
    if (mark %in% toxicMetals){
      message(mark)
      X = as.matrix(ggdf[mark])
      thrsh = quantile(X, probs = 0.90)
      target = ifelse(X > thrsh, T, F)
      tmpdf2  = ggdf[target, c(mark,"sample_id")]
      colnames(tmpdf2) = c("Marker","sample_id")
      mm = match(tmpdf2$sample_id, md2$SampleID)
      tmpdf2$Group = md2$Group[mm]
      tmpdf2$Batch = md2$Batch[mm]
      tmpdf2 = tmpdf2 %>% dplyr::group_by(Group,Batch,sample_id) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
     } else {
      tmpdf = ggdf[,c(mark,"Group","group2","sample_id")]
      colnames(tmpdf)[1] = "Marker"
      tmpdf2 = tmpdf %>% dplyr::group_by(Group,group2,sample_id) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
    }
    tg = length(unique(as.character(tmpdf2$Group)))
    if ((sum(tmpdf2$Marker) > 0.1) & (tg == 2)) {
      gg = ggplot(tmpdf2, aes(x = Group, y = Marker, color = Group)) + 
        geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
        geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
        scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
        scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
        xlab(label = "" ) + ylab(label = "MSI") +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 14), 
              legend.position = "none")  +
        labs(title = clus, subtitle = mark) 
        BTplot[[paste(clus, mark)]]= gg
    }
  }
  return(BTplot)
}
subSample <- function(sample_ids,ncount){ 
  inds <- split(1:length(sample_ids), sample_ids) 
  ncells <- pmin(table(sample_ids), ncount) ## number of cells to downsample per FCS for clsuetering 
  fsom_inds <- lapply(names(inds),
                      function(i){
                        s <- sample(inds[[i]], ncells[i], replace = FALSE) ## randomly choosing 10000 cells from each sample
                      })
  fsom_inds <- unlist(fsom_inds) ## 50000 cells from each sample randomly removing rest of the cells 
  return(fsom_inds) 
}
subSampleByBatch <- function(tmpdf,ncount=3000, md2){
  tmpdf = dplyr::filter(tmpdf, sample_id %in% md2$SampleID)
  sampled_indices <- unlist(lapply(split(seq_len(nrow(tmpdf)), tmpdf$sample_id), function(indices) {
    if (length(indices) <= ncount) {
      indices  # Take all rows if less than 1000
    } else {
      sample(indices, ncount)  # Randomly sample 1000 rows
    }
  }))
  sampled_indices <- unlist(sampled_indices)
  return(sampled_indices)
}
drawUMAP <- function(df, n_neighbors=15, min_dist=0.5, n_components=2, metric = "euclidean",t="Umap",lineage_markers,v="Cell types") {
  title=paste(t,v,sep=":")
  custom.config = umap.defaults
  custom.config$random_state = 123
  custom.config$n_neighbors = n_neighbors ## usually 15
  custom.config$n_components = n_components ## dimension 
  custom.config$metric = metric
  custom.config$min_dist = min_dist
  #run UMAP
  embedding <- umap(df[,lineage_markers],custom.config)
  return(list(embedding))
}

## Total CD8+ T cell analysis ####
## CD3 cell count per sample ##
colnames(cd3) = c("Filename", "count")
mm = match(cd3$Filename, md$FileName)
cd3$SID = md$SampleID[mm]

prop = as.data.frame(table(PIDs))
colnames(prop) = c("SID","count")
mm = match(prop$SID, cd3$SID)
prop$Total = cd3$count[mm]
prop$Prop = (prop$count/prop$Total) * 100
mm = match(prop$SID, md$SampleID)
prop$Group = md$Group[mm]
prop$Ignore = md$Ignore[mm]
prop$PPID = md$PPID[mm]
prop$Batch = md$Batch[mm]
prop  = dplyr::filter(prop, Ignore == "No" & PPID %in%  finalSamples$V1)

prop$Group = factor(ifelse(prop$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
mod = glm(Group~Prop + Batch, data = prop, family = binomial)
out = summary(mod)
Pvalue = as.numeric(out$coefficients["Prop",4])

gg = ggplot(prop, aes(x = Group, y = Prop, color = Group)) + 
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
  scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
  scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
  xlab(label = "" ) + ylab(label = "% CD3 cells") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14), 
        legend.position = "none")  +
  labs(subtitle = paste("CD8+ Tcells | P-value", Pvalue) )

pdf("Supp.Fig.12a.pdf", width = 1.8, height = 4)
print(gg)
dev.off()

### memory vs Naive handgated CD8-T cells ###
prop0 = as.data.frame(table(NaiveCD8_Exp_matrix1$SampleID))
colnames(prop0) = c("SID","count")
prop0$CellType = "Naive"
prop1 = as.data.frame(table(MemoryCD8_Exp_matrix1$SampleID))
colnames(prop1) = c("SID","count")
prop1$CellType = "Memory"
prop = rbind(prop0,prop1)
mm = match(prop$SID, cd3$SID)
prop$Total = cd3$count[mm]
prop$Prop = (prop$count/prop$Total) * 100
mm = match(prop$SID, md$SampleID)
prop$Group = md$Group[mm]
prop$Ignore = md$Ignore[mm]
prop$PPID = md$PPID[mm]
prop$Batch = md$Batch[mm]
prop  = dplyr::filter(prop, Ignore == "No" & PPID %in%  finalSamples)

## P-value calculation of individual cell type ##
for (ct in unique(prop$CellType)){
  tmpdf = dplyr::filter(prop, CellType == ct)
  tmpdf$Group = factor(ifelse(tmpdf$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
  mod = glm(Group~Prop + Batch, data = tmpdf, family = binomial)
  out = summary(mod)
  Pvalue = as.numeric(out$coefficients["Prop",4])
  message(Pvalue)
}

prop$Group = factor(ifelse(prop$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
mod = glm(Group~Prop * CellType, data = prop, family = binomial)
out = summary(mod)
Pvalue_interaction = as.numeric(out$coefficients[4,4])

gg = ggplot(prop, aes(x = Group, y = Prop, color = Group)) + 
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
  scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
  scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
  xlab(label = "" ) + ylab(label = "% CD3 cells") +
  facet_grid(~CellType) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14), 
        legend.position = "none")  +
  labs(subtitle = paste("CD8+ Tcells | P-value", Pvalue_interaction) )

pdf("Supp.Fig.12b.pdf", width = 3.0, height = 4)
print(gg)
dev.off()

## FlowSOM clustering ##
lm = c("PD.1","CXCR3","CCR4","CXCR5","CCR6","CCR5","CRTH2", "HLA.DR", "CD62L","CD127","CD27", "CD28","CCR7","CD45RA")
md2 = dplyr::filter(md, Ignore == "No" & PPID %in% finalSamples)

## FlowSOM output is already stored in variable 'data' ##
## uncomment following lines to re-run FlowSOM ## Warning: Due to stochastic nature of FlowSOM results may show minor differences in each run) ##
# idx =  which(PIDs %in% md2$SampleID)
# data = ExpressionMatrix[idx,]
# sampleIDs = as.character(PIDs[idx])
# lineage_markers = intersect(lm, colnames(data))
# message("creating a flowset..")
# fcs.input <- new("flowFrame", exprs=as.matrix(data[,lineage_markers]))
# message("Running FlowSOM..")
# flowSOM.res <- FlowSOM(fcs.input, colsToUse = lineage_markers, scale = FALSE ,nClus=15) ## 15 clusters
# message("generating Meta-Clusters..")
# labels <- GetMetaclusters(flowSOM.res)
# data = as.data.frame(data)
# data$cluster = labels
# data$sample_id = sampleIDs
pdf7 = as.data.frame(data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarise_all(median))
nam = data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarize(count = n(), proportion = (n() / nrow(data))*100)
rownames(pdf7) = paste(rownames(pdf7), " (", format(nam$proportion, digits=2) , "%)", sep="")
clusOrder = c(1,3:8, 10,11,9,2,12:15) ##manually elected order to show CCR7-CD45RA+ cells above naive cells 
mm = match(clusOrder, pdf7$cluster)
pdf7$cluster = factor(pdf7$cluster , levels = clusOrder)

pdf("Supp.Figure.12e.pdf", width = 8, height = 5)
pheatmap(pdf7[mm,lm], cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white","lightblue", "navy", "orange", "firebrick3"))(50))
dev.off()


## UMAP: Enable this code to run UMAP again ##
#idx = subSampleByBatch(data,1000, md2)
#subsample.expr = data[idx,]
#SampleID = idx
#lineage_markers = intersect(lm, colnames(subsample.expr))
#list.out = drawUMAP(df=subsample.expr[,lineage_markers], n_neighbors=15, min_dist=0.5, n_components=2, metric = "euclidean",t="Umap",lineage_markers,v="Cell types") 

colrs = c( '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
           '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', 
           '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

ggdf = as.data.frame(list.out[[1]]$layout)
colnames(ggdf) = c("UMAP1","UMAP2")
ggdf$cluster = subsample.expr$cluster

g <- ggplot(ggdf,  aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) + theme_bw() + 
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=16),
        panel.grid = element_blank())+
  scale_color_manual(values =colrs )

pdf("Supp.Fig.12d.pdf", width = 5.5, height = 5 )
print(g)
centroids <- ggdf %>% group_by(cluster) %>% summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
g1 = g + geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = cluster), color = "black", size = 5, fontface = "bold")
print(g1)
dev.off()

##  MSI per marker in memory CD8 T cells ###
fm = c("CXCR3","CD27","CD154", "CD294","CD185", "CD62L", "CD183","CD194", "CD196", "CD197",
       "CD195","CD279","CXCR5","CCR4","CCR6","CCR7","CCR5","PD.1", "CRTH2")
functional_markers = fm
MemoryCD8_Exp_matrix0$cluster = "a"
outMSI = MSI_analysis(MemoryCD8_Exp_matrix0, clus='a', fm, md2, toxicMetals)

pdf("Supp.Figure.12c.pdf", width = 1.6, height = 3.5)
plt = outMSI[[1]]
print(plt)
tmpdf2 = outMSI[[1]]$data
tt = glm(Group2~Marker + Batch, data = tmpdf2, family = binomial)
out = summary(tt)
Pvalue = as.numeric(out$coefficients["Marker",4])
print(paste("CXCR3 P-value:", Pvalue))
dev.off()

## cells proportion with hi-marker (e.g. CCR7-hi cells) ##
## Following cods will produce box-plots; Only printing ##
## which we presented in study ##
clus_used = c(1,5,7,9,6)
pdf("Supp.Figure.12f-j.pdf", width = 1.6, height = 3.5)
for (clus in clus_used){
  if(clus %in% c(1)){markers_ = c("CXCR3","PD.1","CD154")}
  if(clus %in% c(5)){markers_ = c("CCR7")}
  if(clus %in% c(6)){markers_ = c("CRTH2","CCR6")}
  if(clus %in% c(7)){markers_ = c("CCR7")}
  if(clus %in% c(9)){markers_ = c("CCR6")}
  out = extract_key_cells(data, keymarkers = markers_ , clus = clus, md2, toxicMetals) ## retrieving cells with "hi" expression of a given marker
  tmpdf = data.frame()
  for(tmp in out){
    gg = tmp + labs(subtitle = paste("Cluster :", clus))
    tmp$data$Cluster = clus
    tmpdf = rbind(tmpdf, tmp$data)
    print(gg)
    }
}
dev.off()

##seurat_output: In each cluster, identify marker (high) cells and compute their proportion SE vs nonSE  
for(clus in unique(data$cluster)){ ## only these are non-naive clusters in which atleast one marker was found to be significant 
  output_1 = data.frame()
  message(clus)
  fm2 = intersect(fm, colnames(data))
  out_fdr = extract_key_cells_seuret(data, fm2, clus, md2, toxicMetals )
  out_fdr$Cluster = clus
  out_fdr$Marker = rownames(out_fdr)
  output_1 = rbind(output_1, out_fdr)
}
print(output_01)

maincluster = c(5,6,7,8,9,10,2,12,13)
pdf("Supp.Figure.12i-j.pdf")
for(clus in maincluster){
  out = MSI_analysis(data, clus, functional_markers , md2, toxicMetals)
  print(out[[1]])
}
dev.off()

### Toxic metal analysis ###
load("Files/CD8_T_HandGated_subsets.Rdata")
out1 = extract_key_cells(data, "110Cd" , "MemoryCD8_2", md2, toxicMetals)
out2 = extract_key_cells(data, "110Cd" , "NaiveCD8_2", md2, toxicMetals)
pdf("Figure4c.pdf", width = 1.6, height = 3.5)
print(out1[[1]])
print(out2[[1]])
dev.off()



