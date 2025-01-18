setwd("/Users/abk347/Library/CloudStorage/OneDrive-HarvardUniversity/FF_integrative/CyTOF/Version-3")
library(ggplot2)
library(dplyr)
library(reshape2)
library(FlowSOM)
library(flowCore)
library(flowVS)
library(flowStats)
library(pheatmap)
library(bayestestR)
library(BayesFactor)
library(rstanarm)
library(patchwork)
library(mclust)
library(Seurat)
library(MAST)
library(umap)
library(Biobase)
set.seed(123)

md = read.csv("FileInfoBcells.csv")
cd3 = read.csv("Data/HandGated/CD3neg_cell_count.csv") 
finalSamples = read.csv("../../Demographics_Table/Final_Samples.csv", header = F)
length(intersect(md$PPID, finalSamples$V1)) == 60 ## if TRUE; you are good to go 
#md = dplyr::filter(md, PPID %in% finalSamples$V1)

############################
########## B cells #########
############################
md$Fpath = paste("Data/HandGated/B_cell/", md$Filename, sep="")
color_conditions <- c('#4363d8','#f58231', 'black')
names(color_conditions) <- levels(factor(md$Group))
# Loading FCS files (You can also load individual file)

fcs_raw <- read.flowSet(as.character(md$Fpath[12]), transformation = FALSE, truncate_max_range = FALSE)
panel_fcs <- pData(parameters(fcs_raw[[1]])) ## getting the meta-data from FCS files object


cofactor = 5
ExpMatrix = list() ## rbinding in for loop can be slow
IDs = list()

flowframes= list() ##<<---

for (j in 1:nrow(md)){
  message(j)
  fcs_raw <- read.flowSet(as.character(md$Fpath[j]), transformation = FALSE, truncate_max_range = FALSE)
  lm = c("HLA.DR", "CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA",
         "TCRgd", "CD66b", "CD25", "CD56", "CD8a", "CD28","CD14", "CD127","CD19")
  fm  = c("CD27","CD154", "CD294","CD185", "CD62L", "CD183", "CD194", "CD196", "CD197", "CD195","CD279",
          "CXCR5","CXCR3","CCR4","CCR6","CCR7","CCR5","PD.1", "CRTH2",
          "79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W",
          "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd") ## includes both proteisn and toxic metals

  aliases <- c(
    "CD294" = "CRTH2",
    "CD185" = "CXCR5",
    "CD183" = "CXCR3",
    "CD194" = "CCR4",
    "CD196" = "CCR6",
    "CD195" = "CCR5",
    "CD279" = "PD.1",
    "CD197" = "CCR7"
  )

  panel_fcs <- pData(parameters(fcs_raw[[1]])) ## getting the meta-data from FCS files object
  panel_fcs$desc = gsub("_HM","",panel_fcs$desc)

  for (i in 1:nrow(panel_fcs)) {
    if (is.na(panel_fcs[i,]$desc))
    {
      panel_fcs[i,]$desc = panel_fcs[i,]$name ## removing all NAs for discription and adding channel name to it
    }
  }
  panel_fcs$desc = gsub(pattern = ".*_", replacement = "", panel_fcs$desc)## to clean the marker name e.g. I want to have marker name as CD8 or CD45
  panel_fcs$desc = gsub(pattern = "-", replacement = ".", panel_fcs$desc)
  panel_fcs <- panel_fcs %>% mutate(desc = ifelse(desc %in% names(aliases), aliases[desc], desc))

  panel_fcs$desc[panel_fcs$name=="Lu175Di"] <- "CD14"
  panel_fcs$desc = gsub("TGRgd","TCRgd",panel_fcs$desc)
  panel_fcs$desc = gsub("CD8A","CD8a",panel_fcs$desc)
  #panel_fcs$name = gsub("Lu176Di","Yb176Di",panel_fcs$name)
  #panel_fcs$name = gsub("Nd142Di","Ce142Di",panel_fcs$name)
  #check if all has been named correctly or not
  lineage_markers = intersect(lm,panel_fcs$desc)
  functional_markers = intersect(panel_fcs$desc, fm)
  all(lineage_markers %in% panel_fcs$desc) ## if TRUE then all are there
  all(functional_markers %in% panel_fcs$desc) ## if TRUE then all are there
  markers <- panel_fcs$name ## this should match with column name of FCS flowset; i.e. for each column name we have a coreesponding marker name
  names(markers) = as.character(panel_fcs$desc)
  markers = data.frame(Marker = names(markers), Metal = as.character(markers))
  markers = rbind(markers, data.frame(Marker = "CD14", Metal = 'Lu175Di'))

  df <- asinh(fsApply(fcs_raw, exprs)/ cofactor)
  mm = match(colnames(df), markers$Metal)
  colnames(df) = markers$Marker[mm]
  df = df[,c(lineage_markers, functional_markers)]
  IDs[[paste(j)]] = factor(rep(md$SampleID[j], nrow(df))) ## unique sample ID )
  #ExpMatrix[[paste(j)]] = df
  fcs.input <- new("flowFrame", exprs=as.matrix(df))
  flowframes[[paste(j)]] = fcs.input
}
#ExpressionMatrix = do.call(rbind, ExpMatrix)
#CD4lm = c("CXCR3","CCR4","CXCR5","CCR6","CCR5","CD127","CD27", "CCR7","CD45RA")
flow_set_1 <- flowSet(flowframes)
#fs.read.warped <- warpSet(flow_set_1, CD4lm)
ExpressionMatrix = fsApply(flow_set_1, exprs)
#ExpressionMatrix = fsApply(fs.read.warped, exprs)
PIDs = do.call(c, IDs)
save(ExpressionMatrix, lineage_markers, functional_markers, md, PIDs, finalSamples, cd3, file = "Expmatrix_Bcell.Rdata")
load("Expmatrix_Bcell.Rdata")


## CD3 cell count per sample ##
colnames(cd3) = c("RegEx", "count")
tmpdf = data.frame()
for(pattern in cd3$RegEx){
  indx = grep(pattern, md$Filename, ignore.case = TRUE)
  fn = md$FileName[indx]
  row = md[indx,]
  tmp = data.frame(pattern, row)
  tmpdf = rbind(tmp,tmpdf)
}
mm = match(cd3$RegEx,tmpdf$pattern)
cd3 = cbind(cd3, tmpdf[mm,])
cd3$SID = cd3$SampleID

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
  xlab(label = "" ) + ylab(label = "% CD3- cells") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14), 
        legend.position = "none")  +
  labs(subtitle = paste("B cells | P-value", Pvalue) )

pdf("Bcells_Boxplot_prop_total_CD3neg.pdf", width = 1.8, height = 4)
print(gg)
dev.off()

##################################################
############ Flow Som clustering   ###############
#################################################
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


md2 = dplyr::filter(md, Ignore == "No" & PPID %in% finalSamples$V1 )
idx =  which(PIDs %in% md2$SampleID)
data = ExpressionMatrix[idx,]
sampleIDs = as.character(PIDs[idx])

lm = c("CCR7", "CCR6", "CXCR5","HLA.DR","CD27", "CD19", "CD62L")

lineage_markers = intersect(lm, colnames(data))
message("creating a flowset..")
fcs.input <- new("flowFrame", exprs=as.matrix(data[,lineage_markers]))
message("Running FlowSOM..")
flowSOM.res <- FlowSOM(fcs.input, colsToUse = lineage_markers, scale = TRUE ,nClus=8) ## 15 clusters
message("generating Meta-Clusters..")
labels <- GetMetaclusters(flowSOM.res)
data = as.data.frame(data)
data$cluster = labels
data$sample_id = sampleIDs

load("Bcell_UMAP.Rdata")

### Heatmap ###
pdf7 = as.data.frame(data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarise_all(median))
nam = data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarize(count = n(), proportion = (n() / nrow(data))*100)
rownames(pdf7) = paste(rownames(pdf7), " (", format(nam$proportion, digits=2) , "%)", sep="")
pdf("Bcell_FlowSOM_clusters_Heatmap.pdf", width = 8, height = 5)
pheatmap(pdf7[,lm], cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white","lightblue", "navy", "orange", "firebrick3"))(50))
dev.off()

########################
######### UMAP #########
########################
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

colrs = c( '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
           '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', 
           '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

idx = subSampleByBatch(data, 1000, md2)
subsample.expr = data[idx,]
SampleID = data$sample_id[idx]
cluster_labels = data$cluster[idx]

lm = c("CCR7", "CCR6", "CXCR5","HLA.DR","CD27", "CD19", "CD62L")

lineage_markers = intersect(lm, colnames(subsample.expr))
list.out = drawUMAP(df=subsample.expr[,lineage_markers], n_neighbors=15, min_dist=0.5, n_components=2, metric = "euclidean",t="Umap",lineage_markers,v="Cell types") 

##save(list.out, data, lineage_markers, subsample.expr, idx, PIDs, file="Bcell_UMAP.Rdata")

ggdf = as.data.frame(list.out[[1]]$layout)
colnames(ggdf) = c("UMAP1","UMAP2")
ggdf$cluster = cluster_labels

g <- ggplot(ggdf,  aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) + theme_bw() + 
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=16),
        panel.grid = element_blank())+
  scale_color_manual(values =colrs )

pdf("Bcell_UMAP_clusters.pdf", width = 5.5, height = 5 )
print(g)
centroids <- ggdf %>% group_by(cluster) %>% summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
g1 = g + geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = cluster), color = "black", size = 5, fontface = "bold")
print(g1)
dev.off()

############################################
### Proportion of cells in each cluster ####
############################################

tcells = as.data.frame(table(data$sample_id))
prop = as.data.frame(table(data$sample_id, data$cluster))
colnames(prop) = c("PID", "ClusterNumber", "count")
mm = match(prop$PID, tcells$Var1)
prop$total = tcells$Freq[mm]
prop$Prop = (prop$count/prop$total) * 100

mm = match(prop$PID, md2$SampleID)
prop$Group = md2$Group[mm]
prop$Batch = md2$Batch[mm]

out_prop = data.frame()
bplots = list()
for(clus in unique(prop$ClusterNumber)){
  ggdf = dplyr::filter(prop,ClusterNumber == clus)
  ggdf$group2 = ifelse(ggdf$Group == "nonSE", 0 , 1)
  ggdf$Group = factor(ifelse(ggdf$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
  
  tt = glm(Group~Prop + Batch, data = ggdf, family = binomial)
  out = summary(tt)
  Pvalue = as.numeric(out$coefficients["Prop",4])
  out = data.frame(pval = Pvalue, cluster = clus)
  #tt = t.test(Prop~Group, data = ggdf)
  #out = data.frame(pval = tt$p.value, cluster = clus)
  out_prop = rbind(out, out_prop)
  
  gg = ggplot(ggdf, aes(x = Group, y = Prop, color = Group)) + 
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
    scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
    scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
    xlab(label = "" ) + ylab(label = "% Total live cells") +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 14), 
          legend.position = "none")  +
    labs(subtitle = clus) 
  
  bplots[[clus]]= gg
}

out_prop$FDR = p.adjust(out_prop$pval, method = "fdr" )

pdf("Bcells_Cluster_cell_Prop_BoxPlot.pdf", width = 10, height = 11)
wrap_plots(plotlist = bplots, ncol = 6, nrow = 3)
dev.off()

write.csv(out_prop, file = "Bcells_Cluster_cell_Prop_BoxPlot.pvalue.csv", quote = F)
############################################
###### MSI per marker for each cluster #####
############################################
MSI_analysis <- function(data, clus, functional_markers, md2){ 
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  ggdf = dplyr::filter(data, cluster == clus)
  mm = match(ggdf$sample_id, md2$SampleID)
  ggdf$Group = md2$Group[mm]
  ggdf$Batch = md2$Batch[mm]
  ggdf$group2 = ifelse(ggdf$Group == "nonSE", 0 , 1)
  ggdf$Group = factor(ifelse(ggdf$Group == "nonSE","nonSE", "SE"), levels = c("nonSE","SE"))
  fm = intersect(functional_markers, colnames(ggdf))
  ## methodology-1 : Non seurat : MSI analysis 
  tmpout = data.frame()
  BTplot = list()
  for (mark in fm){
    tmpdf = ggdf[,c(mark,"Group","group2","Batch","sample_id")]
    colnames(tmpdf)[1] = "Marker"
    tmpdf2 = tmpdf %>% dplyr::group_by(Group,group2,Batch,sample_id) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
    if (sum(tmpdf2$Marker) > 0.1){
      tt = glm(Group~Marker + Batch, data = tmpdf2, family = binomial)
      out = summary(tt)
      Pvalue = as.numeric(out$coefficients["Marker",4])
      out = data.frame(pval = Pvalue, cluster = clus, Marker = mark)
      tmpout = rbind(tmpout, out)
      
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
      tt = lm(Marker~Group, data = tmpdf2)
      out = summary(tt)
      pvalue0 = as.data.frame(out$coefficients)[2,4]
      if ((pvalue0 < 0.06) | (Pvalue < 0.06)){
        BTplot[[paste(clus, mark)]]= gg
      }
    }
  }
  tmpout$FDR = p.adjust(tmpout$pval, method = "fdr" )
  return(list(tmpout,BTplot))
}
## using only markers expressed by b-CELLS 
functional_markers = c("CXCR5","CXCR3","CCR5","PD.1","CCR7","79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W", 
                       "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")

out_MSI = data.frame()
bxplots = list()
for(clus in unique(prop$ClusterNumber)){
  out = MSI_analysis(data, clus, functional_markers, md2)
  out_MSI = rbind(out[[1]], out_MSI)
  bxplots = append(bxplots,out[[2]])
}

pdf("Bcell_cluster_MSI_BoxPlot.pdf", width =2, height = 4.2)
for (p in bxplots){
  print(p)
}
dev.off()

dplyr::filter(out_MSI, pval < 0.05)
write.csv(out_MSI, file = "Bcell_cluster_MSI_BoxPlot.pvalue.csv", quote = F)


############################################
### key Markers + cells for each cluster ###
############################################
## In each cluster, identify marker (high) cells and compute their proportion SE vs nonSE  

toxicMetals= c("79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W", 
               "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")

## Method:1 = GLM based analysis ##
extract_key_cells <- function(data, keymarkers, clus, md2){
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
    prop$Batch = md2$Batch[mm]
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
      tt = glm(Group~Prop + Batch, data = ggdf, family = binomial)
      out = summary(tt)
      Pvalue = as.numeric(out$coefficients["Prop",4])
      out = data.frame(pval = Pvalue, cluster = clus, Marker = mark)
      tempout = rbind(tempout, out)
      
      tt = t.test(Prop~Group, data = ggdf)
      pval = tt$p.value
      
      gg = ggplot(ggdf, aes(x = Group, y = Prop, color = Group)) +
        geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
        geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
        scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
        scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
        xlab(label = "" ) + ylab(label = "% cluster cells") +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 14),
              legend.position = "none")  +  facet_wrap(~Marker) + labs(subtitle = paste(pval,"|",clus))
      
      if(pval <0.05){
        mplots[[paste(mark)]]= gg
      }
    }
  }
  tempout$fdr = p.adjust(tempout$pval, method = "fdr")
  return(list(mplots, tempout))
}

mainClusters = c(1:8)  ### only these are clusters in which atleast one marker was found to be significant 
ggdfg = data.frame()
for (clus in mainClusters){ 
  message(clus)
  fm = intersect(functional_markers, colnames(data))
  out = extract_key_cells(data, fm , clus, md2)
  out_fdr = out[[2]]
  filename = paste("Bcells_", clus, "_prop.pdf", sep="")
  pdf(filename, width = 1.6, height = 3.5)
  for (plt in out[[1]]){
    print(plt)
  }
  dev.off()
  ggdfg = rbind(out_fdr, ggdfg)
}

write.csv(ggdfg, file = "Bcells_MSI_P_value_per_marker_per_clus.csv", quote = F, row.names = F)

## Method:2 = seurat based analysis ##
extract_key_cells_seuret <- function(data, keymarkers, clus, md2){
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  hcless = dplyr::filter(data, cluster %in% clus)
  keymarkers = intersect(colnames(hcless), functional_markers)
  
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

mainClusters = c(1:8) 

seurat_output = data.frame()
## In each cluster, identify marker (high) cells and compute their proportion SE vs nonSE  
for (clus in mainClusters){ ## only these are non-naive clusters in which atleast one marker was found to be significant 
  message(clus)
  fm = intersect(functional_markers, colnames(data))
  out_fdr = extract_key_cells_seuret(data, fm, clus, md2)
  out_fdr$Cluster = clus
  out_fdr$Marker = rownames(out_fdr)
  seurat_output = rbind(seurat_output, out_fdr)
}

write.csv(seurat_output, file = "Bcells_Seuret_P_value_per_marker_per_clus.csv", quote = F, row.names = F)


#################################
##### B memory vs Naive #######
#################################

getExpressionMatrix <- function(mdp){
  cofactor = 5
  ExpMatrix = list() ## rbinding in for loop can be slow
  IDs = list()
  for (j in 1:nrow(mdp)){
    message(j)
    fcs_raw <- read.flowSet(as.character(mdp$Fpath[j]), transformation = FALSE, truncate_max_range = FALSE)
    lm = c("HLA.DR", "CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA",
           "TCRgd", "CD66b", "CD25", "CD56", "CD8a", "CD28","CD14", "CD127","CD19")
    fm  = c("CD27","CD154", "CD294","CD185", "CD62L", "CD183", "CD194", "CD196", "CD197", "CD195","CD279",
            "CXCR5","CXCR3","CCR4","CCR6","CCR7","CCR5","PD.1", "CRTH2",
            "79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W",
            "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd") ## includes both proteisn and toxic metals
    
    aliases <- c(
      "CD294" = "CRTH2",
      "CD185" = "CXCR5",
      "CD183" = "CXCR3",
      "CD194" = "CCR4",
      "CD196" = "CCR6",
      "CD195" = "CCR5",
      "CD279" = "PD.1",
      "CD197" = "CCR7"
    )
    
    panel_fcs <- pData(parameters(fcs_raw[[1]])) ## getting the meta-data from FCS files object
    panel_fcs$desc = gsub("_HM","",panel_fcs$desc)
    
    for (i in 1:nrow(panel_fcs)) {
      if (is.na(panel_fcs[i,]$desc))
      {
        panel_fcs[i,]$desc = panel_fcs[i,]$name ## removing all NAs for discription and adding channel name to it
      }
    }
    panel_fcs$desc = gsub(pattern = ".*_", replacement = "", panel_fcs$desc)## to clean the marker name e.g. I want to have marker name as CD8 or CD45
    panel_fcs$desc = gsub(pattern = "-", replacement = ".", panel_fcs$desc)
    panel_fcs <- panel_fcs %>% mutate(desc = ifelse(desc %in% names(aliases), aliases[desc], desc))
    
    panel_fcs$desc[panel_fcs$name=="Lu175Di"] <- "CD14"
    panel_fcs$desc = gsub("TGRgd","TCRgd",panel_fcs$desc)
    panel_fcs$desc = gsub("CD8A","CD8a",panel_fcs$desc)
    #check if all has been named correctly or not
    lineage_markers = intersect(lm,panel_fcs$desc)
    functional_markers = intersect(panel_fcs$desc, fm)
    all(lineage_markers %in% panel_fcs$desc) ## if TRUE then all are there
    all(functional_markers %in% panel_fcs$desc) ## if TRUE then all are there
    markers <- panel_fcs$name ## this should match with column name of FCS flowset; i.e. for each column name we have a coreesponding marker name
    names(markers) = as.character(panel_fcs$desc)
    markers = data.frame(Marker = names(markers), Metal = as.character(markers))
    markers = rbind(markers, data.frame(Marker = "CD14", Metal = 'Lu175Di'))
    
    df <- asinh(fsApply(fcs_raw, exprs)/ cofactor)
    mm = match(colnames(df), markers$Metal)
    colnames(df) = markers$Marker[mm]
    df = as.data.frame(df[,c(lineage_markers, functional_markers)])
    
    df$SampleID = mdp$SampleID[j] ## unique sample ID )
    ExpMatrix[[paste(j)]] = df
    
  }
  ExpressionMatrix = do.call(rbind, ExpMatrix)
  return(ExpressionMatrix)
}
#################################
load("Bcells_subsets.Rdata")

# # NAIVE B #
# foldername = "Data/HandGated/B_cell_CD27_Pos_Neg/B_CD27neg/"
# fcsfiles = list.files(path = foldername, pattern = "fcs", ignore.case = T)
# matched_indices <- sapply(md$RegEx, function(pattern) {
#   grep(pattern, fcsfiles, ignore.case = TRUE)
# })
# md00 = as.data.frame(as.matrix(unlist(matched_indices)))
# nrow(md00) == 85 ## TRUE good to go 
# md00$Ind = rownames(md00)
# md00$FN = fcsfiles[md00$V1]
# mm = match(md00$Ind, md$RegEx)
# md00 = cbind(md00, md[mm,])
# md00$V1=NULL
# md00$Ind= NULL
# md00$Fpath =  paste(foldername, md00$FN, sep="")
# NaiveB_Exp_matrix0 = getExpressionMatrix(mdp=md00)
# NaiveB_Exp_matrix0 = dplyr::filter(NaiveB_Exp_matrix0, SampleID %in% md2$SampleID)
# 
# # Memory B #
# foldername = "Data/HandGated/B_cell_CD27_Pos_Neg/B_CD27pos/"
# fcsfiles = list.files(path = foldername, pattern = "fcs", ignore.case = T)
# matched_indices <- sapply(md$RegEx, function(pattern) {
#   grep(pattern, fcsfiles, ignore.case = TRUE)
# })
# md01 = as.data.frame(as.matrix(unlist(matched_indices)))
# nrow(md01) == 85 ## TRUE good to go 
# md01$Ind = rownames(md01)
# md01$FN = fcsfiles[md01$V1]
# mm = match(md01$Ind, md$RegEx)
# md01 = cbind(md01, md[mm,])
# md01$V1=NULL
# md01$Ind= NULL
# md01$Fpath =  paste(foldername, md01$FN, sep="")
# MemoryB_Exp_matrix0 = getExpressionMatrix(mdp=md01)
# MemoryB_Exp_matrix0 = dplyr::filter(MemoryB_Exp_matrix0, SampleID %in% md2$SampleID)


## CD3 cell count per sample ##
prop0 = as.data.frame(table(NaiveB_Exp_matrix0$SampleID))
colnames(prop0) = c("SID","count")
prop0$CellType = "Naive"
prop1 = as.data.frame(table(MemoryB_Exp_matrix0$SampleID))
colnames(prop1) = c("SID","count")
prop1$CellType = "Memory"

prop = rbind(prop0,prop1)
mm = match(prop$SID, cd3$SID) ## see above 
prop$Total = cd3$count[mm]
prop$Prop = (prop$count/prop$Total) * 100
mm = match(prop$SID, md$SampleID)
prop$Group = md$Group[mm]
prop$Ignore = md$Ignore[mm]
prop$PPID = md$PPID[mm]
prop$Batch = md$Batch[mm]
prop  = dplyr::filter(prop, Ignore == "No" & PPID %in%  finalSamples$V1)


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
  labs(subtitle = paste("B cells | P-value", Pvalue_interaction) )

pdf("Bcells_Boxplot_prop_Memory_vs_naive.pdf", width = 3.0, height = 4)
print(gg)
dev.off()


#save(NaiveB_Exp_matrix0, MemoryB_Exp_matrix0, md00, md01, file = "Bcells_subsets.Rdata")

#################################################
##### Marker+ cells in gated sub-population #####
#################################################
fm = intersect(functional_markers, colnames(MemoryB_Exp_matrix0))
MemoryB_Exp_matrix0$cluster = "a"
out0 = extract_key_cells(MemoryB_Exp_matrix0, fm , clus='a', md2)
out_surt0 = extract_key_cells_seuret(MemoryB_Exp_matrix0, fm, clus='a', md2)
outMSI0 = MSI_analysis(MemoryB_Exp_matrix0, clus='a', fm, md2)

NaiveB_Exp_matrix0$cluster = "a"
out1 = extract_key_cells(NaiveB_Exp_matrix0, fm , clus="a", md2)
out_surt1 = extract_key_cells_seuret(NaiveB_Exp_matrix0, fm, clus='a', md2)
outMSI1 = MSI_analysis(NaiveB_Exp_matrix0, clus='a', fm, md2)



write.csv(out0[[2]], file = "Bcells_memory_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(out1[[2]], file = "Bcells_naive_P_value_per_marker.csv", quote = F, row.names = F)

write.csv(out_surt0, file = "Bcells_memory_seuret_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(out_surt1, file = "Bcells_naive_seuret_P_value_per_marker.csv", quote = F, row.names = F)

write.csv(outMSI0[[1]], file = "Bcells_memory_MSI_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI1[[1]], file = "Bcells_naive_MSI_P_value_per_marker.csv", quote = F, row.names = F)



pdf("Bcells_memory.pdf", width = 1.6, height = 3.5)
for (plt in out0[[1]]){print(plt)}
dev.off()
pdf("Bcells_naive.pdf", width = 1.6, height = 3.5)
for (plt in out1[[1]]){print(plt)}
dev.off()

pdf("Bcells_memoryMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI0[[2]]){print(plt)}
dev.off()
pdf("Bcells_naiveMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI1[[2]]){print(plt)}
dev.off()


