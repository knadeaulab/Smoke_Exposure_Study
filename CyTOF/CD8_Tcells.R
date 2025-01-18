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
set.seed(123)
md = read.csv("FileInfoCD8.csv")
cd3 = read.csv("Data/HandGated/CD3Pos_cell_counts.csv")
finalSamples = read.csv("../../Demographics_Table/Final_Samples.csv", header = F)
length(intersect(md$PPID, finalSamples$V1)) == 60 ## if TRUE; you are good to go 
#md = dplyr::filter(md, PPID %in% finalSamples$V1)

#################################
########## CD8+ T cells #########
#################################
md$Fpath = paste("Data/HandGated/CD8/", md$FileName, sep="")
color_conditions <- c('#4363d8','#f58231', 'black')
names(color_conditions) <- levels(factor(md$Group))
# Loading FCS files (You can also load individual file)

fcs_raw <- read.flowSet(as.character(md$Fpath[12]), transformation = FALSE, truncate_max_range = FALSE)
panel_fcs <- pData(parameters(fcs_raw[[1]])) ## getting the meta-data from FCS files object


cofactor = 5
ExpMatrix = list() ## rbinding in for loop can be slow
IDs = list()
flowframes= list() 
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

flow_set_1 <- flowSet(flowframes)
ExpressionMatrix = fsApply(flow_set_1, exprs)

#CD4lm = c("CXCR3","CCR4","CXCR5","CCR6","CCR5","CD127","CD27", "CCR7","CD45RA")
#fs.read.warped <- warpSet(flow_set_1, CD4lm)
#ExpressionMatrix = fsApply(fs.read.warped, exprs)

PIDs = do.call(c, IDs)
#save(ExpressionMatrix, lineage_markers, functional_markers, md, PIDs, finalSamples, cd3, file = "Expmatrix_CD8.Rdata")
load("Expmatrix_CD8.Rdata")


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

pdf("CD8_Boxplot_prop_total_CD3.pdf", width = 1.8, height = 4)
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

lm = c("PD.1","CXCR3","CCR4","CXCR5","CCR6","CCR5","CRTH2",
       "HLA.DR", "CD62L","CD127","CD27", "CD28","CCR7","CD45RA")

lineage_markers = intersect(lm, colnames(data))
message("creating a flowset..")
fcs.input <- new("flowFrame", exprs=as.matrix(data[,lineage_markers]))
message("Running FlowSOM..")
flowSOM.res <- FlowSOM(fcs.input, colsToUse = lineage_markers, scale = FALSE ,nClus=15) ## 15 clusters
message("generating Meta-Clusters..")
labels <- GetMetaclusters(flowSOM.res)
data = as.data.frame(data)
data$cluster = labels
data$sample_id = sampleIDs

load("CD8_cell_UMAP.Rdata")
### Heatmap ###
pdf7 = as.data.frame(data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarise_all(median))
nam = data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarize(count = n(), proportion = (n() / nrow(data))*100)
rownames(pdf7) = paste(rownames(pdf7), " (", format(nam$proportion, digits=2) , "%)", sep="")

clusOrder = c(1,3:8, 10,11,9,2,12:15) ##manually elected order to show CCR7-CD45RA+ cells above naive cells 
mm = match(clusOrder, pdf7$cluster)
pdf7$cluster = factor(pdf7$cluster , levels = clusOrder)


pdf("CD8_cell_FlowSOM_clusters_Heatmap.pdf", width = 8, height = 5)
pheatmap(pdf7[mm,lm], cluster_rows = F, cluster_cols = F,
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

idx = subSampleByBatch(data,1000, md2)
subsample.expr = data[idx,]
SampleID = idx
cluster_labels = data$cluster[idx]
lm = c("PD.1","CXCR3","CCR4","CXCR5","CCR6","CCR5","CRTH2",
       "HLA.DR", "CD62L","CD127","CD27", "CD28","CCR7","CD45RA")
lineage_markers = intersect(lm, colnames(subsample.expr))

list.out = drawUMAP(df=subsample.expr[,lineage_markers], n_neighbors=15, min_dist=0.5, n_components=2, metric = "euclidean",t="Umap",lineage_markers,v="Cell types") 

## save(list.out, data, lineage_markers, subsample.expr, idx, PIDs, file="CD8_cell_UMAP.Rdata")

ggdf = as.data.frame(list.out[[1]]$layout)
colnames(ggdf) = c("UMAP1","UMAP2")
ggdf$cluster = cluster_labels

g <- ggplot(ggdf,  aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) + theme_bw() + 
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=16),
        panel.grid = element_blank())+
  scale_color_manual(values =colrs )

pdf("CD8_UMAP_clusters.pdf", width = 5.5, height = 5 )
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

pdf("CD8_Cluster_cell_Prop_BoxPlot.pdf", width = 10, height = 11)
wrap_plots(plotlist = bplots, ncol = 6, nrow = 3)
dev.off()
write.csv(out_prop, file="CD8_Cluster_cell_Prop_BoxPlot.pval.csv", quote = F, row.names = F)

############################################
###### MSI per marker for each cluster #####
############################################
MSI_analysis <- function(data, clus, functional_markers, md2, toxicMetals){ 
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

      # 
      # model <- Mclust(X, G = 2)
      # clusters <- model$classification  # Get the cluster assignments
      # a = median(X[clusters == 1])
      # b = median(X[clusters != 1])
      # if(a > b){target = 1} else{target = 2}
      # tmpdf2 = ggdf[,c(mark,"sample_id","Group","group2","Batch")]
      # tmpdf00 = tmpdf2[clusters==target,]
      # colnames(tmpdf00) = c("Marker","sample_id","Group","group2","Batch")
      # tmpdf2 = tmpdf00 %>% dplyr::group_by(Group,group2,Batch,sample_id) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
    } else {
      tmpdf = ggdf[,c(mark,"Group","group2","Batch","sample_id")]
      colnames(tmpdf)[1] = "Marker"
      tmpdf2 = tmpdf %>% dplyr::group_by(Group,group2,Batch,sample_id) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
    }
    tg = length(unique(as.character(tmpdf2$Group)))
    if ((sum(tmpdf2$Marker) > 0.1) & (tg == 2)) {
      tmpdf2$Group2 = ifelse(tmpdf2$Group=="SE", 1, 0)
      tt = glm(Group2~Marker + Batch, data = tmpdf2, family = binomial)
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
      
      if(is.na(pvalue0)){pvalue0 = 1}
      if(is.na(Pvalue)){pvalue0 = 1}
      
      
      if ((pvalue0 < 0.80)|(Pvalue < 0.80)){
        BTplot[[paste(clus, mark)]]= gg
      }
    }
  }
  tmpout$FDR = p.adjust(tmpout$pval, method = "fdr" )
  return(list(tmpout,BTplot))
}

functional_markers00 = c("CD27","CD28","CD154", "CD294","CD185", "CD62L", "CD183", "CD194", "CD196", "CD197", "CD195","CD279",
                         "CXCR5","CXCR3","CCR4","CCR6","CCR7","CCR5","PD.1", "CRTH2","79Br")

out_MSI = data.frame()
bxplots = list()
for(clus in unique(prop$ClusterNumber)){
  out = MSI_analysis(data, clus, functional_markers00, md2, toxicMetals)
  out_MSI = rbind(out[[1]], out_MSI)
  bxplots = append(bxplots,out[[2]])
}

pdf("CD8_cluster_MSI_BoxPlot.pdf", width =2, height = 4.2)
for (p in bxplots){
  print(p)
}
dev.off()

write.csv(out_MSI, file = "CD8_cluster_MSI_BoxPlot.pvalue.csv", quote = F)


toxicMetals= c("79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W", 
               "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")

out_MSI_metals = data.frame()
bxplotsmetals  = list()
maincluster = c(5,6,7,8,9,10,2,12,13)
for(clus in maincluster){
  out = MSI_analysis(data, clus, toxicMetals, md2, toxicMetals)
  out_MSI_metals = rbind(out[[1]], out_MSI_metals)
  bxplotsmetals  = append(bxplotsmetals ,out[[2]])
}

pdf("CD8_cluster_Metals_MSI_BoxPlot.pdf", width =2, height = 4.2)
for (p in bxplotsmetals){
  print(p)
}
dev.off()


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
      if((mark == "CCR7") & (clus == 5)){
        ggdf = dplyr::filter(ggdf, Prop < 80)
      }
      gg = ggplot(ggdf, aes(x = Group, y = Prop, color = Group)) +
        geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
        geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
        scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
        scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
        xlab(label = "" ) + ylab(label = "% cluster cells") +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 14),
              legend.position = "none")  +  facet_wrap(~Marker) + labs(subtitle = paste(pval,"|",clus))
      
      if(pval <0.30){
        mplots[[paste(mark)]]= gg
      }
    }
  }
  tempout$fdr = p.adjust(tempout$pval, method = "fdr")
  return(list(mplots, tempout))
}

mainClusters = c(1,2,5,6,7,8,9,10,12,13,15) ### only these are clusters in which atleast one marker was found to be significant 
ggdfg = data.frame()
for (clus in mainClusters){ 
  message(clus)
  fm = intersect(functional_markers, colnames(data))
  out = extract_key_cells(data, fm , clus, md2)
  out_fdr = out[[2]]
  filename = paste("CD8_", clus, "_prop.pdf", sep="")
  pdf(filename, width = 1.6, height = 3.5)
  for (plt in out[[1]]){
    print(plt)
  }
  dev.off()
  ggdfg = rbind(out_fdr, ggdfg)
}

write.csv(ggdfg, file = "CD8_Prop_P_value_per_marker_per_clus.csv", quote = F, row.names = F)

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

mainClusters = c(1,2,5,6,7,8,9,10,12,13,15)

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

write.csv(seurat_output, file = "CD8_Seuret_P_value_per_marker_per_clus.csv", quote = F, row.names = F)

#################################
##### CD8 memory vs Naive #######
#################################
getExpressionMatrix <- function(mdp){
  cofactor = 5
  ExpMatrix = list() ## rbinding in for loop can be slow
  IDs = list()
  flowframes = list()
  PIDs = list()
  for (j in 1:nrow(mdp)){
    message(j)
    fcs_raw <- read.flowSet(as.character(mdp$Fpath[j]), transformation = FALSE, truncate_max_range = FALSE)
    lm = c("Ir191", "Ir193", "HLA.DR", "CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA",
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
      "CD197" = "CCR7",
      "Intercalator" = "Ir191", 
      "Intercalator.Ir" = "Ir193",
      "DNA1" = "Ir191", 
      "DNA2" = "Ir193",
      "191Ir" = "Ir191",
      "193Ir" = "Ir193"
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

md = read.csv("FileInfoCD8.csv")
finalSamples = read.csv("../../Demographics_Table/Final_Samples.csv", header = F)
md$Fpath = paste("Data/HandGated/CD8/", md$FileName, sep="")

load("CD8cells_subsets.Rdata")

## "Intercalator", "Intercalator.Ir"
# NAIVE CD8 #
# foldername = "Data/HandGated/CD8_Memory_Naive/Naive_CD8/"
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
# NaiveCD8_Exp_matrix0 = getExpressionMatrix(mdp=md00)
# NaiveCD8_Exp_matrix0 = dplyr::filter(NaiveCD8_Exp_matrix0, SampleID %in% md2$SampleID)
# NaiveCD8_Exp_matrix0$cluster = "NaiveCD8_1"
# 
# ### Memory CD8 ####
# foldername = "Data/HandGated/CD8_Memory_Naive/Memory_CD8//"
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
# MemoryCD8_Exp_matrix0 = getExpressionMatrix(mdp=md01)
# MemoryCD8_Exp_matrix0 = dplyr::filter(MemoryCD8_Exp_matrix0, SampleID %in% md2$SampleID)
# MemoryCD8_Exp_matrix0$cluster = "MemoryCD8_1"


## CD3 cell count per sample ##
cd3 = read.csv("Data/HandGated/CD3Pos_cell_counts.csv")
colnames(cd3) = c("Filename", "count")
mm = match(cd3$Filename, md$FileName)
cd3$SID = md$SampleID[mm]

prop0 = as.data.frame(table(NaiveCD8_Exp_matrix0$SampleID))
colnames(prop0) = c("SID","count")
prop0$CellType = "Naive"
prop1 = as.data.frame(table(MemoryCD8_Exp_matrix0$SampleID))
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
  labs(subtitle = paste("CD8+ Tcells | P-value", Pvalue_interaction) )

pdf("CD8_Boxplot_prop_Memory_vs_naive.pdf", width = 3.0, height = 4)
print(gg)
dev.off()


# # # Excluded NAIVE CD8 #
# foldername = "Data/HandGated/CD8_Memory_Naive/Excluded_naive/"
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
# NaiveCD8_Exp_matrix1 = getExpressionMatrix(mdp=md00)
# NaiveCD8_Exp_matrix1 = dplyr::filter(NaiveCD8_Exp_matrix1, SampleID %in% md2$SampleID)
# NaiveCD8_Exp_matrix1$cluster = "NaiveCD8_2"
# 
# 
# # Excluded Memory CD8 #
# foldername = "Data/HandGated/CD8_Memory_Naive/Excluded_memory/"
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
# MemoryCD8_Exp_matrix1 = getExpressionMatrix(mdp=md01)
# MemoryCD8_Exp_matrix1 = dplyr::filter(MemoryCD8_Exp_matrix1, SampleID %in% md2$SampleID)
# MemoryCD8_Exp_matrix1$cluster = "MemoryCD8_2"

## CD3 cell count per sample ##
cd3 = read.csv("Data/HandGated/CD3Pos_cell_counts.csv")
colnames(cd3) = c("Filename", "count")
mm = match(cd3$Filename, md$FileName)
cd3$SID = md$SampleID[mm]

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
  labs(subtitle = paste("CD8+ Tcells | P-value", Pvalue_interaction) )

pdf("CD8_Boxplot_prop_Excluded_Memory_vs_naive.pdf", width = 3.0, height = 4)
print(gg)
dev.off()

#save(NaiveCD8_Exp_matrix0, MemoryCD8_Exp_matrix0, NaiveCD8_Exp_matrix1, MemoryCD8_Exp_matrix1,  md00, md01, file = "CD8cells_subsets.Rdata")

#################################################
##### Marker+ cells in gated sub-population #####
#################################################
fm = intersect(functional_markers, colnames(MemoryCD8_Exp_matrix0))
MemoryCD8_Exp_matrix0$cluster = "a"
NaiveCD8_Exp_matrix0$cluster = "a"
MemoryCD8_Exp_matrix1$cluster = "a"
NaiveCD8_Exp_matrix1$cluster = "a"

out0 = extract_key_cells(MemoryCD8_Exp_matrix0, fm , clus='a', md2)
out2 = extract_key_cells(MemoryCD8_Exp_matrix1, fm , clus="a", md2)
out1 = extract_key_cells(NaiveCD8_Exp_matrix0, fm , clus="a", md2)
out3 = extract_key_cells(NaiveCD8_Exp_matrix1, fm , clus="a", md2)

write.csv(out0[[2]], file = "CD8_memory_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(out1[[2]], file = "CD8_naive_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(out2[[2]], file = "CD8_excluded_memory_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(out3[[2]], file = "CD8_excluded_P_value_per_marker.csv", quote = F, row.names = F)

pdf("CD8_memory.pdf", width = 1.6, height = 3.5)
for (plt in out0[[1]]){print(plt)}
dev.off()
pdf("CD8_naive.pdf", width = 1.6, height = 3.5)
for (plt in out1[[1]]){print(plt)}
dev.off()
pdf("CD8_excluded_memory.pdf", width = 1.6, height = 3.5)
for (plt in out2[[1]]){print(plt)}
dev.off()
pdf("CD8_excluded_naive.pdf", width = 1.6, height = 3.5)
for (plt in out3[[1]]){print(plt)}
dev.off()


outMSI0 = MSI_analysis(MemoryCD8_Exp_matrix0, clus='a', fm, md2)
outMSI1 = MSI_analysis(NaiveCD8_Exp_matrix0, clus='a', fm, md2)
outMSI2 = MSI_analysis(MemoryCD8_Exp_matrix1, clus='a', fm, md2)
outMSI3 = MSI_analysis(NaiveCD8_Exp_matrix1, clus='a', fm, md2)
write.csv(outMSI0[[1]], file = "CD8_memory_MSI_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI1[[1]], file = "CD8_naive_MSI_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI2[[1]], file = "CD8_excluded_MSI_memory_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI3[[1]], file = "CD8_excluded_MSI_P_value_per_marker.csv", quote = F, row.names = F)
pdf("CD8_memoryMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI0[[2]]){print(plt)}
dev.off()
pdf("CD8_naiveMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI1[[2]]){print(plt)}
dev.off()
pdf("CD8_excluded_memoryMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI2[[2]]){print(plt)}
dev.off()
pdf("CD8_excluded_naiveMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI3[[2]]){print(plt)}
dev.off()


out_surt0 = extract_key_cells_seuret(MemoryCD8_Exp_matrix0, fm, clus='a', md2)
write.csv(out_surt0, file = "CD8_memory_seuret_P_value_per_marker.csv", quote = F, row.names = F)

out_surt1 = extract_key_cells_seuret(NaiveCD8_Exp_matrix0, fm, clus='a', md2)
write.csv(out_surt1, file = "CD8_naive_seuret_P_value_per_marker.csv", quote = F, row.names = F)

out_surt2 = extract_key_cells_seuret(MemoryCD8_Exp_matrix1, fm, clus='a', md2)
write.csv(out_surt2, file = "CD8_excluded_seuret_memory_P_value_per_marker.csv", quote = F, row.names = F)

out_surt3 = extract_key_cells_seuret(NaiveCD8_Exp_matrix1, fm, clus='a', md2)
write.csv(out_surt3, file = "CD8_excluded_seuret_P_value_per_marker.csv", quote = F, row.names = F)


#########################################
########### Toxin metals ################
#########################################
MemoryCD8_Exp_matrix0$cluster = "a"
NaiveCD8_Exp_matrix0$cluster = "a"
MemoryCD8_Exp_matrix1$cluster = "a"
NaiveCD8_Exp_matrix1$cluster = "a"
fm= c("79Br","81Br","209Bi","204Pb", "206Pb", "208Pb","196Hg", "198Hg", "202Hg","182W","184W", 
      "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")
toxicMetals = fm
outMSI0 = MSI_analysis(MemoryCD8_Exp_matrix0, clus='a', fm, md2,toxicMetals)
outMSI1 = MSI_analysis(NaiveCD8_Exp_matrix0, clus='a', fm, md2, toxicMetals)
outMSI2 = MSI_analysis(MemoryCD8_Exp_matrix1, clus='a', fm, md2, toxicMetals)
outMSI3 = MSI_analysis(NaiveCD8_Exp_matrix1, clus='a', fm, md2, toxicMetals)
write.csv(outMSI0[[1]], file = "CD8_memory_Metals_MSI_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI1[[1]], file = "CD8_naive_MSI_Metal_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI2[[1]], file = "CD8_excluded_MSI_Metal_memory_P_value_per_marker.csv", quote = F, row.names = F)
write.csv(outMSI3[[1]], file = "CD8_excluded_MSI_Metal_P_value_per_marker.csv", quote = F, row.names = F)
pdf("CD8_Metal_memoryMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI0[[2]]){print(plt)}
dev.off()
pdf("CD8_Metal_naiveMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI1[[2]]){print(plt)}
dev.off()
pdf("CD8_Metal_excluded_memoryMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI2[[2]]){print(plt)}
dev.off()
pdf("CD8_Metal_excluded_naiveMSI.pdf", width = 1.6, height = 3.5)
for (plt in outMSI3[[2]]){print(plt)}
dev.off()




