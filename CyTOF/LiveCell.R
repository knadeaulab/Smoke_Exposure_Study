library(patchwork)
library(Seurat)
library(MAST)
library(FlowSOM)
library(flowCore)
library(flowViz)
library(flowStats)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mclust)
library(Seurat)
library(umap)
library(Biobase)
library(preprocessCore)

set.seed(123)
dir.create("Data")

#########################################################################################################
## Download Live cell FCS files and R-session from Harvard Dataverse Repository in the 'Data' Directory #
#########################################################################################################
## NOTE: To reproduce the results. We will use R session file which contains exactly the same 
## randomly selected subset of cells used in the study. However, I have also included R codes 
## I used to produce the variables available in Rsession file. 
###########################################################################################

############################# Session-1 ##############################
### Section-1 Preprocess live cells before unsupervised clustering ###
## Note: Skip to session-2 if you want to work on pre-processed data #
######################################################################

md = read.csv("Files/FileInfo.csv") 
md$Fpath = paste("Data/", md$FileName, sep="")
finalSamples = c("1294","1804","2765","2766","5023","5024","5865","5866","FF-024","WF046","WF135","FF-013",
  "WF034","FF-035","FF-022","WF008","WF055","WF045","WF128","WF110","WF081","FF-052","FF-023",
  "WF015","WF057","FF-045","WF023","WF002","WF005","WF042","WF051","WF094","WF102","WF120",
  "FF-042","FF-026","FF-027","FF-028","FF-029","FF-031","FF-033","FF-041","FF-043","FF-016",
  "FF-018","FF-020","FF-021","WF116","FF-008","FF-014","FF-019","FF-038","W0043B","WF066","W0050A",
  "WF006","W0053A","WF069","WF025","WF138")
md = dplyr::filter(md, PPID %in% finalSamples & Internal_control!= "IC") ## IC: Internal Controls; Set to == "IC", if you want to create density plot for internal controls only (Fig. S18a)
color_conditions <- c('#4363d8','#f58231')
names(color_conditions) <- levels(factor(md$Group))

### Reading FCS files ###
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
          "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd") ## includes both protein and toxic metals

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
flow_set_1 <- flowSet(flowframes) ## can be used to plot supp.Fig. 19a
ExpressionMatrix = fsApply(flow_set_1, exprs)
PIDs = do.call(c, IDs)

### FlowSom clustering   ###
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
idx =  which(PIDs %in% md$SampleID)
expr2 = ExpressionMatrix[idx,]
rownames(expr2) = 1:length(idx)
sampleIDs00 = as.character(PIDs[idx])
idx = subSample(sampleIDs00,80000)
data = expr2[idx,]
sampleIDs = as.character(sampleIDs00[idx])
lm = c("CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA","TCRgd", "CD66b", "CD25", "CD56", "CD8a", "CD28","CD14", "CD127","CD19")
lineage_markers = intersect(lm, colnames(data))
message("creating a flowset..")
fcs.input <- new("flowFrame", exprs=as.matrix(data[,lineage_markers]))
message("Running FlowSOM..")
flowSOM.res <- FlowSOM(fcs.input, colsToUse = lineage_markers, scale = FALSE ,nClus=10) ## 10 clusters
message("generating Meta-Clusters..")
labels <- GetMetaclusters(flowSOM.res)
data = as.data.frame(data)
data$cluster = labels
data$sample_id = sampleIDs

### UMAP of live cells from internal controls ###
subSampleByBatch <- function(sample_ids,ncount=4000, md2){ 
  mm = match(sample_ids, md2$SampleID)
  batches = md2$Batch[mm]
  df = data.frame(Index = 1:length(sample_ids), IDs = sample_ids, Batch = batches )
  ##df_unique <- df %>% distinct(Batch, IDs, .keep_all = TRUE)
  fsom_inds <- df %>%
    group_by(IDs, Batch) %>%
    sample_n(ncount)
  return(as.data.frame(fsom_inds))
}
md2 = dplyr::filter(md, Internal_control == "IC" )
idx =  which(PIDs %in% md2$SampleID)
expr2 = ExpressionMatrix[idx,]
rownames(expr2) = 1:length(idx)
idx = subSampleByBatch(sample_ids=PIDs[idx],5000, md2)
subsample.expr = expr2[idx$Index,]
SampleID = idx$IDs
lm = c("CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA","TCRgd", "CD66b", "CD25", "CD56", "CD8a", "CD28","CD14", "CD127","CD19")
lineage_markers = intersect(lm, colnames(subsample.expr))
list.out = drawUMAP(df=subsample.expr[,lineage_markers], n_neighbors=15, min_dist=0.5, n_components=2, metric = "euclidean",t="Umap",lineage_markers,v="Cell types") 


############################# Session-2 #############################
### Section-2 Plot clustering data with UMAP, Heatmap & box-plots ###
######################################################################

load("Data/Live_cells.RData")

### Heatmap ###
lm = c("CD123", "CD11c","CD16", "CD3", "CD4", "CD45RA","TCRgd", "CD66b", "CD25", "CD56", "CD8a", "CD28","CD14", "CD127","CD19")
pdf7 = as.data.frame(data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarise_all(median))
nam = data[,c(lm,"cluster")] %>% dplyr::group_by(cluster) %>% summarize(count = n(), proportion = (n() / nrow(data))*100)
rownames(pdf7) = paste(rownames(pdf7), " (", format(nam$proportion, digits=2) , "%)", sep="")

pdf("Supp.Fig8a.pdf", width = 8, height = 5)
pheatmap(pdf7[,lm], cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "#434CC3", "navy", "orange", "firebrick3"))(50))
dev.off()

### Proportion of live cells in each cluster ####
md2 = dplyr::filter(md, Ignore == "No" )
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

pdf("Supp.Fig8b.pdf", width = 10, height = 11)
wrap_plots(plotlist = bplots, ncol = 6, nrow = 3)
dev.off()

### UMAP of internal controls ###
ggdf = as.data.frame(list.out[[1]]$layout)
colnames(ggdf) = c("UMAP1","UMAP2")
ggdf$SID = SampleID
mm = match(ggdf$SID, md2$SampleID)
ggdf$Batch =md2$Batch[mm]
ggdf = dplyr::filter(ggdf, UMAP1 < 40 & UMAP2 > -10)
g <- ggplot(ggdf,  aes(x = UMAP1, y = UMAP2, color = Batch)) +
  geom_point(size = 0.9, alpha = 0.35) + theme_bw() + 
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=16),
        panel.grid = element_blank())+
  scale_color_manual(values = c("#1F78B4", "#FF7F0E", "#33A02C", "#CDC50E", "#6A3D9A", "#A52A2A", "green") )

pdf("Supp.Fig19B.pdf", height = 5.2, width = 7.5)
print(g)
dev.off()




######### Session-3 #############
### Toxic metal analysis ########
#################################

toxicMetals= c( "196Hg", "198Hg", "202Hg","182W","184W", "79Br","81Br","204Pb","209Bi","206Pb","208Pb",
                "186W","180Ta", "181Ta", "121Sb", "123Sb","114Cd", "112Cd","110Cd","106Cd", "108Cd")

### Density Plot ###
non_zeros_grouped <- data2[,c(toxicMetals, "sample_id")] %>% group_by(sample_id) %>% summarise(across(everything(), ~ sum(. != 0), .names = "{col}"))
ggdf = melt(non_zeros_grouped)
colnames(ggdf)[2] = "Metal"
ggdf$logCell_count = log(ggdf$value)
ggdf$Metal = factor(ggdf$Metal, levels = rev(c("79Br","81Br","209Bi","204Pb", "206Pb","208Pb","196Hg","198Hg","202Hg","182W","184W",
                                               "186W","180Ta","181Ta","121Sb","123Sb","106Cd","108Cd","110Cd","112Cd","114Cd")))
pdf("Figure4a.pdf", height = 3.8, width = 3.9)
ggdf00 = ggdf[grepl("^[0-9]", ggdf$Metal),]
ggplot(ggdf,aes(x = logCell_count , y = Metal, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "logCell_count", option = "C") +
  labs(title = 'Heavy Metal (live cells)') +
  theme_bw()
dev.off()



## Box-plots ##
extract_key_cells <- function(data, keymarkers=fm, clus="Live", md){
  idx = grep(pattern = "sample", colnames(data), ignore.case = T)
  colnames(data)[idx] = "sample_id"
  hcless = dplyr::filter(data, cluster %in% clus)
  mplots = list()
  BTplot = list()
  
  msioutput = data.frame()
  tempout = data.frame()
  tmpout = data.frame()
  for (mark in keymarkers){
    message(mark)
    X = as.matrix(hcless[mark])
    quan = 0.95
    tmpdf1 = data.frame(Value = X, sample_id = hcless$sample_id)
    colnames(tmpdf1) = c("Value", "sample_id")
    SID = tmpdf1$sample_id
    quantiles_by_sample <- tmpdf1 %>% group_by(sample_id) %>% 
      summarise(quantile_90 = quantile(Value, quan, na.rm = TRUE))
    tmpdf1$sample_id = SID
    # Filter for values that exceed the quantile for each sample id sepratelt 
    prop <- tmpdf1 %>% left_join(quantiles_by_sample, by = "sample_id") %>%
      group_by(sample_id) %>% summarise(count_above_90th_quantile = sum(Value > quan, na.rm = TRUE) )
    prop$Batch=NULL
    
    # Step 2: Filter rows where the 'marker' value exceeds the quantile for its batch_id
    tmpdf2 <- tmpdf1 %>% left_join(quantiles_by_sample, by = "sample_id") %>% dplyr::filter(Value > quantile_90)
    tmpdf2$quantile_90=NULL
    tmpdf2$Batch=NULL
    colnames(prop) = c("sample_id", "Freq")
    Total =  as.data.frame(table(hcless$sample_id))
    mm = match(prop$sample_id, Total$Var1)
    prop$Total = Total$Freq[mm]
    prop$Prop = (prop$Freq/prop$Total)*100
    colnames(prop)[1] = c("sample_id")
    mm = match(prop$sample_id, md$SampleID)
    prop$Group = md$Group[mm]
    prop$Group = factor(prop$Group, levels = c("nonSE", "SE"))
    prop$Batch = md$Batch[mm]
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
  return(list(mplots))
}
fm = intersect(toxicMetals, colnames(data))
out = extract_key_cells(data, fm , clus = "Live", md)
#MOI =  c("106Cd", "108Cd", "110Cd", "79Br","180Ta") ##
pdf(file ="Figure4B.pdf", width = 1.6, height = 3.5) ## Metals of interest was selected for figure ###
for (plt in out[[1]]){
  print(plt)}
dev.off()


## seurat significance analysis ##
extract_key_cells_seuret <- function(data, keymarkers, clus, md){
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
  mm = match(hcless$sample_id, md$SampleID)
  hcless$Group = md$Group[mm]
  hcless$Batch = md$Batch[mm]
  metadata <- data.frame(
    CellID = rownames(hcless),
    Group =  hcless$Group,
    Batch = hcless$Batch)
  seurat_object <- AddMetaData(seurat_object, metadata)
  Idents(seurat_object) <- seurat_object@meta.data$Group
  seurat_object <- SetIdent(seurat_object, value = "Group") 
  
  # Example of differential expression analysis using MAST
  differential_markers<-FindMarkers(seurat_object,
    ident.1 = "nonSE",
    ident.2 = "SE",
    test.use = "MAST",
    latent.vars = "Batch"
  )
  return(differential_marker)
}
## Warning: CPU intensive ##
out_fdr = extract_key_cells_seuret(data, fm, clus="Live", md)



