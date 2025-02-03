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

dir.create("Data")

############################################################################################
## Download NKcells_Handgated_subsets.Rdata file from Harvard Dataverse Repository #########
## (URL: https://doi.org/10.7910/DVN/5PWZJM) in the 'Data' Directory #
###########################################################################################

md = read.csv("Files/FileInfoNKcells.csv") ## Sample annotation 
cd3 = read.csv("Files/HandGated/CD3neg_cell_count.csv") ## CD3- cell count per sample 
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

load("Data/NKcells_Handgated_subsets.Rdata")
###############################
######### Functions ###########
###############################

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
              legend.position = "none")  +  facet_wrap(~Marker) + labs(subtitle = paste(pval,"|",clus))
      mplots[[paste(mark)]]= gg
    }
  }
  return(mplots)
}
extract_key_cells_seuret <- function(data, keymarkers, clus, md2, toxicMetals){
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

prop0 = as.data.frame(table(NKcd56_Exp_matrix0$SampleID))
colnames(prop0) = c("SID","count")
prop0$CellType = "CD56bright"
prop1 = as.data.frame(table(NKcd16_Exp_matrix0$SampleID))
colnames(prop1) = c("SID","count")
prop1$CellType = "CD56dimCD16pos"

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
  xlab(label = "" ) + ylab(label = "% CD3- cells") +
  facet_wrap(~CellType, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14), 
        legend.position = "none")  +
  labs(subtitle = paste("NK cells | P-value", Pvalue_interaction) )

pdf("NKcells_Boxplot_prop_CD56bright_vs_CD56dimCD16pos.pdf", width = 3.0, height = 4)
print(gg)
dev.off()

#################################################
##### Marker+ cells in gated sub-population #####
#################################################
fm = intersect(functional_markers, colnames(NKcd16_Exp_matrix0))

NKcd16_Exp_matrix0$cluster = "NK_C1" ##
out0 = extract_key_cells(NKcd16_Exp_matrix0, fm , clus='NK_C1', md2, toxicMetals)
out_surt0 = extract_key_cells_seuret(NKcd16_Exp_matrix0, fm, clus='a', md2, toxicMetals)

NKcd56_Exp_matrix0$cluster = "NK_C2"
out1 = extract_key_cells(NKcd56_Exp_matrix0, fm , clus="NK_C2", md2)
out_surt1 = extract_key_cells_seuret(NKcd56_Exp_matrix0, fm, clus='a', md2)

pdf("CD56dimCD16pos.pdf", width = 1.6, height = 3.5)
for (plt in out0){print(plt)}
dev.off()

pdf("CD56BrightCD16neg.pdf", width = 1.6, height = 3.5)
for (plt in out1){print(plt)}
dev.off()
