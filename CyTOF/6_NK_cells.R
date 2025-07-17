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


md = read.csv("Files/FileInfoNKcells.csv") ## Sample annotation 
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
md2 =  dplyr::filter(md, Ignore == "No" & PPID %in%  finalSamples)
### Functions ###

extract_key_cells <- function(data, keymarkers, clus, md2){
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
      gg = ggplot(ggdf, aes(x = Group, y = Prop, color = Group)) +
        geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
        geom_point(aes(color = Group),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.8) +
        scale_fill_manual(values = c('#4363d8','#f58231',"black")) +
        scale_color_manual(values = c('#4363d8','#f58231')) + theme_bw(base_size = 16) +
        xlab(label = "" ) + ylab(label = "% cluster cells") +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 14),
              legend.position = "none")  +  facet_wrap(~Marker) + labs(subtitle = paste(clus))
      mplots[[paste(mark)]]= gg
    }
  }
  return(list(mplots))
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

### Marker+ cells in gated NK sub-population ###
fm = intersect(functional_markers, colnames(NKcd16_Exp_matrix0))

NKcd16_Exp_matrix0$cluster = "NK_CD16pos" ##
out0 = extract_key_cells(NKcd16_Exp_matrix0, keymarkers = c("CCR6","CD27","CCR5") , clus='NK_CD16pos', md2, toxicMetals)
outMSI = MSI_analysis(NKcd16_Exp_matrix0, clus='NK_CD16pos', functional_markers = "CCR5", md2, toxicMetals)
pdf("Supp.Fig.16a.pdf",width = 1.6, height = 3.5)
print(out0[[1]])
print(out0[[2]])
print(out0[[3]])
print(outMSI[[1]])
dev.off()

NKcd56_Exp_matrix0$cluster = "NK_CD56bright"
out0 = extract_key_cells(NKcd56_Exp_matrix0, keymarkers = c("CCR7","PD.1") , clus='NK_CD56bright', md2, toxicMetals)
pdf("Supp.Fig.16b.pdf",width = 1.6, height = 3.5)
print(out0[[1]])
print(out0[[2]])
dev.off()


out_surt0 = extract_key_cells_seuret(NKcd16_Exp_matrix0, fm, clus='NK_CD16pos', md2, toxicMetals)
out_surt1 = extract_key_cells_seuret(NKcd56_Exp_matrix0, fm, clus='NK_CD56bright', md2, toxicMetals)


### Toxic metal analysis ###
out1 = extract_key_cells(data, "79Br" , "NK_CD16pos", md2)
out2 = extract_key_cells(data, "79Br" , "NK_CD56bright", md2)
pdf("Figure4d.pdf", width = 1.6, height = 3.5)
print(out1[[1]][[1]])
print(out2[[1]][[1]])
dev.off()



