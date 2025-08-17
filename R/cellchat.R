library(CellChat)
library(patchwork)
combined = readRDS("730_combined_seurat.rds")
norm_data = combined@assays$RNA@layers$counts
colnames(norm_data) = rownames(combined@assays$RNA@cells)
rownames(norm_data) = rownames(combined@assays$RNA@features)
meta = combined@meta.data
kidney_metadata = meta[which(meta$tissue == "kidney"),]
liver_metadata = meta[which(meta$tissue == "liver"),]
shared_celltype = intersect(kidney_metadata$cell_type, liver_metadata$cell_type)
meta[which(meta$cell_type %in% shared_celltype),"cell_type"] = 
  paste(meta[which(meta$cell_type %in% shared_celltype),"tissue"], 
        meta[which(meta$cell_type %in% shared_celltype),"cell_type"], sep=" ")

combined_cellchat = createCellChat(object=norm_data, meta=meta, group.by = "cell_type")
combined_cellchat@DB = CellChatDB.human

#################### 加载PPI到cellChat数据库 #############
# creat object
PPI = read.csv("modified_hepatokine_PPI_0.8.csv")
db = combined_cellchat@DB
interaction = data.frame(
  interaction_name = paste(PPI$protein1, PPI$protein2, sep="_"), 
  ligand = PPI$protein1, 
  receptor = PPI$protein2
)
interaction = unique(interaction)
rownames(interaction) = interaction$interaction_name

# combine
merged_interaction = interaction[which(interaction$interaction_name %in% db$interaction$interaction_name == F),]
merged_interaction = merge(merged_interaction, db$interaction, by=c("interaction_name","ligand","receptor"), all=T)

# PPI and CellChatDB.human
library(dplyr)
interaction_info = merged_interaction[,1:4]
rownames(interaction_info) = interaction_info$interaction_name
interaction_info[interaction$interaction_name,"info"] = "PPI"
interaction_info[db$interaction$interaction_name,"info"] = "CellChat DB"
interaction_info[which(interaction_info$ligand %in% intersect(interaction$ligand,db$interaction$ligand)),"info"] = "both"
write.csv(interaction_info, file="/HDD_01/jhy/yumo/interaction_info.csv")

# 
combined_cellchat@DB$interaction = merged_interaction
saveRDS(combined_cellchat,"/HDD_01/jhy/yumo/cellchat.rds")

################### Pre-processing ###############
combined_cellchat = subsetData(combined_cellchat)
combined_cellchat = identifyOverExpressedGenes(combined_cellchat)
combined_cellchat = identifyOverExpressedInteractions(combined_cellchat)

################### cell-cell communication ################
# ligand-receptor
combined_cellchat = computeCommunProb(combined_cellchat, raw.use=T, population.size=T)   # type = "truncatedMean", trim=0.1

saveRDS(combined_cellchat,"/HDD_01/jhy/yumo/cellchat.rds")

combined_cellchat = filterCommunication(combined_cellchat, min.cells=10)

# pathway
combined_cellchat = computeCommunProbPathway(combined_cellchat)
combined_cellchat = aggregateNet(combined_cellchat)

all_net = subsetCommunication(combined_cellchat)
write.csv(net, "/HDD_01/jhy/yumo/all_net.csv")

liver_kidney_net = all_net[grep("Liver",all_net$source),]
liver_kidney_net = liver_kidney_net[grep("Kidney", liver_kidney_net$target),]
liver_kidney_net$interaction_name = paste(liver_kidney_net$ligand, liver_kidney_net$receptor, sep="_")
write.csv(liver_kidney_net, "/HDD_01/jhy/yumo/liver_kidney_net.csv")

# 含hepatokine的net
ligand = unique(PPI$protein1)
hepatokine_net = data.frame()
for (i in ligand) {
  hepatokine_net = rbind(hepatokine_net, liver_kidney_net[which(liver_kidney_net$ligand == i),])
}
hepatokine_net = unique(hepatokine_net)
hepatokine_net[,"cell_cell"] = paste(hepatokine_net$source, hepatokine_net$target, sep = "_")
write.csv(hepatokine_net, "/HDD_01/jhy/yumo/hepatokine_net.csv")

saveRDS(combined_cellchat, "/HDD_01/jhy/yumo/combined_cellchat.rds")

sorted_net = target_new_net[which(target_new_net$pval==0),]

sorted_targed = table(sorted_net$cell_cell) %>% as.data.frame()
sorted_targed = sorted_targed[sorted_targed$Freq>3,]
sorted_targed = sorted_targed[order(sorted_targed$Freq, decreasing=T),]

sorted_net = sorted_net[which(sorted_net$cell_cell %in% sorted_targed$Var1),]
sorted_net$interaction_name_2 = sorted_net$interaction_name



####################### Visualization #################
netVisual_circle(combined_cellchat@net$weight, vertex.weight=as.numeric(table(combined_cellchat@idents)),
                 weight.scale=T, label.edge=F, title.name="Interaction strength")

# bubble plot
netVisual_bubble(combined_cellchat, sources.use=1)


sorted_net$source = as.character(sorted_net$source)
sorted_net$target = as.character(sorted_net$target)
for (i in 1:nrow(sorted_net)) {
  sorted_net[i,"source"] = strsplit(sorted_net[i,"source"],"_")[[1]][2]
  sorted_net[i,"target"] = strsplit(sorted_net[i,"target"],"_")[[1]][2]
}
sorted_net$target[which(sorted_net$target=="kidney epithelial cell")] = "PCT"
sorted_net$source = factor(sorted_net$source)
sorted_net$target = factor(sorted_net$target)


sorted_combined_cellchat = combined_cellchat
sorted_combined_cellchat@net = sorted_net
for (i in nrow(sorted_combined_cellchat@meta)) {
  sorted_combined_cellchat@meta[i,"cell_ontology_class"] = 
    strsplit(sorted_combined_cellchat@meta[i,"cell_ontology_class"],"_")[[1]][2]
}


id = sorted_combined_cellchat@idents
id = as.character(id)
for (i in 1:length(id)) {
  id[i] = strsplit(id[i],"_")[[1]][2]
}
id[which(id =="kidney epithelial cell")] = "PCT"
id = as.factor(id)
sorted_combined_cellchat@idents = id

netVisual_bubble(sorted_combined_cellchat)

# Circle plot
sorted_combined_cellchat = aggregateNet(sorted_combined_cellchat)
celltype = c(unique(paste("Liver", sorted_combined_cellchat@net$source, sep="_")), unique(paste("Kidney", sorted_combined_cellchat@net$target, sep="_")), "Kidney_kidney epithelial cell")
sorted_count = combined_cellchat@net$count[which(rownames(combined_cellchat@net$count) %in% celltype),which(colnames(combined_cellchat@net$count) %in% celltype)]
rownames(sorted_count)[which(rownames(sorted_count) == "Kidney_kidney epithelial cell")] = "Kidney_PCT"
colnames(sorted_count)[which(colnames(sorted_count) == "Kidney_kidney epithelial cell")] = "Kidney_PCT"
for (i in 1:nrow(sorted_count)) {
  rownames(sorted_count)[i] = strsplit(rownames(sorted_count)[i],"_")[[1]][2]
}
for (i in 1:ncol(sorted_count)) {
  colnames(sorted_count)[i] = strsplit(colnames(sorted_count)[i],"_")[[1]][2]
}
netVisual_circle(sorted_combined_cellchat@net$weight, vertex.weight=as.numeric(table(sorted_combined_cellchat@idents)),
                 weight.scale=T, title.name="Interaction strength")

###################### target ligand ####################
target_ligand = L_K_net[which(L_K_net$ligand %in% ligand),]
netVisual_bubble(combined_cellchat, signaling = target_ligand$pathway_name)
