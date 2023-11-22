library(CellChat)
library(patchwork)
library(ComplexHeatmap)

cellchat_WT <- readRDS(file = "C:/Users/hyjfo/Documents/WT_cellchat_1000_test.rds")
cellchat_Mga <- readRDS(file = "C:/Users/hyjfo/Documents/Mga_cellchat_1000_test.rds")
object.list <- list(WT = cellchat_WT, Mga = cellchat_Mga)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

simple_colors=c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896')

compareInteractions(object=cellchat, measure = c("count"), color.use =c('Blue', 'Red'), group = c(1,2), show.legend = FALSE)
compareInteractions(object=cellchat, measure = c("weight"), color.use =c('Blue', 'Red'), group = c(1,2), show.legend = FALSE)
netVisual_diffInteraction(object=cellchat, comparison = c(1, 2), measure = c("count"), color.use =simple_colors, color.edge = c("Red", "Blue"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, weight.scale = TRUE)
netVisual_diffInteraction(object=cellchat, comparison = c(1, 2), measure = c("weight"), color.use =simple_colors, color.edge = c("Red", "Blue"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, weight.scale = TRUE)
netVisual_heatmap(object=cellchat, comparison = c(1, 2), measure = c("count"), signaling = NULL, slot.name = c("netP"), color.use =simple_colors, color.heatmap = c("Blue", "Red"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
netVisual_heatmap(object=cellchat, comparison = c(1, 2), measure = c("weight"), signaling = NULL, slot.name = c("netP"), color.use =simple_colors, color.heatmap = c("Blue", "Red"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], color.use =simple_colors, title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

netAnalysis_signalingChanges_scatter(object=cellchat, idents.use = c('Fibroblasts'), color.use =c('Gray','Blue', 'Red'), comparison = c(1, 2), signaling = NULL, slot.name = "netP", font.size = 8)

cellchat <- computeNetSimilarityPairwise(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), nCores = 32)
netVisual_embeddingPairwise(object=cellchat, slot.name = "netP", type = "functional", color.use =simple_colors, comparison = c(1, 2), dot.size = c(4, 8), label.size = 4)
netVisual_embeddingPairwiseZoomIn(object=cellchat, slot.name = "netP", type = "functional", color.use =simple_colors, comparison = c(1, 2), nCol = 2, dot.size = c(4, 8), label.size = 4)

cellchat <- computeNetSimilarityPairwise(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), nCores = 32)
netVisual_embeddingPairwise(object=cellchat, slot.name = "netP", type = "structural", color.use =simple_colors, comparison = c(1, 2), dot.size = c(4, 8), label.size = 4)
netVisual_embeddingPairwiseZoomIn(object=cellchat, slot.name = "netP", type = "structural", color.use =simple_colors, comparison = c(1, 2), nCol = 2, dot.size = c(4, 8), label.size = 4)

rankSimilarity(object=cellchat, slot.name = "netP", type = "functional", comparison2 = c(1, 2), font.size = 8)

rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('Blue', 'Red'), stacked = TRUE, sources.use = NULL, targets.use = NULL, signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)
rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('Blue', 'Red'), stacked = FALSE, sources.use = NULL, targets.use = NULL, signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object=object.list[[i]], signaling = pathway.union, pattern = "outgoing", slot.name = "netP", color.use =simple_colors, color.heatmap = "BuGn", title = names(object.list)[i], font.size = 8, font.size.title = 8)
ht2 = netAnalysis_signalingRole_heatmap(object=object.list[[i+1]], signaling = pathway.union, pattern = "outgoing", slot.name = "netP", color.use =simple_colors, color.heatmap = "BuGn", title = names(object.list)[i+1], font.size = 8, font.size.title = 8)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object=object.list[[i]], signaling = pathway.union, pattern = "incoming", slot.name = "netP", color.use =simple_colors, color.heatmap = "BuGn", title = names(object.list)[i], font.size = 8, font.size.title = 8)
ht2 = netAnalysis_signalingRole_heatmap(object=object.list[[i+1]], signaling = pathway.union, pattern = "incoming", slot.name = "netP", color.use =simple_colors, color.heatmap = "BuGn", title = names(object.list)[i+1], font.size = 8, font.size.title = 8)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = FALSE, color.text.use = FALSE, angle.x=45)
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = c("MHC-II"), pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = FALSE, color.text.use = FALSE, angle.x=45)
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 2, color.text.use = FALSE, title.name = "Increased signaling in Mga", angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = c("MHC-II"), pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 2, color.text.use = FALSE, title.name = "Increased signaling in Mga", font.size = 20, angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 1, color.text.use = FALSE, title.name = "Decreased signaling in Mga", angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = c("MHC-II"), pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 1, color.text.use = FALSE, title.name = "Decreased signaling in Mga", font.size = 20, angle.x = 45)

cellchat <- identifyOverExpressedGenes(object=cellchat, group.by = NULL, group.dataset = "datasets", pos.dataset = "Mga", features.name = "Mga", only.pos = FALSE, thresh.pc = 0, thresh.fc = 1, thresh.p = 0.05)
net <- netMappingDEG(object=cellchat, features.name = "Mga", thresh= 0.05)
net.up <- subsetCommunication(object=cellchat, net = net, slot.name = "net", sources.use = NULL, targets.use = NULL, signaling = NULL, thresh = 0.05, datasets = "Mga", ligand.pvalues = 0.05, ligand.logFC = 1, receptor.pvalues = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(object=cellchat, net = net, slot.name = "net", sources.use = NULL, targets.use = NULL, signaling = NULL, thresh = 0.05, datasets = "WT", ligand.pvalues = 0.05, ligand.logFC = -1, receptor.pvalues = 0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(pairLR=net.up, object = cellchat)
gene.down <- extractGeneSubsetFromPair(pairLR=net.down, object = cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = FALSE]
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = NULL, pairLR.use = pairLR.use.up, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, title.name = paste0("Up-regulated signaling in Mga", names(object.list)[2]), angle.x = 45)
pairLR.use.down = net.down[, "interaction_name", drop = FALSE]
netVisual_bubble(object=cellchat, sources.use = c(1:6), targets.use = c(7:13), signaling = NULL, pairLR.use = pairLR.use.down, color.heatmap = c("viridis"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, title.name = paste0("Down-regulated signaling in Mga", names(object.list)[2]), angle.x = 45)
netVisual_chord_gene(object=object.list[[2]], slot.name = "net", signaling = NULL, pairLR.use = NULL, net = net.up, sources.use = c(1:6), targets.use = c(7:13), lab.cex = 0.5, title.name = paste0("Up-regulated signaling in Mga", names(object.list)[2]), thresh = 0.05)
netVisual_chord_gene(object=object.list[[1]], slot.name = "net", signaling = NULL, pairLR.use = NULL, net = net.down, sources.use = c(1:6), targets.use = c(7:13), lab.cex = 0.5, title.name = paste0("Down-regulated signaling in Mga", names(object.list)[1]), thresh = 0.05)

weight.max <- getMaxWeight(object.list=object.list, slot.name = c("netP"), attribute = c("MHC-II")) # control the edge weights across different datasets

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object=object.list[[i]], signaling= c("MHC-II"), signaling.name = paste(c("MHC-II"), names(object.list)[i]), color.use =simple_colors, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = 1, weight.scale = TRUE, edge.weight.max = weight.max[1], edge.width.max = 10, layout = c("circle"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object=object.list[[i]], comparison = c(1, 2), measure = c("count"), signaling = c("MHC-II"), slot.name = c("netP"), color.use =simple_colors, color.heatmap = c("Reds"), title.name = paste(c("MHC-II"), "signaling ",names(object.list)[i]), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
}
draw(ht[[1]] + ht[[2]], ht_gap = unit(1, "cm"))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "Mga")) # set factor level
plotGeneExpression(object=cellchat, signaling = c("MHC-II"), enriched.only = TRUE, type = c("violin"), split.by = "datasets",colors.ggplot = TRUE)
