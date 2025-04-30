##' input
scRNA <- readRDS("../2024.4_scATAC/data/palate_epi_MAGIC_periderm_MAGIC_harmony_renamed.rds")


dietPalate <- CreateSeuratObject(scRNA@assays$RNA)
dietPalate@meta.data <- scRNA@meta.data
dietPalate@reductions <- scRNA@reductions
dietPalate@neighbors <- scRNA@neighbors
dietPalate@graphs <- scRNA@graphs


##' output
saveRDS(dietPalate,"process/framework/obj/RNA_E12E14_diet.Rds")
