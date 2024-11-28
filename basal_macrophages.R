library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)

metadata <- read.csv("./data/seurat.integrated.5Ht_6Ho.metadata.csv", row.names = 1)
count <- as.matrix(read.csv("./data/seurat.integrated.5Ht_6Ho.counts.csv", row.names = 1))

# Creating the seurat object and putting metadata and the count into the objects. 
seuratObj <- CreateSeuratObject(counts = count)
seuratObj <- AddMetaData(seuratObj, metadata = metadata)
seuratObj <- UpdateSeuratObject(seuratObj)

# Setting the matrices to be our desired data.
# The ligand-target prior model is a matrix describing the potential that a ligand may regulate a target gene, and it is used to run the ligand activity analysis. 
# The ligand-receptor network contains information on potential ligand-receptor bindings, and it is used to identify potential ligands. 
# The weighted ligand-receptor network contains weights representing the potential that a ligand will bind to a receptor, and it is used for visualization.
lr_network <- readRDS("./data/lr_network_mouse_21122021.rds")
ligand_target_matrix <- readRDS("./data/ligand_target_matrix_nsga2r_final_mouse.rds")
weighted_networks <- readRDS("./data/weighted_networks_nsga2r_final_mouse.rds")

# Getting the seurat object that only consists of the cells that we are focused on.
basal_macrophages <- subset(seuratObj, subset = cluster1 %in% c("Mammary epithelial cells-Basal", "Monocytes-macrophages"))

# Normalize the data, not sure if we need to normalize the data.
# But in most cases it will show a better visualization of the data, and less bias downstream since the analysis of nichenet uses normalized data.  
basal_macrophages <- NormalizeData(basal_macrophages, assay = "RNA")

# First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression. 
# Then standardizes the feature values using the observed mean and expected variance. Feature variance is then calculated on the standardized values after clipping to a maximum.
# nFeatures determines how many features we want to look at, could potentially be too many or too few, but can cause overfitting if too large.
basal_macrophages <- FindVariableFeatures(basal_macrophages, assay = "RNA", selection.method = "vst", nfeatures = 2000)


#Need this because get_expressed_genes uses the Idents column to determine the expression of a given cell type.
Idents(seuratObj) <- "cluster1"

# Need to determine sender and reciever cells. Here we are using basal cells as receiver and macrophages as sender.
# because these are both different types of cells does it matter which is which? Can test by switching later.
# Need to find the expressed genes for each cell type.
# basal expressed genes  = 8610
expressed_genes_basal <- get_expressed_genes("Mammary epithelial cells-Basal", seuratObj, pct = 0.05)

# macrophages = 6364 
expressed_genes_macrophages <- get_expressed_genes("Monocytes-macrophages", seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
# Length 315
expressed_receptors_basal <- intersect(all_receptors, expressed_genes_basal)

# Getting all potential ligands that match from the basal cells to the ligand receptor network
# length 871
potential_ligands <- lr_network %>% filter(to %in% expressed_genes_basal) %>% pull(from) %>% unique()

# Use lapply to get the expressed genes of every sender cell type separately here
# Length = 6364
expressed_genes_sender <- expressed_genes_macrophages %>% unlist() %>% unique()

# length = 153
# The intersection of potential ligand interactions that overlap in basal cells and macrophages. 
# Pool of potential ligands that might mediate interactions between basal and macrophages.
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# 8731
de_results <- FindMarkers(
  basal_macrophages,
  ident.1 = "Mammary epithelial cells-Basal",
  ident.2 = "Monocytes-macrophages",
  group.by = "cluster1",
  assay = "RNA",
  slot = "data",
  min.pct = 0.05
) %>% rownames_to_column("gene")

# From this function we can sort the list and see the p_values that exist, and actually see a lot with really high values
#Example: p_val = 0.9996531, Herc4 
de_results_sorted <- de_results %>%
  arrange(desc(p_val))

# 6425 
# Problem with this from tutorial is that we are sorting on the p_val_adj.
# Geneset of interest is from the differentially expressed data. 
# assessing how well these genes explain the expression of the specific genes in the reciever population, so basal cells. 
# This function is only getting rid of ~2000 potential ligands.
geneset_oi <- de_results %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


# background genes are all genes in the receiver cell 
# should have geneset be 5,000-10,000 values larger than the background set
# Problem here: our geneset of interest is almost the same size as our expressed genes in basal cells.
# 6425 : 8610       geneset : Basal Genes
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_basal,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_focused
)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

# Red line indicates the score for the 3rd ligand.
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(3, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity


best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  


