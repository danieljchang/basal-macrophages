#for next time: Tell me what all variables are and explain what they do.

library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)

#nGene, nUMI, orig.ident, aggregate, res.0.6, celltype, nCount_RNA, nFeature_RNA
seuratObj <- readRDS("./data/toy_seurat/seuratObj.rds")

seuratObj@meta.data %>% head()
#Updates to meet latest requirements of seurat.Compatibility check.
seuratObj <- UpdateSeuratObject(seuratObj)
#Converts aliases to official gene symbols.
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# Note that the number of cells of some cell types is very low and should preferably be higher for a real application
seuratObj@meta.data$celltype %>% table() 

DimPlot(seuratObj, reduction = "tsne")

seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141

DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")


###############################################################
# Nichenet 
###############################################################

organism <- "mouse"
#Read database, and assign to resepcted matrix for ligand receptor interactions and ligand ligand interactions.
# Our data is mouse data. Convert human -> Mouse data 
# look for ligand receptor data on mouse
# Double check date of the data of the ligand receptor matrices. 
# "omnipath" "nichenet_verschueren"
if(organism == "human"){
  lr_network <- readRDS(url("./data/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("./data/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("./data/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS("./data/lr_network_mouse_21122021.rds")
  ligand_target_matrix <- readRDS("./data/ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- readRDS("./data/weighted_networks_nsga2r_final_mouse.rds")
} 

#lr_network is the ligand to receptor matrix, telling us the from ligand to receptor.
# a ligand and its specific receptor. 
# also stores where the source is collected from.
# Has 'from', 'to', 'database', 'source'
lr_network <- lr_network %>% distinct(from, to)

head(lr_network)

# target genes in rows, ligands in columns
# ligands with target genes/all possible down stream genes. 
# potential ligand-receptor bindings
ligand_target_matrix[1:5,1:5]

# interactions and weights in the ligand-receptor and signaling network
# Like the ligand receptor matrix but now we have weights given the strength of interactions.
# describing the potential that a ligand may regulate a target gene
head(weighted_networks$lr_sig) 

# weights representing the potential that a ligand will bind to a receptor.
head(weighted_networks$gr) 




###############################################################
# Step 1: Define a set of potential ligands
###############################################################

# The receiver cell population can only consist of one cell type, so in case of multiple receiver populations, 
# will have to rerun the vignette separately for each one. 
# Only looking at CD8 T cell receivers. 
# CD8 T: a type of T cell that play a key role in the immune system's defense against infections and cancer
receiver = "CD8 T"

# Returns the a list of genes that are expressed in a given cell cluster based on the 
# fraction of cells in that cluster that should express the cell.
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

# A list of unique strings of receptor genes, from the ligand receptor network.
all_receptors <- unique(lr_network$to)  

# A list of receptors that are expressed from the fraction of cells in that cluster.
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

# potential ligands as all ligands whose cognate receptors are expressed according to the clustering.
# Length = 475
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# autocrine signaling if we also include CD8 T. All cell types are pooled together, later filtered by which ligand was expressed.
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")

# Use lapply to get the expressed genes of every sender cell type separately here
# a list of lists, Example: list(Sender1 = c(gene1,gene3,gene4...), sender2 = c(gene1,gene2,gene3),...)
# length = 6, but each list has many gene senders,
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)

# Flattens the list of lists, and then only takes unique expressed values. 
# Length = 8492
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

# By taking the intersection of only the sender_celltypes, and the potential ligands, we can find
# correlated ligands, that would be targetted for analysis.
# Length = 122
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 



###############################################################
# Step 2 Define the gene set of interest
###############################################################

#gene set of interest is determine by the most likely to be influenced by the cell-cell communication.
# The conditions come from the aggregate column in metadata.
# unique(seuratObj@meta.data['aggregate'])
# In this dataset, only have LCMV which is the condition of interest and SS the steady state condition.
condition_oi <-  "LCMV"
condition_reference <- "SS"

# Seurat object based on the identity class receiver = CD8 T, subset of only this cell type.
seurat_obj_receiver <- subset(seuratObj, idents = receiver)

# Using Wilcoxon Rank Sum Test
# Compares the gene expression of LCMV and SS, using a 5% minimum expression of genes to be included.
# Identifies significantly differently expressed CD8 T cells between LCMV and SS conditions.
# Columns: gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
# avg_log2FC: Tells us the log ratio of the differences between the two conditions
# pct.1, pct.2: Percentage of cells expressing the gene, LCMV or SS, respectively. 
# p_val_adj: reduce the chance of incorrectly declaring statistical significance, minimize FP
# 166 rows
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "aggregate",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

# Only includes statistically significant and biologically relevant. 
# Length= 322
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
# filter only the ligand-target gene matrix that are also in our gene set of inteterest.
# Length = 241
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Depending on our dataset, we could be seeing a lot of gene interactions that are 
# statistically significant but not in the ligand-receptor network

###############################################################
# Step 3: Define background genes
# All expressed genes in the receiver cell population is defined as the ‘background set’
# Needs to be "significantly" larger than gene set of interest.
###############################################################

# Line 80: expressed_genes_receiver
# Filtering from the ligand_ligand_target matrix,
# length = 3182
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


###############################################################
# Step 4: Nichenet ligand activity analysis
###############################################################


# Predict_ligand_activities determines how well nichenetr can predict the observed transcriptional response
# This prediction is comparing our geneset of interest to the background expressed genes.
# presence of their target genes in the gene set of interest 
# Columns: test_ligand, auroc, aupr, aupr_corrected, pearson
# 475 rows
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

# Prioritized ligands that induced an antiviral response to CD8 T cells.
# Sorting on -aupr_corrected and adding a column that ranks based on the aupr_corrected
# Can sort on other metrics but they found AUPR is the most informative.
# Columns: test_ligand, auroc, aupr, aupr_corrected, pearson, rank
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

# BEWARE: it is possible that for some gene sets, the target gene prediction performance of the top-ranked ligands would not be much better than random prediction.
# But here we can see a correlation that the top 25 can predict reasonably the response
# top_n should be changed based on statistical significance (Arbitrary choice).
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(25, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

# A list of top 25 based on aupr_corrected and sorting them.
best_upstream_ligands <- ligand_activities %>% top_n(25, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#Creating a dataframe with the test ligand names and the aupr_corrected scores, both from ligand_activity df.
# Sorted from lowest to highest
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

# Plotting a heatmap using ggplot of the sorted list of prioritized ligands based on aupr.
(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

###############################################################
# Step 5: Infer target genes and receptors of top-ranked ligands
###############################################################


# Creating a dataframe with columns, ligand, target, weight.
# Rows = 511
# Infer active ligand target links between possible ligands and genes belonging to a gene set of interest:
# consider the intersect between the top n targets of a ligand and the gene set.
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# regulatory potential scores between all pairs of ligands and targets documented in this tibble
# 74 x 25 table, that maps out the potential scores of ligand and target gene.
# preparing for a heatmap, quantile cutoff on the ligand-target scores
# if we view(active_ligand_target_links) we see 0's across a lot of inputs.
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

# Get the ligands in order and target genes in order.
# 24 ligands, columns
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
# 74 target genes, rows
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

#Transposing the matrix now rows represent the ligands, and columns represent the target genes
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

# can't find function online
# Made a 44x3 tibble, columns: from, to, weight
# Using the weighted matrix data from the ligand receptor and signaling network
# and filtering only ligands and receptors that we have deemed to be relavent.
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

# Making a matrix of the ligand receptor links df and our best upstream ligands.
# Hierarchical clustering using both rows and columns
# 7x25 matrix rows = ligands, columns = receptors.
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prioritized interaction potential"))

###############################################################
# Step: 6 Sender focused
###############################################################


# tibble with columns: test_ligand auroc  aupr aupr_corrected pearson  rank
# 475 rows.
ligand_activities_all <- ligand_activities 
# list of 25 upstream ligands
best_upstream_ligands_all <- best_upstream_ligands

#Doing the same thing as step 4. but instead  ligand activities only containing expressed ligands from all populations 
# 475x6 col namess: test_ligand, auroc, aupr, aupr_corrected, pearson, rank
ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

# the ligand of interest sorted from low to high on aupr_corrected 25 ligands, dataframe
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)

# turning the dataframe into a matrix
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

# Generating heatmap based on the matrix.
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Same as step 5 Visualization of ligand and target genes based on the first 100 most relavent/active.
# tibble, 511 x 3 , ligand target weight.
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# 74 x 25 matrix 
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

# Sorted order of ligands of most interest. length 25
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
# length of 74
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

# Transposing the matrix now rows represent the ligands, and columns represent the target genes
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

# plotting the heatmap of the ligand target matrix, based on predicted activity.
p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target


# Receptor plot
# Creates a dataframe with ligands and receptors in addition to the weights of these links
# From, To, Weights, 44 rows.
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
# Rows ligand, columns receptors. 
# 7 x 25, 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

# This actually tells us how many of the upstream ligands are interacting with target genes, focsued on the sender.
# We find many of the senders to be falsly correlated. since IFN genes are not expressed by the sender cell populations, 
# and it was already filtered out during preprocessing for being too lowly expressed.
## FALSE  TRUE 
##    20     5
best_upstream_ligands_all %in% rownames(seuratObj) %>% table()


# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot


# Gives all the different cell types in the seurat object
# "CD8 T" "CD4 T" "Treg"  "B"     "NK"    "Mono"  "DC"   
celltype_order <- levels(Idents(seuratObj)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
# Used to find the Log Fold Change for specific ligands in each cell type that we are considering.
# how the ligands behave between steady state vs disease. SS vs LVMC
# 30x2 tibble, correlating gene to a cell type, 6 tables
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "aggregate",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

# merging the tibbles together so that it is now a 30x6
DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

# Sorting and turning it into a matrix. based on most active ligands.
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])


# heatmap to show log full change in sender vs our prioritized ligands of interest.
p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

# shows the ligand activity and the distribtuion of expressed ligands across all sender ligands.
(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))


# Combining all the plots together into a grid. Allows for visualization of every
# graphical analysis in one place.

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot



active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
## [1] 637
head(active_ligand_target_links_df)


# Upload and share code/plots. Clean up code, switch ligand and receptor sender receiver cells.

