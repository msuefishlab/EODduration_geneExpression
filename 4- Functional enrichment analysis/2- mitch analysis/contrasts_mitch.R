
### code to 1) run the mitch analysis, 2) identify the gene sets of interest, and 3) find the genes common to these gene sets and the DEG genes from edgeR ###

library(plyr)
library(dplyr)
library(mitch)
library(gplots)
library(writexl)
sessionInfo <- sessionInfo()

### Part 1: Mitch analysis 

# mitch v 1.8.0

## load genesets files
geneset.BP <-gmt_import("gmt_GeneSCF/GO_BP_gid.txt")

# import the Bbrach_genome_v0.4 to zebrafish "translator" file
gene_trans <- read.table("../07_Annotation/R/out/geneIDs_Bbrach_Drerio.txt", sep = '\t', header = TRUE, quote = "", fill = FALSE, colClasses = "character")

## import edgeR results to R and join them in a list
T1day_vs_control <- read.csv("../06_DGE/RSEM.gene.counts.matrix.T_1day_vs_control.edgeR.DE_results", sep = '\t')
T1day_vs_T8day <- read.csv("../06_DGE/RSEM.gene.counts.matrix.T_1day_vs_T_8day.edgeR.DE_results", sep = '\t')
control_vs_T8day <- read.csv("../06_DGE/RSEM.gene.counts.matrix.T_8day_vs_control.edgeR.DE_results", sep = '\t')

#rearrange sampleA,B columns and change sign of logFC in the edgeR results that require this, in order to have sampleA be the treatment that got more testosterone-- this will make results more intuitive: upregulated genes in treatment with more testosterone will have positive logFC.  In this case one edgeR result requires this: T1day_vs_T8day
T1day_vs_T8day <- T1day_vs_T8day[, c("sampleB", "sampleA", "logFC",  "logCPM",  "PValue",  "FDR")]
colnames(T1day_vs_T8day) <- c("sampleA", "sampleB", "logFC",  "logCPM",  "PValue",  "FDR")
T1day_vs_T8day$logFC <- -1*T1day_vs_T8day$logFC

##check that the sign of logFC indeed agrees with sample A, B (positive values indicate upregulation in sampleA)

### Multidimensional analysis
# put all results in a list
allresults <-list("ctrl_vs_T8day"=control_vs_T8day, "ctrl_vs_T1day"=T1day_vs_control, "T1day_vs_T8day"=T1day_vs_T8day)

# import to mitch
allmitch <- mitch_import(allresults, DEtype = 'edgeR', joinType = 'inner', geneTable = gene_trans)

## calculate enrichment.
# the report fails if the number of gene sets is less than the specified by resrows in calc (default of 50). Therefore it is necessary to run this, see how many gene sets are kept after filtering (length(res.BP$enrichment_result$set)), and rerun setting resrows accordingly

# calculate
res.BP.20 <- mitch_calc(allmitch, geneset.BP, priority="significance",  cores = 4, minsetsize = 20)

#filter: keep FDR < 0.01, s> 0.1
res.BP.20$enrichment_result <- filter(res.BP.20$enrichment_result, p.adjustMANOVA < 0.01, s.dist > 0.1)
num_sets <- length(res.BP.20$enrichment_result$set)
res.BP.20 <- mitch_calc(allmitch, geneset.BP, priority="significance",  cores = 4, minsetsize = 20, resrows = num_sets)

# Save results
# generate reports and save charts in pdf format
outfile = "output/report_BP_minset20.html"
#outfile2 = "output/charts_BP_minset20.pdf"
mitch_report(res.BP.20, outfile = outfile)
#mitch_plots(res.BP.20, outfile = outfile2)

# save all gene sets. 
# save as tsv
write.table(res.BP.20$enrichment_result, "output/gene_sets_BP.tsv", sep = '\t', row.names = F, quote = FALSE)
# save to Xcel file
outfile = "output/gene_sets_BP.xlsx"
write_xlsx(res.BP.20$enrichment_result, path = outfile, col_names = TRUE, format_headers = FALSE)


# Assign broad group classification and generate heatmap
# the broad group classification was done with the aid of Amigo2 - Visualize

# generate heatmap palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)

# import file with GO terms classified into broad groups
broad_groups <- read.table("Amigo_visualize/out_BP_GOids_to_broad_groups.csv", sep = '\t', header = TRUE, quote = "", fill = FALSE, row.names = 1)

# prepare the mitch result element used in the heatmap (code is modified from the "mitch.Rmd" file in mitch)

# number of comparisons
d=ncol(res.BP.20$input_profile)
# df with the s values
res_hm_raw <- head (res.BP.20$enrichment_result[,4:(4+d-1)], num_sets)
rownames(res_hm_raw) <- head(res.BP.20$enrichment_result$set, num_sets)
colnames(res_hm_raw) <- gsub("^s.","",colnames(res_hm_raw))

# add the broad group classification 
res_hm <- merge(res_hm_raw, broad_groups, by = "row.names", all = TRUE)
# sort (this is how the GO terms will appear in the heatmap)
res_hm <- res_hm[order(res_hm$sorting_index), ]

# The heatmap is generated with the function heatmap.2. Here, cexCol,cexRow control column and row font size. If font size is too large, some labels will not display
# margin controls the size of the margins: it is a numeric vector with 2 values: top, right. The higher the value, more space available for text labels, less for the plot. So for example, to make the heatmap cells wider, decrease the "right" value.
# To the mitch code, I added the arguments Rowv = FALSE, dendrogram = 'column'. This supresses the row dendrogram and arranges the GO terms as they are in the input matrix 

# But before generating the heatmap, the dataframe needs a final arrangement. The heatmap function only uses the columns with the s values, hence all columns without s values must be dropped. The column with the values to be used as cell labels can be saved as a vector and specified with the argument "labRow"

# gather GOterms as labels
label_goterms <- res_hm[, 1]
# gather the broad groups as labels
label_broad_groups <- res_hm[, "broad.group"]

# dataframe with only the columns with s values.
res_hm_svals <- res_hm[, c("ctrl_vs_T8day", "ctrl_vs_T1day", "T1day_vs_T8day")]

# Save heatmap with GO terms as row labels
png(file = "output/plots/heatmap_goterms.png", res = 600, width = 350, height = 500, units = "mm", bg = "transparent")
heatmap.2(as.matrix(res_hm_svals), scale="none", trace="none", margin=c(10, 50), cexCol=1, cexRow=1, col=my_palette, Rowv = FALSE, dendrogram = 'none', labRow=label_goterms, Colv = c(1,2,3))
dev.off()

# Save heatmap with broad groups as row labels
png(file = "output/plots/heatmap_broadgroups.png", res = 600, width = 350, height = 500, units = "mm", bg = "transparent")
heatmap.2(as.matrix(res_hm_svals), scale="none", trace="none", margin=c(10, 50), cexCol=1, cexRow=1, col=my_palette, Rowv = FALSE, dendrogram = 'none', labRow=label_broad_groups, Colv = c(1,2,3))
dev.off()


### Part 2: identify the gene sets of interest

# arrange tables
# subset significant gene sets, keep all columns
res_high_int <- head(res.BP.20$enrichment_result, num_sets)
rownames(res_high_int) <- res_high_int$set
# add the broad group classification 
res_high_int <- merge(res_high_int, broad_groups, by = "row.names", all = TRUE)

# Filters to identify the gene sets of interest
# 1) they belong to the following broad groups/unclassified gene sets: 
broad_groups_interest <- c("cell adhesion", "cell migration", "cytoskeletal and ECM organization", "lipid metabolism", "muscle contraction", "regulation of membrane potential")
res_high_int_1 <- res_high_int[res_high_int$broad.group %in% broad_groups_interest, ]

# 2) the enrichment change between ctrlvsT1 and T1vsT8 of at least 0.1 (this means there is a palpable change in enrichment after day 1)
res_high_int_2 <- res_high_int_1[abs(res_high_int_1$s.T1day_vs_T8day - res_high_int_1$s.ctrl_vs_T1day) >= 0.1, ]

# 3) Significant change in enrichment direction between ctrlvsT1 and T1vsT8. For this, one of these conditions must be met:
# not significant for ctrlvsT1 but significant for T1vsT8
condA <- res_high_int_2$p.ctrl_vs_T1day > 0.01 & res_high_int_2$p.T1day_vs_T8day <= 0.01
# significant for ctrlvsT1 but not significant for T1vsT8
condB <- res_high_int_2$p.ctrl_vs_T1day <= 0.01 & res_high_int_2$p.T1day_vs_T8day > 0.01
# significant for both ctrlvsT1 and T1vsT8, but of opposite sign
condC <- (res_high_int_2$p.ctrl_vs_T1day <= 0.01 & res_high_int_2$p.T1day_vs_T8day <= 0.01) & (res_high_int_2$s.ctrl_vs_T1day * res_high_int_2$s.T1day_vs_T8day < 0)

res_high_int_3 <- res_high_int_2[(condA | condB | condC), ]

# Results: the gene sets of highest interest are:
mysets.names <- res_high_int_3$set
# These are the gene sets:
# "actin cytoskeleton organization", "actin filament organization", "actomyosin structure organization", "cell adhesion", "cell adhesion mediated by integrin", "cell migration", "cell-cell adhesion", "cell-matrix adhesion", "chromatin organization", "extracellular matrix organization", "GPI anchor biosynthetic process", "homophilic cell adhesion via plasma membrane adhesion molecules", "lipid metabolic process", "microtubule cytoskeleton organization", "muscle contraction", "neural crest cell migration", "positive regulation of cell migration", "regulation of cell migration", "regulation of membrane potential" 

# Save them to file
write.table(mysets.names, "output/gene_sets_highest_interest.tsv", sep = '\t', row.names = F, col.names = F, quote = FALSE)


### Part 3: find the genes common to these gene sets and the DEG genes from edgeR 

# get the genes in the gene sets of highest interest
mysets <- res.BP.20$detailed_sets[mysets.names]

#move to a single dataframe, with gene_set as column
df.mysets <- ldply(mysets, data.table::as.data.table, keep.rownames = TRUE)
colnames(df.mysets)[1] <- "gene_set"
colnames(df.mysets)[2] <- "Drerio_Entrez_geneID"

# Add the Bbrach entrez geneIDs
df.mysets.Bbrach <- dplyr::inner_join(df.mysets, gene_trans, by = "Drerio_Entrez_geneID") %>% dplyr::relocate(Bbrach_Entrez_geneID, .after = Drerio_Entrez_geneID)
# notice that this new df has several more rows than df.mysets, this due to several cases of a Drerio_Entrez_geneID matching multiple Bbrach_Entrez_geneIDs

#drop unwanted columns from the genes in gene sets
df.mysets.short.Bbrach <- dplyr::select(df.mysets.Bbrach, c("gene_set", "Drerio_Entrez_geneID", "Bbrach_Entrez_geneID"))


## get genes in common between the DEG from edgeR and the gene sets of interest

# load the annotated, DEG. These are the three T_*.tsv files in the folder 06_DGE
# make a vector with the names of the files to import
DEG.files <- list.files(path = "../06_DGE", pattern = "*tsv", full.names = TRUE)

#extract the name of each comparison
comparisons <- lapply(DEG.files, function(x) {sub(".*DGE/", "", x) %>% sub(".tsv", "", .)})

#Import the annotated DEG for each comparison
DEG.annot <- lapply(DEG.files, function(x) as.data.frame(read.table(x, header = TRUE, row.names = NULL, sep = "\t", quote = "", stringsAsFactors =FALSE)))

#name each data frame with the respective comparison
names(DEG.annot) <- comparisons

# convert the Bbrach_Entrez_geneID column to character (necessary for the join function in next step)
DEG.annot = lapply(DEG.annot, function(x) transform(x, Bbrach_Entrez_geneID = as.character(Bbrach_Entrez_geneID)))

# inner join each list of DE genes with the genes in gene sets
DEG.in.genesets <- lapply(DEG.annot, function(x) {dplyr::inner_join(df.mysets.short.Bbrach, x, by = "Bbrach_Entrez_geneID")})

## Some genes belong to multiple gene sets. Collapse these to one row
# save all column names minus "gene_set"
cols = setdiff(colnames(DEG.in.genesets$T_1day_vs_control), "gene_set")
# collapse repeated gene sets
DEG.in.genesets.collapsed = lapply(DEG.in.genesets, function(x) {plyr::ddply(x, cols, plyr::summarize, gene_set=paste0(sort(gene_set), collapse = ", "))})
# reorder columns: move gene_set back to first column
DEG.in.genesets.collapsed = lapply(DEG.in.genesets.collapsed, function(x) {dplyr::select(x, "gene_set", everything())})
 

## save files and save workspace

# all genes in gene sets of interest
outfile = "output/all_genes_in_key_genesets.xlsx"
write_xlsx(df.mysets.short.Bbrach, path = outfile, col_names = TRUE, format_headers = FALSE)

# DEG in gene sets of interest
outfile = "output/DEG_genes_in_key_genesets.xlsx"
write_xlsx(DEG.in.genesets.collapsed, path = outfile, col_names = TRUE, format_headers = FALSE)

# save the workspace
save.image(file = "output/gene_sets_BP.RData")
