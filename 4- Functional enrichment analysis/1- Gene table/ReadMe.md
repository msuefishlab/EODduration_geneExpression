The file cmd_protein_dbs.txt contains bash commands (and comments) to identify homologous proteins between B. brachyistius and D. rerio. One step in this file calls for the execution of the script sb_07_blastp.sh. 

Then, use the R code in Annotation_wrangling.Rmd. This code uses mygene to match 1) B. brachyistius proteins to B. brachyistius genes, and 2) D. rerio proteins to D. rerio genes. After this, the code matches the Entrez_geneIDs between both species. The output file geneIDs_Bbrach_Drerio.txt is the "translator" gene table mitch uses.
