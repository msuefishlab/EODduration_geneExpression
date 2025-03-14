The first step is to wrangle the transcripts from the NCBI-annotated (release 100) B. brachyistius genome. From this annotation, we start with the file rna.fna. Our file cmd_wrangle_ncbi_rna_file.txt has the bash commands that process the transcripts. It produces 2 outputs: a "map" file that links transcripts to genes (gene-trans-map.txt, provided in this folder), and transcripts.fna, which includes the predicted transcripts minus the rRNAs.

Then, the code in files sb_04* estimates transcript abundance. Files sb_04a* create the RSEM index, sb_04b* align reads to transcripts and quantify gene expression, and sb_04c* build gene expression matrices. Some of the outputs are used in the Differential Gene Expression Analyses.

Finally, the code in files sb_05* explores the gene expression data. One of the outputs is the gene expression correlation matrix.

