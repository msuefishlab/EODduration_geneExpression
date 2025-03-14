We inspected raw and processed reads with FastQC v0.11.7 (Babraham Bioinformatics) and
used Trimmomatic v.0.39 (Bolger et al., 2014) to remove library adaptors, low quality reads, and
filter small reads. 

The succeeding steps were executed using scripts included with Trinity
v2.11.0 (Grabherr et al., 2011; Haas et al., 2013). We aligned reads from each specimen to the
predicted transcripts of the NCBI-annotated (release 100) B. brachyistius genome, with bowtie2
v2.3.4.1 (Langmead and Salzberg, 2012). Expression quantification was estimated at the gene
level using RSEM v1.3.0 (Li and Dewey, 2011), followed by exploration of the data with a gene
expression correlation matrix based on Euclidean distances and Pearson’s correlation
coefficient (for genes with read counts > 10, Trinity’s default parameters).



singularity for trinity, slurm?, wrangle Bbrach transcripts and Danio prots?
