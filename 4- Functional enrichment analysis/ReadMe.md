Mitch needs three inputs: 

a) profiled expression data (generated in the Differential Expression Analyses)
b) a gene set library (generated with the script update_GO.sh from GeneSCF-v1.1-p3, and provided here: GO_BP_gid.txt)
c) a "translator" gene table that relates gene identifiers from the other two input files (folder 1- Gene table). 

Once all required inputs are available, run mitch (folder 2- mitch analysis) 




To create this file, we first identified homologous proteins predicted from the B. brachyistius reference genome and those predicted
from Danio rerio (GRCz11) by blastp (BLAST+v2.11.0, Camacho et al., 2009). For each protein,
the top hit (e-value ≤1e-10) was used for annotation. 

Then, we used mygene v1.32.0 (Wu et al., 2013; Xin et al., 2016) to match the D. rerio proteins to D. rerio genes.


We excluded gene sets with fewer than 20 genes from the mitch analysis (minsetsize option in
mitch_calc function); and filtered the enrichment result to gene sets with FDR-corrected p-value
< 0.01 and the higher dimensional enrichment score (S) > 0.1, to minimize false positives. 

To further simplify the analysis interpretation, we employed the tool Visualize from AmiGO 2
(Carbon et al., 2009) to build GO graphs with the enriched GO terms, and based on their
hierarchy and connectivity, we manually grouped the GO terms into broad functional categories.
