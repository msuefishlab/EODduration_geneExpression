Mitch needs three inputs: profiled expression data (generated in the Differential Expression Analyses), a gene set library (folder 1- Gene set library), and a "translator" gene table that relates gene identifiers on the other two input files (folder 2- Gene table). 

Once all required inputs are available, run mitch (folder 3- mitch analysis) 



For the gene set library, we used a Danio rerio gmt file with gene ontology (GO) (Ashburner et al., 2000; Carbon et al.,
2021) terms from the GO domain Biological Process as gene sets. To generate this file, we
employed the script update_GO.sh from GeneSCF-v1.1-p3 (Subhash and Kanduri, 2016). 

Since
the gene identifiers in the gmt file are different than those of the edgeR results, mitch requires a
third input file that relates these gene identifiers. 

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
