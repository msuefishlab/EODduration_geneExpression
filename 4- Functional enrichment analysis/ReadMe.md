Mitch needs three inputs: 

a) profiled expression data (generated in the Differential Expression Analyses)  
b) a gene set library (generated with the script update_GO.sh from GeneSCF-v1.1-p3, and provided here: GO_BP_gid.txt)  
c) a "translator" gene table that relates gene identifiers from the other two input files (folder 1- Gene table). 

Once all required inputs are available, run mitch (folder 2- mitch analysis) 


To further simplify the analysis interpretation, we employed the tool Visualize from AmiGO 2
(Carbon et al., 2009) to build GO graphs with the enriched GO terms, and based on their
hierarchy and connectivity, we manually grouped the GO terms into broad functional categories.
