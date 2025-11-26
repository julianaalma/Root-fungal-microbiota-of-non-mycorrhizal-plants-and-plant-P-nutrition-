### Screening plant genomes for Plant AM genes

The main script "Main_Slurm_Script.sh" describes all the steps of the screening analysis using reference sequences from Medicago "Mt_selected" and Rice "Os_selected".
Identified sequences were aligned with Muscle and gene trees were produced with iqtree, including gene reference sequences from the literature listed in files "refseqs".
A species megatree was created with script "sp_tree.R" and file "2_taxa_tree.csv"
Outputs are listed in the excel file "Genomes_151123.xlsx". Options for tree plotting with itol are included.
