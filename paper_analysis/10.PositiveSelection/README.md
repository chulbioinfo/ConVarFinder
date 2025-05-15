# Overview
This section to run dN/dS analysis by using codeML in PAML.
## Input Files
  * **6110.pal**: Mutiple codon sequence alignments of _DRD5_
  * **6110.tre**: Avian familiy tree marked with foreground branches of 3 vocal learning clades as '#1' 

## Control files
  * **6110_h0.ctl**: Control file with options for null model (neutral selection on foreground branches, dN/dS =1)
  * **6110_h1.ctl**: Control file with options for alternative model (positive selection on foreground branches, dN/dS >1)

## Output Files
  * **6110_h0_output.txt**: Output file of null model (neutral selection on foreground branches, dN/dS =1)
  * **6110_h1_output.txt**: Output file of alternative model (positive selection on foreground branches, dN/dS >1)
- - -

# Reference
[PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
