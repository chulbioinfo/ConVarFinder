## Overview ##

The script to run 'fixed difference' analysis compares the position of CSAV in the multiple sequence alignment (MSA) file with the known SNPs in the chromosomal position. The script first converts the CSAV position in the MSA file to chromosomal position based on the BLAST output file and the chromosomal positions of CSAVs are then compared with population variant dataset (vcf file downloaded from Ensembl) to check if CSAV is conserved within the population data. The script here is unable to run on different dataset.

## Requirements ##

The script was previously run with linux (Ubuntu 18.04.3) server with RAM: +100GB and CPU: 24 cores. However, the code does not require much resource. <1GB and 1 core would suffice to run the code. The script is based on python (3.7.3) and hence python (preferrably >3.7.3) needs to be installed prior to running this script. This could be done by:

    sudo apt-get update
    sudo apt-get install python3.8

## Running the script ##

To run this script, the input files needed include:

1. MSA files in "Avian_MSA_Sequence/" directory (please unzip the MSA files before running the code),
2. "TAAS_avians.txt" containing the CSAV variant information,
3. "GALGA_from_Avian_MSA_Sequence.blast_out" - BLAST output file
4. "TAEGU_from_Avian_MSA_Sequence.blast_out" - BLAST output file
5. Variant files (can be downloaded from "https://drive.google.com/file/d/1jZMFxDjGmRdN5lRWKwVVLDcGrPgg_GPy/view?usp=sharing")

Once you have all of these input files, the script can be simply run by:

    python Check_FixedDifference.py

This will generate output files:

1. FixedDifference_Chicken.txt
2. FixedDifference_ZebraFinch.txt

which contains one variant for chicken overlaping with a CSAV and none for Zebra Finch.
