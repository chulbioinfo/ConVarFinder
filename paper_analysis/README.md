# Paper Analysis Code Repository

This repository contains the code used for various analyses in the paper, organized by analysis steps.

## Directory Structure

- `0.rawdata/`: Directory containing raw data files
- `01.CSAV/`: Convergent Single Amino acid Variant analysis code (integrated into ConVarFinder)
- `02.PhylogeneticFeatures/`: Phylogenetic features analysis code
- `03.CSCV/`: Convergent Single Codon Variant analysis code (integrated into ConVarFinder)
- `04.CSNV/`: Convergent Single Nucleotide Variant analysis code (integrated into ConVarFinder)
- `05.Intersection_CSAV_CSCV_CSNV/`: Code for analyzing intersections of CSAV, CSCV, and CSNV results
- `06.CorrelationPlots/`: Correlation analysis and visualization code
- `07.CodonLogo/`: Codon logo visualization
- `08.GOanalysis/`: Gene Ontology analysis code
- `09.FixedDifference/`: Fixed difference analysis code
- `10.PositiveSelection/`: Positive selection analysis data sheets
- `11.DifferentialExpression/`: Differential expression analysis data sheets
- `12.B10K_2nd_363birds/`: Analysis code for 363 bird species from the B10K project
- `13.AncSeq/`: Ancestral sequence estimation (to generate inputs of ConVarFinder)

## Usage

Each directory contains the necessary code and scripts for its respective analysis step. Please refer to the README file within each directory before running the analyses.

## Dependencies

Required dependency packages for each analysis step are specified in the requirements.txt or environment configuration files within each directory.

## Notes

- Please ensure that all required data files are in the correct locations before running the analyses.
- Analysis steps may need to be executed sequentially, so it is recommended to follow the directory numbering order.
- Sufficient memory and disk space may be required for processing large datasets. 
