# Ancestral Sequence Estimation

A Python script for estimating ancestral sequences using RAxML with both codon and binary (gap) models.

Created: May 15, 2025

## Description

This script performs ancestral sequence estimation by combining information from both codon and binary (gap) models. It processes aligned sequences to reconstruct ancestral states while properly handling gaps in the alignment.

## Features

- Converts input FASTA to binary sequences (gap: 1, non-gap: 0)
- Runs RAxML for both binary and codon models
- Combines results: replaces codon sequences with gaps where binary state is 1
- Generates final ancestral sequence FASTA and cladogram tree files
- Creates a dedicated working directory for each analysis
- Handles duplicate sequences automatically
- Output files are compatible with ConVarFinder for subsequent analysis

## Requirements

- Python 3.6 or higher
- Biopython
- RAxML (standard-RAxML/raxmlHPC-PTHREADS-SSE3)
- Required Python packages:
  * argparse
  * os
  * re
  * sys
  * subprocess
  * shutil

## Installation

1. Clone the repository
2. Install required Python packages:
   ```bash
   pip install biopython
   ```
3. Ensure RAxML is installed in the `standard-RAxML` directory

## Usage

Basic usage:
```bash
python ancestral_sequence_estimation.py input.fasta input.tre
```

With custom output name:
```bash
python ancestral_sequence_estimation.py input.fasta input.tre --output my_output
```

## Input Requirements

- FASTA file: Must contain codon-aligned sequences
- Tree file: Must be in Newick format and compatible with RAxML

## Output Files

The script creates a directory named `{input_id}_ancestral` containing:

- `{input_id}_anc.fasta`: Combined FASTA file with original and ancestral sequences
- `{input_id}_anc.tre`: Rooted tree file with node labels
- `binary_input.fasta`: Binary sequence file
- `RAxML_*` files: Intermediate RAxML output files
- `*_reduced.fasta`: Files with duplicate sequences removed

The output files can be directly used as input for ConVarFinder for subsequent analysis.

## RAxML Models

- Binary sequences: BINCAT model
- Codon sequences: GTRCAT model

## Error Handling

The script includes error handling for:
- RAxML executable existence
- Input file format validation
- RAxML execution errors
- File I/O operations

## Example

```bash
# Using example files
python ancestral_sequence_estimation.py 10050.fasta 10050.tre
```

This will create a directory `10050_ancestral` containing all output files.

## Notes

- The script automatically handles duplicate sequences by creating reduced alignment files
- All intermediate files are stored in the working directory
- The final output includes both original and ancestral sequences
- The generated files are compatible with ConVarFinder for further analysis

## License

MIT License

Copyright (c) 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Author

Chul Lee

Postdoctoral Associate

Laboratory of Neurogenetics of Language (http://www.jarvislab.net/)

Rockefeller University

Contact: chul.bioinfo@gmail.com  
GitHub: @chulbioinfo 
