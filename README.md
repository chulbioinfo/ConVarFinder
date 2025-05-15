# ConVarFinder

**ConVarFinder** is a Python tool for detecting **convergent evolutionary variants** at the **codon**, **amino acid**, and **nucleotide** levels. It is designed for comparative genomic studies across multiple species using aligned coding sequences and a phylogenetic tree.

---

## Features

- Detects **target-specific substitutions** at:
  - Codon level
  - Amino acid level
  - Nucleotide positions 1, 2, and 3
- Distinguishes between **convergent**, **parallel**, **divergent**, and **plesiomorphic** evolutionary patterns
- Supports **monophyletic clade detection** via tree parsing
- Generates **rich tab-delimited output** summarizing convergent sites and lineage-specific substitutions
- Written in **pure Python 3**, requires **no external libraries**

---

## Quick Start with Example Files

```bash
python ConVarFinder.py -seq=10050_anc.fasta -tree=10050_anc.tre -tl=TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN -no-outgroup-prompt
```

If you do not use '-no-outgroup-prompt', you will be prompted to select outgroup species (optional).

---

## Input

> **Important:**  
> - The **sequence file** (`-seq`) must include ancestral sequences inferred using tools such as **RAxML**, with internal node names representing ancestral nodes.  
> - The **tree file** (`-tree`) must be a **binary cladogram** (no branch lengths) in **Newick format**, with all internal nodes labeled.  

- **FASTA file** (`-seq=`): Coding sequences for all species. Must be aligned and codon-correct.
- **Newick tree file** (`-tree=`): Phylogenetic tree describing species relationships.
- **Target species list** (`-tl=`): Comma-separated terminal node names.
- **Outgroup species list** (`-ol=`): Optional. Used to define substitution polarity.

---

## Output

Tab-delimited `.txt` file in the format:

```
ConVarFinder_TAEGU-GEOFO-CORBR-MELUN-NESNO-CALAN.txt (ConVarFinder_speciesA-speciesB.txt)
```

Each row corresponds to a codon site with target-specific substitutions, with annotations including:
- Gene and site position
- Substitution patterns
- Evolutionary direction (e.g., `PLE`, `CON`, `DIV`, `PAR`, `MIX`)
- Codon and amino acid states for each species
- Transition paths across ancestral nodes


---

## Options

| Option         | Description |
|----------------|-------------|
| `-seq=FILE`    | Aligned CDS file in FASTA format *(required)* |
| `-tree=FILE`   | Phylogenetic tree in Newick format *(required)* |
| `-tl=LIST`     | Comma-separated target species list *(required)* |
| `-ol=LIST`     | Comma-separated outgroup species list *(optional)* |
| `-ip=PATH`     | Input directory path *(default: ./)* |
| `-op=PATH`     | Output directory path *(default: ./)* |
| `-sfmt=.fasta` | FASTA file extension *(default: .fasta)* |
| `-tfmt=.tre`   | Tree file extension *(default: .tre)* |
| `-no-outgroup-prompt` | Skip prompt for outgroup selection |

---

## Evolutionary Pattern Legend

| Code     | Description |
|----------|-------------|
| `CON`    | Convergent substitution to same state from different ancestors |
| `PAR`    | Parallel substitution from same ancestor |
| `PLE`    | Plesiomorphic (no change from ancestral state) |
| `DIV`    | Divergent substitutions among target clades |
| `cCON`   | Convergent among targets (clade-level) |
| `cDIV`   | Divergent among targets (clade-level) |
| `sCON`, `sPAR`, `sDIV` | Sub-clade level patterns |
| `MIX(...)` | Mixed evolutionary signals |

---

## Preparing Input Files (`Make_inputs/`)

The `Make_inputs/` directory contains helper scripts used to prepare sequence and tree files  
required for running ConVarFinder. This includes steps such as:

- Extracting aligned coding sequences from genome annotations
- Performing ancestral state reconstruction using tools like **RAxML**
- Formatting outputs into compatible `.fasta` and `.tre` files

These scripts ensure that input files contain both extant and inferred ancestral sequences, and that tree files are in binary cladogram format.

---

## Dependencies

- Python ≥ 3.6
- No third-party libraries required

---

## Author

- **Chul Lee**, Rockefeller University  
- Contact: [chul.bioinfo@gmail.com](mailto:chul.bioinfo@gmail.com)  
- GitHub: [@chulbioinfo](https://github.com/chulbioinfo/ConVarFinder)

---

## License

MIT License
