# Overview
This script to run 'SCV_SAV' analysis to idenify single codon and amino acid variants specific to a species which you are interested in. Based on this script and demo data, you can find human-specific codon and amino acid variants mutually exclusive to codona and amino acids of non-human primates.
## Input Files
It needs multiple peptide sequence aglinments as input files (as [fasta](https://en.wikipedia.org/wiki/FASTA_format) format). 'orthologs_information' folder includes descriptions to generate multiple sequence alignments of the primate lineage.
## Variables in Script
  1. input path
  2. input file format
  3. whole species list 
  4. target species list # Human 
  5. outgroup species list # no outgroup species
  6. output path
## Output Files
It generates a text file (.txt) as a output with a summary of CSV_SAV analysis and the full codons and amino acids of whole species list at identified SCV and SAV sites.
- - -

# System Requirments
## Hardware requirements
This script requires only a standard computer with enough RAM to support the in-memory operations. It was developed and tested in a standard computer (AMD Ryzen 5 2400G 3.6GHz and 16GB RAM)

## Software requirements
### OS Requirements
This script is supported for *Windows OS* and *Linux*. The script has been tested on the following systems:
* Windows 10 Education
* Linux: Ubuntu 18.04.3 LTS

## Python Dependencies
This script mainly depends on the Python scientific stack.

    numpy
    scipy

## Running time with the demo data
* Under 1 min
- - -

# Running the script
* (option1) If you want to clone the github link:
<pre>
<code>
git clone https://github.com/chulbioinfo/CSAVanalysis/
cd CSAVanalysis/
</code>
</pre>

* (option2) If you use the zip file from the github link:
<pre>
<code>
wget https://github.com/chulbioinfo/CSAVanalysis/archive/master.zip
unzip master.zip
cd CSAVanalysis-master/
</code>
</pre>

* Run the script
<pre>
<code>
cd 12.HumanSAV/Analysis/
python SCV_SAV.py
</code>
</pre>

(Optional step to check inputs before running the script)
  - This sciprt requires following variables:
  1. input path
  2. input file format
  3. whole species list 
  4. target species list
  5. outgroup species list
  6. output path

  - Example variables in the script
<pre>
<code>
    seqPATH = "./"
    seqFORM = ".best.fas"
    sID_LIST = ["PANTR","PANPA","HOMSA","GORGO","PONAB","NOMLE","CERAT","MANLE","PAPAN","MACFA","MACMU","MACNE","CHLSA","RHIRO","RHIBI","COLAN","CEBCA","SAIBO","CALJA","AOTNA","CARSY","MICMU","PROCO","OTOGA"]
    nTargets = 'HOMSA'
    outgroup_LIST = []
    oPATH = "./"
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.

