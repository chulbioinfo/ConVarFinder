# Overview
This script to run 'CSCV' analysis to idenify convergent single codon variants (CSCVs) of target species which share a convergent trait you are interested in. Based on this script and demo data, you can find CSCVs of avian vocal learners mutually exclusive to codons of vocal non-learners.
## Input Files
It needs multiple codon sequence aglinments as input files (as [fasta](https://en.wikipedia.org/wiki/FASTA_format) format).
## Variables in Script
  1. input path
  2. input file format
  3. whole species list
  4. target species list # avian vocal learners
  5. outgroup species list # Rifleman with the uncertainty for vocal learning ability
  6. output path
## Output Files
It generates a text file (.txt) as a output with a summary of CSCV analysis and the full codon and amino acids of whole species list at identified CSCV sites.
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
* 4 min 52.236422 sec
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

* Decompress the demo data of input files (as [fasta](https://en.wikipedia.org/wiki/FASTA_format) format)
<pre>
<code>
cd 0.rawdata/MSA/
unzip MSA.zip
cd ../../
</code>
</pre>

* Run the script
<pre>
<code>
cd 03.CSCV/bin/
python CSCV.py
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
    seqPATH = "../../0.rawdata/MSA/cds/"
    seqFORM = ".sate.default.pep2cds.removed.shortname.fasta"
    sID_LIST = ['TAEGU','GEOFO','CORBR','MELUN','NESNO','CALAN','MANVI','FALPE','CARCR','MERNU','PICPU','BUCRH','APAVI','LEPDI','COLST','TYTAL','HALLE','HALAL','CATAU','PELCR','EGRGA','NIPNI','PHACA','FULGL','PYGAD','APTFO','GAVST','PHALE','EURHE','CHAVO','BALRE','OPHHO','CHAPE','CAPCA','CHLUN','TAUER','CUCCA','MESUN','PTEGU','COLLI','PHORU','PODCR','GALGA','MELGA','ANAPL','TINMA','STRCA','HUMAN','ACACH']
    nTargets = 'TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN'
    outgroup_LIST = ['HUMAN','ACACH']
    oPATH = "../output/"
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.
