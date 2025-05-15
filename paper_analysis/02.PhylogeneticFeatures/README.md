# Overview
This script to run 'PhylogeneticFeatureCalculator' to calculate 4 types of phylogenetic features: product of origin branch lengths (POB), product of termianl branch lengths (PTB), distance betwen terminal nodes (DTN), and distance between terminal branches (DTB).
## Input Files
It does not need any input file.
## Variables in Script
  1. target species list # avian vocal learners
  2. phylogenetic tree with branch lengths
## Output Files
It generates a text file (.txt) as a output with a summary of phylogenetic features of target species.
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
This script does not need any python library except for 'sys' library which is a default module of Python.

## Running time with the demo data
* under 1 sec
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
cd 02.PhylogeneticFeatures/
python PhylogeneticFeatureCalculator.py
</code>
</pre>

(Optional step to check inputs before running the script)
  - This sciprt requires following variables:
  1. target species
  2. phylogenetic tree with branch lengths

  - Example variables in the script
<pre>
<code>
    nTargets ="TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN" # 6 vocal learning birds
    nTree = "(((((((((((((TAEGU:0.042,GEOFO:0.0372):0.0245,CORBR:0.0385):0.044,MANVI:0.0693):0.0479,(MELUN:0.0585,NESNO:0.042):0.0528):0.0021,FALPE:0.0762):0.0013,CARCR:0.055):0.0015,(((((((MERNU:0.103,PICPU:0.1434):0.0039,BUCRH:0.0998):0.0025,APAVI:0.1048):0.0012,LEPDI:0.067):0.0014,COLST:0.1205):0.0006,TYTAL:0.0675):0.0006,((HALLE:0.001,HALAL:0.0011):0.0427,CATAU:0.0301):0.002):0.0008):0.0041,((((((PELCR:0.0417,EGRGA:0.0562):0.0013,NIPNI:0.0398):0.001,PHACA:0.0639):0.0016,(FULGL:0.0342,(PYGAD:0.0135,APTFO:0.0109):0.0214):0.0012):0.0019,GAVST:0.0445):0.0021,(PHALE:0.0591,EURHE:0.0947):0.0032):0.0012):0.0011,((CHAVO:0.0602,BALRE:0.0571):0.0015,OPHHO:0.0756):0.0009):0.001,(((CALAN:0.1126,CHAPE:0.0908):0.0204,CAPCA:0.0697):0.0038,((CHLUN:0.0731,TAUER:0.0708):0.0014,CUCCA:0.1115):0.0011):0.0011):0.0011,(((MESUN:0.0871,PTEGU:0.0717):0.0027,COLLI:0.0963):0.0014,(PHORU:0.0312,PODCR:0.0584):0.013):0.0009):0.032,((GALGA:0.0391,MELGA:0.0449):0.1131,ANAPL:0.1003):0.0394):0.0456,(TINMA:0.1471,STRCA:0.0669):0.0456);"
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.
