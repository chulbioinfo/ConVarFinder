# Overview
This section for gene ontology analysis for genes with molecular convergences.
* **gp_in/iCSAV/TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN.txt**: Input file of _gprofiler_
* **run_gprofiler.py**: Python script to run _gprofiler_ by using _Python_
* **gp_out/iCSAV/TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN.txt**: Output file of _gprofiler_
* **summary_gprofiler_20190323_gene_and_gobp_over1.txt**: Input file for correlations between the number of genes with molecular convergences and both of the number of GO terms and p-values of GO terms
* **summary_gprofiler_20190324.R**: R script to visualize correlations between the number of genes with molecular convergences and both of the number of GO terms and p-values of GO terms
* **Fig6_20190323.pdf**: Visual item of correlations between the number of genes with molecular convergences and both of the number of GO terms and p-values of GO terms

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
This script mainly depends on several Python packages.

    from gprofiler import gprofiler
    import pandas
    import sys
    import re
    import os
    import string
    import glob
    import queue
    import threading
    import time

## R Dependencies
This script mainly depends on the ggplot 2.

    ggplot2

## Running time with the demo data
* under 1 min


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
cd 08.GOanalysis/bin/
python run_gprofiler.py
Rscript summary_gprofiler_20190324.R
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.
