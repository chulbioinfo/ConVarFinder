# Overview
The scripts to visualize correlations among the number of molecular convergences at 3 levels (amino acids, codons, and nucleotides) and the phylogenetic features (POB, PTB, DTN, and DTB).
## Input Files
It needs text files (.txt) with counts of types of molecular convergences as input files.
## Variables in Script
  1. working directory
  2. input files
  3. output names

## Output Files
It generates visualization items and data matrix for correlations (.pdf and .csv) as outputs.
- - -

# System Requirments
## Hardware requirements
This script requires only a standard computer with enough RAM to support the in-memory operations. It was developed and tested in a standard computer (AMD Ryzen 5 2400G 3.6GHz and 16GB RAM)

## Software requirements
### OS Requirements
This script is supported for *Windows OS* and *Linux*. The script has been tested on the following systems:
* Windows 10 Education
* Linux: Ubuntu 18.04.3 LTS

## R Dependencies
This script mainly depends on the following packages in R (v3.5.3).

    install.packages("rlang") # v0.4.5
    install.packages("car") # v3.0-10

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
cd 06.CorrelationPlots/
Rscript Fig1_script_20201107.R
Rscript Fig2_script_20201107.R
Rscript Fig3_script_20201107.R
</code>
</pre>

(Optional step to check inputs before running the script)
  - This sciprt requires following variables:
  1. working directory
  2. input files
  3. output names

  - Example variables in the script
<pre>
<code>
    # Fig1_script_20201107.R
    ## set working directory
    setwd("./")
    ## input files
    data <- read.table("1000ctrl/0.summary_Phy_CSAV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("1000ctrl/0.summary_Phy_CSCV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("1000ctrl/0.summary_Phy_CSNV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSAV_20201011.txt",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSCV_20201011.txt",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSNV_20201011.txt",sep='\t',header = T)
    ## output files
    pdf("./Fig1d_iCSAV_dCSAV_randomctrl_v20201011.pdf",width=6,height=6)
    pdf("./Fig1e_iCSAV_dCSAV_randomctrl_v20201011.pdf",width=6,height=6)
    pdf("./Fig1f_iCSAV_dCSAV_randomctrl_v20201011.pdf",width=6,height=6)
    pdf("./Fig1g_iCSAV_dCSAV_randomctrl_v20201011.pdf",width=6,height=6)

    # Fig2_script_20201107.R
    ## set working directory
    setwd("./")
    ## input files
    fNAME = "1000ctrl/0.summary_Phy_CSAV_20190321.marked_targets.info"
    fNAME = "1000ctrl/0.summary_Phy_CSCV_20190321.marked_targets.info"
    fNAME = "1000ctrl/0.summary_Phy_CSNV_20190321.marked_targets.info"
    fNAME = "61ctrl/0.summary_Phy_CSAV_20201011.txt"
    fNAME = "61ctrl/0.summary_Phy_CSCV_20201011.txt"
    fNAME = "61ctrl/0.summary_Phy_CSNV_20201011.txt"
    ## output files
    pdf("./Fig2_1000randomctrl_v20201011.pdf",width=8,height=6)
    pdf("./FigS4_61corectrl_v20201011.pdf",width=8,height=6)

    # Fig4_script_20201107.R
    ## set working directory
    setwd("./")
    ## input files
    data <- read.table("1000ctrl/0.summary_Phy_CSAV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("1000ctrl/0.summary_Phy_CSCV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("1000ctrl/0.summary_Phy_CSNV_20190321.marked_targets.info",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSAV_20201011.txt",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSCV_20201011.txt",sep='\t',header = T)
    data <- read.table("61ctrl/0.summary_Phy_CSNV_20201011.txt",sep='\t',header = T)
    ## output files
    pdf("./Fig4_1000randomctrl_v20201022.pdf",width=16,height=16)
    write.csv(data_matrix, file = "FigS6a.csv")
    write.csv(data_matrix, file = "FigS6c.csv")
    pdf("./FigS5_61corectrl_v20201022.pdf",width=16,height=16)
    write.csv(data_matrix, file = "FigS6b.csv")
    write.csv(data_matrix, file = "FigS6d.csv")
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.
