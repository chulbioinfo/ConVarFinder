# CSAV program
# Version 1.0 (16.Nov.2020)
# written by Chul Lee (e-mail: chul.bioinfo@gmail.com)
# This code was developed and conducted in Python 3.7.1 (v3.7.1:260ec2c36a, Oct 20 2018, 14:57:15) [MSC v.1915 64 bit (AMD64)]
# OS for development and analysis: Windows 10 Education

# Libraries
import scipy.stats as stats
from datetime import datetime
import os
import glob
import sys
import math
import time


# Global variable
## input
seqPATH = "../../0.rawdata/MSA/cds/"
seqFORM = ".sate.default.pep2cds.removed.shortname.fasta"
### list of species that is analyzed in this code: If species is not included in this list, it would not be included as the other species
sID_LIST = ['TAEGU','GEOFO','CORBR','MELUN','NESNO','CALAN','MANVI','FALPE','CARCR','MERNU','PICPU','BUCRH','APAVI','LEPDI','COLST','TYTAL','HALLE','HALAL','CATAU','PELCR','EGRGA','NIPNI','PHACA','FULGL','PYGAD','APTFO','GAVST','PHALE','EURHE','CHAVO','BALRE','OPHHO','CHAPE','CAPCA','CHLUN','TAUER','CUCCA','MESUN','PTEGU','COLLI','PHORU','PODCR','GALGA','MELGA','ANAPL','TINMA','STRCA','HUMAN','ACACH']
### list of target species with a convergent trait 
nTargets = 'TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN'
### list of outgroup species excluded in the other species but recorded with amino acid of the outgroup species
outgroup_LIST = ['HUMAN','ACACH']
## output PATH
oPATH   = "../output/"
## making directory of output PATH
#os.system("mkdir "+oPATH)

# Automated variables
targetID_LIST = nTargets.split(',')

othersID_LIST  = []
for sID in sID_LIST:
 if not sID in targetID_LIST:
  if not sID in outgroup_LIST:
   othersID_LIST.append(sID)
sID_LIST = []
sID_LIST = targetID_LIST + othersID_LIST + outgroup_LIST

tmpTargets = nTargets.replace(',','_')

CNT_target    = len(targetID_LIST)
CNT_others    = len(othersID_LIST)

## setting start time
startTime = datetime.today()
print("Start time")
print(startTime)

nYear   = str(datetime.today().year)
nMonth  = str(datetime.today().month)
nDay    = str(datetime.today().day)
nHour   = str(datetime.today().hour)
nMin    = str(datetime.today().minute)
if len(nMonth)==1:
    nMonth = "0" + nMonth
if len(nDay)==1:
    nDay = "0" + nDay
if len(nHour)==1:
    nHour = "0" + nHour
if len(nMin)==1:
    nMin = "0" + nMin
nDate   = nYear + nMonth + nDay + "_" + nHour + nMin

## Nucleotides: must not include '.' in the sequences, because '.' is used as a marker of loss of gene
AA_SET        = ["A","T","G","C","N","-",'*',"X"] 

## output NAME
fNAME_SUMMARY = "CSNV_output_" + nTargets + "_" + nDate + ".txt"

# writing output
FPOUT_SUMMARY = open(oPATH + fNAME_SUMMARY,'w')
tmpline = "Gene" + '\t' + "Position" + '\t' + "CSNV_Type" + '\t' + "TargetNUC" + '\t' + "OthersNUC" + '\t' + "Pvalue" + '\t' + "Adj.pvalue" + '\t' + "Total_Entropy" + '\t' + "Target_Entropy" + '\t' + "Others_Entropy" + "\t"+ "\t".join(sID_LIST)+'\n'
FPOUT_SUMMARY.write(tmpline)


################################ Funtions #################################

########
# CSAV #
########

def EnrichmentTEST(targetAA_list, othersAA_list):
    targetAA_set = list(set(targetAA_list))
    targetAA_set.sort()
    othersAA_set = list(set(othersAA_list))
    othersAA_set.sort()
    CNT_target_targetAA     = 0
    CNT_target_othersAA     = 0
    CNT_others_targetAA     = 0
    CNT_others_othersAA     = 0
    for nAA in AA_SET:
        if nAA in targetAA_set:
            targetAA = nAA
            CNT_target_targetAA += targetAA_list.count(targetAA)  
            CNT_others_targetAA += othersAA_list.count(targetAA)
        else:
            othersAA = nAA
            CNT_target_othersAA += targetAA_list.count(othersAA)
            CNT_others_othersAA += othersAA_list.count(othersAA)
    a = CNT_target_targetAA
    b = CNT_target_othersAA
    c = CNT_others_targetAA
    d = CNT_others_othersAA
    contingency_table = [[a,b],[c,d]]
    # Fisher Exact
    oddsratio, pvalue = stats.fisher_exact(contingency_table)
    return ''.join(targetAA_set),''.join(othersAA_set), oddsratio, pvalue


def CSAVanalysis(targetAA_list, othersAA_list):
    targetAA_set = list(set(targetAA_list))
    othersAA_set = list(set(othersAA_list))
    CSAV  = ""
    iFlag = 0
    CNT_targetETC = 0
    for targetAA in targetAA_set:
        if targetAA in AA_SET:
            if targetAA in othersAA_set:
                iFlag = 1
        else:
            CNT_targetETC += 1
    if iFlag == 0:
        CSAV = "MutuallyExclusive"
    return CSAV

	
def EntropyTEST(AA_list):
    CNT_AAs = float(0)
    for nAA in AA_list:
        if nAA in AA_SET:
            CNT_AAs += 1
            
    Entropy = 0
    tmpAA_set = list(set(AA_list))
    for nAA in tmpAA_set:
        if nAA in AA_SET:
            freqAA = AA_list.count(nAA)/CNT_AAs
            if freqAA == 0:
                entropyAA = 0
            else:
                entropyAA = -freqAA * math.log(freqAA,2)
            Entropy += entropyAA
    return Entropy


def MAKE_SEQ_DB(fNAME):
    sID_list = []
    sID_nSeq_dic = {}
    fpin = open(fNAME,'r')
    for line in fpin:
        if line[0] == ">":
            sID = line[1:].strip()
            sID_nSeq_dic.setdefault(sID,'')
            sID_list.append(sID)
        else:
            sID_nSeq_dic[sID]+= line.strip().replace(' ','')
    fpin.close()
    seq_len = len(sID_nSeq_dic[sID])
    tmpline = ''
    for i in range(seq_len):
        tmpline += "."
    sID_nSeq_matrix = []
    for sID in sID_LIST:
        try:
            sID_nSeq_matrix.append(sID_nSeq_dic[sID])
        except:
            sID_nSeq_matrix.append(tmpline)
    return sID_list, sID_nSeq_matrix


def RUN_CSAV_analysis(sID_list, sID_nSeq_matrix, oNAME):
    tmp_nGene_iCNT_CSAV_dic = nGene_iCNT_CSAV_dic

    print_geneNAME = os.path.basename(oNAME).split(seqFORM)[0]
    
    print_Seq_Len  = 0
    print_Entropy  = 0
    print_Num_TAAS = 0
    
    seq_len 		= len(sID_nSeq_matrix[0])
    print_Seq_Len 	= seq_len
    BonferroniALPHA     = 0.05/seq_len
	
    nGene = os.path.basename(oNAME).split("_")[0]
    for iColumn in range(seq_len):
        iPos = str(iColumn + 1)

        # Target
        targetAA_list = []
        for iRow in range(CNT_target):
            targetAA = sID_nSeq_matrix[iRow][iColumn].upper()
            targetAA_list.append(targetAA)

        # Others
        othersAA_list = []
        for iRow in range(CNT_target,CNT_target+CNT_others):
            othersAA = sID_nSeq_matrix[iRow][iColumn].upper()
            othersAA_list.append(othersAA)

        # Outgroup
        outgroupAA_list = []
        for iRow in range(len(sID_LIST)-len(outgroup_LIST),len(sID_LIST)):
            outgroupAA = sID_nSeq_matrix[iRow][iColumn].upper()
            outgroupAA_list.append(outgroupAA)

        ##############################
        # STEP1. Fisher's exact test #
        ##############################
        targetAAs, othersAAs, oddsratio, pvalue = EnrichmentTEST(targetAA_list, othersAA_list)
        adj_pvalue = pvalue*seq_len

        ########################
        # STEP2. CSAV analysis #
        ########################
        CSAV = CSAVanalysis(targetAA_list, othersAA_list)
        
        ######################
        # STEP3. EntropyTest #
        ######################
        Total_Entropy  = EntropyTEST(targetAA_list + othersAA_list)
        Target_Entropy = EntropyTEST(targetAA_list)
        Others_Entropy = EntropyTEST(othersAA_list)
        targetAAs = targetAAs.replace('.','')
        othersAAs = othersAAs.replace('.','')
        if not CSAV == "":
            if Target_Entropy == 0:
                if Others_Entropy == 0:
                    CSAV = "1"
                else: # Others_Entropy > 0
                    CSAV = "2"
            else: # Target_Entropy > 0:
                if Others_Entropy == 0:
                    CSAV = "3"
                else: # Others_Entropy > 0
                    CSAV = "4"
        Result_position = nGene + '\t' + iPos + '\t' + CSAV + '\t' + targetAAs+'\t' + othersAAs + '\t'+str(pvalue)+'\t'+str(adj_pvalue)+ '\t' + str(Total_Entropy)+ '\t' + str(Target_Entropy)+ '\t' + str(Others_Entropy)
        print_Entropy += Total_Entropy
        if not CSAV == "":
                print_Num_TAAS += 1
                totalAA_list = targetAA_list + othersAA_list + outgroupAA_list
                for nAA in totalAA_list:
                    if nAA == ".":
                        Result_position += '\t' +""
                    else:
                        Result_position += '\t' +nAA
                FPOUT_SUMMARY.write(Result_position+'\n')
                

    
################################# Main ####################################


# CSAV analysis
CNT_TEST = 0
flist = glob.glob(seqPATH+"*"+seqFORM)
iCNT_ortho = len(flist)
print("Number of orthologous genes: " + str(iCNT_ortho))
print("")
print("Running CSNV analysis")
iProcess = 0
iPercent = 10
nGene_iCNT_CSAV_dic = {}
for fNAME in flist:
    iProcess += 1
    if iCNT_ortho >= 10:
        if iProcess in range(int(iCNT_ortho/10),int(iCNT_ortho),int(iCNT_ortho/10)):
            print("Processing",iPercent,"%")
            iPercent += 10
    sID_list, sID_nSeq_matrix = MAKE_SEQ_DB(fNAME)
    oNAME = oPATH + os.path.basename(fNAME)
    RUN_CSAV_analysis(sID_list, sID_nSeq_matrix, oNAME)
FPOUT_SUMMARY.close()



endTime = datetime.today()
print("Running time")
print(endTime - startTime)
print(".........................................Finish CSNV analysis")

