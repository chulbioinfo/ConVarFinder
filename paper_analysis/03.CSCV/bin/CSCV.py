# CSCV program to detect convergent single codon variants
# Version 1.0 (14.Nov.2020)
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





print("CSCV analysis"+'\n'+'\n')
 
#from scipy import stats
from datetime import datetime
import os, glob, sys, math, time
import scipy.stats as stats
startTime = datetime.today()
print("Start time")
print(startTime)
################################ Global variable ############################
nTargets = 'TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN'
#nTargets = 'FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU' # birds of prey
#nTargets = 'PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL' # waterbirds

## setting date
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

## Amino acids: must not include '.' in the sequences, because '.' is used as a marker of loss of gene
AA_SET        = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","U","O","-",'*',"X"] 
CODONTABLE_dic  = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'Z','TAG':'Z','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'Z','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',"---":"-",'CTN':'L','GTN':'V','TCN':'S','CCN':'P','ACN':'T','GCN':'A','CGN':'R','GGN':'G'} # Stop codon = 'Z', NNN = 'X'


## output NAME
fNAME_SUMMARY = "CSCV_output_" + nTargets + "_" + nDate + ".txt"
FPOUT_SUMMARY = open(oPATH + fNAME_SUMMARY,'w')
tmpline = "Gene" + '\t' + "Position_Nuc" + '\t' + "CSCV_Type" + '\t' + "TargetCODON" + '\t' + "OthersCODON" + '\t' + "Pvalue_CODON" + '\t' + "Adj.pvalue_CODON" + '\t' + "Total_Entropy_CODON" + '\t' + "Target_Entropy_CODON" + '\t' + "Others_Entropy_CODON" + "\t"+ "\t".join(sID_LIST)+'\t' + "Position_AA" + '\t' + "CSAV" + '\t' + "TargetAA" + '\t' + "OthersAA" + '\t' + "Pvalue_AA" + '\t' + "Adj.pvalue_AA" + '\t' + "Total_Entropy_AA" + '\t' + "Target_Entropy_AA" + '\t' + "Others_Entropy_AA" + "\t"+ "\t".join(sID_LIST)+'\n'
FPOUT_SUMMARY.write(tmpline)


################################ Funtions #################################

########
# CSCV #
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
    #oddsratio, pvalue = "#N/A","#N/A" #stats.fisher_exact(contingency_table)
    return ''.join(targetAA_set),''.join(othersAA_set), oddsratio, pvalue

def EnrichmentTEST_CODON(targetCODON_list, othersCODON_list):
    targetCODON_set = []
    for targetCODON in targetCODON_list:
        if not targetCODON == "...":
            targetCODON_set.append(targetCODON)
    othersCODON_set = []
    for othersCODON in othersCODON_list:
        if not othersCODON == "...":
            othersCODON_set.append(othersCODON)
    targetCODON_set = list(set(targetCODON_set))
    targetCODON_set.sort()
    othersCODON_set = list(set(othersCODON_set))
    othersCODON_set.sort()
    CNT_target_targetCODON     = 0
    CNT_target_othersCODON     = 0
    CNT_others_targetCODON     = 0
    CNT_others_othersCODON     = 0
    for nCODON in CODONTABLE_dic.keys():
        if nCODON in targetCODON_set:
            targetCODON = nCODON
            CNT_target_targetCODON += targetCODON_list.count(targetCODON)  
            CNT_others_targetCODON += othersCODON_list.count(targetCODON)
        else:
            othersCODON = nCODON
            CNT_target_othersCODON += targetCODON_list.count(othersCODON)
            CNT_others_othersCODON += othersCODON_list.count(othersCODON)
    a = CNT_target_targetCODON
    b = CNT_target_othersCODON
    c = CNT_others_targetCODON
    d = CNT_others_othersCODON
    contingency_table = [[a,b],[c,d]]
    # Fisher Exact
    oddsratio, pvalue = stats.fisher_exact(contingency_table)
    #oddsratio, pvalue = "#N/A","#N/A" #stats.fisher_exact(contingency_table)
    return ','.join(targetCODON_set),','.join(othersCODON_set), oddsratio, pvalue


def TAASanalysis(targetAA_list, othersAA_list):
    targetAA_set = list(set(targetAA_list))
    othersAA_set = list(set(othersAA_list))
    TAAS  = ""
    iFlag = 0
    CNT_targetETC = 0
    for targetAA in targetAA_set:
        if targetAA in AA_SET:
            if targetAA in othersAA_set:
                iFlag = 1
        else:
            CNT_targetETC += 1
    if iFlag == 0:
        TAAS = "MutuallyExclusive"
    return TAAS


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

def TAASanalysis_CODON(nGene, iPos, targetCODON_list, othersCODON_list):
    targetCODON_set = list(set(targetCODON_list))
    othersCODON_set = list(set(othersCODON_list))
    TAAS_CODON  = ""
    iFlag = 0
    CNT_targetETC = 0
    for targetCODON in targetCODON_set:
        if targetCODON in CODONTABLE_dic.keys():
            if not targetCODON == "...": 
                if targetCODON in othersCODON_set:
                    iFlag = 1
            else:
                CNT_targetETC += 1
                
    if iFlag == 0:
        TAAS_CODON = "MutuallyExclusive"
    return TAAS_CODON

def EntropyTEST_CODON(CODON_list):
    CNT_CODONs = float(0)
    for nCODON in CODON_list:
        if nCODON in CODONTABLE_dic.keys():
            CNT_CODONs += 1
            
    Entropy = 0
    tmpCODON_set = list(set(CODON_list))
    for nCODON in tmpCODON_set:
        if nCODON in CODONTABLE_dic.keys():
            freqCODON = CODON_list.count(nCODON)/CNT_CODONs
            if freqCODON == 0:
                entropyCODON = 0
            else:
                entropyCODON = -freqCODON * math.log(freqCODON,2)
            
            Entropy += entropyCODON
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

def RUN_TAAS_analysis(sID_list, sID_nSeq_matrix, oNAME):
    
    print_geneNAME = os.path.basename(oNAME).split(seqFORM)[0]
    
    print_Seq_Len  = 0
    print_Entropy  = 0
    print_Num_TAAS = 0
    
    seq_len 		= len(sID_nSeq_matrix[0])
    print_Seq_Len 	= seq_len
    BonferroniALPHA     = 0.05/seq_len
	
    nGene = os.path.basename(oNAME).split("_")[0]
    iPos_AA = 0
    for iColumn in range(0,seq_len,3):
        
        iPos = str(iColumn + 1)
        iPos_AA = int(iPos_AA)
        iPos_AA +=1  

        # Target
        targetCODON_list = []
        targetAA_list = []
        for iRow in range(CNT_target):
            targetCODON = sID_nSeq_matrix[iRow][iColumn:(iColumn+3)]
            targetCODON_list.append(targetCODON)
            try:
                targetAA = CODONTABLE_dic[targetCODON]
                
            except:
                if targetCODON == "...":
                    targetAA = "."
                else:
                    targetAA = "X"
            targetAA_list.append(targetAA)
            

        # Others
        othersCODON_list = []
        othersAA_list = []
        for iRow in range(CNT_target,CNT_target+CNT_others):
            othersCODON = sID_nSeq_matrix[iRow][iColumn:(iColumn+3)]
            othersCODON_list.append(othersCODON)
            try:
                othersAA = CODONTABLE_dic[othersCODON]    
            except:
                if othersCODON == "...":
                    othersAA = "."
                else:
                    othersAA = "X"
            othersAA_list.append(othersAA)
            

        # Outgroup
        outgroupCODON_list = []
        outgroupAA_list = []
        for iRow in range(len(sID_LIST)-len(outgroup_LIST),len(sID_LIST)):
            outgroupCODON = sID_nSeq_matrix[iRow][iColumn:(iColumn+3)]
            outgroupCODON_list.append(outgroupCODON)
            try:
                outgroupAA = CODONTABLE_dic[outgroupCODON]
            except:
                if outgroupCODON == "...":
                    outgroupAA = "."
                else:
                    outgroupAA = "X"
            outgroupAA_list.append(outgroupAA)

        # CSCV analysis    
        #############################################
        # STEP1. Fisher's exact test at codon level #
        #############################################
        targetCODONs, othersCODONs, oddsratio_CODON, pvalue_CODON = EnrichmentTEST_CODON(targetCODON_list, othersCODON_list)
        adj_pvalue_CODON = pvalue_CODON*seq_len
        
        #######################################
        # STEP2. TAAS analysis at codon level #
        #######################################
        TAAS_CODON = TAASanalysis_CODON(nGene, iPos, targetCODON_list, othersCODON_list)
        
        #####################################
        # STEP3. EntropyTest at codon level #
        #####################################
        Total_Entropy_CODON  = EntropyTEST_CODON(targetCODON_list + othersCODON_list)
        Target_Entropy_CODON = EntropyTEST_CODON(targetCODON_list)
        Others_Entropy_CODON = EntropyTEST_CODON(othersCODON_list)
        if not TAAS_CODON == "":
            if Target_Entropy_CODON == 0:
                if Others_Entropy_CODON == 0:
                    TAAS_CODON = "1"
                else: # Others_Entropy_CODON > 0
                    TAAS_CODON = "2"
            else: # Target_Entropy_CODON > 0:
                if Others_Entropy_CODON == 0:
                    TAAS_CODON = "3"
                else: # Others_Entropy > 0
                    TAAS_CODON = "4"


        # CSAV analysis

        ##############################
        # STEP1. Fisher's exact test #
        ##############################
        targetAAs, othersAAs, oddsratio, pvalue = EnrichmentTEST(targetAA_list, othersAA_list)
        adj_pvalue = pvalue*seq_len

        ########################
        # STEP2. TAAS analysis #
        ########################
        TAAS = TAASanalysis(targetAA_list, othersAA_list)
        
        ######################
        # STEP3. EntropyTest #
        ######################
        Total_Entropy  = EntropyTEST(targetAA_list + othersAA_list)
        Target_Entropy = EntropyTEST(targetAA_list)
        Others_Entropy = EntropyTEST(othersAA_list)
        targetAAs = targetAAs.replace('.','')
        othersAAs = othersAAs.replace('.','')
        if not TAAS == "":
            if Target_Entropy == 0:
                if Others_Entropy == 0:
                    TAAS = "1"
                else: # Others_Entropy > 0
                    TAAS = "2"
            else: # Target_Entropy > 0:
                if Others_Entropy == 0:
                    TAAS = "3"
                else: # Others_Entropy > 0
                    TAAS = "4"
                    
        iPos = str(iPos)
        iPos_AA = str(iPos_AA)
        Result_position = nGene + '\t' + iPos + '\t' + TAAS_CODON + '\t' + targetCODONs+'\t' + othersCODONs + '\t'+str(pvalue_CODON)+'\t'+str(adj_pvalue_CODON)+ '\t' + str(Total_Entropy_CODON)+ '\t' + str(Target_Entropy_CODON)+ '\t' + str(Others_Entropy_CODON)
        
        print_Entropy += Total_Entropy
        if not TAAS_CODON == "":
            
                print_Num_TAAS += 1
                totalCODON_list = targetCODON_list + othersCODON_list + outgroupCODON_list
                for nCODON in totalCODON_list:
                    if nCODON == "...":
                        Result_position += '\t' +""
                    else:
                        Result_position += '\t' +nCODON
                Result_position += '\t' + iPos_AA + '\t' + TAAS + '\t' + targetAAs+'\t' + othersAAs + '\t'+str(pvalue)+'\t'+str(adj_pvalue)+ '\t' + str(Total_Entropy)+ '\t' + str(Target_Entropy)+ '\t' + str(Others_Entropy)
                totalAA_list = targetAA_list + othersAA_list + outgroupAA_list
                for nAA in totalAA_list:
                    if nAA == ".":
                        Result_position += '\t' +""
                    else:
                        Result_position += '\t' +nAA
                FPOUT_SUMMARY.write(Result_position+'\n')
                
    
################################# Main ####################################


# CSCV analysis
CNT_TEST = 0
flist = glob.glob(seqPATH+"*"+seqFORM)
iCNT_ortho = len(flist)
print("Number of input genes : " + str(iCNT_ortho))
print("")
print("Running CSCV")
iProcess = 0
iPercent = 10
for fNAME in flist:
    iProcess += 1
    if iCNT_ortho >= 10:
        if iProcess in range(int(iCNT_ortho/10),int(iCNT_ortho),int(iCNT_ortho/10)):
            print("Processing",iPercent,"%")
            iPercent += 10
    sID_list, sID_nSeq_matrix = MAKE_SEQ_DB(fNAME)
    oNAME = oPATH + os.path.basename(fNAME)
    RUN_TAAS_analysis(sID_list, sID_nSeq_matrix, oNAME)
FPOUT_SUMMARY.close()


endTime = datetime.today()
print("Running time")
print(endTime - startTime)
print(".........................................Finish CSCV analysis")
