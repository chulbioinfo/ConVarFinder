import sys
import os
import glob
fNAME_maf_list = glob.glob("363-avian-2020_Taeniopygia_guttata_chr*_*-*.maf")

sID_list = []
fpin = open("sID_list.txt",'r')
for line in fpin:
    sID = line.strip()
    if len(sID)>=1:
        sID_list.append(sID)
fpin.close()


def reorder_maf(fNAME_maf):
    fNAME_maf = os.path.basename(fNAME_maf)
    fpout_maf = open(fNAME_maf+".ordered.maf",'w')
    fpout_fasta = open(fNAME_maf+".ordered.fasta",'w')
    fpout_maf.write("##maf version=1 scoring=N/A\n\na\n")
    

    fpin = open(fNAME_maf,'r')
    iFlag = 0
    for line in fpin:
        if line[0]=="a":
            if iFlag == 0:
                iFlag = 1
                sID_nInfo_dic = {}
            else:
                for sID in sID_list:
                    if sID in sID_nInfo_dic.keys():
                        nInfo = sID_nInfo_dic[sID]
                        nSeq = nInfo.strip().split("\t")[-1]
                        tmpline_maf = nInfo
                        fpout_maf.write(tmpline_maf)
                        tmpline_fasta = ">"+sID+'\n'+nSeq+'\n'
                        fpout_fasta.write(tmpline_fasta)
                sID_nInfo_dic = {}
                
        elif line[0]=="s":
            part = line.strip().split('\t')
            sID = part[1].split('.')[0]
            nInfo = line
            sID_nInfo_dic.setdefault(sID,nInfo)
        else:
            pass
    fpin.close()

    for sID in sID_list:
        if sID in sID_nInfo_dic.keys():
            nInfo = sID_nInfo_dic[sID]
            nSeq = nInfo.strip().split("\t")[-1]
            tmpline_maf = nInfo
            fpout_maf.write(tmpline_maf)
            tmpline_fasta = ">"+sID+'\n'+nSeq+'\n'
            fpout_fasta.write(tmpline_fasta)
    fpout_maf.close()
    fpout_fasta.close()

for fNAME_maf in fNAME_maf_list:
    reorder_maf(fNAME_maf)
