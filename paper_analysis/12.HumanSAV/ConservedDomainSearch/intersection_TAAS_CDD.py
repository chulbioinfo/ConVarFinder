import sys, os, glob

TAAS_dic = {}
fpin = open("AVL_HUMAN_TAAS.txt",'r')
fpin.readline() #Symbol	Ensembl_stable_ID(v93)	Ortho_ID_(primates)	Position_AA_MSA	Position_AA_HOMSA
for line in fpin:
  part = line.strip().split("\t")
  nSymbol = part[0]
  posHUMAN= int(part[4])
  TAAS_dic.setdefault(nSymbol,[])
  TAAS_dic[nSymbol].append(posHUMAN)
fpin.close()


fpout  =  open("interection_TAAS_CDD.txt",'w')

CDD_dic = {}
fpin = open("cdd.txt",'r')
fpin.readline() #Symbol	Ensembl_stable_ID(v93)	Ortho_ID_(primates)	Position_AA_MSA	Position_AA_HOMSA
for line in fpin:
  part = line.strip().split("\t")
  nSymbol = part[0]
  posS = int(part[3])
  posE = int(part[4])
  CDD = part[8]
  if nSymbol in TAAS_dic.keys():
    for posHUMAN in TAAS_dic[nSymbol]:
      if posS <= posHUMAN:
        if posHUMAN <= posE:
          tmpCDD = CDD + "("+str(posS)+"-"+str(posE)+")"
          fpout.write(nSymbol+"_"+str(posHUMAN)+'\t'+tmpCDD+'\n')
fpin.close()
fpout.close()
