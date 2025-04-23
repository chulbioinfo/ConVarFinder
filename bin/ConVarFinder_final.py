# ConVarFinder is a python script to detect convergent variants at codon and amino acid levels
# version v20210721.1
# Written by Chul Lee (c) Seoul National Univ. email: chul.bioinfo@gmail.com


# Library
import sys
import os
import glob

# Global variables
AA_SET          = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","U","O","-",'*',"X"] 
CODONTABLE_dic  = {'NNN':'X','TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'Z','TAG':'Z','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'Z','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',"---":"-",'CTN':'L','GTN':'V','TCN':'S','CCN':'P','ACN':'T','GCN':'A','CGN':'R','GGN':'G'} # Stop codon = 'Z', NNN = 'X'
NUC_SET         = ["A","T","G","C","N","X","-"]


# Functions
def read_seq_file(fNAME_seq):
  sID_nSeq_dic = {}
  fpin = open(fNAME_seq,'r')
  for line in fpin:
    if line[0]==">":
      sID = line[1:].strip('\n')
      sID_nSeq_dic.setdefault(sID,'')
    else:
      sID_nSeq_dic[sID]+=line.strip()  
  fpin.close()
  return(sID_nSeq_dic)


def read_tree_file(fNAME_tree):
  fpin = open(fNAME_tree,'r')
  nTree = ''
  for line in fpin:
    nTree += line.strip()
  fpin.close()
  return(nTree)

    
def read_tre_removing_branchlength(nTree):
  line = nTree.strip(";")
  if ":" in line:
    tmp_node_winfo_list = []
    node_winfo_list = line.split(")")
    for node_winfo in node_winfo_list:
      if ',' in node_winfo:
        nodeID_list = []
        for eachnode_winfo in node_winfo.split(","):
          nodeID = eachnode_winfo.split(":")[0]
          nodeID_list.append(nodeID)
        tmp_node_winfo_list.append(','.join(nodeID_list))
      else:
        nodeID = node_winfo.split(":")[0]
        tmp_node_winfo_list.append(nodeID)
    line = ')'.join(tmp_node_winfo_list)+';'
  return line


def terminal_node_from_tree(nTree):
  ## list of terminal nodes
  terNode_list = []
  for terInfo in nTree.split(","):
    if "(" in terInfo:
      sID = terInfo.split("(")[-1]
    elif ")" in terInfo:
      sID = terInfo.split(")")[0]
    else:
      pass
    terNode_list.append(sID)
  return(terNode_list)


def ancestral_node_from_tree(nTree):
  ## list of ancestral nodes
  ancNode_list = []
  for ancInfo in nTree.split(")"):
    if "(" in ancInfo:
      pass
    else:
      if ',' in ancInfo:
        sID = ancInfo.split(",")[0]
      elif ';' in ancInfo:
        sID = ancInfo.split(";")[0]
      else:
        sID = ancInfo.strip()
      ancNode_list.append(sID)
  return(ancNode_list)


def reordering_targetID_list(terNode_list, targetID_list, MonophyleticPair_list):
  tmp_targetID_list = []
  mono_list = []
  poly_list = []
  for targetID in terNode_list:
    if targetID in targetID_list:
      if targetID in MonophyleticPair_list:
        mono_list.append(targetID)
      else:
        poly_list.append(targetID)
  for targetID in mono_list:
    tmp_targetID_list.append(targetID)
  for targetID in poly_list:
    tmp_targetID_list.append(targetID)
  return(tmp_targetID_list)


def make_mID_sID_dic(partial_nTree, mID_sID_dic): ## dictionary maternal and sister nodes
  # mother ID
  mID = partial_nTree.split(")")[-1]
  mID = mID.replace(";",'')
  
  # set partial_tree excluding root
  tmp_partial_nTree = partial_nTree[1:].split(")"+mID)[0]

  # set sister ID 1 and 2
  if not ")" in tmp_partial_nTree:            # case: "s1 , s2"    (mID = a1)
    sID_1 = tmp_partial_nTree.split(",")[0]
    sID_2 = tmp_partial_nTree.split(",")[1]
    mID_sID_dic.setdefault(mID,[sID_1,sID_2])
    return(mID_sID_dic)
  else: #"( and )" is present in partial tree
    sID_info_list = tmp_partial_nTree.split(",")
    FirstID_info  = sID_info_list[0]
    LastID_info   = sID_info_list[-1]
    if FirstID_info.count("(") > 0:
      if LastID_info.count(")") > 0:          # case "(s1,s2)a1 , (s3,s4)a2"
        p1_tree = FirstID_info
        for i in range(1,len(sID_info_list),1):
          p1_tree += ','+sID_info_list[i]
          if p1_tree.count("(") == p1_tree.count(")"):
            p2_tree = tmp_partial_nTree.split(p1_tree+',')[1]
            if not p2_tree.count("(") == p2_tree.count(")"):
              print("# Errorneous_Topology_Of_2nd_sister_lineage:",tmp_partial_nTree)
              sys.exit()
            sID_1 = p1_tree.split(")")[-1]
            sID_2 = p2_tree.split(")")[-1]
            break
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p1_tree, mID_sID_dic)
        mID_sID_dic = make_mID_sID_dic(p2_tree, mID_sID_dic)
        return(mID_sID_dic)
      else:                                   # case "(s1,s2)a1 , s3"
        sID_2 = sID_info_list[-1]
        p1_tree = tmp_partial_nTree.split(","+sID_2)[0]
        sID_1 = p1_tree.split(")")[-1]
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p1_tree, mID_sID_dic)
        return(mID_sID_dic)
    else:
      if LastID_info.count(")") > 0:          # case "s1 , (s2,s3)a1"
        sID_1 = sID_info_list[0]
        p2_tree = tmp_partial_nTree.split(sID_1+",")[1]
        sID_2 = p2_tree.split(")")[-1]
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p2_tree, mID_sID_dic)
        return(mID_sID_dic)
      else:
        print("Error-parsing_tree:",tmp_partial_nTree)
        sys.exit()


def make_sID_mID_dic(mID_sID_dic):
  sID_mID_dic = {}
  for mID in mID_sID_dic.keys():
    for sID in mID_sID_dic[mID]:
      sID_mID_dic.setdefault(sID,mID)
  return(sID_mID_dic)


def remove_ID_from_list(tmp_targetID_list,tmp_targetID):
  iPos_LastElement  = len(tmp_targetID_list)-1
  iPos_tmp_targetID = tmp_targetID_list.index(tmp_targetID) 
  tmp_targetID_list[iPos_tmp_targetID]  = tmp_targetID_list[iPos_LastElement]
  tmp_targetID_list[iPos_LastElement]   = tmp_targetID
  tmp_targetID_list = tmp_targetID_list[:-1]
  return(tmp_targetID_list)
  
  
def find_monophyletic_clade(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic):
  ## First time of the reference target
  iFlag = 0
  iFlag_End_MARCA_targetID = False
  if not ref_targetID in targetID_mID_dic.keys():
    for tmp_targetID in targetID_list:
      if not ref_targetID == tmp_targetID:
        ### case finding monophyletic pairs within target list
        if sID_mID_dic[ref_targetID] == sID_mID_dic[tmp_targetID]:
          iFlag= 1
          if ref_targetID in targetID_list:
            targetID_list = remove_ID_from_list(targetID_list, ref_targetID)
          targetID_list = remove_ID_from_list(targetID_list, tmp_targetID)
          mID = sID_mID_dic[ref_targetID]
          targetID_mID_dic.setdefault(ref_targetID,mID)
          targetID_mID_dic.setdefault(tmp_targetID,mID)
          return(targetID_list, targetID_mID_dic, iFlag_End_MARCA_targetID)
    ### case without monophyletic pairs within target list
    if iFlag == 0:
      mID = ref_targetID
      targetID_mID_dic.setdefault(ref_targetID,mID)
      targetID_list = remove_ID_from_list(targetID_list, ref_targetID)
      iFlag_End_MARCA_targetID = True
      return(targetID_list, targetID_mID_dic, iFlag_End_MARCA_targetID)
  ## Next time of the reference target
  else:
    ref_mID = targetID_mID_dic[ref_targetID]
    for tmp_targetID in targetID_list:
      ### case finding additional monophyletic species in target list
      if sID_mID_dic[ref_mID] == sID_mID_dic[tmp_targetID]:
        iFlag = 1
        targetID_list = remove_ID_from_list(targetID_list, tmp_targetID)
        mID = sID_mID_dic[ref_mID]
        for targetID in targetID_mID_dic.keys():
          targetID_mID_dic[targetID] = mID
        targetID_mID_dic.setdefault(tmp_targetID,mID)
        return(targetID_list, targetID_mID_dic, iFlag_End_MARCA_targetID)
    ### case without additional monophyletic species within target list
    if iFlag == 0:
      iFlag_End_MARCA_targetID = True
      return(targetID_list, targetID_mID_dic, iFlag_End_MARCA_targetID)

  
def make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list):
  targetID_list, targetID_mID_dic, iFlag_End_MARCA_targetID = find_monophyletic_clade(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic)
  
  if iFlag_End_MARCA_targetID == True:
    for targetID in targetID_mID_dic.keys():
        mID = targetID_mID_dic[targetID]
        mID_targetID_dic.setdefault(mID,[])
        mID_targetID_dic[mID].append(targetID)
    targetID_mID_dic = {}
    if len(targetID_list) > 0:
      targetID_list = reordering_targetID_list(terNode_list, targetID_list ,MonophyleticPair_list)
    return(targetID_list, targetID_mID_dic, mID_targetID_dic)
  else:
    ### Recursive function
    targetID_list, targetID_mID_dic, mID_targetID_dic = make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list)
    targetID_list = reordering_targetID_list(terNode_list, targetID_list,MonophyleticPair_list)
    return(targetID_list, targetID_mID_dic, mID_targetID_dic)


def scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list):
  ref_targetID = targetID_list[0]
  targetID_list, targetID_mID_dic, mID_targetID_dic = make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list)
  if len(targetID_list) > 0:
    ### Recursive function
    targetID_list, targetID_mID_dic , mID_targetID_dic = scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list)
  return(targetID_list, targetID_mID_dic , mID_targetID_dic)


def parsing_MonophyleticPair_list(mID_sID_dic,terNode_list):
  MonophyleticPair_list = []
  for mID in mID_sID_dic.keys():
    iCNT_terNode = 0
    for sID in mID_sID_dic[mID]:
      if sID in terNode_list:
        iCNT_terNode +=1
    if iCNT_terNode == 2:
      for sID in mID_sID_dic[mID]:
        MonophyleticPair_list.append(sID)
  return(MonophyleticPair_list)
    

def Anc_Finder(nTree, nTargets):
  targetID_list = nTargets.split(",")
  nTree = read_tre_removing_branchlength(nTree)
  terNode_list = terminal_node_from_tree(nTree)
  ancNode_list = ancestral_node_from_tree(nTree)
  mID_sID_dic = {}
  mID_sID_dic = make_mID_sID_dic(nTree, mID_sID_dic)
  sID_mID_dic = make_sID_mID_dic(mID_sID_dic)
  MonophyleticPair_list = parsing_MonophyleticPair_list(mID_sID_dic,terNode_list)
  targetID_mID_dic = {}
  mID_targetID_dic = {}
  targetID_list, targetID_mID_dic , mID_targetID_dic = scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list)
  
  mID_mmID_dic = {}
  for mID in mID_targetID_dic.keys():
    if len(mID_targetID_dic[mID]) > 1:
      if "ROOT" == mID:
        mmID = "AncROOT"  
      else:
        mmID = sID_mID_dic[mID]
      mID_mmID_dic.setdefault(mID,mmID)
  return(mID_targetID_dic, mID_mmID_dic)


def parse_target_from_tree(nTargets,nTree):
  target_list = nTargets.split(",")
  terNode_list = terminal_node_from_tree(nTree)
  tmp_target_list = []
  for sID in terNode_list:
    if sID in target_list:
      tmp_target_list.append(sID)
  return(tmp_target_list)

def parse_target_others_from_list(nTargets,sID_list,nOutgroup):
  if ',' in nTargets:
    target_list = nTargets.split(",")
  else:
    target_list = [nTargets]
  if ',' in nOutgroup:
    outgroup_list = nOutgroup.aplit(",")
  else:
    outgroup_list = [nOutgroup]
  tmp_target_list = []
  tmp_others_list = []
  tmp_outgroup_list = []
  for sID in sID_list:
    if sID in outgroup_list:
      tmp_outgroup_list.append(sID)
    else:
      if sID in target_list:
        tmp_target_list.append(sID)
      else:
        tmp_others_list.append(sID)
  return(tmp_target_list, tmp_others_list, tmp_outgroup_list)

def codon2aa(nCodon):
  try:
    nAA = CODONTABLE_dic[nCodon]
  except:
    nAA = 'X'
  return(nAA)
  
def set_options():
  # default options
  #nTargets = "MELUN,NESNO,GALGA,MELGA,TINMA,STRCA"
  nTargets = "TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN"
  nOutgroup = "ACACH"
  templete_nTree = "((((((((((((((TAEGU,GEOFO)62,CORBR)61,MANVI)60,ACACH)59,(MELUN,NESNO)63)58,FALPE)57,CARCR)56,(((((((MERNU,PICPU)70,BUCRH)69,APAVI)68,LEPDI)67,COLST)66,TYTAL)65,((HALLE,HALAL)72,CATAU)71)64)55,((((((PELCR,EGRGA)78,NIPNI)77,PHACA)76,(FULGL,(PYGAD,APTFO)80)79)75,GAVST)74,(PHALE,EURHE)81)73)54,((CHAVO,BALRE)83,OPHHO)82)53,(((CALAN,CHAPE)86,CAPCA)85,((CHLUN,TAUER)88,CUCCA)87)84)52,(((MESUN,PTEGU)91,COLLI)90,(PHORU,PODCR)92)89)51,((GALGA,MELGA)94,ANAPL)93)50,(TINMA,STRCA)49)ROOT;"
  iPATH = "../RAxML_merged/"
  fNAME_seq_format = ".fasta"
  fNAME_tree_format = ".tre"
  oPATH = "./"
  # end of default options

  # imported options
  option_list = ["-tl=",    # input: List of target species
                 "-ol=",   # input: List of outgroup species
                 "-ttf=",  # input: file name of templete tree
                 "-ip=",   # input: input path
                 "-tfmt=", # input: file format of binary tree with labelled ancestral and terminal nodes (output of RAxML: )
                 "-sfmt=", # input: file format of reconstruced sequences with ancestral and terminal nodes (output of RAxML: )
                 "-op=",   # output: output path
                 "-h"]     # help

  description_list = ["START OF LOG FILE: ConVarFinder v210721",
                      "USAGE:  ./ConVarFinder.py [-options]",
                      "|------------------------------ HELP: -------------------------------+",
                      "|-h    help                                                          |",
                      "|-tl=   strings of list of target species as comma seperated (csv)   |",
                      "|-ol=   strings of list of outgroup species as comma seperated (csv) |",
                      "|-ttf=  input file of templete tree                                  |",
                      "|-tfmt=  file format of binary tree with ancestral nodes             |",
                      "|-sfmt=  file format of input with ancestral reconstructions         |",
                      "|-ip=  input path                                                    |",
                      "|-op=  output path                                                   |",
                      "+--------------------------------------------------------------------+"]

  # update options
  if len(sys.argv) > 1:
    for i in range(1,len(sys.argv),1):
      nOpt = sys.argv[i]
      if not nOpt == '-h':
        if not nOpt.split("=")[0]+'=' in option_list:
          print("# Argument error:", nOpt)
          sys.exit()
      if "-tl=" in nOpt:
        nTargets = nOpt.split("-tl=")[1]
      if "-ol=" in nOpt:
        nOutgroup = nOpt.split("-ol=")[1]
      if "-ttf=" in nOpt:
        templete_fNAME_tree = nOpt.split("-ttf=")[1]
        templete_nTree = read_tree_file(templete_fNAME_tree)
        templete_nTree = read_tre_removing_branchlength(templete_nTree)
      if "-sfmt=" in nOpt:
        fNAME_seq_format = nOpt.split("-sfmt=")[1]
      if "-tfmt=" in nOpt:
        fNAME_seq_format = nOpt.split("-tfmt=")[1]
      if "-ip=" in nOpt:
        iPATH = nOpt.split("-ip=")[1]
      if "-op=" in nOpt:
        oPATH = nOpt.split("-op=")[1]
      if "-h" in nOpt:
        for nDescription in description_list:
          print(nDescription)
        sys.exit()

  flist_seq = glob.glob(iPATH+"*"+fNAME_seq_format)
  flist_tree = glob.glob(iPATH+"*"+fNAME_tree_format)
  templete_mID_targetID_dic, templete_mID_mmID_dic = Anc_Finder(templete_nTree, nTargets)
  templete_sID_list = terminal_node_from_tree(templete_nTree)
  target_list, others_list, outgroup_list = parse_target_others_from_list(nTargets, templete_sID_list, nOutgroup)
  templete_sID_list = target_list + others_list

  # check input files
  iFlag_err = 0
  for fNAME_seq in flist_seq:
    nGene = os.path.basename(fNAME_seq).split(fNAME_seq_format)[0]
    fNAME_tree = os.path.dirname(fNAME_seq) +"/"+ nGene + fNAME_tree_format
    if os.path.isfile(fNAME_tree) == False:
      print("Error-Absent tree:",fNAME_tree)
      iFlag_err = 1
  for fNAME_tree in flist_tree:
    nGene = os.path.basename(fNAME_tree).split(fNAME_tree_format)[0]
    fNAME_seq = os.path.dirname(fNAME_tree) +"/"+ nGene + fNAME_seq_format
    if os.path.isfile(fNAME_seq) == False:
      print("Error-Absent sequence:",fNAME_seq)
      iFlag_err = 1
  if iFlag_err == 1:
    sys.exit()
  option_list = [nTargets,                  # 0
                 templete_mID_targetID_dic, # 1
                 iPATH,                     # 2
                 fNAME_seq_format,          # 3
                 fNAME_tree_format,         # 4
                 oPATH,                     # 5
                 templete_sID_list,         # 6
                 outgroup_list]             # 7
  return(option_list)


def scan_ancestors_from_mID_to_targetID(mID, targetID, sID_mID_dic, mID_to_targetID_list):
  tmp_mID = sID_mID_dic[targetID]
  mID_to_targetID_list.append(tmp_mID)
  if not mID == tmp_mID:
    tmp_mID, mID_to_targetID_list = scan_ancestors_from_mID_to_targetID(mID, tmp_mID, sID_mID_dic, mID_to_targetID_list)
  return(tmp_mID, mID_to_targetID_list)



def ConVarFinder(fNAME_seq, terNode_list, mID_targetID_dic, mID_mmID_dic, sID_mID_dic, local_target_list, templete_sID_list, fpout, outgroup_list):
  # terNode_list      : terminal nodes of gene tree of partial family tree
  # templete_sID_list : terminal nodes of raw family tree
  sID_nSeq_dic = read_seq_file(fNAME_seq)
  for sID in sID_nSeq_dic.keys():
    iLen_nSeq = len(sID_nSeq_dic[sID])
    break
  if not iLen_nSeq%3 == 0:
    print("Errorneous sequence - not codon-wise:",fNAME_seq)
    sys.exit()
  #mID_clade_sID_dic = {}
  #for mID in mID_targetID_dic.keys():
    #ancestors_targets_list = []
    #ancestorID_list, ancestors_targets_list = scan_ancestors(mID, local_target_list, sID_mID_dic, ancestors_targets_list)
    #mID_clade_sID_dic.setdefault(mID,ancestors_targets_list)
  
  
    
  iPos_Codon = -2
  iPos_AA = 0
  iPos_Nuc1 = -2
  iPos_Nuc2 = -1
  iPos_Nuc3 = 0
  for i in range(0,iLen_nSeq,3):
    iPos_Codon  += 3
    iPos_AA     += 1
    iPos_Nuc1 += 3
    iPos_Nuc2 += 3
    iPos_Nuc3 += 3
    
    # set codon and aa lists of target species and the others
    tmp_target_codon_list = []
    tmp_others_codon_list = []
    tmp_target_aa_list = []
    tmp_others_aa_list = []
    tmp_target_nuc1_list = []
    tmp_target_nuc2_list = []
    tmp_target_nuc3_list = []
    tmp_others_nuc1_list = []
    tmp_others_nuc2_list = []
    tmp_others_nuc3_list = []

    for sID in terNode_list:
      if not sID in outgroup_list:
        if sID in local_target_list:
          nCodon_target = sID_nSeq_dic[sID][i:i+3]
          tmp_target_codon_list.append(nCodon_target)
          nNuc1_target = sID_nSeq_dic[sID][i]
          nNuc2_target = sID_nSeq_dic[sID][i+1]
          nNuc3_target = sID_nSeq_dic[sID][i+2]
          tmp_target_nuc1_list.append(nNuc1_target)
          tmp_target_nuc2_list.append(nNuc2_target)
          tmp_target_nuc3_list.append(nNuc3_target)
          nAA_target = codon2aa(nCodon_target)
          tmp_target_aa_list.append(nAA_target)
        else:
          nCodon_others = sID_nSeq_dic[sID][i:i+3]
          tmp_others_codon_list.append(nCodon_others)
          nNuc1_others = sID_nSeq_dic[sID][i]
          nNuc2_others = sID_nSeq_dic[sID][i+1]
          nNuc3_others = sID_nSeq_dic[sID][i+2]
          tmp_others_nuc1_list.append(nNuc1_others)
          tmp_others_nuc2_list.append(nNuc2_others)
          tmp_others_nuc3_list.append(nNuc3_others)
          nAA_others = codon2aa(nCodon_others)
          tmp_others_aa_list.append(nAA_others)
            
    # check mutually exclusive codon substitutions
    tmp_target_codon_list = list(set(tmp_target_codon_list))
    tmp_others_codon_list = list(set(tmp_others_codon_list))
    iFlag_Codon_Sub = 0
    for nCodon_target in tmp_target_codon_list:
      if nCodon_target in tmp_others_codon_list:
        iFlag_Codon_Sub = 1

    # check mutually exclusive amino acid substitutions
    tmp_target_aa_list = list(set(tmp_target_aa_list))
    tmp_others_aa_list = list(set(tmp_others_aa_list))
    iFlag_AA_Sub = 0
    for nAA_target in tmp_target_aa_list:
      if nAA_target in tmp_others_aa_list:
        iFlag_AA_Sub = 1

    # check mutually exclusive nucleotide substitutions at 1st position in codon
    tmp_target_nuc1_list = list(set(tmp_target_nuc1_list))
    tmp_others_nuc1_list = list(set(tmp_others_nuc1_list))
    iFlag_Nuc1_Sub = 0
    for nNuc_target in tmp_target_nuc1_list:
      if nNuc_target in tmp_others_nuc1_list:
        iFlag_Nuc1_Sub = 1

    # check mutually exclusive nucleotide substitutions at 2nd position in codon
    tmp_target_nuc2_list = list(set(tmp_target_nuc2_list))
    tmp_others_nuc2_list = list(set(tmp_others_nuc2_list))
    iFlag_Nuc2_Sub = 0
    for nNuc_target in tmp_target_nuc2_list:
      if nNuc_target in tmp_others_nuc2_list:
        iFlag_Nuc2_Sub = 1

    # check mutually exclusive nucleotide substitutions at 3rd position in codon
    tmp_target_nuc3_list = list(set(tmp_target_nuc3_list))
    tmp_others_nuc3_list = list(set(tmp_others_nuc3_list))
    iFlag_Nuc3_Sub = 0
    for nNuc_target in tmp_target_nuc3_list:
      if nNuc_target in tmp_others_nuc3_list:
        iFlag_Nuc3_Sub = 1

    
    if iFlag_Codon_Sub == 0: # if Target-specific codon substitutions
      # nGene                gene name
      # nPhyRel_Amrca2mrca   phylogenetic relationship from Ancestors of MRCAs to MRCAs of each clade
      # nPhyRel_mrca2target  phylogenetic relationship from Ancestors of MRCAs to each target species

      # nPos_codon                      start position of codon substitutions (bp)
      # nType_EvoDir_Amrca2mrca_Codon   evolutionary direction between Ancestors of MRCAs and MRCAs of each clade
      # nType_EvoDir_Amrca2target_Codon  evolutionary direction from Ancestors of MRCAs to each target species
      # nType_CodonSub_target           type of codon substitutions of targets
      # nType_CodonSub_others           type of codon substitutions of others
      # nType_CodonSub                  type of codon substitutions: identical, and different
      # nCodon_Amrca2mrca               codon substitutions from Ancestors of MRCAs to MRCAs of each clade
      # nCodon_mrca2target              codon substitutions from Ancestors of MRCAs to each target species
      # nCodon_targets                  codon set of targets
      # nCodon_others                   codon set of the others
      
      # nPos_AA                      start position of codon substitutions (bp)
      # nType_EvoDir_Amrca2mrca_AA   evolutionary direction between Ancestors of MRCAs and MRCAs of each clade
      # nType_EvoDir_Amrca2target_AA  evolutionary direction from Ancestors of MRCAs to each target species
      # nType_AASub_target           type of amino acid substitutions of targets
      # nType_AASub_others           type of amino acid substitutions of others
      # nType_AASub                  type of amino acid subsitutions: identical, and different
      # nAA_Amrca2mrca               amino acid substitutions from Ancestors of MRCAs to MRCAs of each clade
      # nAA_mrca2target              amino acid substitutions from Ancestors of MRCAs to each target species
      # nAA_targets                  amino acid set of targets
      # nAA_others                   amino acid set of the others

      # nPos_Nuc1                      start position of codon substitutions (bp)
      # nType_EvoDir_Amrca2mrca_Nuc1   evolutionary direction between Ancestors of MRCAs and MRCAs of each clade
      # nType_EvoDir_Amrca2target_Nuc1  evolutionary direction from Ancestors of MRCAs to each target species
      # nType_Nuc1Sub_target           type of amino acid substitutions of targets
      # nType_Nuc1Sub_others           type of amino acid substitutions of others
      # nType_Nuc1Sub                  type of amino acid subsitutions: identical, and different
      # nNuc1_Amrca2mrca               amino acid substitutions from Ancestors of MRCAs to MRCAs of each clade
      # nNuc1_mrca2target              amino acid substitutions from Ancestors of MRCAs to each target species
      # nNuc1_targets                  amino acid set of targets
      # nNuc1_others                   amino acid set of the others


      AmrcaCodon_list    = []
      AmrcaAA_list       = []
      AmrcaNuc1_list     = []
      AmrcaNuc2_list     = []
      AmrcaNuc3_list     = []
      
      mrcaCodon_list     = []
      mrcaAA_list        = []
      mrcaNuc1_list      = []
      mrcaNuc2_list      = []
      mrcaNuc3_list      = []

      Amrca2mrca_ID_list    = []
      Amrca2mrca_Codon_list = []
      Amrca2mrca_AA_list    = []
      Amrca2mrca_Nuc1_list    = []
      Amrca2mrca_Nuc2_list    = []
      Amrca2mrca_Nuc3_list    = []
      
      mrca2target_ID_list     = []
      mrca2target_Codon_list  = []
      mrca2target_AA_list     = []
      mrca2target_Nuc1_list     = []
      mrca2target_Nuc2_list     = []
      mrca2target_Nuc3_list     = []

      iCNT_clade_with_substitutions_codon = 0
      iCNT_clade_with_substitutions_AA    = 0
      iCNT_clade_with_substitutions_nuc1  = 0
      iCNT_clade_with_substitutions_nuc2  = 0
      iCNT_clade_with_substitutions_nuc3  = 0
      
      # generate dictionary of targetID for mrca_2_targetID_list
      targetID_mrca2target_ID_dic = {}
      for mID in mID_targetID_dic.keys():
        for targetID in mID_targetID_dic[mID]:
          tmp_mID_to_targetID_list = []
          if not mID == targetID:# a clade with multiple species (n>=2)
            mID, tmp_mID_to_targetID_list = scan_ancestors_from_mID_to_targetID(mID, targetID, sID_mID_dic, tmp_mID_to_targetID_list)
          tmp_mID_to_targetID_list.reverse()
          tmp_mID_to_targetID_list.append(targetID)
          targetID_mrca2target_ID_dic.setdefault(targetID,tmp_mID_to_targetID_list)
      
      for mrcaID in mID_targetID_dic.keys():
        mrcaCodon = sID_nSeq_dic[mrcaID][i:i+3]
        mrcaAA    = codon2aa(mrcaCodon)
        mrcaNuc1  = sID_nSeq_dic[mrcaID][i]
        mrcaNuc2  = sID_nSeq_dic[mrcaID][i+1]
        mrcaNuc3  = sID_nSeq_dic[mrcaID][i+2]
        mrcaCodon_list.append(mrcaCodon)
        mrcaAA_list.append(mrcaAA)
        mrcaNuc1_list.append(mrcaNuc1)
        mrcaNuc2_list.append(mrcaNuc2)
        mrcaNuc3_list.append(mrcaNuc3)
        
        AmrcaID    = sID_mID_dic[mrcaID]
        AmrcaCodon = sID_nSeq_dic[AmrcaID][i:i+3]
        AmrcaAA    = codon2aa(AmrcaCodon)
        AmrcaNuc1    = sID_nSeq_dic[AmrcaID][i]
        AmrcaNuc2    = sID_nSeq_dic[AmrcaID][i+1]
        AmrcaNuc3    = sID_nSeq_dic[AmrcaID][i+2]
        AmrcaCodon_list.append(AmrcaCodon)
        AmrcaAA_list.append(AmrcaAA)
        AmrcaNuc1_list.append(AmrcaNuc1)
        AmrcaNuc2_list.append(AmrcaNuc2)
        AmrcaNuc3_list.append(AmrcaNuc3)
        
        Amrca2mrca_ID_list.append(    AmrcaID+">"+mrcaID      )
        Amrca2mrca_Codon_list.append( AmrcaCodon+'>'+mrcaCodon)
        Amrca2mrca_AA_list.append(    AmrcaAA+'>'+mrcaAA      )
        Amrca2mrca_Nuc1_list.append(    AmrcaNuc1+'>'+mrcaNuc1      )
        Amrca2mrca_Nuc2_list.append(    AmrcaNuc2+'>'+mrcaNuc2      )
        Amrca2mrca_Nuc3_list.append(    AmrcaNuc3+'>'+mrcaNuc3      )
        

        tmp_targets_id_list     = mID_targetID_dic[mrcaID]
        for targetID in tmp_targets_id_list:
          tmp_mrca_to_target_ID_list = targetID_mrca2target_ID_dic[targetID]
          mrca2target_ID_list.append(   ">".join(tmp_mrca_to_target_ID_list))

          tmp_mrca_to_target_Codon_list = []
          tmp_mrca_to_target_AA_list = []
          tmp_mrca_to_target_Nuc1_list = []
          tmp_mrca_to_target_Nuc2_list = []
          tmp_mrca_to_target_Nuc3_list = []
          for nodeID in tmp_mrca_to_target_ID_list:
            nodeCodon = sID_nSeq_dic[nodeID][i:i+3]
            nodeAA = codon2aa(nodeCodon)
            nodeNuc1 = sID_nSeq_dic[nodeID][i]
            nodeNuc2 = sID_nSeq_dic[nodeID][i+1]
            nodeNuc3 = sID_nSeq_dic[nodeID][i+2]
            
            tmp_mrca_to_target_Codon_list.append(nodeCodon)
            tmp_mrca_to_target_AA_list.append(nodeAA)
            tmp_mrca_to_target_Nuc1_list.append(nodeNuc1)
            tmp_mrca_to_target_Nuc2_list.append(nodeNuc2)
            tmp_mrca_to_target_Nuc3_list.append(nodeNuc3)
            
          mrca2target_Codon_list.append(">".join(tmp_mrca_to_target_Codon_list))
          mrca2target_AA_list.append(   ">".join(tmp_mrca_to_target_AA_list))
          mrca2target_Nuc1_list.append(">".join(tmp_mrca_to_target_Nuc1_list))
          mrca2target_Nuc2_list.append(">".join(tmp_mrca_to_target_Nuc2_list))
          mrca2target_Nuc3_list.append(">".join(tmp_mrca_to_target_Nuc3_list))
          

          # Check synapormorphy
          tmp_mrca_to_target_Codon_list = list(set(tmp_mrca_to_target_Codon_list))
          if len(tmp_mrca_to_target_Codon_list) >= 2:
            iCNT_clade_with_substitutions_codon += 1

          tmp_mrca_to_target_AA_list = list(set(tmp_mrca_to_target_AA_list))
          if len(tmp_mrca_to_target_AA_list) >= 2:
            iCNT_clade_with_substitutions_AA += 1

          tmp_mrca_to_target_Nuc1_list = list(set(tmp_mrca_to_target_Nuc1_list))
          if len(tmp_mrca_to_target_Nuc1_list) >= 2:
            iCNT_clade_with_substitutions_nuc1 += 1
            
          tmp_mrca_to_target_Nuc2_list = list(set(tmp_mrca_to_target_Nuc2_list))
          if len(tmp_mrca_to_target_Nuc2_list) >= 2:
            iCNT_clade_with_substitutions_nuc2 += 1
            
          tmp_mrca_to_target_Nuc3_list = list(set(tmp_mrca_to_target_Nuc3_list))
          if len(tmp_mrca_to_target_Nuc3_list) >= 2:
            iCNT_clade_with_substitutions_nuc3 += 1
          


            
      mrcaCodon_list    = list(set(mrcaCodon_list))
      mrcaAA_list       = list(set(mrcaAA_list))
      mrcaNuc1_list     = list(set(mrcaNuc1_list))
      mrcaNuc2_list     = list(set(mrcaNuc2_list))
      mrcaNuc3_list     = list(set(mrcaNuc3_list))

      AmrcaCodon_list   = list(set(AmrcaCodon_list))
      AmrcaAA_list      = list(set(AmrcaAA_list))
      AmrcaNuc1_list    = list(set(AmrcaNuc1_list))
      AmrcaNuc2_list    = list(set(AmrcaNuc2_list))
      AmrcaNuc3_list    = list(set(AmrcaNuc3_list))

      targetCodon_list  = tmp_target_codon_list
      targetAA_list     = tmp_target_aa_list
      targetNuc1_list     = tmp_target_nuc1_list
      targetNuc2_list     = tmp_target_nuc2_list
      targetNuc3_list     = tmp_target_nuc3_list
      
      # set local variables
      ## for gene and species information
      nGene      = os.path.basename(fNAME_seq).split(".")[0]
      nPhyRel_Amrca2mrca    = '/'.join(Amrca2mrca_ID_list)
      nPhyRel_mrca2target   = '/'.join(mrca2target_ID_list)

      ## for codon information
      nPos_codon = str(iPos_Codon)
      
      #nType_EvoDir_Amrca2mrca_Codon
      if len(mrcaCodon_list) == 1:
        mrcaCodon = mrcaCodon_list[0]
        if len(AmrcaCodon_list) ==1:
          AmrcaCodon = AmrcaCodon_list[0]
          if mrcaCodon == AmrcaCodon:
            nType_EvoDir_Amrca2mrca_Codon = "PLE" # Plesiomorphic: AAA > AAA / AAA > AAA
          else:
            nType_EvoDir_Amrca2mrca_Codon = "PAR" # Parallel:    AAA > TTT / AAA > TTT
        else:
          iFlag_ple = 0
          for AmrcaCodon in AmrcaCodon_list:
            if mrcaCodon == AmrcaCodon:
              iFlag_ple += 1
          if iFlag_ple == 0:
            nType_EvoDir_Amrca2mrca_Codon = "CON" # Convergent:  AAA > TTT / GGG > TTT
          else:
            nType_EvoDir_Amrca2mrca_Codon = "MIX(Con+Ple)" # Convergent+Ple
      else:
        iFlag_ple = 0
        for mrcaCodon in mrcaCodon_list:
          for AmrcaCodon in AmrcaCodon_list:
            if mrcaCodon == AmrcaCodon:
              iFlag_ple += 1
        if iFlag_ple == 0:
          nType_EvoDir_Amrca2mrca_Codon = "DIV" # Divergent:   AAA > TTT / AAA > GGG
        else:
          nType_EvoDir_Amrca2mrca_Codon = "MIX(Div+Ple)" # Divergent+Ple
         
      #nType_CodonSub of targets
      if len(tmp_target_codon_list) == 1:
        nType_CodonSub_target = "i" # Identical
      else:
        nType_CodonSub_target = "d" # Different
        #nType_CodonSub of others
      if len(tmp_others_codon_list) == 1:
        nType_CodonSub_others = "i" # Identical
      else:
        nType_CodonSub_others = "d" # Different

      #nType_EvoDir_Codon_from_Amrca2target
      if nType_CodonSub_target == "i":
        if nType_EvoDir_Amrca2mrca_Codon[:3]  == "PLE":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "cCON"
          else:
            nType_EvoDir_Amrca2target_Codon = "cCON"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "PAR":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "sPAR"
          else:
            nType_EvoDir_Amrca2target_Codon = "cCON"        
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "CON":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "sCON"
          else:
            nType_EvoDir_Amrca2target_Codon = "cCON"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "DIV":
          if iCNT_clade_with_substitutions_codon == 0:
            print("# Error - divergent codon sub at MRCA without substitutions b/w MRCA and targets but identical sub among targets")
            print("#", nGene, nPos_codon)
            print("#", Amrca2mrca_Codon_list)
            print("#", mrca2target_Codon_list)
            print("# End of report")
            sys.exit()
          else:
            nType_EvoDir_Amrca2target_Codon = "cCON"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "MIX":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "cCon"
          else:
            nType_EvoDir_Amrca2target_Codon = "cCON"
      if nType_CodonSub_target == "d":
        if nType_EvoDir_Amrca2mrca_Codon[:3]  == "PLE":
          if iCNT_clade_with_substitutions_codon == 0:
            print("# Error - Plesiomorphic codon sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
            print("#", nGene, nPos_codon)
            print("#", Amrca2mrca_Codon_list)
            print("#", mrca2target_Codon_list)
            print("# End of report")
            sys.exit()
          else:
            nType_EvoDir_Amrca2target_Codon = "cDIV"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "PAR":
          if iCNT_clade_with_substitutions_codon == 0:
            print("# Error - Parallel codon sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
            print("#", nGene, nPos_codon)
            print("#", Amrca2mrca_Codon_list)
            print("#", mrca2target_Codon_list)
            print("# End of report")
            sys.exit()
          else:
            nType_EvoDir_Amrca2target_Codon = "cDIV"        
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "CON":
          if iCNT_clade_with_substitutions_codon == 0:
            print("# Error - Convergent codon sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
            print("#", nGene, nPos_codon)
            print("#", Amrca2mrca_Codon_list)
            print("#", mrca2target_Codon_list)
            print("# End of report")
            sys.exit()
          else:
            nType_EvoDir_Amrca2target_Codon = "cDIV"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "DIV":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "sDIV"
          else:
            nType_EvoDir_Amrca2target_Codon = "cDIV"
        elif nType_EvoDir_Amrca2mrca_Codon[:3]  == "MIX":
          if iCNT_clade_with_substitutions_codon == 0:
            nType_EvoDir_Amrca2target_Codon = "cDIV"
          else:
            nType_EvoDir_Amrca2target_Codon = "cDIV"

      #nType_CodonSub
      if nType_CodonSub_target == "i":
        if nType_CodonSub_others == "i":
          nType_CodonSub = "1"
        else:
          nType_CodonSub = "2"
      else:
        if nType_CodonSub_others == "i":
          nType_CodonSub = "3"
        else:
          nType_CodonSub = "4"
      
      nCodon_targets      = ','.join(tmp_target_codon_list)
      nCodon_others       = ','.join(tmp_others_codon_list)
      nCodon_Amrca2mrca   = '/'.join(Amrca2mrca_Codon_list)
      nCodon_mrca2target  = '/'.join(mrca2target_Codon_list)
      
      ## for amino acid information
      if iFlag_AA_Sub == 0: # if Target-specific amino acid substitutions
        nPos_AA            = str(iPos_AA)
      
        #nType_EvoDir_Amrca2mrca_AA
        if len(mrcaAA_list) == 1:
          mrcaAA = mrcaAA_list[0]
          if len(AmrcaAA_list) ==1:
            AmrcaAA = AmrcaAA_list[0]
            if mrcaAA == AmrcaAA:
              nType_EvoDir_Amrca2mrca_AA = "PLE" # Plesiomorphic: AAA > AAA / AAA > AAA
            else:
              nType_EvoDir_Amrca2mrca_AA = "PAR" # Parallel:    AAA > TTT / AAA > TTT
          else:
            iFlag_ple = 0
            for AmrcaAA in AmrcaAA_list:
              if mrcaAA == AmrcaAA:
                iFlag_ple += 1
            if iFlag_ple == 0:
              nType_EvoDir_Amrca2mrca_AA = "CON" # Convergent:  AAA > TTT / GGG > TTT
            else:
              nType_EvoDir_Amrca2mrca_AA = "MIX(Con+Ple)" # Convergent+Ple
        else:
          iFlag_ple = 0
          for mrcaAA in mrcaAA_list:
            for AmrcaAA in AmrcaAA_list:
              if mrcaAA == AmrcaAA:
                iFlag_ple += 1
          if iFlag_ple == 0:
            nType_EvoDir_Amrca2mrca_AA = "DIV" # Divergent:   AAA > TTT / AAA > GGG
          else:
            nType_EvoDir_Amrca2mrca_AA = "MIX(Div+Ple)" # Divergent+Ple
            
        #nType_AASub of targets
        if len(tmp_target_aa_list) == 1:
          nType_AASub_target = "i" # Identical
        else:
          nType_AASub_target = "d" # Different
        #nType_AASub of others
        if len(tmp_others_aa_list) == 1:
          nType_AASub_others = "i" # Identical
        else:
          nType_AASub_others = "d" # Different

        #nType_EvoDir_AA_from_Amrca2target
        if nType_AASub_target == "i":
          if nType_EvoDir_Amrca2mrca_AA[:3]  == "PLE":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "cCON"
            else:
              nType_EvoDir_Amrca2target_AA = "cCON"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "PAR":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "sPAR"
            else:
              nType_EvoDir_Amrca2target_AA = "cCON"        
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "CON":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "sCON"
            else:
              nType_EvoDir_Amrca2target_AA = "cCON"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "DIV":
            if iCNT_clade_with_substitutions_AA == 0:
              print("# Error - divergent AA sub at MRCA without substitutions b/w MRCA and targets but identical sub among targets")
              print("#", nGene, nPos_AA)
              print("#", Amrca2mrca_AA_list)
              print("#", mrca2target_AA_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_AA = "cCON"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "MIX":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "cCON"
            else:
              nType_EvoDir_Amrca2target_AA = "cCON"
        if nType_AASub_target == "d":
          if nType_EvoDir_Amrca2mrca_AA[:3]  == "PLE":
            if iCNT_clade_with_substitutions_AA == 0:
              print("# Error - Plesiomorphic AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_AA)
              print("#", Amrca2mrca_AA_list)
              print("#", mrca2target_AA_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_AA = "cDIV"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "PAR":
            if iCNT_clade_with_substitutions_AA == 0:
              print("# Error - Parallel AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_AA)
              print("#", Amrca2mrca_AA_list)
              print("#", mrca2target_AA_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_AA = "cDIV"        
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "CON":
            if iCNT_clade_with_substitutions_AA == 0:
              print("# Error - Convergent AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_AA)
              print("#", Amrca2mrca_AA_list)
              print("#", mrca2target_AA_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_AA = "cDIV"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "DIV":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "sDIV"
            else:
              nType_EvoDir_Amrca2target_AA = "cDIV"
          elif nType_EvoDir_Amrca2mrca_AA[:3]  == "MIX":
            if iCNT_clade_with_substitutions_AA == 0:
              nType_EvoDir_Amrca2target_AA = "cDIV"
            else:
              nType_EvoDir_Amrca2target_AA = "cDIV"

        #nType_AASub
        if nType_AASub_target == "i":
          if nType_AASub_others == "i":
            nType_AASub = "1"
          else:
            nType_AASub = "2"
        else:
          if nType_AASub_others == "i":
            nType_AASub = "3"
          else:
            nType_AASub = "4"
        
        nAA_targets      = ','.join(tmp_target_aa_list)
        nAA_others       = ','.join(tmp_others_aa_list)
        nAA_Amrca2mrca   = '/'.join(Amrca2mrca_AA_list)
        nAA_mrca2target  = '/'.join(mrca2target_AA_list)

      else: # if Target-specific amino acid substitutions
        nPos_AA            = ''
        nType_EvoDir_Amrca2mrca_AA   = ''
        nType_EvoDir_Amrca2target_AA = ''
        nType_AASub_target = ''
        nType_AASub_others = ''
        nType_AASub        = ''
        nAA_Amrca2mrca     = ''
        nAA_mrca2target    = ''
        nAA_targets        = ''
        nAA_others         = ''





      ## for nucleotide 1st pos information
      if iFlag_Nuc1_Sub == 0: # if Target-specific nucleotide substitutions
        nPos_Nuc1            = str(iPos_Nuc1)
      
        #nType_EvoDir_Amrca2mrca_AA
        if len(mrcaNuc1_list) == 1:
          mrcaNuc1 = mrcaNuc1_list[0]
          if len(AmrcaNuc1_list) ==1:
            AmrcaNuc1 = AmrcaNuc1_list[0]
            if mrcaNuc1 == AmrcaNuc1:
              nType_EvoDir_Amrca2mrca_Nuc1 = "PLE" # Plesiomorphic: AAA > AAA / AAA > AAA
            else:
              nType_EvoDir_Amrca2mrca_Nuc1 = "PAR" # Parallel:    AAA > TTT / AAA > TTT
          else:
            iFlag_ple = 0
            for AmrcaNuc1 in AmrcaNuc1_list:
              if mrcaNuc1 == AmrcaNuc1:
                iFlag_ple += 1
            if iFlag_ple == 0:
              nType_EvoDir_Amrca2mrca_Nuc1 = "CON" # Convergent:  AAA > TTT / GGG > TTT
            else:
              nType_EvoDir_Amrca2mrca_Nuc1 = "MIX(Con+Ple)" # Convergent+Ple
        else:
          iFlag_ple = 0
          for mrcaNuc1 in mrcaNuc1_list:
            for AmrcaNuc1 in AmrcaNuc1_list:
              if mrcaNuc1 == AmrcaNuc1:
                iFlag_ple += 1
          if iFlag_ple == 0:
            nType_EvoDir_Amrca2mrca_Nuc1 = "DIV" # Divergent:   AAA > TTT / AAA > GGG
          else:
            nType_EvoDir_Amrca2mrca_Nuc1 = "MIX(Div+Ple)" # Divergent+Ple
            
        #nType_Nuc1Sub of targets
        if len(tmp_target_nuc1_list) == 1:
          nType_Nuc1Sub_target = "i" # Identical
        else:
          nType_Nuc1Sub_target = "d" # Different
        #nType_Nuc1Sub of others
        if len(tmp_others_nuc1_list) == 1:
          nType_Nuc1Sub_others = "i" # Identical
        else:
          nType_Nuc1Sub_others = "d" # Different

        #nType_EvoDir_AA_from_Amrca2target
        if nType_Nuc1Sub_target == "i":
          if nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "sPAR"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"        
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "sCON"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              print("# Error - divergent nucleotide sub at MRCA without substitutions b/w MRCA and targets but identical sub among targets")
              print("#", nGene, nPos_Nuc1)
              print("#", Amrca2mrca_Nuc1_list)
              print("#", mrca2target_Nuc1_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cCON"
        if nType_Nuc1Sub_target == "d":
          if nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              print("# Error - Plesiomorphic nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc1)
              print("#", Amrca2mrca_Nuc1_list)
              print("#", mrca2target_Nuc1_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              print("# Error - Parallel nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc1)
              print("#", Amrca2mrca_Nuc1_list)
              print("#", mrca2target_Nuc1_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"        
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              print("# Error - Convergent AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc1)
              print("#", Amrca2mrca_Nuc1_list)
              print("#", mrca2target_Nuc1_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "sDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc1[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc1 == 0:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc1 = "cDIV"

        #nType_AASub
        if nType_Nuc1Sub_target == "i":
          if nType_Nuc1Sub_others == "i":
            nType_Nuc1Sub = "1"
          else:
            nType_Nuc1Sub = "2"
        else:
          if nType_Nuc1Sub_others == "i":
            nType_Nuc1Sub = "3"
          else:
            nType_Nuc1Sub = "4"
        
        nNuc1_targets      = ','.join(tmp_target_nuc1_list)
        nNuc1_others       = ','.join(tmp_others_nuc1_list)
        nNuc1_Amrca2mrca   = '/'.join(Amrca2mrca_Nuc1_list)
        nNuc1_mrca2target  = '/'.join(mrca2target_Nuc1_list)

      else: # if Target-specific nucleotide substitutions
        nPos_Nuc1            = ''
        nType_EvoDir_Amrca2mrca_Nuc1   = ''
        nType_EvoDir_Amrca2target_Nuc1 = ''
        nType_Nuc1Sub_target = ''
        nType_Nuc1Sub_others = ''
        nType_Nuc1Sub        = ''
        nNuc1_Amrca2mrca     = ''
        nNuc1_mrca2target    = ''
        nNuc1_targets        = ''
        nNuc1_others         = ''




      ## for nucleotide 2nd pos information
      if iFlag_Nuc2_Sub == 0: # if Target-specific nucleotide substitutions
        nPos_Nuc2            = str(iPos_Nuc2)
      
        #nType_EvoDir_Amrca2mrca_AA
        if len(mrcaNuc2_list) == 1:
          mrcaNuc2 = mrcaNuc2_list[0]
          if len(AmrcaNuc2_list) ==1:
            AmrcaNuc2 = AmrcaNuc2_list[0]
            if mrcaNuc2 == AmrcaNuc2:
              nType_EvoDir_Amrca2mrca_Nuc2 = "PLE" # Plesiomorphic: AAA > AAA / AAA > AAA
            else:
              nType_EvoDir_Amrca2mrca_Nuc2 = "PAR" # Parallel:    AAA > TTT / AAA > TTT
          else:
            iFlag_ple = 0
            for AmrcaNuc2 in AmrcaNuc2_list:
              if mrcaNuc2 == AmrcaNuc2:
                iFlag_ple += 1
            if iFlag_ple == 0:
              nType_EvoDir_Amrca2mrca_Nuc2 = "CON" # Convergent:  AAA > TTT / GGG > TTT
            else:
              nType_EvoDir_Amrca2mrca_Nuc2 = "MIX(Con+Ple)" # Convergent+Ple
        else:
          iFlag_ple = 0
          for mrcaNuc2 in mrcaNuc2_list:
            for AmrcaNuc2 in AmrcaNuc2_list:
              if mrcaNuc2 == AmrcaNuc2:
                iFlag_ple += 1
          if iFlag_ple == 0:
            nType_EvoDir_Amrca2mrca_Nuc2 = "DIV" # Divergent:   AAA > TTT / AAA > GGG
          else:
            nType_EvoDir_Amrca2mrca_Nuc2 = "MIX(Div+Ple)" # Divergent+Ple
            
        #nType_Nuc1Sub of targets
        if len(tmp_target_nuc2_list) == 1:
          nType_Nuc2Sub_target = "i" # Identical
        else:
          nType_Nuc2Sub_target = "d" # Different
        #nType_Nuc1Sub of others
        if len(tmp_others_nuc2_list) == 1:
          nType_Nuc2Sub_others = "i" # Identical
        else:
          nType_Nuc2Sub_others = "d" # Different

        #nType_EvoDir_AA_from_Amrca2target
        if nType_Nuc2Sub_target == "i":
          if nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "sPAR"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"        
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "sCON"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              print("# Error - divergent nucleotide sub at MRCA without substitutions b/w MRCA and targets but identical sub among targets")
              print("#", nGene, nPos_Nuc2)
              print("#", Amrca2mrca_Nuc2_list)
              print("#", mrca2target_Nuc2_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cCON"
        if nType_Nuc2Sub_target == "d":
          if nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              print("# Error - Plesiomorphic nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc2)
              print("#", Amrca2mrca_Nuc2_list)
              print("#", mrca2target_Nuc2_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              print("# Error - Parallel nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc2)
              print("#", Amrca2mrca_Nuc2_list)
              print("#", mrca2target_Nuc2_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"        
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              print("# Error - Convergent AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc2)
              print("#", Amrca2mrca_Nuc2_list)
              print("#", mrca2target_Nuc2_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "sDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc2[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc2 == 0:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc2 = "cDIV"

        #nType_AASub
        if nType_Nuc2Sub_target == "i":
          if nType_Nuc2Sub_others == "i":
            nType_Nuc2Sub = "1"
          else:
            nType_Nuc2Sub = "2"
        else:
          if nType_Nuc2Sub_others == "i":
            nType_Nuc2Sub = "3"
          else:
            nType_Nuc2Sub = "4"
        
        nNuc2_targets      = ','.join(tmp_target_nuc2_list)
        nNuc2_others       = ','.join(tmp_others_nuc2_list)
        nNuc2_Amrca2mrca   = '/'.join(Amrca2mrca_Nuc2_list)
        nNuc2_mrca2target  = '/'.join(mrca2target_Nuc2_list)

      else: # if Target-specific nucleotide substitutions
        nPos_Nuc2            = ''
        nType_EvoDir_Amrca2mrca_Nuc2   = ''
        nType_EvoDir_Amrca2target_Nuc2 = ''
        nType_Nuc2Sub_target = ''
        nType_Nuc2Sub_others = ''
        nType_Nuc2Sub        = ''
        nNuc2_Amrca2mrca     = ''
        nNuc2_mrca2target    = ''
        nNuc2_targets        = ''
        nNuc2_others         = ''



      ## for nucleotide 3rd pos information
      if iFlag_Nuc3_Sub == 0: # if Target-specific nucleotide substitutions
        nPos_Nuc3            = str(iPos_Nuc3)
      
        #nType_EvoDir_Amrca2mrca_AA
        if len(mrcaNuc3_list) == 1:
          mrcaNuc3 = mrcaNuc3_list[0]
          if len(AmrcaNuc3_list) ==1:
            AmrcaNuc3 = AmrcaNuc3_list[0]
            if mrcaNuc3 == AmrcaNuc3:
              nType_EvoDir_Amrca2mrca_Nuc3 = "PLE" # Plesiomorphic: AAA > AAA / AAA > AAA
            else:
              nType_EvoDir_Amrca2mrca_Nuc3 = "PAR" # Parallel:    AAA > TTT / AAA > TTT
          else:
            iFlag_ple = 0
            for AmrcaNuc3 in AmrcaNuc3_list:
              if mrcaNuc3 == AmrcaNuc3:
                iFlag_ple += 1
            if iFlag_ple == 0:
              nType_EvoDir_Amrca2mrca_Nuc3 = "CON" # Convergent:  AAA > TTT / GGG > TTT
            else:
              nType_EvoDir_Amrca2mrca_Nuc3 = "MIX(Con+Ple)" # Convergent+Ple
        else:
          iFlag_ple = 0
          for mrcaNuc3 in mrcaNuc3_list:
            for AmrcaNuc3 in AmrcaNuc3_list:
              if mrcaNuc3 == AmrcaNuc3:
                iFlag_ple += 1
          if iFlag_ple == 0:
            nType_EvoDir_Amrca2mrca_Nuc3 = "DIV" # Divergent:   AAA > TTT / AAA > GGG
          else:
            nType_EvoDir_Amrca2mrca_Nuc3 = "MIX(Div+Ple)" # Divergent+Ple
            
        #nType_Nuc1Sub of targets
        if len(tmp_target_nuc3_list) == 1:
          nType_Nuc3Sub_target = "i" # Identical
        else:
          nType_Nuc3Sub_target = "d" # Different
        #nType_Nuc1Sub of others
        if len(tmp_others_nuc3_list) == 1:
          nType_Nuc3Sub_others = "i" # Identical
        else:
          nType_Nuc3Sub_others = "d" # Different

        #nType_EvoDir_AA_from_Amrca2target
        if nType_Nuc3Sub_target == "i":
          if nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "sPAR"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"        
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "sCON"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              print("# Error - divergent nucleotide sub at MRCA without substitutions b/w MRCA and targets but identical sub among targets")
              print("#", nGene, nPos_Nuc3)
              print("#", Amrca2mrca_Nuc3_list)
              print("#", mrca2target_Nuc3_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cCON"
        if nType_Nuc3Sub_target == "d":
          if nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "PLE":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              print("# Error - Plesiomorphic nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc3)
              print("#", Amrca2mrca_Nuc3_list)
              print("#", mrca2target_Nuc3_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "PAR":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              print("# Error - Parallel nucleotide sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc3)
              print("#", Amrca2mrca_Nuc3_list)
              print("#", mrca2target_Nuc3_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"        
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "CON":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              print("# Error - Convergent AA sub at MRCA without substitutions b/w MRCA and targets but different sub among targets")
              print("#", nGene, nPos_Nuc3)
              print("#", Amrca2mrca_Nuc3_list)
              print("#", mrca2target_Nuc3_list)
              print("# End of report")
              sys.exit()
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "DIV":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "sDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"
          elif nType_EvoDir_Amrca2mrca_Nuc3[:3]  == "MIX":
            if iCNT_clade_with_substitutions_nuc3 == 0:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"
            else:
              nType_EvoDir_Amrca2target_Nuc3 = "cDIV"

        #nType_AASub
        if nType_Nuc3Sub_target == "i":
          if nType_Nuc3Sub_others == "i":
            nType_Nuc3Sub = "1"
          else:
            nType_Nuc3Sub = "2"
        else:
          if nType_Nuc3Sub_others == "i":
            nType_Nuc3Sub = "3"
          else:
            nType_Nuc3Sub = "4"
        
        nNuc3_targets      = ','.join(tmp_target_nuc3_list)
        nNuc3_others       = ','.join(tmp_others_nuc3_list)
        nNuc3_Amrca2mrca   = '/'.join(Amrca2mrca_Nuc3_list)
        nNuc3_mrca2target  = '/'.join(mrca2target_Nuc3_list)

      else: # if Target-specific nucleotide substitutions
        nPos_Nuc3            = ''
        nType_EvoDir_Amrca2mrca_Nuc3   = ''
        nType_EvoDir_Amrca2target_Nuc3 = ''
        nType_Nuc3Sub_target = ''
        nType_Nuc3Sub_others = ''
        nType_Nuc3Sub        = ''
        nNuc3_Amrca2mrca     = ''
        nNuc3_mrca2target    = ''
        nNuc3_targets        = ''
        nNuc3_others         = ''




      
      # site-wise codon and amino acid sequence information of each species
      nCodon_line = ''
      nAA_line    = ''
      for sID in templete_sID_list:
        if sID in sID_nSeq_dic.keys():
          nCodon = sID_nSeq_dic[sID][i:i+3]
          nAA = codon2aa(nCodon)
        else:
          nCodon = ''
          nAA    = ''
        nCodon_line += '\t' + nCodon
        nAA_line    += '\t' + nAA

      # write results
      tmpline = nGene +'\t'+ nPhyRel_Amrca2mrca +'\t'+ nPhyRel_mrca2target
      tmpline += '\t'+ nPos_codon +'\t'+ nType_EvoDir_Amrca2mrca_Codon +'\t'+ nType_EvoDir_Amrca2target_Codon +'\t'+ nType_CodonSub_target +'\t'+ nType_CodonSub_others +'\t'+ nType_CodonSub +'\t'+ nCodon_Amrca2mrca +'\t'+ nCodon_mrca2target +'\t'+ nCodon_targets +'\t'+ nCodon_others
      tmpline += '\t'+ nPos_AA    +'\t'+ nType_EvoDir_Amrca2mrca_AA    +'\t'+ nType_EvoDir_Amrca2target_AA    +'\t'+ nType_AASub_target    +'\t'+ nType_AASub_others    +'\t'+ nType_AASub    +'\t'+ nAA_Amrca2mrca    +'\t'+ nAA_mrca2target    +'\t'+ nAA_targets    +'\t'+ nAA_others
      tmpline += '\t'+ nPos_Nuc1  +'\t'+ nType_EvoDir_Amrca2mrca_Nuc1  +'\t'+ nType_EvoDir_Amrca2target_Nuc1  +'\t'+ nType_Nuc1Sub_target  +'\t'+ nType_Nuc1Sub_others  +'\t'+ nType_Nuc1Sub  +'\t'+ nNuc1_Amrca2mrca  +'\t'+ nNuc1_mrca2target  +'\t'+ nNuc1_targets  +'\t'+ nNuc1_others
      tmpline += '\t'+ nPos_Nuc2  +'\t'+ nType_EvoDir_Amrca2mrca_Nuc2  +'\t'+ nType_EvoDir_Amrca2target_Nuc2  +'\t'+ nType_Nuc2Sub_target  +'\t'+ nType_Nuc2Sub_others  +'\t'+ nType_Nuc2Sub  +'\t'+ nNuc2_Amrca2mrca  +'\t'+ nNuc2_mrca2target  +'\t'+ nNuc2_targets  +'\t'+ nNuc2_others
      tmpline += '\t'+ nPos_Nuc3  +'\t'+ nType_EvoDir_Amrca2mrca_Nuc3  +'\t'+ nType_EvoDir_Amrca2target_Nuc3  +'\t'+ nType_Nuc3Sub_target  +'\t'+ nType_Nuc3Sub_others  +'\t'+ nType_Nuc3Sub  +'\t'+ nNuc3_Amrca2mrca  +'\t'+ nNuc3_mrca2target  +'\t'+ nNuc3_targets  +'\t'+ nNuc3_others
      tmpline += nCodon_line + nAA_line +'\n'
      fpout.write(tmpline)
    else: # if not Target-specific codon substitutions
      pass
  return(fpout)
      
      
      

def main(option_list, nGene, fpout):
  # local vaiables for each gene
  nGene                     = nGene
  nTargets                  = option_list[0]
  templete_mID_targetID_dic = option_list[1]
  iPATH                     = option_list[2]
  fNAME_seq_format          = option_list[3]
  fNAME_tree_format         = option_list[4]
  templete_sID_list         = option_list[6]
  outgroup_list             = option_list[7]
  
  fNAME_seq   = iPATH + nGene + fNAME_seq_format
  fNAME_tree  = iPATH + nGene + fNAME_tree_format

  #targetID_list       = nTargets.split(',') 
  nGene_nTree         = read_tree_file(fNAME_tree)
  nGene_nTree         = read_tre_removing_branchlength(nGene_nTree)
  nGene_terNode_list  = terminal_node_from_tree(nGene_nTree)
  nGene_target_list   = parse_target_from_tree(nTargets,nGene_nTree)

  # check_clade
  iCNT_clade = 0
  for mID in templete_mID_targetID_dic.keys():
    iFlag_clade = 0
    for targetID in templete_mID_targetID_dic[mID]:
      if targetID in nGene_target_list:
        iFlag_clade += 1
    if iFlag_clade > 0:
      iCNT_clade += 1
  if not iCNT_clade == len(templete_mID_targetID_dic.keys()):
    return(fpout)
  else:
    # Phylogenetic relationships
    nGene_mID_targetID_dic, nGene_mID_mmID_dic = Anc_Finder(nGene_nTree, ','.join(nGene_target_list))
    nGene_mID_sID_dic = {}
    nGene_mID_sID_dic = make_mID_sID_dic(nGene_nTree, nGene_mID_sID_dic)
    nGene_sID_mID_dic = make_sID_mID_dic(nGene_mID_sID_dic)
    # Convergent variant finder
    fpout = ConVarFinder(fNAME_seq, nGene_terNode_list, nGene_mID_targetID_dic, nGene_mID_mmID_dic, nGene_sID_mID_dic, nGene_target_list, templete_sID_list, fpout, outgroup_list)
  return(fpout)


if __name__ == "__main__":
  # Input
  option_list       = set_options()
  nTargets          = option_list[0]
  oPATH             = option_list[5]
  iPATH             = option_list[2]
  fNAME_seq_format  = option_list[3]
  templete_sID_list = option_list[6]
  outgroup_list     = option_list[7]
  oNAME             = oPATH + "ConVarFinder_"+nTargets + '.txt'


  # analysis for each gene
  flist_seq = glob.glob(iPATH+"*"+fNAME_seq_format)
  fpout = open(oNAME,'w')
  # nGene                gene name
  # nPhyRel_Amrca2mrca   phylogenetic relationship from Ancestors of MRCAs to MRCAs of each clade
  # nPhyRel_mrca2target  phylogenetic relationship from Ancestors of MRCAs to each target species

  # nPos_codon                      start position of codon substitutions (bp)
  # nType_EvoDir_Amrca2mrca_Codon   evolutionary direction between Ancestors of MRCAs and MRCAs of each clade
  # nType_EvoDir_Amrca2target_Codon  evolutionary direction from Ancestors of MRCAs to each target species
  # nType_CodonSub_target           type of codon substitutions of targets
  # nType_CodonSub_others           type of codon substitutions of others
  # nType_CodonSub                  type of codon substitutions: identical, and different
  # nCodon_Amrca2mrca               codon substitutions from Ancestors of MRCAs to MRCAs of each clade
  # nCodon_mrca2target              codon substitutions from Ancestors of MRCAs to each target species
  # nCodon_targets                  codon set of targets
  # nCodon_others                   codon set of the others
  
  # nPos_AA                      start position of codon substitutions (bp)
  # nType_EvoDir_Amrca2mrca_AA   evolutionary direction between Ancestors of MRCAs and MRCAs of each clade
  # nType_EvoDir_Amrca2target_AA  evolutionary direction from Ancestors of MRCAs to each target species
  # nType_AASub_target           type of amino acid substitutions of targets
  # nType_AASub_others           type of amino acid substitutions of others
  # nType_AASub                  type of amino acid subsitutions: identical, and different
  # nAA_Amrca2mrca               amino acid substitutions from Ancestors of MRCAs to MRCAs of each clade
  # nAA_mrca2target              amino acid substitutions from Ancestors of MRCAs to each target species
  # nAA_targets                  amino acid set of targets
  # nAA_others                   amino acid set of the others
  # "\t".join(templete_sID_list)\n"   list of codon
  # '\t".join(templete_sID_list)+'\n' list of aa
  
  nHeader = "Gene" +'\t'+ "Phy.Rel._Amrca2mrca" +'\t'+ "Phy.Rel._mrca2target"
  nHeader += '\t'+ "Pos_codon"+'\t'+ "EvoDir_Amrca2mrca_Codon" +'\t'+ "EvoDir_Amrca2target_Codon" +'\t'+ 'Type_CodonSub_target' +'\t'+ 'Type_CodonSub_others' +'\t'+ 'Type_CodonSub' +'\t'+ "Codon_Amrca2mrca" +'\t'+ "Codon_mrca2target" +'\t'+ "Codon_targets" +'\t'+ "Codon_others"
  nHeader += '\t'+ "Pos_AA"   +'\t'+ "EvoDir_Amrca2mrca_AA"    +'\t'+ "EvoDir_Amrca2target_AA"    +'\t'+ 'Type_AASub_target'    +'\t'+ 'Type_AASub_others'    +'\t'+ 'Type_AASub'    +'\t'+ "AA_Amrca2mrca"    +'\t'+ "AA_mrca2target"    +'\t'+ "AA_targets"    +'\t'+ "AA_others"
  nHeader += '\t'+ "Pos_Nuc1" +'\t'+ "EvoDir_Amrca2mrca_Nuc1"  +'\t'+ "EvoDir_Amrca2target_Nuc1"  +'\t'+ 'Type_Nuc1Sub_target'  +'\t'+ 'Type_Nuc1Sub_others'  +'\t'+ 'Type_Nuc1Sub'  +'\t'+ "Nuc1_Amrca2mrca"  +'\t'+ "Nuc1_mrca2target"  +'\t'+ "Nuc1_targets"  +'\t'+ "Nuc1_others"
  nHeader += '\t'+ "Pos_Nuc2" +'\t'+ "EvoDir_Amrca2mrca_Nuc2"  +'\t'+ "EvoDir_Amrca2target_Nuc2"  +'\t'+ 'Type_Nuc2Sub_target'  +'\t'+ 'Type_Nuc2Sub_others'  +'\t'+ 'Type_Nuc2Sub'  +'\t'+ "Nuc2_Amrca2mrca"  +'\t'+ "Nuc2_mrca2target"  +'\t'+ "Nuc2_targets"  +'\t'+ "Nuc2_others"
  nHeader += '\t'+ "Pos_Nuc3" +'\t'+ "EvoDir_Amrca2mrca_Nuc3"  +'\t'+ "EvoDir_Amrca2target_Nuc3"  +'\t'+ 'Type_Nuc3Sub_target'  +'\t'+ 'Type_Nuc3Sub_others'  +'\t'+ 'Type_Nuc3Sub'  +'\t'+ "Nuc3_Amrca2mrca"  +'\t'+ "Nuc3_mrca2target"  +'\t'+ "Nuc3_targets"  +'\t'+ "Nuc3_others"

  nHeader += '\t'+ "\t".join(templete_sID_list) +'\t'+ "\t".join(templete_sID_list) +'\n'
  fpout.write(nHeader)
  for fNAME_seq in flist_seq:
    nGene = os.path.basename(fNAME_seq).split(fNAME_seq_format)[0]
    fpout = main(option_list, nGene, fpout)
  fpout.close()
