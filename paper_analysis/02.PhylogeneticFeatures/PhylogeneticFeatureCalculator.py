# Calculator of phylogenetic features
# Version 1.0 (14.Nov.2020)
# written by Chul Lee (e-mail: chul.bioinfo@gmail.com)
# This code was developed and conducted in Python 3.7.1 (v3.7.1:260ec2c36a, Oct 20 2018, 14:57:15) [MSC v.1915 64 bit (AMD64)]
# OS for development and analysis: Windows 10 Education

# Libraries
import sys


# Global variable
## Input

## Input of target species list
nTargets ="TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN" # 6 vocal learning birds
## nTargets = "PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL" # 15 waterbirds
## nTargets = "FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU" # 6 birds of prey
## Explanation about target list
## - The above list of targets is default in this program.
## - Have to be seperated by comma(,).
## - If you want to change targets, you can choice one of two ways.
##
## (1) Change nTargets in the first line in this code
##     nTargets = "speciesA,speciesB,...,speciesZ"
##
## (2) Insert argument in terminal (recommended)
##     python ProOriBL.py speciesA,speciesB,...,speciesZ
##

## Input of phylogenetic tree
nTree = "(((((((((((((TAEGU:0.042,GEOFO:0.0372):0.0245,CORBR:0.0385):0.044,MANVI:0.0693):0.0479,(MELUN:0.0585,NESNO:0.042):0.0528):0.0021,FALPE:0.0762):0.0013,CARCR:0.055):0.0015,(((((((MERNU:0.103,PICPU:0.1434):0.0039,BUCRH:0.0998):0.0025,APAVI:0.1048):0.0012,LEPDI:0.067):0.0014,COLST:0.1205):0.0006,TYTAL:0.0675):0.0006,((HALLE:0.001,HALAL:0.0011):0.0427,CATAU:0.0301):0.002):0.0008):0.0041,((((((PELCR:0.0417,EGRGA:0.0562):0.0013,NIPNI:0.0398):0.001,PHACA:0.0639):0.0016,(FULGL:0.0342,(PYGAD:0.0135,APTFO:0.0109):0.0214):0.0012):0.0019,GAVST:0.0445):0.0021,(PHALE:0.0591,EURHE:0.0947):0.0032):0.0012):0.0011,((CHAVO:0.0602,BALRE:0.0571):0.0015,OPHHO:0.0756):0.0009):0.001,(((CALAN:0.1126,CHAPE:0.0908):0.0204,CAPCA:0.0697):0.0038,((CHLUN:0.0731,TAUER:0.0708):0.0014,CUCCA:0.1115):0.0011):0.0011):0.0011,(((MESUN:0.0871,PTEGU:0.0717):0.0027,COLLI:0.0963):0.0014,(PHORU:0.0312,PODCR:0.0584):0.013):0.0009):0.032,((GALGA:0.0391,MELGA:0.0449):0.1131,ANAPL:0.1003):0.0394):0.0456,(TINMA:0.1471,STRCA:0.0669):0.0456);"

## Explanation about tree
## - Tree from Fig S7 in Erich et al. (2014) Science.
## - have to be bifurcating rooted tree
## - have to contain branch lengths
## - have to remove bootstrap values
## 24 primates
## nTree = "(((((((((PANTR:0.00246757,PANPA:0.00307243):0.00430055,HOMSA:0.00660945):0.00175688,GORGO:0.00867627):0.00836254,PONAB:0.01726309):0.00282649,NOMLE:0.01961136):0.01077990,(((((CERAT:0.00580441,MANLE:0.00664559):0.00064477,PAPAN:0.00665523):0.00119510,((MACFA:0.00215180,MACMU:0.00232820):0.00094595,MACNE:0.00418405):0.00432688):0.00401008,CHLSA:0.01169715):0.00475309,((RHIRO:0.00209680,RHIBI:0.00301320):0.00939164,COLAN:0.01232336):0.00504790):0.01387906):0.01760134,(((CEBCA:0.02593806,SAIBO:0.02639194):0.00278398,CALJA:0.02999602):0.00010000,AOTNA:0.02501050):0.02356746):0.02752902,CARSY:0.07990043):0.00098575,((MICMU:0.04120270,PROCO:0.04239730):0.02776095,OTOGA:0.07574905):0.00730385);"


## Input of decimal places
iDP = 8 # Calculate down to 'iDP' decimal points of phylogenetic features




######################
# Preparing variables

## make info_tree: topology, branch lengths
info_tree = nTree.split(',')


## make sID_list
sID_list = []
part = nTree.replace("(",'').split(',')
for sIDwInfo in part:
    sID = sIDwInfo.split(':')[0]
    sID_list.append(sID)

sID_LIST = []
part = nTree.replace("(",'').split(',')
for sIDwInfo in part:
    sID = sIDwInfo.split(':')[0]
    sID_LIST.append(sID)


## make target_list
if len(sys.argv)== 2:
    target_list = sys.argv[1].strip().split(',')
elif len(sys.argv)== 1:
    target_list = nTargets.strip().split(',')
else:
    print("Input error",sys.argv)
    sys.exit()
    
## make others_list
others_list = []
for sID in sID_list:
    if not sID in target_list:
        others_list.append(sID)
        
#print(target_list)
#sys.exit()
## check if targetID in sID_list?
iCNT_target = 0
for targetID in target_list:
    if targetID in sID_list:
        iCNT_target += 1
if not len(target_list)==iCNT_target:
    print("ERROR: targets -"+','.join(target_list))
    sys.exit()

## Sort target_list by the order of sID_list
tmp_target_list = []
for sID in sID_list:
    if sID in target_list:
        tmp_target_list.append(sID)
target_list = tmp_target_list

## Output
oNAME = "PhylogeneticFeatures_output_"+','.join(target_list)+'.txt'
fpout = open(oNAME,'w')
fpout.write("# Targets: "+','.join(target_list)+'\n')


#############################
# Local function

def ordering_sID(sTarget_list):
    tmp_sTarget_list = []
    for sID in sID_list:
        if sID in sTarget_list:
            tmp_sTarget_list.append(sID)
    sTarget_list = tmp_sTarget_list
    return sTarget_list


####################################
# Print variables
print("### Input")
#print("# Tree: "+ nTree)
print("# Targets: "+','.join(target_list))
print("#")


####################################
# Function
def find_originbranchlength(monoTarget_list, info_tree, sID_list):
    iCNT_R_parenthesis = 0
    iCNT_L_parenthesis = 0
    for sID in monoTarget_list:
        i = sID_list.index(sID)
        iCNT_R_parenthesis += info_tree[i].count(")") 
        iCNT_L_parenthesis += info_tree[i].count("(")
    diff_RminusL = iCNT_R_parenthesis - iCNT_L_parenthesis
    if iCNT_R_parenthesis < iCNT_L_parenthesis:
        iLen_originbranch = float(info_tree[sID_list.index(monoTarget_list[-1])].split(')')[-1].split(':')[-1])
    elif iCNT_R_parenthesis > iCNT_L_parenthesis:
        iLen_originbranch = float(info_tree[i].split(")")[-diff_RminusL-1].split(":")[-1])
    else:#iCNT_R_parenthesis == iCNT_L_parenthesis
        iLen_originbranch = 0
    return iLen_originbranch


def calculate_sum_branchlength(monoTarget_list, info_tree, sID_list):
    if len(monoTarget_list)==1:
        #print(info_tree[sID_list.index(monoTarget_list[0])])
        iSum_branch = float(info_tree[sID_list.index(monoTarget_list[0])].split(')')[0].split(":")[-1])

    else:
        iCNT_R_parenthesis = 0
        iCNT_L_parenthesis = 0
        for sID in monoTarget_list:
            i = sID_list.index(sID)
            iCNT_R_parenthesis += info_tree[i].count(")") 
            iCNT_L_parenthesis += info_tree[i].count("(")
        diff_RminusL = iCNT_R_parenthesis - iCNT_L_parenthesis

        iSum_branch = 0
        if iCNT_R_parenthesis < iCNT_L_parenthesis:
            for sID in monoTarget_list:
                BL_info = info_tree[sID_list.index(sID)]
                BL_info = BL_info.replace(");","")
                #print(BL_info)
                if ")" in BL_info:
                    BL_info_list = BL_info.split(")")
                    for sID_nBL in BL_info_list:
                        iLen_branch = float(sID_nBL.split(":")[-1])
                        iSum_branch += iLen_branch
                else:
                    sID_nBL = BL_info
                    iLen_branch = float(sID_nBL.split(":")[-1])
                    iSum_branch += iLen_branch
            
        elif iCNT_R_parenthesis > iCNT_L_parenthesis:
            for sID in monoTarget_list[:-1]:
                BL_info = info_tree[sID_list.index(sID)]
                BL_info = BL_info.replace(");","")
                #print(BL_info)
                if ")" in BL_info:
                    BL_info_list = BL_info.split(")")
                    for sID_nBL in BL_info_list:
                        iLen_branch = float(sID_nBL.split(":")[-1])
                        iSum_branch += iLen_branch
                else:
                    sID_nBL = BL_info
                    iLen_branch = float(sID_nBL.split(":")[-1])
                    iSum_branch += iLen_branch
            sID = monoTarget_list[-1]
            BL_info_list = info_tree[sID_list.index(sID)].split(")")[:-diff_RminusL]
            for sID_nBL in BL_info_list:
                iLen_branch = float(sID_nBL.split(":")[-1])
                iSum_branch += iLen_branch
            #iLen_OriBranch = float(info_tree[i].split(")")[-diff_RminusL-1].split(":")[-1])

    
        else:#iCNT_R_parenthesis == iCNT_L_parenthesis
            # Most common ancestral branch
            SUM_OF_WHOLE_BL = 0
            for sID in sID_list:
                BL_info = info_tree[sID_list.index(sID)]
                #print(BL_info)
                #print(BL_info)
                if ")" in BL_info:
                    BL_info_list = BL_info.split(")")
                    for sID_nBL in BL_info_list:
                        iLen_branch = float(sID_nBL.split(":")[-1])
                        SUM_OF_WHOLE_BL += iLen_branch
                else:
                    sID_nBL = BL_info
                    iLen_branch = float(sID_nBL.split(":")[-1])
                    SUM_OF_WHOLE_BL += iLen_branch
            iLen_MRCA_branch = find_originbranchlength(sID_list, info_tree, sID_list)
            SUM_OF_WHOLE_BL -= iLen_MRCA_branch
            SUM_OF_WHOLE_BL = round(SUM_OF_WHOLE_BL,iDP)
            iSum_branch = SUM_OF_WHOLE_BL
            
    #print(iSum_branch)
    return iSum_branch


def scan_right_parenthesis(iPos_sTarget, iCNT_L_parenthesis, info_tree, sID_list):
    iCNT_L_parenthesis += info_tree[iPos_sTarget+1].count('(')
    iCNT_R_parenthesis = info_tree[iPos_sTarget+1].count(')')
    if iCNT_L_parenthesis <= iCNT_R_parenthesis:
        nInternalBranch = info_tree[iPos_sTarget+1].split(')')[iCNT_L_parenthesis]
        if ":" in nInternalBranch:
            iLen_InternalBranch_target = float(nInternalBranch.split(':')[-1])
        else:
            iLen_InternalBranch_target = "Branch of most common ancestor"
        endTarget = sID_list[iPos_sTarget+1]
        return endTarget, iLen_InternalBranch_target
    else:
        iCNT_L_parenthesis -= iCNT_R_parenthesis
        endTarget, iLen_InternalBranch_target = scan_right_parenthesis(iPos_sTarget+1, iCNT_L_parenthesis, info_tree, sID_list)
        return endTarget, iLen_InternalBranch_target


def scan_left_parenthesis(iPos_sTarget, iCNT_R_parenthesis, info_tree, sID_list):
    iCNT_R_parenthesis += info_tree[iPos_sTarget-1].count(')')
    iCNT_L_parenthesis = info_tree[iPos_sTarget-1].count('(')
    if iCNT_R_parenthesis <= iCNT_L_parenthesis:
        startTarget = sID_list[iPos_sTarget-1]
        return startTarget
    else:
        iCNT_R_parenthesis -= iCNT_L_parenthesis
        startTarget = scan_left_parenthesis(iPos_sTarget-1, iCNT_R_parenthesis, info_tree, sID_list)
        return startTarget


def find_pairedgroup(pairedID, direction_parenthesis, info_tree, sID_list):# Right parenthesis
    if direction_parenthesis == ")":
        if info_tree[sID_list.index(pairedID)].count(')')>0:
            pairedID_list = [pairedID]
        else:
            iCNT_R_parenthesis = 0
            iCNT_L_parenthesis = 0
            iRange_list = range(sID_list.index(pairedID),len(sID_list),1)
            for i in iRange_list:
                sID = sID_list[i]
                if i == 0:
                    iCNT_L_parenthesis -= 1
                iCNT_R_parenthesis += info_tree[i].count(")") 
                iCNT_L_parenthesis += info_tree[i].count("(")
                if iCNT_R_parenthesis >= iCNT_L_parenthesis:
                    pairedID_list = sID_list[sID_list.index(pairedID):sID_list.index(sID)+1]
                    break
                else:pass
    else:# direction_parenthesis == "("
        if info_tree[sID_list.index(pairedID)].count('(')>0:
            pairedID_list = [pairedID]
        else:
            iCNT_R_parenthesis = 0
            iCNT_L_parenthesis = 0
            iRange_list = range(sID_list.index(pairedID),-1,-1)
            for i in iRange_list:
                sID = sID_list[i]
                if i == len(sID_list)-1:
                    iCNT_R_parenthesis -= 1
                iCNT_R_parenthesis += info_tree[i].count(")") 
                iCNT_L_parenthesis += info_tree[i].count("(")
                if iCNT_L_parenthesis >= iCNT_R_parenthesis: 
                    pairedID_list = sID_list[sID_list.index(sID):sID_list.index(pairedID)+1]
                    break
                else:pass
    return pairedID_list


def find_upper_clade(sTarget_list,iSum_OriginBranch, target_list, info_tree, sID_list):#sTarget_list = target species in a monophyletic clade 
    startID = sTarget_list[0]
    endID = sTarget_list[-1]
    info_tree_of_sTarget_list = info_tree[sID_list.index(startID):sID_list.index(endID)+1]
    iCNT_L_parenthesis = 0 # (
    iCNT_R_parenthesis = 0 # )
    for info_tree_sTarget in info_tree_of_sTarget_list:
        iCNT_L_parenthesis += info_tree_sTarget.count('(')
        iCNT_R_parenthesis += info_tree_sTarget.count(')')
        
    if iCNT_L_parenthesis > iCNT_R_parenthesis:
        pairedID              = sID_list[sID_list.index(endID)+1]
        direction_parenthesis = ")"
        pairedID_list         = find_pairedgroup(pairedID, direction_parenthesis, info_tree, sID_list)
        iCNT_target = 0
        for pairedID in pairedID_list:
            if pairedID in target_list:
                iCNT_target +=1
        if iCNT_target == len(pairedID_list):
            ORIGIN_sTarget_list = sTarget_list + pairedID_list
            ORIGIN_sTarget_list = ordering_sID(ORIGIN_sTarget_list)#, info_tree, sID_list
            ORIGIN_iSum_OriginBranch = calculate_sum_branchlength(ORIGIN_sTarget_list, info_tree, sID_list)
        else:
            iSum_OriginBranch = calculate_sum_branchlength(sTarget_list, info_tree, sID_list)
            ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch = sTarget_list, iSum_OriginBranch
            
    elif iCNT_L_parenthesis < iCNT_R_parenthesis:
        pairedID              = sID_list[sID_list.index(startID)-1]
        direction_parenthesis = "("
        pairedID_list         = find_pairedgroup(pairedID, direction_parenthesis, info_tree, sID_list)
        iCNT_target = 0
        for pairedID in pairedID_list:
            if pairedID in target_list:
                iCNT_target +=1
        if iCNT_target == len(pairedID_list):
            ORIGIN_sTarget_list = sTarget_list + pairedID_list
            ORIGIN_sTarget_list = ordering_sID(ORIGIN_sTarget_list)
            ORIGIN_iSum_OriginBranch = calculate_sum_branchlength(ORIGIN_sTarget_list, info_tree, sID_list)
        else:
            iSum_OriginBranch = calculate_sum_branchlength(sTarget_list, info_tree, sID_list)
            ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch = sTarget_list, iSum_OriginBranch
            
    else:#iCNT_L_parenthesis == iCNT_R_parenthesis
        iSum_OriginBranch = calculate_sum_branchlength(sTarget_list, info_tree, sID_list)
        ORIGIN_sTarget_list = sTarget_list
        ORIGIN_iSum_OriginBranch = iSum_OriginBranch
    
    return ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch 
    

def scan_outer_internalbranch(sTarget_list,iSum_OriginBranch, target_list, info_tree, sID_list):
    ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch = find_upper_clade(sTarget_list,iSum_OriginBranch, target_list, info_tree, sID_list)
    if not ORIGIN_sTarget_list == sTarget_list:
        sTarget_list = ORIGIN_sTarget_list
        iSum_OriginBranch = ORIGIN_iSum_OriginBranch
        ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch = scan_outer_internalbranch(sTarget_list,iSum_OriginBranch, target_list, info_tree, sID_list)
    return ORIGIN_sTarget_list, ORIGIN_iSum_OriginBranch

    
def trim_unselected_target_list(sTarget_list,unselected_target_list):
    nTarget_list = []
    for nTarget in unselected_target_list:
        if not nTarget in sTarget_list:
            nTarget_list.append(nTarget)
    unselected_target_list = nTarget_list
    return unselected_target_list


def scan_OriBranch(sTarget, unselected_target_list, target_list, info_tree, sID_list, nMarker_target):
    sTarget_list = [sTarget]
    
    # 1. Set Origin branch length = terminal branch length
    iSum_TerminalBranch_target = calculate_sum_branchlength(sTarget_list, info_tree, sID_list)
    iSum_OriginBranch = iSum_TerminalBranch_target

    #sys.exit()

    # 2.find ancestral branch
    iPos_sTarget = sID_list.index(sTarget)
    sID_in_clade_list = [sTarget]
    if '(' in info_tree[iPos_sTarget]:
        iCNT_L_parenthesis = 1
        endTarget, iLen_InternalBranch_target = scan_right_parenthesis(iPos_sTarget, iCNT_L_parenthesis, info_tree, sID_list)
        sID_inInternalBranch_list = sID_list[iPos_sTarget:sID_list.index(endTarget)+1]
        
    else:
        nInternalBranch = info_tree[iPos_sTarget].split(')')[1]
        if ":" in nInternalBranch:
            iLen_InternalBranch_target = float(nInternalBranch.split(':')[-1])
        else:
            iLen_InternalBranch_target = 0 #"Branch of most common ancestor"

        iCNT_R_parenthesis = 1
        startTarget = scan_left_parenthesis(iPos_sTarget, iCNT_R_parenthesis, info_tree, sID_list)
        sID_inInternalBranch_list = sID_list[sID_list.index(startTarget):iPos_sTarget+1]

    sID_AncestralBranch_list    = sID_inInternalBranch_list
    iLen_AncestralBranch        = calculate_sum_branchlength(sID_AncestralBranch_list, info_tree, sID_list)
    
    ###
    # targets == species under ancestral branch?
    # Yes: iFlag = 1 / No: iFlag = 0
    iFlag = 0
    iCNT_eqID = 0
    for sID in sID_AncestralBranch_list:
        if sID in target_list:
            iCNT_eqID +=1
    if iCNT_eqID == len(sID_AncestralBranch_list):
        sTarget_list = sID_AncestralBranch_list
        iSum_OriginBranch = iLen_AncestralBranch
        iFlag = 1

    ###
    #
    if iFlag == 0:
        unselected_target_list = trim_unselected_target_list(sTarget_list,unselected_target_list)
        print_line = "# "+ nMarker_target+" group: "+','.join(sTarget_list)+' / SumBL: '+str(round(iSum_OriginBranch,iDP))
        print(print_line)
        fpout.write(print_line+'\n')
        return iSum_OriginBranch, unselected_target_list

    else:#iFlag == 1
        sTarget_list, iSum_OriginBranch = scan_outer_internalbranch(sTarget_list, iSum_OriginBranch, target_list, info_tree, sID_list)
        unselected_target_list = trim_unselected_target_list(sTarget_list,unselected_target_list)
        print_line = "# "+ nMarker_target+" group: "+','.join(sTarget_list)+' / SumBL: '+str(round(iSum_OriginBranch,iDP))
        print(print_line)
        fpout.write(print_line+'\n')
        return iSum_OriginBranch, unselected_target_list



def find_upper_monophyleticlineage(selected_sID_list, info_tree, sID_list):#sTarget_list = target species in a monophyletic clade
    # Description for iFlag_mono: if selected target list (sTarget_list) does not include all of targets: 0, else: 1
    startID = selected_sID_list[0]
    endID = selected_sID_list[-1]
    info_tree_of_selected_sID_list = info_tree[sID_list.index(startID):sID_list.index(endID)+1]
    iCNT_L_parenthesis = 0 # (
    iCNT_R_parenthesis = 0 # )
    for info_tree_selected_sID in info_tree_of_selected_sID_list:
        iCNT_L_parenthesis += info_tree_selected_sID.count('(')
        iCNT_R_parenthesis += info_tree_selected_sID.count(')')
        
    if iCNT_L_parenthesis > iCNT_R_parenthesis:
        pairedID              = sID_list[sID_list.index(endID)+1]
        direction_parenthesis = ")"
        pairedID_list         = find_pairedgroup(pairedID, direction_parenthesis, info_tree, sID_list)
        monophyletic_sID_list   = selected_sID_list + pairedID_list
        monophyletic_sID_list   = ordering_sID(monophyletic_sID_list)
        
    elif iCNT_L_parenthesis < iCNT_R_parenthesis:
        pairedID              = sID_list[sID_list.index(startID)-1]
        direction_parenthesis = "("
        pairedID_list         = find_pairedgroup(pairedID, direction_parenthesis, info_tree, sID_list)
        monophyletic_sID_list   = selected_sID_list + pairedID_list
        monophyletic_sID_list   = ordering_sID(monophyletic_sID_list)
            
    else:#iCNT_L_parenthesis == iCNT_R_parenthesis
        monophyletic_sID_list = sID_list
        
    return monophyletic_sID_list


def scan_outer_internalbranch_to_find_monophyleticlineage_of_MRCA_of_targets(selected_sID_list, target_list, info_tree, sID_list):
    selected_sID_list = find_upper_monophyleticlineage(selected_sID_list, info_tree, sID_list)
    iCNT_target_in_monophyclade = 0
    for sID in target_list:
        if sID in selected_sID_list:
            iCNT_target_in_monophyclade += 1
    if not iCNT_target_in_monophyclade == len(target_list):
        selected_sID_list = scan_outer_internalbranch_to_find_monophyleticlineage_of_MRCA_of_targets(selected_sID_list, target_list, info_tree, sID_list)
    return selected_sID_list


def find_MRCA_of_target(target_list, info_tree, sID_list):#sTarget_list = target species in a monophyletic clade
    selected_sID_list = []
    selected_sID_list.append(target_list[0])
    selected_sID_list = find_upper_monophyleticlineage(selected_sID_list, info_tree, sID_list)
    iCNT_target_in_monophyclade = 0
    for sID in target_list:
        if sID in selected_sID_list:
            iCNT_target_in_monophyclade += 1
    if not iCNT_target_in_monophyclade == len(target_list):
        selected_sID_list = scan_outer_internalbranch_to_find_monophyleticlineage_of_MRCA_of_targets(selected_sID_list, target_list, info_tree, sID_list)
    return selected_sID_list


def make_trimmed_tree_selected_speceis(selected_sID_list):
    iCNT_selected_species  = len(selected_sID_list)
    tmp_selected_info_tree_list = nTree.replace(";",":0").split(',')[sID_list.index(selected_sID_list[0]):sID_list.index(selected_sID_list[-1])+1]
    iCNT_right_parenthesis = 0
    for sID_info in tmp_selected_info_tree_list[:-1]:
        iCNT_right_parenthesis += sID_info.count(")")
    
    iCNT_right_parenthesis = iCNT_selected_species - iCNT_right_parenthesis - 1
    lastID_info_list = tmp_selected_info_tree_list[-1].strip().split(')')[:iCNT_right_parenthesis]
    lastID_info = ')'.join(lastID_info_list)+'):0'
    selected_info_tree = tmp_selected_info_tree_list[:-1] + [lastID_info]
    return selected_info_tree


def Run_SumOfBL(sID_list, info_tree, target_list):
    ## calculate sum of whole branch lengths
    SUM_OF_WHOLE_BL = 0
    for sID in sID_list:
        BL_info = info_tree[sID_list.index(sID)]
        #print(BL_info)
        #print(BL_info)
        if ")" in BL_info:
            BL_info_list = BL_info.split(")")
            for sID_nBL in BL_info_list:
                iLen_branch = float(sID_nBL.split(":")[-1])
                SUM_OF_WHOLE_BL += iLen_branch
        else:
            sID_nBL = BL_info
            iLen_branch = float(sID_nBL.split(":")[-1])
            SUM_OF_WHOLE_BL += iLen_branch
    iLen_MRCA_branch = find_originbranchlength(sID_list, info_tree, sID_list)#monotarget, info_tree, sID_list
    SUM_OF_WHOLE_BL -= iLen_MRCA_branch
    SUM_OF_WHOLE_BL = round(SUM_OF_WHOLE_BL,iDP)
    print("#SUM_OF_WHOLE_BL_under_MRCAofTargets:",SUM_OF_WHOLE_BL)
    print("### Result")


    SUM_OF_TERMINAL_BL = 0
    for sTarget in target_list:
        iLen_termianl_BL = find_originbranchlength([sTarget], info_tree, sID_list)
        SUM_OF_TERMINAL_BL += iLen_termianl_BL
    SUM_OF_TERMINAL_BL = round(SUM_OF_TERMINAL_BL,iDP)

    if sID_list == target_list:# Target species are in a monophletic lineage
        print("# Target group: " + ','.join(target_list) + ' / DistanceTN: '+ str(SUM_OF_WHOLE_BL)+' / DistanceTB: '+str(SUM_OF_WHOLE_BL-SUM_OF_TERMINAL_BL))
        print("# Others group: None")
        target_clade_list = [target_list]
        DistanceTN = SUM_OF_WHOLE_BL
        DistanceTB = SUM_OF_WHOLE_BL - SUM_OF_TERMINAL_BL
        


    else:# Target species are in independent lineages (=polyphyletic lineages)
        ordered_target_list = target_list
        ordered_target_list = ordering_sID(ordered_target_list)
        unselected_target_list = ordered_target_list
        iCNT_origins_targets   = 0
        target_clade_list      = []
        nMarker_target = "Target"
        iCNT_origins_targets = 0
        for sTarget in ordered_target_list:
            if sTarget in unselected_target_list:
                iSum_OriginBranch, tmp_unselected_target_list = scan_OriBranch(sTarget, unselected_target_list, target_list, info_tree, sID_list, nMarker_target)

                selected_target_list = []
                for target_ID in unselected_target_list:
                    if not target_ID in tmp_unselected_target_list:
                        selected_target_list.append(target_ID)
                target_clade_list.append(selected_target_list)
                unselected_target_list = tmp_unselected_target_list
            
                if len(unselected_target_list)==0:
                    break

        others_list = []
        for sID in sID_list:
            if not sID in target_list:
                others_list.append(sID)
                
        ordered_others_list = others_list
        ordered_others_list = ordering_sID(ordered_others_list)
        unselected_others_list = ordered_others_list
        iSum_OriginBranch_others_list = [] # branches under most common ancestor of monophyletic others (including MRCA branch)
        others_clade_list      = []
        nMarker_target = "Others"
        for sOthers in ordered_others_list:
            if sOthers in unselected_others_list:
                iSum_OriginBranch, tmp_unselected_others_list = scan_OriBranch(sOthers, unselected_others_list, others_list, info_tree, sID_list, nMarker_target)
                iSum_OriginBranch_others_list.append(iSum_OriginBranch)

                selected_others_list = []
                for others_ID in unselected_others_list:
                    if not others_ID in tmp_unselected_others_list:
                        selected_others_list.append(others_ID)
                others_clade_list.append(selected_others_list)
                unselected_others_list = tmp_unselected_others_list
            
                if len(unselected_others_list)==0:
                    break

        DistanceTN = SUM_OF_WHOLE_BL
        for iSum_BranchesUnderOriginBranch in iSum_OriginBranch_others_list:
            DistanceTN -= iSum_BranchesUnderOriginBranch
        DistanceTB = DistanceTN - SUM_OF_TERMINAL_BL

    DistanceTN = round(DistanceTN,iDP)
    DistanceTB = round(DistanceTB,iDP)

    ProOriBL = 1
    SumOriBL = 0
    iLen_origin_branch_list = []
    for monophyletic_target_list in target_clade_list:
        iLen_origin_branch = find_originbranchlength(monophyletic_target_list, nTree.split(','), sID_LIST)#monotarget, info_tree, sID_list
        iLen_origin_branch_list.append(iLen_origin_branch)
        ProOriBL *= round(iLen_origin_branch,iDP)
        SumOriBL += round(iLen_origin_branch,iDP)
    
    ProTerBL = 1
    SumTerBL = 0
    iLen_terminal_branch_list = []
    for targetID in target_list:
        iLen_terminal_branch = find_originbranchlength([targetID], nTree.split(','), sID_LIST)
        iLen_terminal_branch_list.append(iLen_terminal_branch)
        ProTerBL *= round(iLen_terminal_branch,iDP)
        SumTerBL += round(iLen_terminal_branch,iDP)

    

    print("#")
    print("### Summary")
    print("# Number of independent origins of targets: " + str(len(target_clade_list)))    
    print("# Distance between Termianl Nodes (DTN): " + str(DistanceTN))
    print("# Distance between Terminal Branches (DTB): " + str(DistanceTB))
    print("# Product of origin branch length (POB): " + str(ProOriBL))
    print("# Sum of origin branch length (SOB): " + str(SumOriBL))
    print("# Product of terminal branch length (PTB): " + str(ProTerBL))
    print("# Sum of terminal branch length (STB): " + str(SumTerBL))
    print("# DistanceTN,DistanceTB,ProOriBL,SumOriBL,ProTerBL,SumTerBL")
    print(str(DistanceTN)+','+str(DistanceTB)+','+str(ProOriBL)+','+str(SumOriBL)+','+str(ProTerBL)+','+str(SumTerBL)+'\n')
    
    fpout.write("# Number of independent origins of targets: " + str(len(target_clade_list))+'\n')
    fpout.write("# Distance between Termianl Nodes (DTN): " + str(DistanceTN)+'\n')
    fpout.write("# Distance between Terminal Branches (DTB): " + str(DistanceTB)+'\n')
    fpout.write("# Product of origin branch length (POB): " + str(ProOriBL)+'\n')
    fpout.write("# Sum of origin branch length (SOB): " + str(SumOriBL)+'\n')
    fpout.write("# Product of terminal branch length (PTB): " + str(ProTerBL)+'\n')
    fpout.write("# Sum of terminal branch length (STB): " + str(SumTerBL)+'\n')
    fpout.write("# DistanceTN,DistanceTB,ProOriBL,SumOriBL,ProTerBL,SumTerBL" +'\n')
    fpout.write(str(DistanceTN)+','+str(DistanceTB)+','+str(ProOriBL)+','+str(SumOriBL)+','+str(ProTerBL)+','+str(SumTerBL)+'\n')
    fpout.close()


# main
selected_sID_list = find_MRCA_of_target(target_list, info_tree, sID_list)
selected_info_tree = make_trimmed_tree_selected_speceis(selected_sID_list)
print("# Tree information under MRCA of targets")
print("#",selected_info_tree)
Run_SumOfBL(selected_sID_list, selected_info_tree, target_list)

####################### End of Code, written by Chul Lee
