# ConVarFinder is a python script to detect convergent variants at codon and amino acid levels
# version v1 (05-14-2025)
# Originally written by Chul Lee (c) Rockefeller Univ. email: chul.bioinfo@gmail.com
# GitHub: https://github.com/chulbioinfo/ConVarFinder

"""
ConVarFinder: A tool for detecting convergent variants at codon and amino acid levels

This script analyzes genetic sequences from multiple species to identify convergent
evolutionary variants. It can detect convergence at three levels:
1. Codon level
2. Amino acid level
3. Nucleotide level (positions 1, 2, and 3 within codons)

Usage:
    python ConVarFinder.py [options]

Options:

"""

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
  """Reads a FASTA sequence file and returns a dictionary of sequences."""
  sID_nSeq_dic = {}
  try:
    with open(fNAME_seq, 'r') as fpin:
      sID = ''
      for line in fpin:
        if line.startswith(">"):
          sID = line[1:].strip()
          sID_nSeq_dic[sID] = ''
        elif sID: # Ensure sID is set before adding sequence
          sID_nSeq_dic[sID] += line.strip()
  except FileNotFoundError:
    print(f"Error: Sequence file not found: {fNAME_seq}")
    sys.exit()
  return sID_nSeq_dic


def read_tree_file(fNAME_tree):
  """Reads a Newick tree file and returns the tree string."""
  nTree = ''
  try:
    with open(fNAME_tree, 'r') as fpin:
      for line in fpin:
        nTree += line.strip()
  except FileNotFoundError:
    print(f"Error: Tree file not found: {fNAME_tree}")
    sys.exit()
  return nTree

    
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
    target_list = [nTargets] if nTargets else [] # Handle empty nTargets

  if nOutgroup and isinstance(nOutgroup, str) and ',' in nOutgroup:
    outgroup_list = nOutgroup.split(",")
  elif nOutgroup and isinstance(nOutgroup, str): # Single outgroup species
    outgroup_list = [nOutgroup]
  else: # No outgroup or nOutgroup is not a string (e.g. already a list or None)
    outgroup_list = []

  tmp_target_list = []
  tmp_others_list = []
  tmp_outgroup_list = []

  for sID in sID_list:
    if sID in outgroup_list:
      tmp_outgroup_list.append(sID)
    elif sID in target_list:
      tmp_target_list.append(sID)
    else:
      tmp_others_list.append(sID)
  return(tmp_target_list, tmp_others_list, tmp_outgroup_list)

def codon2aa(nCodon, codon_table):
  """Translates a codon to an amino acid using the provided codon table."""
  try:
    # Ensure codon is uppercase to match keys in CODONTABLE_dic
    nAA = codon_table[nCodon.upper()]
  except KeyError:
    nAA = 'X' # Default for unknown codons
  return nAA

def _validate_required_options(options):
  """Validates if required options are provided."""
  missing_required = []
  if not options.get('seq_file'):
    missing_required.append("-seq (Sequence file)")
  if not options.get('tree_file'):
    missing_required.append("-tree (Tree file)")

  if missing_required:
    print("Error: The following required options are missing:")
    for missing in missing_required:
      print(f"  {missing}")
    sys.exit(1) # Exit with a non-zero code to indicate error

def _get_available_species(tree_file, seq_file):
  """Gets available species present in both tree and sequence files."""
  templete_nTree = read_tree_file(tree_file)
  templete_nTree = read_tre_removing_branchlength(templete_nTree)
  terminal_nodes = terminal_node_from_tree(templete_nTree)
  sID_nSeq_dic = read_seq_file(seq_file)
  
  available_species = [species for species in terminal_nodes if species in sID_nSeq_dic]
  
  if not available_species:
    print("Error: No matching species found between tree and sequence files.")
    sys.exit(1)
  return available_species, templete_nTree, sID_nSeq_dic # Return sID_nSeq_dic as well

def _prompt_species_selection(prompt_message, species_list, allow_empty=False, item_name="species"):
    """Prompts user to select species from a list and returns selected species as a comma-separated string."""
    if not species_list:
        if not allow_empty:
             print(f"Error: No {item_name} available for selection.")
             return "" # Or sys.exit(1) if this is a critical error
        else:
            return ""

    print(f"\n{prompt_message}")
    for i, species in enumerate(species_list, 1):
        print(f"{i}. {species}")
    
    while True:
        try:
            selection_prompt = f"\nEnter the numbers of {item_name} (comma-separated"
            if allow_empty:
                selection_prompt += ", or press Enter for none): "
            else:
                selection_prompt += "): "
            
            selection_input = input(selection_prompt)
            if not selection_input.strip() and allow_empty:
                return ""
            
            selected_indices = [int(x.strip()) for x in selection_input.split(',')]
            selected_species_list = [species_list[i-1] for i in selected_indices if 0 < i <= len(species_list)] 

            if not selected_species_list and not allow_empty:
                print(f"Error: At least one {item_name} must be selected.")
                continue
            
            return ','.join(selected_species_list)
        except (ValueError, IndexError):
            print(f"Error: Please enter valid numbers (1-{len(species_list)}) separated by commas.")
        except KeyboardInterrupt: 
            print("\nSelection cancelled by user.")
            sys.exit(0)

def set_options():
  """Sets program options based on command-line arguments and user prompts. Returns a dictionary of options."""
  # Default options
  options = {
    'nTargets': "",
    'nOutgroup': "",
    'iPATH': "./",
    'fNAME_seq_format': ".fasta",
    'fNAME_tree_format': ".tre",
    'oPATH': "./",
    'seq_file': "",
    'tree_file': "",
    'skip_outgroup_prompt': False,
    'templete_mID_targetID_dic': {},
    'templete_mID_mmID_dic': {},
    'templete_sID_list_for_analysis': [],
    'sID_nSeq_dic': {} # To store sequence data after reading it once
  }

  # Option mapping: command-line argument -> options key
  option_mapping = {
    '-tl=': 'nTargets',
    '-ol=': 'nOutgroup',
    '-ip=': 'iPATH',
    '-tfmt=': 'fNAME_tree_format',
    '-sfmt=': 'fNAME_seq_format',
    '-op=': 'oPATH',
    '-seq=': 'seq_file',
    '-tree=': 'tree_file',
    '-no-outgroup-prompt': 'skip_outgroup_prompt' 
  }

  help_message = [
    "ConVarFinder v1 - Convergent Variant Finder", 
    "Usage: python ConVarFinder.py [options]",
    "",
    "Dependencies:",
    "  - Python 3.x",
    "  - No additional packages required (uses only standard library)",
    "",
    "Required options:",
    "  -seq=FILE          Sequence file in FASTA format (e.g., -seq=10050.fasta)",
    "  -tree=FILE         Tree file in Newick format (e.g., -tree=10050.tre)",
    "",
    "Optional options:",
    "  -tl=SPECIES_LIST    List of target species (comma-separated).",
    "                      If not specified, prompts selection from available terminal nodes.",
    "  -ol=SPECIES_LIST    List of outgroup species (comma-separated).",
    "                      If not specified, prompts selection from remaining species.",
    "                      Press Enter to skip outgroup selection if prompting.",
    "  -no-outgroup-prompt Skip the outgroup selection prompt.",
    "  -sfmt=FORMAT        Sequence file format (Default: .fasta)",
    "  -tfmt=FORMAT        Tree file format (Default: .tre)",
    "  -ip=PATH            Input directory path (Default: ./)",
    "  -op=PATH            Output directory path (Default: ./)",
    "",
    "Example:",
    "  python ConVarFinder.py -seq=sequences.fasta -tree=tree.tre -tl=speciesA,speciesB",
    "",
    "Note:",
    "  - Options can be specified in any order.",
    "  - Output file: 'ConVarFinder_[TARGET_SPECIES].txt' in the output directory."
  ]

  args = sys.argv[1:]
  for arg in args:
    if arg == '-h' or arg == '--help':
      for line in help_message:
        print(line)
      sys.exit(0)
    
    processed_arg = False
    for opt_prefix, opt_key in option_mapping.items():
      if opt_key == 'skip_outgroup_prompt': 
          if arg == '-no-outgroup-prompt':
              options[opt_key] = True
              processed_arg = True
              break
      elif arg.startswith(opt_prefix):
        value = arg[len(opt_prefix):]
        if opt_key in ['fNAME_tree_format', 'fNAME_seq_format'] and not value.startswith('.'):
          value = '.' + value
        options[opt_key] = value
        processed_arg = True
        break
    if not processed_arg and arg not in ('-no-outgroup-prompt'): 
        print(f"Warning: Unknown or malformed argument: {arg}")

  _validate_required_options(options)

  available_species, templete_nTree, sID_nSeq_dic_from_file = _get_available_species(options['tree_file'], options['seq_file'])
  options['sID_nSeq_dic'] = sID_nSeq_dic_from_file # Store sequence data in options
  
  if not options['nTargets']:
    options['nTargets'] = _prompt_species_selection(
        "Available terminal nodes in tree (for target species selection):",
        available_species,
        allow_empty=False,
        item_name="target species"
    )
    if not options['nTargets']: 
        print("Error: No target species selected. Exiting.")
        sys.exit(1)
    print(f"\nSelected target species: {options['nTargets']}")

  selected_targets_list = options['nTargets'].split(',') if options['nTargets'] else []
  remaining_species = [sp for sp in available_species if sp not in selected_targets_list]
  
  if not options['nOutgroup'] and not options['skip_outgroup_prompt']:
    if remaining_species:
      options['nOutgroup'] = _prompt_species_selection(
          "Remaining species available for outgroup:",
          remaining_species,
          allow_empty=True,
          item_name="outgroup species"
      )
      if options['nOutgroup']:
        print(f"\nSelected outgroup species: {options['nOutgroup']}")
      else:
        print("\nNo outgroup species selected (skipped by user).")
    else: 
      print("\nNo remaining species available for outgroup selection. No outgroup will be used.")
      options['nOutgroup'] = "" 
  elif options['skip_outgroup_prompt']:
      print("\nOutgroup selection prompt skipped by -no-outgroup-prompt option. No outgroup will be used unless specified with -ol.")
      if not options['nOutgroup']: 
          options['nOutgroup'] = ""

  templete_mID_targetID_dic, templete_mID_mmID_dic = Anc_Finder(templete_nTree, options['nTargets'])
  options['templete_mID_targetID_dic'] = templete_mID_targetID_dic
  options['templete_mID_mmID_dic'] = templete_mID_mmID_dic
  
  tree_terminal_nodes = terminal_node_from_tree(templete_nTree) 
  
  _, _, outgroup_list_parsed = parse_target_others_from_list(
      options['nTargets'], 
      tree_terminal_nodes, # Use all terminal nodes from the tree for this parsing
      options['nOutgroup']
  )
  options['outgroup_list'] = outgroup_list_parsed # Store the parsed list of outgroups

  # This list is for the species columns in the output. It should contain all species analyzed (targets + others).
  options['templete_sID_list_for_analysis'] = [s for s in tree_terminal_nodes if s not in outgroup_list_parsed]

  return options


def scan_ancestors_from_mID_to_targetID(mID, targetID, sID_mID_dic, mID_to_targetID_list):
  tmp_mID = sID_mID_dic[targetID]
  mID_to_targetID_list.append(tmp_mID)
  if not mID == tmp_mID:
    tmp_mID, mID_to_targetID_list = scan_ancestors_from_mID_to_targetID(mID, tmp_mID, sID_mID_dic, mID_to_targetID_list)
  return(tmp_mID, mID_to_targetID_list)


def _extract_site_specific_data(sID_nSeq_dic, terNode_list, local_target_list, outgroup_list, site_idx, codon_table):
    """Extracts codon, AA, and nucleotide data for a specific site from sequences."""
    site_data = {
        'target_codons': [], 'others_codons': [],
        'target_aas': [], 'others_aas': [],
        'target_nuc1s': [], 'others_nuc1s': [],
        'target_nuc2s': [], 'others_nuc2s': [],
        'target_nuc3s': [], 'others_nuc3s': []
    }
    for sID in terNode_list:
        if sID not in outgroup_list and sID in sID_nSeq_dic: # Ensure sID is in sequence data
            seq = sID_nSeq_dic[sID]
            if site_idx + 3 <= len(seq): # Boundary check for codon extraction
                codon = seq[site_idx : site_idx+3]
                nuc1, nuc2, nuc3 = codon[0], codon[1], codon[2]
                aa = codon2aa(codon, codon_table)

                if sID in local_target_list:
                    site_data['target_codons'].append(codon)
                    site_data['target_aas'].append(aa)
                    site_data['target_nuc1s'].append(nuc1)
                    site_data['target_nuc2s'].append(nuc2)
                    site_data['target_nuc3s'].append(nuc3)
                else:
                    site_data['others_codons'].append(codon)
                    site_data['others_aas'].append(aa)
                    site_data['others_nuc1s'].append(nuc1)
                    site_data['others_nuc2s'].append(nuc2)
                    site_data['others_nuc3s'].append(nuc3)
            else:
                # Handle cases where sequence is too short for the site_idx
                # print(f"Warning: Sequence for {sID} is too short for site index {site_idx}.")
                pass # Or append placeholder values like '---' or 'X' if needed for alignment
    return site_data

def _check_mutually_exclusive_substitutions(site_data):
    """Checks for mutually exclusive substitutions between target and other groups."""
    target_codons_set = set(site_data['target_codons'])
    others_codons_set = set(site_data['others_codons'])
    target_aas_set = set(site_data['target_aas'])
    others_aas_set = set(site_data['others_aas'])
    target_nuc1s_set = set(site_data['target_nuc1s'])
    others_nuc1s_set = set(site_data['others_nuc1s'])
    target_nuc2s_set = set(site_data['target_nuc2s'])
    others_nuc2s_set = set(site_data['others_nuc2s'])
    target_nuc3s_set = set(site_data['target_nuc3s'])
    others_nuc3s_set = set(site_data['others_nuc3s'])

    # If no target codons/AAs/Nucs, there can't be an exclusive substitution.
    # This also handles cases where target sets might be empty due to no valid data for targets at a site.
    if not target_codons_set and not target_aas_set and not target_nuc1s_set and not target_nuc2s_set and not target_nuc3s_set:
        return {
            'codon': True, 'aa': True, 
            'nuc1': True, 'nuc2': True, 'nuc3': True
        } # All True means no exclusive substitution (overlap or empty target set)

    flags = {
        # True if overlap (NOT exclusive), False if exclusive (no overlap)
        'codon': bool(target_codons_set.intersection(others_codons_set)) if target_codons_set else True,
        'aa':    bool(target_aas_set.intersection(others_aas_set))    if target_aas_set    else True,
        'nuc1':  bool(target_nuc1s_set.intersection(others_nuc1s_set))  if target_nuc1s_set  else True,
        'nuc2':  bool(target_nuc2s_set.intersection(others_nuc2s_set))  if target_nuc2s_set  else True,
        'nuc3':  bool(target_nuc3s_set.intersection(others_nuc3s_set))  if target_nuc3s_set  else True 
    }
    return flags

def _collect_ancestral_path_data(mID_targetID_dic, sID_mID_dic, sID_nSeq_dic, site_idx, codon_table):
    """Collects data along ancestral paths for MRCA, AMRCA and targets."""
    path_data = {
        'Amrca_states': {'codon': [], 'aa': [], 'nuc1': [], 'nuc2': [], 'nuc3': []},
        'mrca_states': {'codon': [], 'aa': [], 'nuc1': [], 'nuc2': [], 'nuc3': []},
        'Amrca2mrca_IDs': [], 'Amrca2mrca_paths': {'codon': [], 'aa': [], 'nuc1': [], 'nuc2': [], 'nuc3': []},
        'mrca2target_IDs': [], 'mrca2target_paths': {'codon': [], 'aa': [], 'nuc1': [], 'nuc2': [], 'nuc3': []},
        'clade_substitutions': {'codon': 0, 'aa': 0, 'nuc1': 0, 'nuc2': 0, 'nuc3': 0}
    }

    targetID_mrca2target_ID_paths_dic = {}
    for mID in mID_targetID_dic.keys():
        for targetID in mID_targetID_dic[mID]:
            current_path_IDs = []
            if not mID == targetID: # A clade with multiple species (n>=2)
                # scan_ancestors_from_mID_to_targetID is defined globally
                # It returns (last_mID_reached, path_list_excluding_targetID)
                _, path_nodes = scan_ancestors_from_mID_to_targetID(mID, targetID, sID_mID_dic, []) # Pass empty list for path
                current_path_IDs.extend(path_nodes)
            current_path_IDs.reverse() # Path was from target to mID, reverse to be mID to target_ancestor
            current_path_IDs.append(targetID) # Add the target itself to the end of the path
            targetID_mrca2target_ID_paths_dic.setdefault(targetID, current_path_IDs)

    for mrcaID in mID_targetID_dic.keys():
        mrca_seq = sID_nSeq_dic.get(mrcaID)
        if not mrca_seq or site_idx + 3 > len(mrca_seq):
            continue # Skip if MRCA not in seq or seq too short

        mrcaCodon = mrca_seq[site_idx : site_idx+3]
        mrcaAA    = codon2aa(mrcaCodon, codon_table)
        mrcaNuc1, mrcaNuc2, mrcaNuc3 = mrcaCodon[0], mrcaCodon[1], mrcaCodon[2]

        path_data['mrca_states']['codon'].append(mrcaCodon)
        path_data['mrca_states']['aa'].append(mrcaAA)
        path_data['mrca_states']['nuc1'].append(mrcaNuc1)
        path_data['mrca_states']['nuc2'].append(mrcaNuc2)
        path_data['mrca_states']['nuc3'].append(mrcaNuc3)
        
        AmrcaID = sID_mID_dic.get(mrcaID)
        if not AmrcaID:
            # This MRCA might be the root of the (sub)tree being analyzed, or an error.
            # For ConVarFinder, typically we expect an ancestor unless it's the ultimate root.
            # print(f"Warning: MRCA {mrcaID} has no parent (AmrcaID) in sID_mID_dic.")
            continue 
            
        Amrca_seq = sID_nSeq_dic.get(AmrcaID)
        if not Amrca_seq or site_idx + 3 > len(Amrca_seq):
            continue # Skip if AmrcaID invalid or seq too short

        AmrcaCodon = Amrca_seq[site_idx : site_idx+3]
        AmrcaAA    = codon2aa(AmrcaCodon, codon_table)
        AmrcaNuc1, AmrcaNuc2, AmrcaNuc3 = AmrcaCodon[0], AmrcaCodon[1], AmrcaCodon[2]

        path_data['Amrca_states']['codon'].append(AmrcaCodon)
        path_data['Amrca_states']['aa'].append(AmrcaAA)
        path_data['Amrca_states']['nuc1'].append(AmrcaNuc1)
        path_data['Amrca_states']['nuc2'].append(AmrcaNuc2)
        path_data['Amrca_states']['nuc3'].append(AmrcaNuc3)
        
        path_data['Amrca2mrca_IDs'].append(f"{AmrcaID}>{mrcaID}")
        path_data['Amrca2mrca_paths']['codon'].append(f"{AmrcaCodon}>{mrcaCodon}")
        path_data['Amrca2mrca_paths']['aa'].append(f"{AmrcaAA}>{mrcaAA}")
        path_data['Amrca2mrca_paths']['nuc1'].append(f"{AmrcaNuc1}>{mrcaNuc1}")
        path_data['Amrca2mrca_paths']['nuc2'].append(f"{AmrcaNuc2}>{mrcaNuc2}")
        path_data['Amrca2mrca_paths']['nuc3'].append(f"{AmrcaNuc3}>{mrcaNuc3}")
        
        current_mrca_targets = mID_targetID_dic[mrcaID]
        for targetID in current_mrca_targets:
            if targetID not in targetID_mrca2target_ID_paths_dic:
                continue 
            
            id_path_list = targetID_mrca2target_ID_paths_dic[targetID]
            path_data['mrca2target_IDs'].append(">".join(id_path_list))

            states_on_path = {'codon': [], 'aa': [], 'nuc1': [], 'nuc2': [], 'nuc3': []}
            valid_path_for_target = True
            for nodeID_on_path in id_path_list:
                node_seq = sID_nSeq_dic.get(nodeID_on_path)
                if not node_seq or site_idx + 3 > len(node_seq):
                    valid_path_for_target = False; break
                
                nodeCodon = node_seq[site_idx : site_idx+3]
                states_on_path['codon'].append(nodeCodon)
                states_on_path['aa'].append(codon2aa(nodeCodon, codon_table))
                states_on_path['nuc1'].append(nodeCodon[0])
                states_on_path['nuc2'].append(nodeCodon[1])
                states_on_path['nuc3'].append(nodeCodon[2])
            
            if not valid_path_for_target:
                # Append placeholders if path is invalid to maintain structure for join, or skip
                for key in path_data['mrca2target_paths']:
                    path_data['mrca2target_paths'][key].append("PATH_ERR") # Or some other indicator
                continue

            for key in states_on_path:
                path_data['mrca2target_paths'][key].append(">".join(states_on_path[key]))
                if len(set(states_on_path[key])) >= 2:
                    path_data['clade_substitutions'][key] += 1
    return path_data

def _determine_evolutionary_pattern(mrca_states_unique, amrca_states_unique, target_states_unique, has_clade_substitutions):
    """Determines evolutionary pattern (Amrca2mrca and Amrca2target) for a single level (codon, aa, or nuc)."""
    evo_dir_amrca2mrca = "UNDEF"
    evo_dir_amrca2target = "UNDEF"

    # Determine EvoDir_Amrca2mrca
    if not mrca_states_unique or not amrca_states_unique: # Not enough data to determine
        evo_dir_amrca2mrca = "MIXED_DATA_DEFICIENCY"
    elif len(mrca_states_unique) == 1:
        mrca_state = list(mrca_states_unique)[0]
        if len(amrca_states_unique) == 1:
            amrca_state = list(amrca_states_unique)[0]
            if mrca_state == amrca_state:
                evo_dir_amrca2mrca = "PLE"  # Plesiomorphic: AAA > AAA / AAA > AAA
            else:
                evo_dir_amrca2mrca = "PAR"  # Parallel:    AAA > TTT / AAA > TTT (all MRCAs changed to same state from same AMrca state)
        else: # Multiple Amrca states, one Mrca state
            if mrca_state in amrca_states_unique:
                evo_dir_amrca2mrca = "MIX(Con+Ple)" # Convergent to Mrca + Plesiomorphic from one Amrca
            else:
                evo_dir_amrca2mrca = "CON"  # Convergent:  AAA > TTT / GGG > TTT (all MRCAs changed to same state from different AMrca states)
    else: # Multiple Mrca states (len(mrca_states_unique) > 1)
        is_ple_component_present = any(m_state in amrca_states_unique for m_state in mrca_states_unique)
        if is_ple_component_present:
            evo_dir_amrca2mrca = "MIX(Div+Ple)" # Divergent from Amrca + Plesiomorphic for some
        else:
            evo_dir_amrca2mrca = "DIV"  # Divergent:   AAA > TTT / AAA > GGG (MRCAs differ, and all differ from AMrcas)

    # Determine EvoDir_Amrca2target based on Amrca2mrca and target uniformity + clade substitutions
    is_target_uniform = len(target_states_unique) == 1

    if is_target_uniform:
        if evo_dir_amrca2mrca.startswith("PLE"):
            evo_dir_amrca2target = "cCON" # Assumed convergent at target level if plesiomorphic at MRCA
        elif evo_dir_amrca2mrca.startswith("PAR"):
            evo_dir_amrca2target = "sPAR" if not has_clade_substitutions else "cCON"
        elif evo_dir_amrca2mrca.startswith("CON"):
            evo_dir_amrca2target = "sCON" if not has_clade_substitutions else "cCON"
        elif evo_dir_amrca2mrca.startswith("DIV"):
            # This case was an error in the original: identical targets from divergent MRCAs implies convergence at target level
            evo_dir_amrca2target = "cCON" 
        elif evo_dir_amrca2mrca.startswith("MIX"):
             evo_dir_amrca2target = "cCON" # If MRCAs are mixed but targets uniform, it implies convergence to that state
        else: # UNDEF or MIXED_DATA_DEFICIENCY
            evo_dir_amrca2target = "cCON_UNSURE_MRCA" if target_states_unique else "UNDEF_TARGET_STATE"
    else: # Targets are NOT uniform (diverse)
        if evo_dir_amrca2mrca.startswith("DIV"):
            evo_dir_amrca2target = "sDIV" if not has_clade_substitutions else "cDIV"
        elif evo_dir_amrca2mrca.startswith("PLE") or evo_dir_amrca2mrca.startswith("PAR") or evo_dir_amrca2mrca.startswith("CON"):
            # If MRCAs were uniform (or parallel/convergent to one state) but targets diverse, implies divergence from MRCAs
            evo_dir_amrca2target = "cDIV"
        elif evo_dir_amrca2mrca.startswith("MIX"):
            evo_dir_amrca2target = "cDIV"
        else: # UNDEF or MIXED_DATA_DEFICIENCY
            evo_dir_amrca2target = "cDIV_UNSURE_MRCA"

    return evo_dir_amrca2mrca, evo_dir_amrca2target


def ConVarFinder(fNAME_seq, terNode_list, mID_targetID_dic_gene, mID_mmID_dic_gene, sID_mID_dic_gene, local_target_list_gene, templete_sID_list_cols, fpout, outgroup_list, codon_table, sID_nSeq_dic_main):

  iLen_nSeq = 0
  if sID_nSeq_dic_main: 
      # Attempt to get a sequence to determine length; robustly handle if empty or first key has no seq
      try:
          first_sID = next(iter(sID_nSeq_dic_main))
          if first_sID in sID_nSeq_dic_main and sID_nSeq_dic_main[first_sID]:
              iLen_nSeq = len(sID_nSeq_dic_main[first_sID])
      except StopIteration: # sID_nSeq_dic_main is empty
          pass

  if iLen_nSeq == 0:
    # print(f"Warning: No sequences found or sequence length is 0 in {fNAME_seq} for analysis. Skipping.")
    return fpout 
  
  if not iLen_nSeq % 3 == 0:
    print(f"Error: Sequence length ({iLen_nSeq}) in {fNAME_seq} is not a multiple of 3. Skipping.")
    return fpout
  
  gene_name = os.path.basename(fNAME_seq)
  if '.' in gene_name: # Basic way to get name before first dot
      gene_name = gene_name.split('.')[0]

  iPos_Codon, iPos_AA, iPos_Nuc1, iPos_Nuc2, iPos_Nuc3 = -2, 0, -2, -1, 0

  for i in range(0, iLen_nSeq, 3):
    iPos_Codon, iPos_AA = iPos_Codon + 3, iPos_AA + 1
    iPos_Nuc1, iPos_Nuc2, iPos_Nuc3 = iPos_Nuc1 + 3, iPos_Nuc2 + 3, iPos_Nuc3 + 3
    
    site_data = _extract_site_specific_data(sID_nSeq_dic_main, terNode_list, local_target_list_gene, outgroup_list, i, codon_table)
    substitution_flags = _check_mutually_exclusive_substitutions(site_data)

    if not substitution_flags['codon']: # If codon substitution IS target-specific (exclusive)
      # Collect ancestral and path data for this specific site
      # mID_targetID_dic_gene and sID_mID_dic_gene are for the current gene tree
      ancestral_path_info = _collect_ancestral_path_data(mID_targetID_dic_gene, sID_mID_dic_gene, sID_nSeq_dic_main, i, codon_table)

      # Unique states for determining evolutionary patterns
      amrca_codons_unique = set(ancestral_path_info['Amrca_states']['codon'])
      mrca_codons_unique = set(ancestral_path_info['mrca_states']['codon'])
      target_codons_unique = set(site_data['target_codons'])
      
      evo_dir_amrca2mrca_codon, evo_dir_amrca2target_codon = _determine_evolutionary_pattern(
          mrca_codons_unique, amrca_codons_unique, target_codons_unique, 
          ancestral_path_info['clade_substitutions']['codon'] > 0
      )

      type_codon_sub_target = "i" if len(target_codons_unique) == 1 else "d"
      type_codon_sub_others = "i" if len(set(site_data['others_codons'])) == 1 else "d"
      if type_codon_sub_target == "i": type_codon_sub = "1" if type_codon_sub_others == "i" else "2"
      else: type_codon_sub = "3" if type_codon_sub_others == "i" else "4"

      # Initialize AA and Nuc result fields to empty or placeholder
      results_aa = {k: '' for k in ['pos', 'evo_amrca2mrca', 'evo_amrca2target', 'type_target', 'type_others', 'type_sub', 'path_amrca2mrca', 'path_mrca2target', 'states_target', 'states_others']}
      results_nuc1 = results_aa.copy()
      results_nuc2 = results_aa.copy()
      results_nuc3 = results_aa.copy()

      if not substitution_flags['aa']: # If AA substitution is also target-specific
          amrca_aas_unique = set(ancestral_path_info['Amrca_states']['aa'])
          mrca_aas_unique = set(ancestral_path_info['mrca_states']['aa'])
          target_aas_unique = set(site_data['target_aas'])
          evo_dir_amrca2mrca_aa, evo_dir_amrca2target_aa = _determine_evolutionary_pattern(
              mrca_aas_unique, amrca_aas_unique, target_aas_unique, 
              ancestral_path_info['clade_substitutions']['aa'] > 0
          )
          results_aa['pos'] = str(iPos_AA)
          results_aa['evo_amrca2mrca'] = evo_dir_amrca2mrca_aa
          results_aa['evo_amrca2target'] = evo_dir_amrca2target_aa
          results_aa['type_target'] = "i" if len(target_aas_unique) == 1 else "d"
          results_aa['type_others'] = "i" if len(set(site_data['others_aas'])) == 1 else "d"
          if results_aa['type_target'] == "i": results_aa['type_sub'] = "1" if results_aa['type_others'] == "i" else "2"
          else: results_aa['type_sub'] = "3" if results_aa['type_others'] == "i" else "4"
          results_aa['path_amrca2mrca'] = '/'.join(ancestral_path_info['Amrca2mrca_paths']['aa'])
          results_aa['path_mrca2target'] = '/'.join(ancestral_path_info['mrca2target_paths']['aa'])
          results_aa['states_target'] = ','.join(target_aas_unique)
          results_aa['states_others'] = ','.join(set(site_data['others_aas']))

      # Nuc1 processing
      if not substitution_flags['nuc1']:
          amrca_nuc1s_unique = set(ancestral_path_info['Amrca_states']['nuc1'])
          mrca_nuc1s_unique = set(ancestral_path_info['mrca_states']['nuc1'])
          target_nuc1s_unique = set(site_data['target_nuc1s'])
          evo_dir_amrca2mrca_nuc1, evo_dir_amrca2target_nuc1 = _determine_evolutionary_pattern(
              mrca_nuc1s_unique, amrca_nuc1s_unique, target_nuc1s_unique, 
              ancestral_path_info['clade_substitutions']['nuc1'] > 0
          )
          results_nuc1['pos'] = str(iPos_Nuc1)
          results_nuc1['evo_amrca2mrca'] = evo_dir_amrca2mrca_nuc1
          results_nuc1['evo_amrca2target'] = evo_dir_amrca2target_nuc1
          results_nuc1['type_target'] = "i" if len(target_nuc1s_unique) == 1 else "d"
          results_nuc1['type_others'] = "i" if len(set(site_data['others_nuc1s'])) == 1 else "d"
          if results_nuc1['type_target'] == "i": results_nuc1['type_sub'] = "1" if results_nuc1['type_others'] == "i" else "2"
          else: results_nuc1['type_sub'] = "3" if results_nuc1['type_others'] == "i" else "4"
          results_nuc1['path_amrca2mrca'] = '/'.join(ancestral_path_info['Amrca2mrca_paths']['nuc1'])
          results_nuc1['path_mrca2target'] = '/'.join(ancestral_path_info['mrca2target_paths']['nuc1'])
          results_nuc1['states_target'] = ','.join(target_nuc1s_unique)
          results_nuc1['states_others'] = ','.join(set(site_data['others_nuc1s']))

      # Nuc2 processing
      if not substitution_flags['nuc2']:
          amrca_nuc2s_unique = set(ancestral_path_info['Amrca_states']['nuc2'])
          mrca_nuc2s_unique = set(ancestral_path_info['mrca_states']['nuc2'])
          target_nuc2s_unique = set(site_data['target_nuc2s'])
          evo_dir_amrca2mrca_nuc2, evo_dir_amrca2target_nuc2 = _determine_evolutionary_pattern(
              mrca_nuc2s_unique, amrca_nuc2s_unique, target_nuc2s_unique, 
              ancestral_path_info['clade_substitutions']['nuc2'] > 0
          )
          results_nuc2['pos'] = str(iPos_Nuc2)
          results_nuc2['evo_amrca2mrca'] = evo_dir_amrca2mrca_nuc2
          results_nuc2['evo_amrca2target'] = evo_dir_amrca2target_nuc2
          results_nuc2['type_target'] = "i" if len(target_nuc2s_unique) == 1 else "d"
          results_nuc2['type_others'] = "i" if len(set(site_data['others_nuc2s'])) == 1 else "d"
          if results_nuc2['type_target'] == "i": results_nuc2['type_sub'] = "1" if results_nuc2['type_others'] == "i" else "2"
          else: results_nuc2['type_sub'] = "3" if results_nuc2['type_others'] == "i" else "4"
          results_nuc2['path_amrca2mrca'] = '/'.join(ancestral_path_info['Amrca2mrca_paths']['nuc2'])
          results_nuc2['path_mrca2target'] = '/'.join(ancestral_path_info['mrca2target_paths']['nuc2'])
          results_nuc2['states_target'] = ','.join(target_nuc2s_unique)
          results_nuc2['states_others'] = ','.join(set(site_data['others_nuc2s']))

      # Nuc3 processing
      if not substitution_flags['nuc3']:
          amrca_nuc3s_unique = set(ancestral_path_info['Amrca_states']['nuc3'])
          mrca_nuc3s_unique = set(ancestral_path_info['mrca_states']['nuc3'])
          target_nuc3s_unique = set(site_data['target_nuc3s'])
          evo_dir_amrca2mrca_nuc3, evo_dir_amrca2target_nuc3 = _determine_evolutionary_pattern(
              mrca_nuc3s_unique, amrca_nuc3s_unique, target_nuc3s_unique, 
              ancestral_path_info['clade_substitutions']['nuc3'] > 0
          )
          results_nuc3['pos'] = str(iPos_Nuc3)
          results_nuc3['evo_amrca2mrca'] = evo_dir_amrca2mrca_nuc3
          results_nuc3['evo_amrca2target'] = evo_dir_amrca2target_nuc3
          results_nuc3['type_target'] = "i" if len(target_nuc3s_unique) == 1 else "d"
          results_nuc3['type_others'] = "i" if len(set(site_data['others_nuc3s'])) == 1 else "d"
          if results_nuc3['type_target'] == "i": results_nuc3['type_sub'] = "1" if results_nuc3['type_others'] == "i" else "2"
          else: results_nuc3['type_sub'] = "3" if results_nuc3['type_others'] == "i" else "4"
          results_nuc3['path_amrca2mrca'] = '/'.join(ancestral_path_info['Amrca2mrca_paths']['nuc3'])
          results_nuc3['path_mrca2target'] = '/'.join(ancestral_path_info['mrca2target_paths']['nuc3'])
          results_nuc3['states_target'] = ','.join(target_nuc3s_unique)
          results_nuc3['states_others'] = ','.join(set(site_data['others_nuc3s']))

      # Prepare site-wise codon and AA sequence info for all species in templete_sID_list_cols
      site_codon_line_parts = []
      site_aa_line_parts = []
      for sID_col in templete_sID_list_cols:
        seq_for_sp = sID_nSeq_dic_main.get(sID_col)
        if seq_for_sp and i + 3 <= len(seq_for_sp):
          codon_val = seq_for_sp[i:i+3]
          aa_val = codon2aa(codon_val, codon_table)
          site_codon_line_parts.append(codon_val)
          site_aa_line_parts.append(aa_val)
        else:
          site_codon_line_parts.append('---')
          site_aa_line_parts.append('-')
      site_codon_line_str = '\t'.join(site_codon_line_parts)
      site_aa_line_str    = '\t'.join(site_aa_line_parts)

      output_values = [
          gene_name, 
          '/'.join(ancestral_path_info['Amrca2mrca_IDs']), 
          '/'.join(ancestral_path_info['mrca2target_IDs']),
          str(iPos_Codon), evo_dir_amrca2mrca_codon, evo_dir_amrca2target_codon, 
          type_codon_sub_target, type_codon_sub_others, type_codon_sub, 
          '/'.join(ancestral_path_info['Amrca2mrca_paths']['codon']), 
          '/'.join(ancestral_path_info['mrca2target_paths']['codon']), 
          ','.join(target_codons_unique), ','.join(set(site_data['others_codons']))
      ]
      for res_level in [results_aa, results_nuc1, results_nuc2, results_nuc3]:
          output_values.extend([
              res_level['pos'], res_level['evo_amrca2mrca'], res_level['evo_amrca2target'],
              res_level['type_target'], res_level['type_others'], res_level['type_sub'],
              res_level['path_amrca2mrca'], res_level['path_mrca2target'],
              res_level['states_target'], res_level['states_others']
          ])
      
      tmpline = '\t'.join(map(str, output_values))
      if templete_sID_list_cols: # Add species data only if there are species columns
          tmpline += '\t' + site_codon_line_str + '\t' + site_aa_line_str
      tmpline += '\n'
      fpout.write(tmpline)
  return fpout

def main(options, nGene_from_filename, fpout, codon_table_main):
  """Main processing logic for a single gene, using options dictionary."""
  # Extract necessary info from the options dictionary
  nTargets                  = options['nTargets']
  templete_mID_targetID_dic = options['templete_mID_targetID_dic']
  templete_mID_mmID_dic     = options['templete_mID_mmID_dic']
  templete_sID_list_for_analysis = options['templete_sID_list_for_analysis'] # Species for output columns
  outgroup_list             = options['outgroup_list'] # Parsed list of outgroup species
  seq_file                  = options['seq_file']
  tree_file                 = options['tree_file']
  sID_nSeq_dic_main         = options['sID_nSeq_dic'] # Sequence data, read once

  # The nGene_from_filename is passed from the __main__ block
  nGene = nGene_from_filename 

  # Current gene's tree processing (original logic used fNAME_tree which was derived from iPATH and gene name)
  # Now, tree_file is the direct path to the tree provided by the user.
  # This script assumes one tree file is used for all genes if multiple seq files were intended to be processed
  # in a loop (which is not the current structure). If each gene has its own tree, fNAME_tree logic would need to be per-gene.
  # For now, we use the globally provided tree_file for this specific gene's analysis.
  nGene_nTree = read_tree_file(tree_file) # This is the user-provided tree file
  nGene_nTree_no_branch = read_tre_removing_branchlength(nGene_nTree)
  nGene_terNode_list = terminal_node_from_tree(nGene_nTree_no_branch) # Terminal nodes from this specific tree
  
  # Parse target species specific to THIS gene tree (if they exist in it)
  # nTargets is the list of desired targets from user input/prompt against the MASTER tree.
  # We need to see which of these are present in the nGene_terNode_list.
  nGene_target_list = [t for t in nTargets.split(',') if t in nGene_terNode_list and t] if nTargets else []

  if not nGene_target_list:
    # print(f"Warning: None of the specified target species ({nTargets}) are found in the tree for gene {nGene}. Skipping this gene.")
    return fpout # Skip if no target species for this gene are in its tree

  # Check clade consistency - original logic
  # This compares mIDs from the master tree (templete_mID_targetID_dic) with targets found in the current gene tree.
  iCNT_clade = 0
  if templete_mID_targetID_dic: # Ensure it's not empty
    for mID in templete_mID_targetID_dic.keys():
      iFlag_clade = 0
      for targetID_in_mID_clade in templete_mID_targetID_dic[mID]:
        if targetID_in_mID_clade in nGene_target_list: # Check against targets present in current gene tree
          iFlag_clade += 1
      if iFlag_clade > 0:
        iCNT_clade += 1
    if not iCNT_clade == len(templete_mID_targetID_dic.keys()):
      # print(f"Warning: Clade consistency check failed for gene {nGene}. Skipping.")
      return fpout
  elif nTargets: # If there are targets, but templete_mID_targetID_dic is empty (e.g. all targets are singletons)
      pass # Allow to proceed, Anc_Finder for nGene_nTree will handle it.
  else: # No targets specified, no specific clades to check from master tree processing
      pass 

  # Phylogenetic relationships for the current gene tree and its specific targets
  # Use nGene_target_list which contains only targets present in the current nGene_nTree
  current_gene_mID_targetID_dic, current_gene_mID_mmID_dic = Anc_Finder(nGene_nTree_no_branch, ','.join(nGene_target_list))
  
  # We need sID_mID_dic for the current gene tree to pass to ConVarFinder
  current_gene_mID_sID_dic = {}
  current_gene_mID_sID_dic = make_mID_sID_dic(nGene_nTree_no_branch, current_gene_mID_sID_dic)
  current_gene_sID_mID_dic = make_sID_mID_dic(current_gene_mID_sID_dic)

  # Call Convergent variant finder
  fpout = ConVarFinder(
      seq_file,      # Use the original seq_file path for fNAME_seq in ConVarFinder
      nGene_terNode_list, 
      current_gene_mID_targetID_dic, 
      current_gene_mID_mmID_dic, 
      current_gene_sID_mID_dic, 
      nGene_target_list, 
      options['templete_sID_list_for_analysis'], # Use the correctly scoped list for columns
      fpout, 
      options['outgroup_list'], 
      codon_table_main, 
      options['sID_nSeq_dic'] # Pass the main sequence dictionary
  )
  return fpout


if __name__ == "__main__":
  # Global constants are defined at the top of the file.
  # CODONTABLE_dic is passed to main, then to ConVarFinder.
  
  current_options = set_options() # Returns a dictionary
  
  nTargets          = current_options['nTargets']
  oPATH             = current_options['oPATH']
  header_sID_list   = current_options['templete_sID_list_for_analysis'] 
  seq_file          = current_options['seq_file']
  tree_file         = current_options['tree_file']
  fNAME_seq_format  = current_options['fNAME_seq_format']
  outgroup_list_final = current_options['outgroup_list'] # Parsed outgroup list
  
  safe_targets = nTargets.replace(",", "-") if nTargets else "all_selected_targets"
  oNAME = os.path.join(oPATH, f"ConVarFinder_{safe_targets}.txt") 
  
  if not os.path.isdir(oPATH):
    try:
      os.makedirs(oPATH, exist_ok=True)
      print(f"Output directory created: {oPATH}")
    except OSError as e:
      print(f"Error creating output directory {oPATH}: {e}")
      sys.exit(1)

  print(f"Target species: {nTargets if nTargets else 'Prompted or all available if not specified'}")
  if outgroup_list_final:
      print(f"Outgroup species: {','.join(outgroup_list_final)}")
  else:
      print("No outgroup species selected or specified.")
  print(f"Output file path: {oNAME}")

  try:
    with open(oNAME, 'w') as fpout:
      nHeader_parts = [
          "Gene", "Phy.Rel._Amrca2mrca", "Phy.Rel._mrca2target",
          "Pos_codon", "EvoDir_Amrca2mrca_Codon", "EvoDir_Amrca2target_Codon", 'Type_CodonSub_target', 'Type_CodonSub_others', 'Type_CodonSub', "Codon_Amrca2mrca", "Codon_mrca2target", "Codon_targets", "Codon_others",
          "Pos_AA", "EvoDir_Amrca2mrca_AA", "EvoDir_Amrca2target_AA", 'Type_AASub_target', 'Type_AASub_others', 'Type_AASub', "AA_Amrca2mrca", "AA_mrca2target", "AA_targets", "AA_others",
          "Pos_Nuc1", "EvoDir_Amrca2mrca_Nuc1", "EvoDir_Amrca2target_Nuc1", 'Type_Nuc1Sub_target', 'Type_Nuc1Sub_others', 'Type_Nuc1Sub', "Nuc1_Amrca2mrca", "Nuc1_mrca2target", "Nuc1_targets", "Nuc1_others",
          "Pos_Nuc2", "EvoDir_Amrca2mrca_Nuc2", "EvoDir_Amrca2target_Nuc2", 'Type_Nuc2Sub_target', 'Type_Nuc2Sub_others', 'Type_Nuc2Sub', "Nuc2_Amrca2mrca", "Nuc2_mrca2target", "Nuc2_targets", "Nuc2_others",
          "Pos_Nuc3", "EvoDir_Amrca2mrca_Nuc3", "EvoDir_Amrca2target_Nuc3", 'Type_Nuc3Sub_target', 'Type_Nuc3Sub_others', 'Type_Nuc3Sub', "Nuc3_Amrca2mrca", "Nuc3_mrca2target", "Nuc3_targets", "Nuc3_others"
      ]
      if header_sID_list: # Only add species columns if list is not empty
        nHeader_parts.extend(header_sID_list) # Codon columns
        nHeader_parts.extend(header_sID_list) # AA columns
      fpout.write('\t'.join(nHeader_parts) + '\n')

      nGene_from_filename = os.path.basename(seq_file)
      if fNAME_seq_format and nGene_from_filename.endswith(fNAME_seq_format):
          nGene_from_filename = nGene_from_filename[:-len(fNAME_seq_format)]
      
      print(f"Analysis started: file {seq_file}, {tree_file}, nGene={nGene_from_filename}")
      main(current_options, nGene_from_filename, fpout, CODONTABLE_dic) # Pass CODONTABLE_dic to main
      
      print(f"Analysis completed: Results saved to {oNAME}")

  except IOError as e:
    print(f"Error writing to output file {oNAME}: {e}")
    sys.exit(1)
  except Exception as e: 
    print(f"An unexpected error occurred: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
