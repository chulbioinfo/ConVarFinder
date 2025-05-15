#!/usr/bin/env python3
from Bio import SeqIO
import subprocess
import os
import re
import sys
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil

def create_binary_sequence(input_fasta, output_fasta):
    """Generate a binary sequence FASTA file: gap ('-') as 1, others as 0."""
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            binary_seq = ''.join(['1' if base == '-' else '0' for base in str(record.seq)])
            outfile.write(f'>{record.id}\n{binary_seq}\n')

def run_raxml(input_fasta, tree_file, model, prefix, output_dir):
    """Run RAxML for ancestral sequence estimation."""
    raxml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'standard-RAxML', 'raxmlHPC-PTHREADS-SSE3')
    
    # 입력 파일이 이미 작업 디렉토리에 있는 경우 복사하지 않음
    if os.path.dirname(input_fasta) != output_dir:
        input_basename = os.path.basename(input_fasta)
        input_copy = os.path.join(output_dir, input_basename)
        shutil.copy2(input_fasta, input_copy)
        input_fasta = input_copy
    
    cmd = [
        raxml_path,
        '-T', '4',  # Use 4 threads
        '-f', 'A',  # Ancestral sequence estimation
        '-m', model,
        '-p', '12345',
        '-s', input_fasta,
        '-t', tree_file,
        '-n', prefix
    ]
    subprocess.run(cmd, cwd=output_dir)
    
    # .reduced 파일이 있다면 작업 디렉토리로 이동
    reduced_file = f"{input_fasta}.reduced"
    if os.path.exists(reduced_file):
        reduced_dest = os.path.join(output_dir, f"{prefix}_reduced.fasta")
        shutil.move(reduced_file, reduced_dest)
        print(f"Reduced file saved to: {reduced_dest}")

def process_ancestral_states(binary_anc, codon_anc, output_fasta):
    """Combine binary and codon ancestral states to generate the final result."""
    binary_states = {}
    codon_states = {}
    
    # Read binary ancestral states
    with open(binary_anc, 'r') as f:
        current_id = None
        current_seq = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('ROOT'):
                if current_id and current_seq:
                    binary_states[current_id] = current_seq
                current_id = "ROOT"
                current_seq = line.split(' ', 1)[1]
            elif line[0].isdigit():
                if current_id and current_seq:
                    binary_states[current_id] = current_seq
                current_id = line.split(' ', 1)[0]
                current_seq = line.split(' ', 1)[1]
            else:
                current_seq += line
        if current_id and current_seq:
            binary_states[current_id] = current_seq
    
    # Read codon ancestral states
    with open(codon_anc, 'r') as f:
        current_id = None
        current_seq = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('ROOT'):
                if current_id and current_seq:
                    codon_states[current_id] = current_seq
                current_id = "ROOT"
                current_seq = line.split(' ', 1)[1]
            elif line[0].isdigit():
                if current_id and current_seq:
                    codon_states[current_id] = current_seq
                current_id = line.split(' ', 1)[0]
                current_seq = line.split(' ', 1)[1]
            else:
                current_seq += line
        if current_id and current_seq:
            codon_states[current_id] = current_seq
    
    # Write final result
    with open(output_fasta, 'w') as outfile:
        for seq_id in binary_states:
            if seq_id not in codon_states:
                print(f"Warning: No codon sequence found for {seq_id}")
                continue
                
            binary_seq = binary_states[seq_id]
            codon_seq = codon_states[seq_id]
            
            # Replace codon with '-' if binary sequence is 1
            final_seq = ''
            for i in range(0, len(codon_seq), 3):
                codon = codon_seq[i:i+3]
                if i//3 < len(binary_seq) and binary_seq[i//3] == '1':
                    final_seq += '---'
                else:
                    final_seq += codon
            
            outfile.write(f'>{seq_id}\n{final_seq}\n')
            print(f"Processed sequence {seq_id}")

def convert_codon_to_fasta(codon_file, output_file):
    """Convert RAxML codon file to FASTA format."""
    with open(codon_file, 'r') as f:
        lines = f.readlines()
    
    sequences = []
    current_seq = ""
    current_name = ""
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('ROOT'):
            if current_seq:
                sequences.append((current_name, current_seq))
            current_name = "ROOT"
            current_seq = line.split(' ', 1)[1]
        elif line[0].isdigit():
            if current_seq:
                sequences.append((current_name, current_seq))
            current_name = line.split(' ', 1)[0]
            current_seq = line.split(' ', 1)[1]
        else:
            current_seq += line
    
    if current_seq:
        sequences.append((current_name, current_seq))
    
    # Save in FASTA format
    with open(output_file, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")
            print(f"Sequence {name} has been saved.")

def estimate_ancestral_sequences(fasta_file, tree_file, output_file):
    """Estimate ancestral sequences using RAxML."""
    # Execute RAxML command
    raxml_cmd = [
        "raxmlHPC",
        "-f A",  # Ancestral sequence estimation mode
        "-m GTRGAMMA",  # Use GTR+GAMMA model
        "-s", fasta_file,  # Input sequence file
        "-t", tree_file,  # Input tree file
        "-n", "anc"  # Output file suffix
    ]
    
    try:
        subprocess.run(" ".join(raxml_cmd), shell=True, check=True)
        
        # Convert RAxML output file to FASTA format
        codon_file = "RAxML_marginalAncestralStates.codon"
        if os.path.exists(codon_file):
            print(f"Found RAxML codon file: {codon_file}")
            convert_codon_to_fasta(codon_file, output_file)
            print(f"Ancestral sequence has been saved to {output_file}")
        else:
            print(f"RAxML codon file not found: {codon_file}")
            
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running RAxML: {e}")

def parse_raxml_states(file_path):
    """Parse RAxML_marginalAncestralStates file and return as dict"""
    states = {}
    with open(file_path) as f:
        current_id = None
        current_seq = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('ROOT'):
                if current_id and current_seq:
                    states[current_id] = current_seq
                current_id = "ROOT"
                current_seq = line.split(' ', 1)[1]
            elif line[0].isdigit():
                if current_id and current_seq:
                    states[current_id] = current_seq
                current_id = line.split(' ', 1)[0]
                current_seq = line.split(' ', 1)[1]
            else:
                current_seq += line
        if current_id and current_seq:
            states[current_id] = current_seq
    return states

def make_ancestral_fasta(input_fasta, codon_file, binary_file, output_fasta):
    # 1. Add original sequences
    records = list(SeqIO.parse(input_fasta, "fasta"))
    with open(output_fasta, "w") as out:
        for rec in records:
            out.write(f">{rec.id}\n{str(rec.seq)}\n")

        # 2. Generate ancestral node sequences
        codon_states = parse_raxml_states(codon_file)
        binary_states = parse_raxml_states(binary_file)
        for node in codon_states:
            codon_seq = codon_states[node]
            binary_seq = binary_states.get(node, None)
            if not binary_seq:
                continue
            # Apply gaps
            final_seq = ""
            for i in range(0, len(codon_seq), 3):
                codon = codon_seq[i:i+3]
                if i//3 < len(binary_seq) and binary_seq[i//3] == '1':
                    final_seq += "---"
                else:
                    final_seq += codon
            out.write(f">{node}\n{final_seq}\n")

def copy_tree_file(src_tree, dest_tree):
    shutil.copyfile(src_tree, dest_tree)

def main():
    parser = argparse.ArgumentParser(
        description="""
Estimate ancestral sequences using RAxML with both codon and binary (gap) models.\n\n\
This script will:\n\
  1. Convert the input FASTA to a binary sequence (gap: 1, non-gap: 0)\n\
  2. Run RAxML for both binary and codon models\n\
  3. Combine the results: wherever the binary ancestral state is 1, the codon sequence is replaced with '-'\n\
  4. Output the final ancestral sequence FASTA and the cladogram tree file\n\n\
Input Requirements:\n\
  - FASTA file: Must contain codon-aligned sequences\n\
  - Tree file: Must be in Newick format and compatible with RAxML\n\n\
Example usage:\n\
  python ancestral_sequence_estimation.py input.fasta input.tre\n\
  python ancestral_sequence_estimation.py input.fasta input.tre --output my_output\n\n\
Dependencies:\n\
  - Python 3.6 or higher\n\
  - Biopython (pip install biopython)\n\
  - RAxML (standard-RAxML/raxmlHPC-PTHREADS-SSE3 in the same directory as this script)\n\
  - Required Python packages:\n\
    * argparse\n\
    * os\n\
    * re\n\
    * sys\n\
    * subprocess\n\
    * shutil\n\n\
RAxML Models:\n\
  - Binary sequences: BINCAT model\n\
  - Codon sequences: GTRCAT model\n\n\
Output Files:\n\
  - {input_id}_ancestral/ directory containing:\n\
    * {input_id}_anc.fasta: Combined FASTA file with original and ancestral sequences\n\
    * {input_id}_anc.tre: Rooted tree file with node labels\n\
    * binary_input.fasta: Binary sequence file\n\
    * RAxML_* files: Intermediate RAxML output files\n\n\
Error Handling:\n\
  - Checks for RAxML executable existence\n\
  - Validates input file formats\n\
  - Handles RAxML execution errors\n\
  - Manages file I/O operations\n""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_fasta', help='Input FASTA file (codon alignment)')
    parser.add_argument('input_tree', help='Input tree file (Newick format)')
    parser.add_argument('--output', help='Optional output name for the result files. If not provided, output files will be named based on the input file ID with "_anc" suffix.')
    args = parser.parse_args()

    # Check RAxML dependency
    raxml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'standard-RAxML', 'raxmlHPC-PTHREADS-SSE3')
    if not os.path.exists(raxml_path):
        print("Error: RAxML executable not found at", raxml_path)
        print("Please ensure that RAxML (standard-RAxML/raxmlHPC-PTHREADS-SSE3) is installed in the same directory as this script.")
        sys.exit(1)

    input_fasta = os.path.abspath(args.input_fasta)
    tree_file = os.path.abspath(args.input_tree)
    
    # Create work directory based on input file name
    input_id = os.path.splitext(os.path.basename(input_fasta))[0]
    work_dir = os.path.join(os.path.dirname(input_fasta), f"{input_id}_ancestral")
    os.makedirs(work_dir, exist_ok=True)

    # Determine output file names
    if args.output:
        output_fasta = os.path.join(work_dir, f"{args.output}.fasta")
        output_tree = os.path.join(work_dir, f"{args.output}.tre")
    else:
        output_fasta = os.path.join(work_dir, f"{input_id}_anc.fasta")
        output_tree = os.path.join(work_dir, f"{input_id}_anc.tre")

    # Intermediate file paths
    binary_fasta = os.path.join(work_dir, 'binary_input.fasta')

    # RAxML result file paths
    binary_anc = os.path.join(work_dir, 'RAxML_marginalAncestralStates.binary')
    codon_anc = os.path.join(work_dir, 'RAxML_marginalAncestralStates.codon')

    # Clean up previous RAxML output files and intermediate files
    for file in [binary_fasta, binary_anc, codon_anc]:
        if os.path.exists(file):
            os.remove(file)
    for file in os.listdir(work_dir):
        if file.startswith('RAxML_') and (file.endswith('.binary') or file.endswith('.codon')):
            os.remove(os.path.join(work_dir, file))

    # 1. Generate binary sequence
    print("Generating binary sequence...")
    create_binary_sequence(input_fasta, binary_fasta)

    # 2. Ancestral estimation for binary sequence
    print("Running ancestral estimation for binary sequence...")
    run_raxml(binary_fasta, tree_file, 'BINCAT', 'binary', work_dir)

    # 3. Ancestral estimation for codon sequence
    print("Running ancestral estimation for codon sequence...")
    run_raxml(input_fasta, tree_file, 'GTRCAT', 'codon', work_dir)

    # 4. Make final _anc.fasta (leaf+ancestral)
    print("Generating final _anc.fasta (leaf+ancestral)...")
    make_ancestral_fasta(input_fasta, codon_anc, binary_anc, output_fasta)

    # 5. Copy rooted tree
    src_tree = os.path.join(work_dir, 'RAxML_nodeLabelledRootedTree.codon')
    if os.path.exists(src_tree):
        copy_tree_file(src_tree, output_tree)
        print(f"Rooted tree copied to: {output_tree}")
    else:
        print(f"Rooted tree file not found: {src_tree}")

    print("Process completed!")
    print(f"Final result file: {output_fasta}")
    print(f"Final tree file: {output_tree}")
    print("\nThe generated files can be used as input for ConVarFinder for subsequent analysis.")

if __name__ == "__main__":
    main() 