#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert muscle output format and provide alignment statistics
"""

import sys
import re
import argparse
from collections import defaultdict
from math import gcd
from functools import reduce

# Degenerate base codes
CODE = {
    '-': 0,
    '.': 0,
    'A': 2,
    'T': 3,
    'C': 5,
    'G': 7,
    'R': 14,  # A|G
    'Y': 15,  # C|T
    'M': 20,  # A|C
    'K': 21,  # G|T
    'S': 35,  # G|C
    'W': 24,  # A|T
    'H': 30,  # A|T|C
    'B': 105, # G|T|C
    'V': 70,  # G|A|C
    'D': 42,  # G|A|T
    'N': 210  # A|G|C|T
}

def lcm(a, b):
    return a * b // gcd(a, b)

def judge(bases):
    """Determine consensus symbol for alignment column"""
    if not bases:
        return ' '
        
    # Remove gaps
    bases = [b for b in bases if b not in ['-', '.']]
    if not bases:
        return '-'
    
    # Check if all bases are the same
    first_base = bases[0].upper()
    if all(b.upper() == first_base for b in bases):
        return '*' if first_base in ['A','T','C','G'] else ' '
    
    # Calculate product of codes for all bases
    product = reduce(lambda x, y: x * y, [CODE.get(b.upper(), 1) for b in bases])
    
    # Determine consensus symbol based on product
    if product == 0:
        return ' '
    elif product % 210 == 0:
        return 'N'
    elif product % 42 == 0:
        return 'D'
    elif product % 70 == 0:
        return 'V'
    elif product % 105 == 0:
        return 'B'
    elif product % 30 == 0:
        return 'H'
    elif product % 24 == 0:
        return 'W'
    elif product % 35 == 0:
        return 'S'
    elif product % 21 == 0:
        return 'K'
    elif product % 20 == 0:
        return 'M'
    elif product % 15 == 0:
        return 'Y'
    elif product % 14 == 0:
        return 'R'
    elif product % 7 == 0:
        return 'G'
    elif product % 5 == 0:
        return 'C'
    elif product % 3 == 0:
        return 'T'
    elif product % 2 == 0:
        return 'A'
    else:
        return ' '

def trim_sequence(seq, trim_len):
    """Trim sequence from both ends"""
    return seq[trim_len:-trim_len] if trim_len > 0 else seq

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_file', help='Input file in CLUSTALW or FASTA format')
    parser.add_argument('--trim', action='store_true', help='Whether to trim sequences')
    parser.add_argument('--del', type=int, default=50, help='Number of bases to trim')
    parser.add_argument('--tag', help='Sequence ID to use as trimming standard')
    parser.add_argument('--primer', action='store_true', help='Include primer information')
    parser.add_argument('--sanger', action='store_true', help='Include degenerate base alignment info')
    parser.add_argument('-o', '--out', default=None, help='Output file name')
    
    args = parser.parse_args()
    
    if args.out is None:
        args.out = f"{args.input_file}.link"
    
    # Read input file and determine format
    with open(args.input_file) as f:
        first_line = f.readline().strip()
        
    if first_line.startswith('>'):
        file_format = 'fasta'
    elif first_line.startswith('MUSCLE'):
        file_format = 'clustal'
    else:
        sys.stderr.write("Cannot recognize input file format!\n")
        sys.exit(1)
    
    # Process sequences
    sequences = defaultdict(str)
    id_lengths = {}
    max_id_len = 0
    
    with open(args.input_file) as f:
        if file_format == 'fasta':
            current_id = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]
                    id_lengths[current_id] = len(current_id)
                    max_id_len = max(max_id_len, len(current_id))
                elif current_id:
                    sequences[current_id] += line.upper()
        else:
            for line in f:
                line = line.strip()
                if line and not line.startswith('MUSCLE') and not line.startswith('*'):
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq = parts[-1]
                        id_lengths[seq_id] = len(seq_id)
                        max_id_len = max(max_id_len, len(seq_id))
                        sequences[seq_id] += seq.upper()
    
    # Generate consensus and statistics
    seq_list = list(sequences.values())
    length = len(seq_list[0]) if seq_list else 0
    
    with open(args.out, 'w') as out:
        # Write sequences with proper spacing
        for seq_id, seq in sequences.items():
            spacing = ' ' * (max_id_len - len(seq_id) + 8)
            out.write(f"{seq_id}{spacing}{seq}\n")
        
        # Write consensus line
        out.write(">>consensus" + ' ' * (max_id_len + 8 - 11))
        
        # Calculate consensus symbols
        consensus = []
        for i in range(length):
            column_bases = [seq[i] for seq in seq_list]
            consensus.append(judge(column_bases))
        out.write(''.join(consensus) + '\n')
        
        # Calculate statistics
        consensus_str = ''.join(consensus)
        # In the main function, modify the statistics calculation:
        mismatch = consensus_str.count('|')
        conserved = consensus_str.count('*')
        gap_whole = consensus_str.count('-')
        rate = conserved / (length - gap_whole) if (length - gap_whole) > 0 else 0
        
        # Output statistics
        out.write("\n----BASIC INFO----\n")
        out.write(f"whole_alignment_length\t\t{length}\n")
        out.write(f"Conserved\t{conserved}\n")
        out.write(f"gap_whole\t{gap_whole}\n")
        out.write(f"mismatch_whole\t{mismatch}\n")
        out.write("conserved_rate_whole:")
        out.write(f"{rate:.3f}\n")
        
        # Additional trimmed statistics if requested
        if args.trim:
            if args.tag:
                standard_seq = None
                for seq_id in sequences:
                    if args.tag in seq_id:
                        standard_seq = sequences[seq_id]
                        break
                
                if standard_seq:
                    # Implement trimming logic here
                    pass

if __name__ == "__main__":
    main()