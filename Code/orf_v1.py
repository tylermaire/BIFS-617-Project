# Tyler Maire
# --- Get user input ---
print("Enter FASTA file path (or '-' for manual input):")
fasta = input("File path: ").strip().strip('"').strip("'")
minlen = input("Minimum ORF length (default 50): ").strip()
minlen = int(minlen) if minlen else 50

def parse_fasta(lines):
    """Parse FASTA input into (header, sequence) pairs."""
    seqs = []
    header = None
    seq = ''
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header:
                seqs.append((header, seq))
            header = line[1:]
            seq = ''
        else:
            seq += ''.join(line.split())
    if header:
        seqs.append((header, seq))
    return seqs

def revcomp(seq):
    """Return reverse complement of DNA sequence."""
    comp = ''
    for b in seq:
        if b == 'A': comp += 'T'
        elif b == 'T': comp += 'A'
        elif b == 'C': comp += 'G'
        elif b == 'G': comp += 'C'
        elif b == 'a': comp += 't'
        elif b == 't': comp += 'a'
        elif b == 'c': comp += 'g'
        elif b == 'g': comp += 'c'
        else: comp += b
    return comp[::-1]

# Nick Lawson
# 13JUL2025
# ORF Finder - Version 1

# List of valid stop codons
STOP_CODONS = ['TAA', 'TAG', 'TGA']

# Read the input FASTA file path that was entered earlier
with open(fasta, 'r') as f:
    lines = f.readlines()

# Parse the FASTA records using the existing parse_fasta function
records = parse_fasta(lines)

# Create six reading frames from a given sequence
def get_frames(seq):
    seq = seq.upper()
    rc = revcomp(seq)  # Get the reverse complement
    return [
        seq,       # Frame 1
        seq[1:],   # Frame 2
        seq[2:],   # Frame 3
        rc,        # Frame 4 (reverse)
        rc[1:],    # Frame 5 (reverse)
        rc[2:]     # Frame 6 (reverse)
    ]

# Find all ORFs in a single reading frame
def find_orfs(seq, frame_num, is_rev):
    i = 0
    results = []
    while i <= len(seq) - 3:
        if seq[i:i+3] == 'ATG':  # Look for start codon
            j = i + 3
            while j <= len(seq) - 3:
                if seq[j:j+3] in STOP_CODONS:  # Look for stop codon
                    orf_seq = seq[i:j+3]
                    if len(orf_seq) >= minlen:  # Check minimum length
                        pos = -(len(seq)-i) if is_rev else i+1
                        results.append((frame_num, pos, len(orf_seq), orf_seq))
                    break
                j += 3
            i = j
        else:
            i += 3
    return results

# Process all parsed sequences to find ORFs in all frames
all_orfs = []
for header, sequence in records:
    frames = get_frames(sequence)
    for n in range(6):
        rev = n >= 3  # Reverse frames are 4, 5, 6
        found = find_orfs(frames[n], n+1, rev)
        for item in found:
            all_orfs.append((header, *item))

#Naziha James
#defining function to format the sequence into multiple lines
def format_fasta_sequence(seq, codons_per_line=15):
    formatted_lines = []#empty list to store the lines
    for i in range(0, len(seq), codons_per_line *3):#loops through the sequence based on the number of codons
        chunk = seq[i:i + codons_per_line * 3]#adding a space and *3 accounts for the three nucleotides per codon 
        codon_spaced = ' '.join([chunk[j:j+3] for j in range(0, len(chunk), 3)])#appending formatted chunk
        formatted_lines.append(codon_spaced)#joining formatted lines together
    return '\n'.join(formatted_lines)
    
if all_orfs:#checking if all_orfs is empty
    print("\nFound ORFs:")#if ORFS have been found 
    for header, frame, start, length, orf_seq in all_orfs:#iterating through each ORF
        direction = "Reverse" if frame > 3 else "Forward"#determining direction based on the frame number
        print(f">{header} | Frame {frame} ({direction}) | Start: {start} | Length: {length}")#format
        print(orf_seq)
else:#if no ORFS are found, else is executed
    print("No ORFs found.")


