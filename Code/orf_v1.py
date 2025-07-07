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



if __name__ == '__main__':
    main()

