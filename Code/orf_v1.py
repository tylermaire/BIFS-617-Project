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

def find_orfs(seq, frame, minlen, rc=False):
    """Find ORFs in a specific frame (1–6)."""
    orfs = []
    n = len(seq)
    i = (frame - 1) % 3
    while i + 3 <= n:
        if seq[i:i+3] == 'ATG':
            j = i + 3
            while j + 3 <= n:
                codon = seq[j:j+3]
                if codon in ('TAA', 'TAG', 'TGA'):
                    length = j + 3 - i
                    if length >= minlen:
                        pos = i + 1 if not rc else -(n - i)
                        orfs.append((frame, pos, length, seq[i:j+3]))
                    break
                j += 3
        i += 3
    return orfs

def format_orf(header, frame, pos, length, seq):
    """Format ORF info in FASTA with codons spaced, 15 per line."""
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    lines = [' '.join(codons[i:i+15]) for i in range(0, len(codons), 15)]
    return f"> {header} | FRAME = {frame} POS = {pos} LEN = {length}\n" + '\n'.join(lines)

def main():
    total = 0
    found = []
    per_seq = []

    if fasta == '-':
        print("Paste your FASTA sequence (end with blank line):")
        lines = []
        while True:
            line = input()
            if line.strip() == '':
                break
            lines.append(line)
    else:
        with open(fasta, 'r') as f:
            lines = f.readlines()

    for header, seq in parse_fasta(lines):
        seq = seq.upper()
        count = 0
        for f in range(1, 4):
            for orf in find_orfs(seq, f, minlen):
                print(format_orf(header, *orf))
                found.append((header, orf[0], orf[1], orf[2]))
                count += 1
        rc = revcomp(seq)
        for f in range(4, 7):
            for orf in find_orfs(rc, f, minlen, rc=True):
                print(format_orf(header, *orf))
                found.append((header, orf[0], orf[1], orf[2]))
                count += 1
        per_seq.append((header, count))
        total += count

    # Summary
    print("\n--- ORF SUMMARY ---")
    print(f"Total ORFs found: {total}")
    for h, c in per_seq:
        print(f"Sequence '{h}': {c} ORFs")
    if total:
        print("\nDetails:")
        for i, (h, f, p, l) in enumerate(found, 1):
            print(f"{i}. {h}, Frame: {f}, Pos: {p}, Length: {l}")
    else:
        print("No ORFs found.")

if __name__ == '__main__':
    main()

