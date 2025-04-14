import streamlit as st

# Define helper functions for the analysis
def calculate_length(seq):
    return len(seq)

def count_nucleotides(seq):
    # Make sure input is in uppercase for consistent counting.
    seq = seq.upper()
    nucleotides = ['A', 'T', 'G', 'C']
    return {nucleotide: seq.count(nucleotide) for nucleotide in nucleotides}

def calculate_gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100

def reverse_complement(seq):
    # Define a mapping for nucleotide complement.
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = seq.upper()
    # Reverse the sequence and obtain complement
    rev_compl = ''.join([complement.get(nuc, nuc) for nuc in reversed(seq)])
    return rev_compl

def find_motifs(seq, motif):
    """Find all starting indices (0-indexed) of motif appearances in the sequence."""
    seq = seq.upper()
    motif = motif.upper()
    indices = []
    start = 0
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        indices.append(idx)
        # Move one position forward to check for overlapping occurrences.
        start = idx + 1
    return indices

def transcribe_dna_to_rna(seq):
    seq = seq.upper()
    return seq.replace("T", "U")

def translate_rna_to_protein(rna_seq):
    # A simple codon table for translation
    codon_table = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C', 'UGC': 'C',
        'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        # Stop codons are represented by an underscore
        'UAA': '_', 'UAG': '_', 'UGA': '_'
    }
    
    protein = ""
    rna_seq = rna_seq.upper()
    # Process the RNA sequence in chunks (codons) of 3
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        protein += codon_table.get(codon, 'X')  # 'X' for any unknown codon
    return protein

# Streamlit app interface
st.title("DNA Sequence Analysis App")

# Input DNA sequence (multiline text area)
dna_sequence = st.text_area("Enter DNA sequence:", value="", height=150)

if dna_sequence:
    # Provide a select box to choose the analysis option.
    analysis_choice = st.selectbox("Select Analysis Option:", 
                                   [
                                    "Calculate Length", 
                                    "Count Nucleotides", 
                                    "Calculate GC Content", 
                                    "Find Reverse Complement", 
                                    "Find Motifs", 
                                    "Transcribe to RNA", 
                                    "Translate to Protein (Basic)"
                                   ])

    # Make sure the sequence is stripped of any whitespace/newlines.
    dna_sequence = dna_sequence.replace("\n", "").replace(" ", "")

    if analysis_choice == "Calculate Length":
        length = calculate_length(dna_sequence)
        st.write(f"**Length:** {length}")

    elif analysis_choice == "Count Nucleotides":
        nucleotide_counts = count_nucleotides(dna_sequence)
        st.write("**Nucleotide Counts:**")
        st.json(nucleotide_counts)

    elif analysis_choice == "Calculate GC Content":
        gc_content = calculate_gc_content(dna_sequence)
        st.write(f"**GC Content:** {gc_content:.2f}%")

    elif analysis_choice == "Find Reverse Complement":
        rev_compl = reverse_complement(dna_sequence)
        st.write("**Reverse Complement:**")
        st.code(rev_compl)

    elif analysis_choice == "Find Motifs":
        # Provide an additional input for the motif
        motif = st.text_input("Enter a motif to search for:").strip().upper()
        if motif:
            indices = find_motifs(dna_sequence, motif)
            if indices:
                st.write(f"Motif '**{motif}**' found at indices:", indices)
            else:
                st.write(f"Motif '**{motif}**' not found.")

    elif analysis_choice == "Transcribe to RNA":
        rna_sequence = transcribe_dna_to_rna(dna_sequence)
        st.write("**RNA Sequence:**")
        st.code(rna_sequence)

    elif analysis_choice == "Translate to Protein (Basic)":
        rna_sequence = transcribe_dna_to_rna(dna_sequence)
        protein_sequence = translate_rna_to_protein(rna_sequence)
        st.write("**Protein Sequence (Basic):**")
        st.code(protein_sequence)
