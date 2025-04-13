import streamlit as st
import time

def calculate_length(sequence):
    return len(sequence)

def count_nucleotides(sequence):
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for nucleotide in sequence:
        if nucleotide in counts:
            counts[nucleotide] += 1
    return counts

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    if total_length > 0:
        gc_content = (gc_count / total_length) * 100
        return f"{gc_content:.2f}%"
    else:
        return "0%"

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = sequence[::-1]
    try:
        complement_seq = ''.join(complement[base] for base in reverse_seq)
        return complement_seq
    except KeyError:
        return "Invalid characters in sequence."

# ---------------- Streamlit App ----------------
st.set_page_config(page_title="DNA Sequence Analyzer", layout="centered")
st.title("🧬 DNA Sequence Analyzer")

with st.expander("ℹ️ Input Criteria", expanded=False):
    st.markdown("""
    - Only enter DNA bases: `A`, `T`, `C`, `G` (case-insensitive).
    - No spaces, numbers, or special characters allowed.
    - Example: `ATGCGTACGTTAGC`
    """)

# Input and Button
sequence = st.text_area("Enter a DNA sequence:", height=150).upper()
analyze_button = st.button("🔍 Analyze")

# Process only when button is clicked
if analyze_button:
    valid_chars = {'A', 'T', 'C', 'G'}

    if all(char in valid_chars for char in sequence) and sequence:
        with st.spinner("Analyzing your DNA sequence..."):
            time.sleep(1.5)  # simulate loading

            length = calculate_length(sequence)
            counts = count_nucleotides(sequence)
            gc_content = calculate_gc_content(sequence)
            reverse_comp = reverse_complement(sequence)

        st.success("✅ Analysis Complete!")

        st.subheader("🔍 Analysis Results")
        st.write(f"**Length:** {length} bases")
        st.write("**Nucleotide Counts:**")
        st.write(f"- A: {counts['A']}")
        st.write(f"- T: {counts['T']}")
        st.write(f"- C: {counts['C']}")
        st.write(f"- G: {counts['G']}")
        st.write(f"**GC Content:** {gc_content}")
        st.write(f"**Reverse Complement:** `{reverse_comp}`")

    else:
        st.error("❌ Invalid input! Please enter a valid DNA sequence with only A, T, C, G.")
