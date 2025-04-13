# def get_dna_sequence():
#     """Gets the DNA sequence from the user."""
#     sequence = input("Enter the DNA sequence: ").upper()  # Convert to uppercase
#     return sequence

# if __name__ == "__main__":
#     dna = get_dna_sequence()
#     print(f"You entered: {dna}")
    
import streamlit as st

def get_dna_sequence():
    """Gets the DNA sequence from user input through Streamlit."""
    sequence = st.text_input("Enter the DNA sequence:").upper()
    return sequence

def main():
    st.title("DNA Sequence Input App")

    dna = get_dna_sequence()

    if dna:
        st.success(f"You entered: {dna}")

if __name__ == "__main__":
    main()
