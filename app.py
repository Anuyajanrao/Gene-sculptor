# -------------------- All your imports --------------------
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import kaleido
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import tempfile
import subprocess
import logomaker
from collections import Counter
import numpy as np
import altair as alt
import base64
import mimetypes
from PIL import Image
import os
import re
from io import StringIO, BytesIO
import gzip
from fpdf import FPDF
import io

# -------------------- Streamlit Config --------------------

st.set_page_config(page_title="Genome Sculptor", layout="wide", page_icon="🧬")

def set_bg_hack_url():
    '''
    A function to set a dimmed background image using linear gradient overlay.
    '''
    st.markdown(
        f"""
        <style>
        .stApp {{
            background: linear-gradient(rgba(255, 255, 255, 0.6), rgba(255, 255, 255, 0.6)),
                        url(https://www.shutterstock.com/image-illustration/3d-illustration-x-chromosome-medicine-600nw-2163174573.jpg);
            background-size: cover;
            background-attachment: fixed;
            background-position: center;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

# Call the function to set the dimmed background image
set_bg_hack_url()



# -------------------- Initialize session state --------------------
if 'current_section' not in st.session_state:
    st.session_state.current_section = "Home"
if 'fastq_data' not in st.session_state:
    st.session_state.fastq_data = None
if 'qual_sums' not in st.session_state:
    st.session_state.qual_sums = {} # Initialize quality sums in state

# Function to navigate by updating session state
def navigate_to(section):
    st.session_state.current_section = section
    st.rerun() # Rerun to show the new section immediately

# --- Custom CSS for Cards ---
st.markdown("""
<style>
/* Style for buttons that look like cards */
div.stButton > button:first-child {
    background-color: rgba(255, 255, 255, 0.9);
    border-radius: 12px;
    padding: 20px;
    margin-top: 15px;
    margin-bottom: 15px;
    cursor: pointer;
    transition: transform 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
    text-align: center;
    color: black;
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
    width: 100%; /* Make button take full width of its container */
    border: none; /* Remove default button border */
    font-size: 1.1rem; /* Adjust font size */
    height: 100%; /* Make cards in a row have similar height */
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    white-space: normal; /* Allow text to wrap */
    line-height: 1.5;
}
div.stButton > button:first-child:hover {
    transform: translateY(-5px);
    box-shadow: 0 8px 16px rgba(0, 0, 0, 0.2);
}

/* Ensure button text is black */
div.stButton > button:first-child > p {
    color: black !important;
    margin: 0; /* Remove default paragraph margin */
    word-break: break-word; /* Break long words */
}

/* Adjust padding for smaller screens if needed */
@media (max-width: 768px) {
    div.stButton > button:first-child {
        padding: 15px;
        font-size: 1rem;
    }
}

/* Specific styles for the big cards on the Home page */
.home-card-button > button:first-child {
    height: 250px; /* Make the card taller */
    font-size: 1.5rem; /* Larger font size */
    font-weight: bold;
    background-color: rgba(255, 255, 255, 0.95); /* Slightly less transparent */
    box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
}

.home-card-button > button:first-child:hover {
    transform: translateY(-8px); /* More pronounced hover effect */
    box-shadow: 0 12px 24px rgba(0, 0, 0, 0.3);
}

/* Style for the icons/emojis within the big cards */
.home-card-button > button:first-child .card-icon {
    font-size: 3em; /* Large icon size */
    margin-bottom: 10px;
}

/* Custom CSS for About section */
.about-section-container {
    background-color: rgba(255, 255, 255, 0.9);
    padding: 20px;
    border-radius: 12px;
    color: black;
    margin-bottom: 20px;
}

.about-section-container h3,
.about-section-container h4 {
    color: #000000; /* Darker color for headings */
}

.about-section-container p {
     color: #333333; /* Darker color for text */
     font-size: 1rem; /* Standard paragraph size */
}

.about-section-container a {
     color: #0066cc; /* Link color */
     text-decoration: none;
}

.about-section-container a:hover {
     text-decoration: underline;
}

.profile-image {
    border-radius: 50%;
    border: 3px solid #6A5ACD;
    width: 120px;
    height: 120px;
    object-fit: cover; /* Ensure image covers the circular area */
}

</style> 
""", unsafe_allow_html=True)

## -------------------- Custom Styling --------------------
st.markdown("""
<style>
/* Global font smoothing */
html, body, [class*="css"] {
    -webkit-font-smoothing: antialiased;
    font-family: 'Arial', sans-serif;
}

/* ---------- Sidebar Styling ---------- */
[data-testid="stSidebar"] {
    background: linear-gradient(135deg, #dfe9f3, #ffffff);
    padding-top: 30px;
    border-right: 2px solid #d1d9e6;
}

/* Sidebar Title */
.sidebar-title {
    font-size: 28px;
    font-weight: bold;
    color: #264653;
    text-align: center;
    margin-bottom: 30px;
}

/* Sidebar Buttons */
div.stButton > button {
    width: 100%;
    padding: 12px 0;
    margin-bottom: 15px;
    font-size: 18px;
    font-weight: 600;
    color: #264653;
    background: linear-gradient(135deg, #d0f0f4, #e5f7fb);
    border-radius: 12px;
    border: none;
    box-shadow: 0 3px 6px rgba(0,0,0,0.1);
    transition: all 0.3s ease-in-out;
}

div.stButton > button:hover {
    background: linear-gradient(135deg, #d8f3dc, #a8dadc);
    color: #1d3557;
    transform: translateY(-2px);
}

/* ---------- Home Page Styling ---------- */
.home-title {
    font-size: 3.5rem;
    font-weight: 800;
    text-align: center;
    margin-bottom: 25px;
    letter-spacing: 2px;
    text-transform: uppercase;
    background: linear-gradient(135deg, #000000, #434343); /* Black gradient */
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    text-decoration: none; /* Remove underline */
    padding: 10px;
}

/* White background box behind the title */
.home-title-container {
    background-color: white; /* White background for the title box */
    display: inline-block; /* Makes the background box just around the text */
    padding: 20px 40px; /* Add some padding around the text */
    border-radius: 15px; /* Rounded corners for the background */
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1); /* Subtle shadow for depth */
}



/* Welcome Box */
.home-welcome-box {
    background: #ffffff;  /* Pure white background */
    padding: 30px;
    border-radius: 15px;
    box-shadow: 0 6px 12px rgba(0,0,0,0.1);
    margin-bottom: 30px;
}

/* Welcome Box Text */
.home-welcome-box h2 {
    color: #000000; /* Black color for text */
    text-align: center;
    font-size: 26px;
    font-weight: bold;
}

.home-welcome-box p {
    font-size: 19px;
    color: #000000; /* Black color for text */
    text-align: center;
}



/* Home Buttons */
.home-card-button > button:first-child {
    height: 220px;
    font-size: 22px !important;
    font-weight: bold;
    color: white !important;
    background: linear-gradient(135deg, #a18cd1, #fbc2eb);
    border-radius: 15px;
    box-shadow: 0 6px 12px rgba(0,0,0,0.1);
    transition: transform 0.3s ease, box-shadow 0.3s ease, background 0.3s ease;
    padding: 25px;
    border: none;
}

.home-card-button > button:first-child:hover {
    transform: translateY(-6px);
    box-shadow: 0 12px 24px rgba(0,0,0,0.15);
    background: linear-gradient(135deg, #fbc2eb, #a18cd1);
}

/* Responsive Adjustments */
@media (max-width: 768px) {
    .home-card-button > button:first-child {
        height: 180px;
        font-size: 18px !important;
    }
}
</style>
""", unsafe_allow_html=True)

# -------------------- Sidebar Navigation --------------------
with st.sidebar:
    st.markdown('<div class="sidebar-title">🧭 Navigation</div>', unsafe_allow_html=True)
    if st.button("🏠 Home"):
        st.session_state.current_section = "Home"
    if st.button("📖 User Guide"):
        st.session_state.current_section = "User Guide"
    if st.button("ℹ️ About"):
        st.session_state.current_section = "About Creator"
    if st.button("❓ FAQ"):
        st.session_state.current_section = "FAQ"

# -------------------- Display Section Content --------------------
if st.session_state.current_section == "Home":
    st.markdown(f"<div class='home-title'>🧬 <u>Genome Sculptor</u></div>", unsafe_allow_html=True)
    st.markdown("""
    <div class='home-welcome-box'>
        <h2><strong>🔬 Welcome!</strong></h2>
        <p>This interactive platform lets you <b>analyze</b> and <b>visualize</b> <strong>DNA</strong>, <strong>RNA</strong>, and <strong>Protein sequences</strong> with ease.</p>
        <p>👇 <b>Select an analysis type below</b> to get started:</p>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown('<div class="home-card-button">', unsafe_allow_html=True)
        if st.button("🧬 FASTA Analysis", key="fasta_analysis_card_btn", help="Click to analyze FASTA sequences", use_container_width=True):
            navigate_to("FASTA Analysis")
        st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="home-card-button">', unsafe_allow_html=True)
        if st.button("🔬 FASTQ Analysis", key="fastq_qc_card_btn", help="Click to perform QC on FASTQ files", use_container_width=True):
            navigate_to("FASTQ QC")
        st.markdown('</div>', unsafe_allow_html=True)

# --- FASTA Analysis Section ---
elif st.session_state.current_section == "FASTA Analysis":
    st.title("📂 FASTA Analysis Tools")
    st.markdown("""
    <div style='background-color: rgba(255, 255, 255, 0.9); padding: 20px; border-radius: 12px; color: black;'>
        <p>Select a tool from the cards below to perform analysis on your sequences.</p>
    </div>
    """, unsafe_allow_html=True)

    # Back button
    if st.button("⬅️ Back to Home", key="fasta_back_home"):
        navigate_to("Home") # Correct navigation back to Home

    st.markdown("---")

    # Grid of FASTA tool cards
    fasta_tools = [
        ("GC Content Plot", "🟢 Analyze GC content over sequence windows."),
        ("Nucleotide Composition", "🧩 Visualize the percentage of A, T, G, C bases."),
        ("Sequence Logo Generator", "🧬 Create a sequence logo from aligned sequences."),
        ("Amino Acid Frequency", "🔡 View the frequency of amino acids in a protein sequence."),
        ("Gene Expression Heatmap", "🌡️ Visualize gene expression matrix data."),
        ("Motif Search", "🔍 Find specific patterns (motifs) in a sequence."),
        ("Repeat Analysis", "🔄 Identify and count repetitive k-mers."),
        ("Primer Design", "🧪 Suggest basic forward and reverse primers."),
        ("Circular Genome Plot", "⭕ Plot features on a circular genome map.")
    ]

    cols = st.columns(3) # Use 3 columns for the grid

    for i, (tool_name, tool_desc) in enumerate(fasta_tools):
        with cols[i % 3]: # Cycle through columns
             # Use f-string with emojis and bold text for card content
             if st.button(f"**{tool_name}**\n\n{tool_desc}", key=f"fasta_tool_card_{tool_name.replace(' ', '_').lower()}", help=tool_desc, use_container_width=True):
                 navigate_to(tool_name)


# --- FASTQ QC Section ---
elif st.session_state.current_section == "FASTQ QC":
    st.title("📊 FASTQ Quality Control")

    if st.button("⬅️ Back to Home", key="fastq_back_home"):
        navigate_to("Home")

    st.markdown("---")

    st.markdown("""
    <div style='background-color: rgba(255, 255, 255, 0.9); padding: 20px; border-radius: 12px; color: black;'>
        <p>Select a report section from the cards below to analyze your FASTQ file.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")

    if st.session_state.current_section == "FASTQ QC":
        st.title("🧪 FASTQ Quality Control Dashboard")
        st.markdown("Click a QC check below to view details and run the analysis.")

    card_sections = [
        ("📊 Summary Metrics", "summary"),
        ("📈 Read Length Distribution", "read_len"),
        ("🧬 GC Content Distribution", "gc_dist"),
        ("🔍 Per-Base Quality Scores", "pb_qual"),
        ("🔄 Sequence Duplication Levels", "dup_levels"),
        ("🔗 Adapter Content", "adapter"),
        ("📌 Overrepresented Sequences", "overrep"),
        ("🧩 K-mer Content", "kmer"),
        ("❓ N Content per Base", "n_content"),
        ("🎯 Base Composition Bias Heatmap", "base_bias"),
    ]

    # Layout cards in rows of 3
    for i in range(0, len(card_sections), 3):
        cols = st.columns(3)
        for col, (title, section_key) in zip(cols, card_sections[i:i+3]):
            with col:
                with st.container(border=True):
                    if st.button(title, key=f"card_{section_key}"):
                        st.session_state.current_section = title.replace("📊 ", "").replace("📈 ", "").replace("🧬 ", "")\
                                                                 .replace("🔍 ", "").replace("🔄 ", "").replace("🔗 ", "")\
                                                                 .replace("📌 ", "").replace("🧩 ", "").replace("❓ ", "")\
                                                                 .replace("🎯 ", "")
                        st.rerun()


# --- Individual FASTA Tool Pages ---

elif st.session_state.current_section == "GC Content Plot":
    st.title("🟢 GC Content Plot")
    if st.button("⬅️ Back to FASTA Tools", key="gc_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    uploaded_file = st.file_uploader("Upload FASTA file", type="fasta", key="gc_fasta_uploader")
    sequence_input = st.text_area("Or paste a DNA sequence", key="gc_seq_input")
    window_size = st.slider("Window size", min_value=10, max_value=200, value=100, step=10, key="gc_window_size")
    if st.button("Generate Output", key="gc_generate_btn"):
        sequences = []
        if uploaded_file:
            try:
                records = list(SeqIO.parse(uploaded_file, "fasta"))
                sequences.extend(records)
            except Exception as e:
                st.error(f"Error reading FASTA file: {e}")
                sequences = [] # Clear sequences on error
        elif sequence_input:
            sequences.append(SeqIO.SeqRecord(seq=sequence_input.upper(), id="User_Input"))

        if sequences:
            for record in sequences:
                seq = str(record.seq) # Ensure it's a string
                if not all(base in 'ATGCatgc' for base in seq):
                     st.warning(f"Skipping {record.id}: Contains non-DNA characters.")
                     continue
                if len(seq) < window_size:
                     st.warning(f"Skipping {record.id}: Sequence length ({len(seq)}) is less than window size ({window_size}).")
                     continue

                gc_content = [
                    (seq[i:i+window_size].upper().count("G") + seq[i:i+window_size].upper().count("C")) / window_size
                    for i in range(len(seq) - window_size + 1)
                ]
                x = np.arange(len(gc_content))
                plt.figure(figsize=(10, 5))
                plt.plot(x, gc_content, color="green")
                plt.fill_between(x, gc_content, color="lightgreen", alpha=0.5)
                plt.title(f"GC Content Plot for {record.id}")
                plt.xlabel("Position")
                plt.ylabel("GC Content Fraction")
                st.pyplot(plt)
        else:
             st.info("Please upload a FASTA file or paste a DNA sequence.")


elif st.session_state.current_section == "Nucleotide Composition":
    st.title("🧩 Nucleotide Composition")
    if st.button("⬅️ Back to FASTA Tools", key="nuc_comp_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    seq = st.text_area("Enter DNA sequence", key="nuc_comp_seq_input")
    if st.button("Generate Output", key="nuc_comp_generate_btn"):
        if seq:
            counts = Counter(seq.upper())
            valid_counts = {base: count for base, count in counts.items() if base in 'ATGC'}
            if not valid_counts:
                st.warning("No valid DNA nucleotides (A, T, G, C) found in the sequence.")
            else:
                fig = px.pie(names=list(valid_counts.keys()), values=list(valid_counts.values()), title="Nucleotide Composition")
                st.plotly_chart(fig)
        else:
            st.info("Please enter a DNA sequence.")


elif st.session_state.current_section == "Sequence Logo Generator":
    st.title("🧬 Sequence Logo Generator (with Automatic Alignment)")
    if st.button("⬅️ Back to FASTA Tools", key="logo_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    st.info("Paste multiple **unaligned** DNA sequences in FASTA format or upload a FASTA file. The tool will align sequences and generate a high-quality logo.")
    uploaded_file = st.file_uploader("📁 Upload FASTA file", type="fasta", key="logo_fasta_uploader")
    msa_input = st.text_area("🧾 Or paste raw FASTA sequences below", height=200, key="logo_msa_input")

    if st.button("Generate Sequence Logo", key="logo_generate_btn"):
        records = []
        if uploaded_file:
            try:
                records = list(SeqIO.parse(uploaded_file, "fasta"))
            except Exception as e:
                st.error(f"Error reading FASTA file: {e}")
                records = []
        elif msa_input:
            msa_handle = StringIO(msa_input)
            try:
                records = list(SeqIO.parse(msa_handle, "fasta"))
            except Exception as e:
                 st.error(f"Error parsing pasted FASTA sequence: {e}")
                 records = []


        if not records or len(records) < 2:
            st.error("❌ Please provide at least two sequences.")
        else:
            if not all(all(base in 'ATGCatgc\n\r-' for base in str(rec.seq)) for rec in records):
                 st.error("❌ Please ensure sequences contain only valid DNA characters (A, T, G, C) and alignment gaps (-) if pasting aligned sequences.")
                 st.stop()

            input_fasta_path = None
            aligned_output_path = None

            try:
                with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as input_fasta:
                    SeqIO.write(records, input_fasta, "fasta")
                    input_fasta_path = input_fasta.name

                aligned_output_path = input_fasta_path + "_aligned.fasta"

                try:
                    # Check if clustalo is available
                    subprocess.run(["clustalo", "--version"], check=True, capture_output=True)
                except (subprocess.CalledProcessError, FileNotFoundError):
                    st.error("❌ Clustal Omega is not found or not in your system's PATH. Please install it to use this feature.")
                    st.stop()

                process = subprocess.run([
                    "clustalo", "-i", input_fasta_path, "-o", aligned_output_path,
                    "--force", "--outfmt=fasta"
                ], check=True, capture_output=True, text=True)

                aligned_records = list(SeqIO.parse(aligned_output_path, "fasta"))
                aligned_seqs = [str(rec.seq).upper() for rec in aligned_records]

                if not aligned_seqs:
                     st.error("❌ Alignment produced no sequences. Please check your input.")
                     st.stop()

                if len(set(len(seq) for seq in aligned_seqs)) != 1:
                    st.error("❌ Alignment failed. Sequences must be same length after alignment.")
                    st.stop()

                if any(len(seq) == 0 for seq in aligned_seqs):
                    st.error("❌ Alignment resulted in empty sequences. Please check your input.")
                    st.stop()

                matrix = logomaker.alignment_to_matrix(aligned_seqs, to_type='information')

                fig, ax = plt.subplots(figsize=(len(aligned_seqs[0]) * 0.5, 6), dpi=150)
                logo = logomaker.Logo(matrix, ax=ax,
                                      color_scheme={'A': 'green', 'T': 'red', 'G': 'orange', 'C': 'blue'})
                ax.set_title("Sequence Logo", fontsize=20, weight='bold')
                ax.set_xlabel("Position (5′ → 3′)", fontsize=16)
                ax.set_ylabel("Information (bits)", fontsize=16)
                ax.set_ylim(0, 0.5) # Adjusted y-limit for better visualization
                ax.axhline(0, color='black', linewidth=0.5)
                ax.tick_params(labelsize=12)

                img_buffer = BytesIO()
                fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                img_buffer.seek(0)
                encoded = base64.b64encode(img_buffer.read()).decode()

                st.components.v1.html(
                    f"""
                    <div style="overflow-x: scroll; border:1px solid #ccc; padding:10px;">
                        <img src="data:image/png;base64,{encoded}" style="height:300px">
                    </div>
                    """,
                    height=350,
                )

                st.download_button(
                    label="📥 Download Sequence Logo as PNG",
                    data=BytesIO(base64.b64decode(encoded)),
                    file_name="sequence_logo.png",
                    mime="image/png"
                )

            except subprocess.CalledProcessError as e:
                st.error(f"❌ Clustal Omega execution failed. Check your input or Clustal Omega installation.")
                st.error(f"Stderr: {e.stderr.decode()}")
                st.error(f"Stdout: {e.stdout.decode()}")
            except FileNotFoundError:
                 st.error("❌ Clustal Omega command not found. Please ensure it is installed and in your system's PATH.")
            except Exception as e:
                 st.error(f"An unexpected error occurred during logo generation: {e}")
                 st.exception(e) # Show traceback for other errors
            finally:
                if input_fasta_path and os.path.exists(input_fasta_path):
                     os.remove(input_fasta_path)
                if aligned_output_path and os.path.exists(aligned_output_path):
                     os.remove(aligned_output_path)


elif st.session_state.current_section == "Amino Acid Frequency":
    st.title("🧬 Amino Acid Frequency")
    if st.button("⬅️ Back to FASTA Tools", key="aa_freq_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    protein_seq = st.text_area("Enter protein sequence", key="aa_freq_seq_input")
    if st.button("Generate Output", key="aa_freq_generate_btn"):
        if protein_seq:
            counts = Counter(protein_seq.upper())
            valid_counts = {aa: count for aa, count in counts.items() if aa in "ACDEFGHIKLMNPQRSTVWY"}
            if not valid_counts:
                 st.warning("No valid amino acid characters found in the sequence. Please enter a protein sequence.")
            else:
                df = pd.DataFrame(list(valid_counts.items()), columns=["Amino Acid", "Count"])
                fig = px.bar(df, x="Amino Acid", y="Count", title="Amino Acid Frequency")
                st.plotly_chart(fig)
        else:
             st.info("Please enter a protein sequence.")


elif st.session_state.current_section == "Gene Expression Heatmap":
    st.title("🌡️ Gene Expression Heatmap")
    if st.button("⬅️ Back to FASTA Tools", key="heatmap_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    uploaded_file = st.file_uploader("Upload gene expression matrix (CSV)", type="csv", key="heatmap_csv_uploader")
    if uploaded_file and st.button("Generate Output", key="heatmap_generate_btn"):
        try:
            df = pd.read_csv(uploaded_file, index_col=0)
            df_numeric = df.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
            if df_numeric.empty:
                st.error("❌ The CSV file does not contain valid numeric data for the heatmap.")
            else:
                fig, ax = plt.subplots(figsize=(12, max(6, len(df_numeric.index) * 0.3))) # Adjust figure height
                sns.heatmap(df_numeric, cmap="viridis", ax=ax)
                plt.title("Gene Expression Heatmap")
                st.pyplot(fig)
        except Exception as e:
             st.error(f"❌ Error reading or processing the CSV file: {e}")
             st.exception(e)
    elif not uploaded_file:
         st.info("Please upload a CSV file.")


elif st.session_state.current_section == "Motif Search":
    st.title("🔍 Motif Search")
    if st.button("⬅️ Back to FASTA Tools", key="motif_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    sequence = st.text_area("Enter DNA sequence", key="motif_seq_input")
    motif = st.text_input("Enter motif (regex)", key="motif_regex_input")
    if st.button("Generate Output", key="motif_generate_btn"):
        if sequence and motif:
            try:
                matches = list(re.finditer(motif, sequence))
                if matches:
                    st.success(f"Found {len(matches)} matches")
                    for match in matches:
                        st.write(f"Position: {match.start()} | Match: {match.group()}")
                else:
                    st.warning("No matches found.")
            except re.error as e:
                 st.error(f"❌ Invalid regular expression: {e}")
        else:
            st.info("Please enter both a DNA sequence and a motif to search.")


elif st.session_state.current_section == "Repeat Analysis":
    st.title("🔄 Repeat Analysis")
    if st.button("⬅️ Back to FASTA Tools", key="repeat_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    sequence = st.text_area("Enter DNA sequence", key="repeat_seq_input")
    k = st.slider("Select k-mer size", 2, 10, 3, key="repeat_kmer_size") # Increased max k-mer size
    if sequence and st.button("Generate Output", key="repeat_generate_btn"):
        seq = sequence.upper() # Convert to upper case for consistency
        if not all(base in 'ATGC' for base in seq):
            st.warning("Please ensure the sequence contains only DNA characters (A, T, G, C).")
        elif len(seq) < k:
             st.warning(f"Sequence length ({len(seq)}) is less than k-mer size ({k}). Cannot perform analysis.")
        else:
            kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
            counts = Counter(kmers)
            if counts:
                df = pd.DataFrame(counts.items(), columns=["k-mer", "Count"]).sort_values(by="Count", ascending=False)
                st.write(df)
            else:
                 st.info("No k-mers found in the sequence.")
    elif not sequence:
         st.info("Please enter a DNA sequence to perform repeat analysis.")


elif st.session_state.current_section == "Primer Design":
    st.title("🧪 Primer Design")
    if st.button("⬅️ Back to FASTA Tools", key="primer_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    sequence = st.text_area("Enter DNA sequence", key="primer_seq_input")
    primer_len = st.slider("Suggested Primer Length", 18, 25, 20, key="primer_length") # Allow adjusting length
    if sequence and st.button("Generate Output", key="primer_generate_btn"):
        seq = sequence.upper()
        if not all(base in 'ATGC' for base in seq):
            st.warning("Please ensure the sequence contains only DNA characters (A, T, G, C).")
        elif len(seq) < primer_len * 2:
             st.warning(f"Sequence is too short for primer design with length {primer_len} (minimum {primer_len * 2} bases recommended).")
        else:
            st.info(f"Suggested primer pairs (basic example - length {primer_len})")
            st.write(f"Forward Primer: {seq[:primer_len]}")
            # Simple reverse complement
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            reverse_complement = "".join([complement[base] for base in seq[-primer_len:][::-1]])
            st.write(f"Reverse Primer: {reverse_complement}")
            st.caption("Note: This is a very basic primer design suggestion based on length. For complex or critical applications (e.g., specific melting temperatures, avoiding hairpins/dimers), use dedicated primer design software.")
    elif not sequence:
        st.info("Please enter a DNA sequence for primer design.")


elif st.session_state.current_section == "Circular Genome Plot":
    st.title("⭕ Circular Genome Plot")
    if st.button("⬅️ Back to FASTA Tools", key="circular_back_fasta"):
        navigate_to("FASTA Analysis")
    st.markdown("---")
    st.info("Upload a CSV file with columns: `start`, `end`, `label`. Optionally include `color`.")
    uploaded_file = st.file_uploader("Upload genome feature file (CSV)", type="csv", key="circular_csv_uploader")
    if uploaded_file and st.button("Generate Output", key="circular_generate_btn"):
        try:
            df = pd.read_csv(uploaded_file)
            required_cols = ['start', 'end', 'label']
            if not all(col in df.columns for col in required_cols):
                 st.error(f"❌ CSV must contain columns: {', '.join(required_cols)}. Optional: 'color'.")
            else:
                if not pd.api.types.is_numeric_dtype(df['start']) or not pd.api.types.is_numeric_dtype(df['end']):
                     st.error("❌ 'start' and 'end' columns must contain numeric values.")
                elif (df['start'] < 0).any() or (df['end'] < 0).any() or (df['start'] > df['end']).any():
                    st.error("❌ Invalid start/end positions. Start must be <= End and both must be non-negative.")
                else:
                    # Determine genome length for better scaling if possible (e.g., from max end)
                    # For a true circular plot, you might need the total genome length from the user
                    # Here, we'll just scale based on the max end position found
                    max_pos = df['end'].max() if not df.empty else 1
                    if max_pos == 0: max_pos = 1 # Avoid division by zero

                    fig = go.Figure()

                    for _, row in df.iterrows():
                        start = row['start']
                        end = row['end']
                        label = row['label']
                        # Use color column if available, otherwise default
                        color = row['color'] if 'color' in df.columns else None

                        # Calculate the angle range for the feature (scaled to 0-360 degrees)
                        # This is a simplified representation; a proper circular plot library is better for complex genomes
                        start_angle = (start / max_pos) * 360
                        end_angle = (end / max_pos) * 360

                        # Add a sector or bar polar trace
                        fig.add_trace(go.Barpolar(
                            r=[1], # Radial position (making it a circle)
                            theta=[start_angle], # Start angle
                            width=(end_angle - start_angle), # Angle width
                            marker_color=color if color else px.colors.qualitative.Plotly[0], # Use default Plotly colors if no color column
                            name=label, # Label for legend/hover
                            # text=label, # Hover text
                            # hoverinfo='text'
                        ))

                    fig.update_layout(
                        title="Circular Genome Plot",
                        polar=dict(
                            radialaxis=dict(visible=False), # Hide radial axis
                            angularaxis=dict(showticklabels=False, ticks='') # Hide angular axis ticks/labels
                        ),
                        showlegend=True # Show legend for labels
                    )
                    st.plotly_chart(fig)
        except Exception as e:
             st.error(f"❌ Error reading or processing the CSV file: {e}")
             st.exception(e) # Display traceback for debugging
    elif not uploaded_file:
         st.info("Please upload a CSV file with 'start', 'end', and 'label' columns.")


# --- Individual FASTQ QC Report Sections (Individual Pages) ---

# Centralized FASTQ file uploader and processing
def process_fastq_file_and_store(uploaded_file):
    if uploaded_file is None:
        st.session_state.fastq_data = None
        st.session_state.qual_sums = {}
        return None

    file_name = uploaded_file.name
    try:
        if file_name.endswith(".gz"):
            file_handle = gzip.open(BytesIO(uploaded_file.getvalue()), "rt")
        else:
            file_handle = io.TextIOWrapper(BytesIO(uploaded_file.getvalue()))

    except Exception as e:
        st.error(f"Error opening file: {e}")
        st.session_state.fastq_data = None
        st.session_state.qual_sums = {}
        st.exception(e)
        return None

    st.info("⚙️ Processing FASTQ file (please wait)...")

    sequences = []
    qualities = [] # Store quality lists per read
    gc_counts = []
    seq_counter = Counter()
    total_reads = 0
    sum_read_len = 0
    sum_gc = 0
    q20_bases = 0 # Count bases with Q >= 20
    q30_bases = 0 # Count bases with Q >= 30
    total_bases = 0
    lengths_hist = Counter()
    gc_hist = Counter()
    adapter_hits = Counter()
    adapters = ["AGATCGGAAGAGC", "GATCGGAAGAGC"] # Common Illumina adapters
    max_read_len = 0
    base_counts = {"A": Counter(), "T": Counter(), "G": Counter(), "C": Counter(), "N": Counter()}

    # Initialize qual_sums for per-base quality calculation
    st.session_state.qual_sums = {}


    try:
        for title, seq, qual_str in FastqGeneralIterator(file_handle):
            seq = seq.upper()
            qual = [ord(ch) - 33 for ch in qual_str]

            read_len = len(seq)
            if read_len == 0:
                continue

            max_read_len = max(max_read_len, read_len)
            total_reads += 1
            sum_read_len += read_len
            total_bases += read_len

            gc_count = seq.count('G') + seq.count('C')
            sum_gc += gc_count

            # Count bases with Q20/Q30
            q20_bases += sum(1 for q in qual if q >= 20)
            q30_bases += sum(1 for q in qual if q >= 30)


            lengths_hist[read_len] += 1
            gc_pct = int((gc_count / read_len) * 100) if read_len > 0 else 0
            gc_hist[gc_pct] += 1

            seq_counter[seq] += 1

            for adapter in adapters:
                if adapter in seq:
                    adapter_hits[adapter] += 1

            # Per-base nucleotide counts
            for i, base in enumerate(seq):
                if base in base_counts:
                    base_counts[base][i] += 1
                else:
                    # Handle bases not A, T, G, C (e.g., N)
                    if 'N' not in base_counts:
                         base_counts['N'] = Counter()
                    base_counts["N"][i] += 1


            # Per-base quality sums - update session state directly
            for i, q in enumerate(qual):
                 if i not in st.session_state.qual_sums:
                    st.session_state.qual_sums[i] = 0
                 st.session_state.qual_sums[i] += q

            # Store the quality list for this read (needed for avg read quality calculation later)
            qualities.append(qual)


        if total_reads == 0:
            st.warning("No valid reads found in the file.")
            st.session_state.fastq_data = None
            st.session_state.qual_sums = {}
            return None

        avg_read_len = sum_read_len / total_reads if total_reads > 0 else 0
        avg_gc_content = (sum_gc / sum_read_len) * 100 if sum_read_len > 0 else 0
        pct_q20_bases = (q20_bases / total_bases) * 100 if total_bases > 0 else 0
        pct_q30_bases = (q30_bases / total_bases) * 100 if total_bases > 0 else 0

        # Store processed data in session state
        st.session_state.fastq_data = {
            "total_reads": total_reads,
            "avg_read_len": avg_read_len,
            "avg_gc_content": avg_gc_content,
            "pct_q20_bases": pct_q20_bases,
            "pct_q30_bases": pct_q30_bases,
            "lengths_hist": lengths_hist,
            "gc_hist": gc_hist,
            "seq_counter": seq_counter,
            "adapter_hits": adapter_hits,
            "max_read_len": max_read_len,
            "base_counts": base_counts,
            "total_bases": total_bases,
            # qual_sums is already in session state
        }
        st.success("✅ FASTQ file processed successfully!")
        return st.session_state.fastq_data

    except Exception as e:
        st.error(f"An error occurred during file processing: {e}")
        st.session_state.fastq_data = None
        st.session_state.qual_sums = {}
        st.exception(e) # Display traceback
        return None


# Function to display the file uploader and processing button
def show_fastq_uploader_and_process(section_key):
    uploaded_file = st.file_uploader("Upload your FASTQ file (supports FASTQ, FASTQ.gz)", type=["fastq", "fq", "gz"], key=f"fastq_uploader_{section_key}")
    if uploaded_file and st.button("Process File and Generate Output", key=f"process_fastq_btn_{section_key}"):
        # Clear previous data before processing a new file
        st.session_state.fastq_data = None
        st.session_state.qual_sums = {} # Important: Reset quality sums
        process_fastq_file_and_store(uploaded_file)
        # No explicit rerun needed here, as the button click and state update trigger it


# --- FASTQ QC Report Sections (Individual Pages) ---

if st.session_state.current_section == "Summary Metrics":
    st.title("📊 FASTQ Quality Control Summary")
    if st.button("⬅️ Back to FASTQ QC", key="summary_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("summary")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("📊 Summary Metrics")
        st.success(f"✅ Processed {data['total_reads']:,} reads totaling {data['total_bases']:,} bases successfully!")
        col1, col2, col3 = st.columns(3)
        col1.metric("📊 Total Reads", f"{data['total_reads']:,}")
        col2.metric("🧬 Avg Read Length", f"{data['avg_read_len']:.2f} bp")
        col3.metric("🧪 Avg GC Content", f"{data['avg_gc_content']:.2f}%")
        col4, col5 = st.columns(2)
        col4.metric("Bases Q ≥ 20", f"{data['pct_q20_bases']:.2f}%")
        col5.metric("Bases Q ≥ 30", f"{data['pct_q30_bases']:.2f}%")


elif st.session_state.current_section == "Read Length Distribution":
    st.title("📈 Read Length Distribution")
    if st.button("⬅️ Back to FASTQ QC", key="read_len_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("read_len")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("📈 Read Length Distribution")
        if data['lengths_hist']:
            fig1, ax1 = plt.subplots()
            sns.histplot(list(data['lengths_hist'].elements()), bins=min(30, data['max_read_len']+1), kde=False, ax=ax1, color='skyblue') # Use max_read_len for bins
            ax1.set_xlabel("Read Length (bp)")
            ax1.set_ylabel("Frequency")
            ax1.set_title("Distribution of Read Lengths")
            st.pyplot(fig1)
            buf = BytesIO()
            fig1.savefig(buf, format="png")
            st.download_button(
                label="Download plot as PNG",
                data=buf,
                file_name="read_length_distribution.png",
                mime="image/png"
            )
        else:
             st.info("No read length data to display.")


elif st.session_state.current_section == "GC Content Distribution":
    st.title("🧬 GC Content Distribution")
    if st.button("⬅️ Back to FASTQ QC", key="gc_dist_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("gc_dist")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("🧬 GC Content Distribution")
        if data['gc_hist']:
            # Use a more robust plotting method like Plotly for better interactivity and handling of sparse data
            gc_df = pd.DataFrame(list(data['gc_hist'].items()), columns=['GC_Percentage', 'Count']).sort_values(by='GC_Percentage')
            fig2 = px.bar(gc_df, x='GC_Percentage', y='Count', title='Distribution of GC Content per Read')
            fig2.update_layout(xaxis_title="GC Content (%)", yaxis_title="Number of Reads")
            st.plotly_chart(fig2, use_container_width=True)

            buf = BytesIO()
            fig2.write_image(buf, format="png")
            st.download_button(
                label="Download plot as PNG",
                data=buf,
                file_name="gc_content_distribution.png",
                mime="image/png"
            )
        else:
             st.info("No GC content data to display.")


elif st.session_state.current_section == "Per-Base Quality Scores":
    st.title("🔍 Per-Base Quality Scores")
    if st.button("⬅️ Back to FASTQ QC", key="pb_qual_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("pb_qual")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("🔍 Per-Base Quality Scores")
        if st.session_state.get('qual_sums') and data['total_reads'] > 0:
            # Convert qual_sums to a list ordered by position
            # Use data['max_read_len'] to ensure all positions are included up to the max read length
            per_base_mean_qual = [st.session_state.qual_sums.get(i, 0) / data['total_reads'] for i in range(data['max_read_len'])]

            if per_base_mean_qual:
                fig3 = px.line(y=per_base_mean_qual, labels={"x": "Base Position", "y": "Mean Quality Score"}, title="Mean Quality Score per Base Position")
                # Add quality thresholds
                fig3.add_shape(type="line", x0=0, y0=20, x1=len(per_base_mean_qual)-1, y1=20, line=dict(color="orange", width=2, dash="dash"), name="Q20")
                fig3.add_shape(type="line", x0=0, y0=30, x1=len(per_base_mean_qual)-1, y1=30, line=dict(color="red", width=2, dash="dash"), name="Q30")
                fig3.update_layout(showlegend=True)
                st.plotly_chart(fig3, use_container_width=True)

                buf = BytesIO()
                fig3.write_image(buf, format="png")
                st.download_button(
                    label="Download plot as PNG",
                    data=buf,
                    file_name="per_base_quality_scores.png",
                    mime="image/png"
                )
            else:
                 st.info("No per-base quality data to display.")
        else:
            st.warning("Quality score data not available or no reads processed.")

        st.subheader("🧩 Per-Base Nucleotide Composition")
        if data['base_counts']:
            # Find the maximum position across all bases
            max_pos = 0
            for base_counter in data['base_counts'].values():
                if base_counter:
                     max_pos = max(max_pos, max(base_counter.keys()))
            max_pos += 1 # Adjust to include the last position

            if max_pos > 0:
                # Convert base_counts counters to a DataFrame
                base_df_data = {base: [data['base_counts'][base].get(i, 0) for i in range(max_pos)] for base in data['base_counts']}
                base_df = pd.DataFrame(base_df_data)

                # Calculate percentage composition, handling rows with zero total counts
                row_sums = base_df.sum(axis=1)
                base_df_pct = base_df.div(row_sums, axis=0) * 100
                base_df_pct = base_df_pct.fillna(0) # Replace NaN (from division by zero) with 0

                base_df_pct = base_df_pct.reset_index().rename(columns={'index': 'Base Position'})
                base_df_long = base_df_pct.melt(id_vars='Base Position', var_name='Nucleotide', value_name='Percentage')

                fig_comp = px.line(base_df_long, x="Base Position", y="Percentage", color='Nucleotide', title="Per-Base Nucleotide Composition")
                st.plotly_chart(fig_comp, use_container_width=True)

                buf = BytesIO()
                fig_comp.write_image(buf, format="png")
                st.download_button(
                    label="Download composition plot as PNG",
                    data=buf,
                    file_name="per_base_nucleotide_composition.png",
                    mime="image/png"
                )
            else:
                 st.info("No per-base nucleotide composition data to display.")
        else:
             st.info("Nucleotide count data not available.")


elif st.session_state.current_section == "Sequence Duplication Levels":
    st.title("🔄 Sequence Duplication Levels")
    if st.button("⬅️ Back to FASTQ QC", key="dup_levels_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("dup_levels")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("🔄 Sequence Duplication Levels")
        dup_counts = pd.Series(list(data['seq_counter'].values())).value_counts().sort_index()
        if not dup_counts.empty:
            # Limit the number of bars for better visualization if many duplication levels exist
            max_bars = 50
            if len(dup_counts) > max_bars:
                 st.warning(f"Showing only the first {max_bars} duplication levels for clarity.")
                 dup_counts = dup_counts.head(max_bars)

            fig4 = px.bar(x=dup_counts.index, y=dup_counts.values, labels={'x': 'Number of Duplicates', 'y': 'Number of Sequences'}, title="Duplication Levels")
            st.plotly_chart(fig4, use_container_width=True)
            buf = BytesIO()
            fig4.write_image(buf, format="png")
            st.download_button(
                label="Download plot as PNG",
                data=buf,
                file_name="sequence_duplication_levels.png",
                mime="image/png"
            )

            st.subheader("🔥 Top 10 Most Frequent Sequences") # Increased to top 10
            top10 = data['seq_counter'].most_common(10)
            if top10:
                df_top10 = pd.DataFrame(top10, columns=["Sequence", "Count"])
                st.dataframe(df_top10)
            else:
                st.info("No frequent sequences found.")
        else:
            st.info("No sequence duplication data to display.")


elif st.session_state.current_section == "Adapter Content":
    st.title("🔗 Adapter Content (simple check)")
    if st.button("⬅️ Back to FASTQ QC", key="adapter_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("adapter")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("🔗 Adapter Content (simple check)")
        adapters = ["AGATCGGAAGAGC", "GATCGGAAGAGC"] # Common Illumina adapters
        adapter_counts_list = [(a, data['adapter_hits'].get(a, 0), f"{(data['adapter_hits'].get(a, 0) / data['total_reads'] * 100 if data['total_reads'] > 0 else 0):.2f}%") for a in adapters]
        df_adapt = pd.DataFrame(adapter_counts_list, columns=["Adapter Sequence", "Count", "% Reads Containing"])
        st.dataframe(df_adapt)
        st.caption("Note: This is a basic check for the presence of common Illumina adapters. For comprehensive adapter trimming and analysis, use dedicated tools like Cutadapt.")

# Additional FASTQ QC Sections

elif st.session_state.current_section == "Overrepresented Sequences":
    st.title("📌 Overrepresented Sequences")
    if st.button("⬅️ Back to FASTQ QC", key="overrep_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("overrep")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("📌 Overrepresented Sequences")
        top20 = data['seq_counter'].most_common(20)
        if top20:
            df_top20 = pd.DataFrame(top20, columns=["Sequence", "Count"])
            df_top20["% of Total Reads"] = df_top20["Count"] / data["total_reads"] * 100
            st.dataframe(df_top20)
        else:
            st.info("No overrepresented sequences found.")


elif st.session_state.current_section == "K-mer Content":
    st.title("🧩 K-mer Content")
    if st.button("⬅️ Back to FASTQ QC", key="kmer_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("kmer")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        k = 6
        st.subheader(f"🧩 K-mer Content (k={k})")
        kmer_counter = Counter()
        for seq in data['seq_counter']:
            count = data['seq_counter'][seq]
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_counter[kmer] += count

        top_kmers = kmer_counter.most_common(20)
        if top_kmers:
            df_kmer = pd.DataFrame(top_kmers, columns=["K-mer", "Count"])
            st.dataframe(df_kmer)
        else:
            st.info("No significant k-mer content found.")


elif st.session_state.current_section == "N Content per Base":
    st.title("❓ N Content per Base")
    if st.button("⬅️ Back to FASTQ QC", key="n_content_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("n_content")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("❓ N Content per Base")
        n_counts = data['base_counts'].get('N', Counter())
        if n_counts:
            max_pos = max(n_counts.keys()) + 1
            n_data = [n_counts.get(i, 0) / data['total_reads'] * 100 for i in range(max_pos)]
            fig_n = px.line(y=n_data, labels={"x": "Base Position", "y": "% N Bases"}, title="N Content per Base Position")
            st.plotly_chart(fig_n, use_container_width=True)
        else:
            st.info("No 'N' base content found across positions.")


elif st.session_state.current_section == "Base Composition Bias Heatmap":
    st.title("🎯 Base Composition Bias Heatmap")
    if st.button("⬅️ Back to FASTQ QC", key="base_bias_back_fastq"):
        navigate_to("FASTQ QC")
    st.markdown("---")
    show_fastq_uploader_and_process("base_bias")
    if st.session_state.fastq_data:
        data = st.session_state.fastq_data
        st.subheader("🎯 Base Composition Bias Heatmap")
        base_counts = data['base_counts']
        max_pos = max(max(v.keys()) for v in base_counts.values() if v) + 1

        base_df = pd.DataFrame({
            base: [base_counts[base].get(i, 0) for i in range(max_pos)]
            for base in base_counts
        })

        row_sums = base_df.sum(axis=1)
        base_pct_df = base_df.div(row_sums, axis=0).fillna(0) * 100
        fig_heat = px.imshow(base_pct_df.T, labels=dict(x="Base Position", y="Nucleotide", color="% Composition"),
                             aspect="auto", color_continuous_scale="RdBu")
        st.plotly_chart(fig_heat, use_container_width=True)


# --- User Guide ---
elif st.session_state.current_section == "User Guide":
    st.title("📖 User Guide")
    if st.button("⬅️ Back to Home", key="guide_back_home"):
        navigate_to("Home")
    st.markdown("---")
    st.markdown("""
    <div style="background-color:#e6f7ff; padding:25px; border-radius:12px; color: black;">

    <h2 style="color:#000000; margin-bottom:15px;">🛠️ <b>User Guide</b></h2>

    <p style="font-size:16px; color:#333333; line-height:1.6;">
      This interactive platform provides a variety of bioinformatics tools to help you analyze and visualize DNA, RNA, and protein sequences, as well as perform Quality Control on FASTQ files.
    </p>

    <h3 style="color:#000000; margin-top:20px;">📂 FASTA Analysis Tools</h3>
    <p style="font-size:16px; color:#333333;">Access these tools by clicking the "FASTA Analysis Tools" card on the Home page.</p>
    <ul style="font-size:16px; color:#333333; line-height:1.8; padding-left:20px; margin-top:10px;">
      <li>🧬 <b>GC Content Plot:</b> Upload FASTA or paste sequence. Adjust window size. Plot shows GC % along sequence.</li>
      <li>🧩 <b>Nucleotide Composition:</b> Paste DNA sequence. Pie chart shows A, T, G, C distribution.</li>
      <li>🎨 <b>Sequence Logo Generator:</b> Upload FASTA (unaligned or aligned) or paste sequences. Requires external Clustal Omega. Generates logo showing conservation.</li>
      <li>🔡 <b>Amino Acid Frequency:</b> Paste protein sequence. Bar chart shows amino acid counts.</li>
      <li>🌡️ <b>Gene Expression Heatmap:</b> Upload CSV (rows=genes, cols=samples, values=expression). Visualizes expression patterns.</li>
      <li>🔍 <b>Motif Search:</b> Paste DNA sequence and regex motif. Finds and lists motif occurrences.</li>
      <li>🔁 <b>Repeat Analysis:</b> Paste DNA sequence and select k-mer size. Lists frequent k-mers.</li>
      <li>🧪 <b>Primer Design:</b> Paste DNA sequence. Provides basic forward/reverse primer suggestions.</li>
      <li>🌀 <b>Circular Genome Plot:</b> Upload CSV (start, end, label). Plots features on a circle (requires genome length for accurate scale).</li>
    </ul>

    <h3 style="color:#000000; margin-top:20px;">📊 FASTQ Quality Control</h3>
     <p style="font-size:16px; color:#333333;">Access these reports by clicking the "FASTQ Quality Control" card on the Home page. For each report section, upload your FASTQ file (supports .fastq and .fastq.gz) and click "Process File and Generate Output". The app will remember the last processed file for the current session.</p>
    <ul style="font-size:16px; color:#333333; line-height:1.8; padding-left:20px; margin-top:10px;">
      <li>📊 <b>Summary Metrics:</b> Displays total reads, total bases, average length, average GC, and % bases with high quality scores (Q20, Q30).</li>
      <li>📈 <b>Read Length Distribution:</b> Histogram showing the frequency of different read lengths.</li>
      <li>🧬 <b>GC Content Distribution:</b> Histogram showing the distribution of GC percentage across reads.</li>
      <li>🔍 <b>Per-Base Quality Scores:</b> Line plot showing the average quality score at each position in the reads. Also includes per-base nucleotide composition.</li>
      <li>🔄 <b>Sequence Duplication Levels:</b> Bar chart showing how many sequences appear multiple times. Also lists the top 10 most frequent sequences.</li>
      <li>🔗 <b>Adapter Content:</b> Simple check for the presence of common Illumina adapter sequences.</li>
    </ul>

    <p style="font-size:16px; color:#333333; margin-top:20px;">
      <b>General Tips:</b>
      <ul>
          <li>Ensure your input data is in the correct format (FASTA, CSV, FASTQ).</li>
          <li>For large files, processing may take some time.</li>
          <li>Plots can often be downloaded directly from the plot area or via a download button.</li>
      </ul>
    </p>

    </div>
    """, unsafe_allow_html=True)

# --- About Creator ---
elif st.session_state.current_section == "About Creator":
    st.title("👥 About")
    if st.button("⬅️ Back to Home", key="about_back_home"):
        navigate_to("Home")
    st.markdown("---")

    # Author Section
    st.markdown("""
    <div class="about-section-container">
        <h3>👩‍🔬 About the Author</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src=https://media.licdn.com/dms/image/v2/D4D03AQFlS6DAPY4sKA/profile-displayphoto-shrink_400_400/B4DZV3ObW5HwAg-/0/1741462027428?e=1752105600&v=beta&t=JbaMb6m2xkknXjxJVlv0XV_Dq2ixrzl25riU2FYAgR0
                 class="profile-image">
            <div>
                <h4>Anuya Janrao </h4>
                <p>M.Sc. Bioinformatics Student at DES Pune University</p>
                <p>TMSc Bioinformatics candidate passionate about computational biology, rare disease research, and data visualization. I apply bioinformatics, machine learning, and network analysis to uncover gene-disease associations and therapeutic targets, with a focus on rare diseases. Skilled in analyzing high-throughput sequencing and multi-omics data. Experienced in building interactive dashboards using Power BI, Tableau, Streamlit, and Python/R. I aim to bridge biological insights with data-driven solutions for precision medicine and drug discovery.

</p>
                 <p>🔗 <a href="https://www.linkedin.com/in/arati-joshi-a24939264" target="_blank">Connect on LinkedIn</a></p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Mentorship Section
    st.markdown("""
    <div class="about-section-container">
        <h3>👨‍🏫 Mentorship</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src="https://media.licdn.com/dms/image/v2/D5603AQF9gsU7YBjWVg/profile-displayphoto-shrink_400_400/B56ZZI.WrdH0Ag-/0/1744981029051?e=1752105600&v=beta&t=F4QBDSEgjUvnBS00xPkKqPTLI0jQaMpYefaOzARY1Yg"
                 class="profile-image">
            <div>
                <h4>Dr. Kushagra Kashyap</h4>
                <p>Assistant Professor (Bioinformatics), Department of Life Sciences, School of Science and Mathematics, DES Pune University</p>
                <p>This project was developed under the guidance of Dr. Kashyap, who provided valuable insights and mentorship
                throughout the development process. His expertise in bioinformatics and computational biology was instrumental
                in shaping this project.</p>
                <p>🔗 <a href="https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb" target="_blank">Connect on LinkedIn</a></p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

 # Web Server Info Section
    st.markdown("""
    <div style="background-color: #f0f2f6; padding: 30px; border-radius: 15px; color: #333; box-shadow: 0 6px 12px rgba(0, 0, 0, 0.1);">
        <h3 style="color: #0056b3; margin-bottom: 20px; text-align: center;">
            <span style="font-size: 1.5em; margin-right: 10px;">✨</span> Web Server Information & Purpose
        </h3>
        <div style="display: flex; flex-direction: column; gap: 15px;">
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">Server Name:</span> Genome Sculptor
            </p>
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">Framework:</span> Streamlit
            </p>
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">Purpose:</span> Genome Sculptor is designed as an interactive and user-friendly platform for fundamental bioinformatics analysis and visualization. Its primary goal is to make common tasks in sequence analysis and quality control accessible to students, researchers, and anyone interested in exploring genomic data without requiring extensive command-line expertise.
            </p>
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">It provides tools for:</span>
            </p>
            <ul style="font-size: 1.1em; line-height: 1.8; padding-left: 30px;">
                <li>🧬 Analyzing and visualizing DNA, RNA, and protein sequences (FASTA analysis).</li>
                <li>🔬 Performing essential Quality Control checks on High-Throughput Sequencing data (FASTQ analysis).</li>
                <li>📊 Generating informative plots and reports for better understanding of biological data.</li>
            </ul>
            <p style="font-size: 1.1em; line-height: 1.6;">
                This server serves as an educational tool and a convenient resource for quick analyses in computational biology.
            </p>
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">Libraries Used:</span> Biopython, pandas, matplotlib, seaborn, plotly, logomaker, numpy, altair, fpdf, etc.
            </p>
            <p style="font-size: 1.1em; line-height: 1.6;">
                <span style="font-weight: bold; color: #007bff;">Source Code:</span> Currently available on GitHub (https://github.com/Anuyajanrao/Gene-sculptor/blob/main/app.py).
            </p>
        </div>
    </div>
    """, unsafe_allow_html=True)

# --- FAQ ---
elif st.session_state.current_section == "FAQ":
    st.title("❓ FAQ")
    if st.button("⬅️ Back to Home", key="faq_back_home"):
        navigate_to("Home")
    st.markdown("---")

    st.markdown("### 💡 Frequently Asked Questions")

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q1: How do I upload a file?</p>
        <p style="margin-left: 20px;">Use the <b>file uploader widget</b> available in each analysis section after selecting a tool or QC report.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q2: Is my data stored?</p>
        <p style="margin-left: 20px;"><b>No</b>. All processing happens in-memory during your session. Data is not stored on the server after you close the tab.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q3: Can I analyze protein sequences?</p>
        <p style="margin-left: 20px;">Yes, use the <b>Amino Acid Frequency</b> tool for protein composition. The <b>Sequence Logo Generator</b> can also work with protein alignments.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q4: What file formats are supported?</p>
        <p style="margin-left: 20px;">The app supports <b>FASTA</b> (.fasta, .fa, .fna), <b>FASTQ</b> (.fastq, .fq, .fastq.gz), and <b>CSV</b> (.csv) for specific tools like Gene Expression Heatmap and Circular Genome Plot.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q5: The Sequence Logo Generator gives an error about Clustal Omega. What does that mean?</p>
        <p style="margin-left: 20px;">This tool requires the external program <a href="http://www.clustal.org/omega/" target="_blank" style="color: #007bff; text-decoration: none;">Clustal Omega</a> to perform multiple sequence alignment. It needs to be installed on the system running the Streamlit app and accessible in the system's PATH. If you are running this locally, you might need to install it yourself.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q6: Why isn't the FASTQ analysis working for my file?</p>
        <p style="margin-left: 20px;">Ensure your file is a valid FASTQ format (plain text or gzip compressed). Very large files might take a long time to process, or could exceed memory limits depending on the server environment. The app currently supports .fastq, .fq, and .fastq.gz extensions.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa; margin-bottom: 20px;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q7: The PDF report is missing plots. Why?</p>
        <p style="margin-left: 20px;">Generating plots directly within the PDF is complex. The current PDF report provides a text summary of key metrics and top sequences. To get plots, you should download them as PNGs directly from the individual report sections.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div style="padding: 15px; border-radius: 10px; background-color: #f8f9fa;">
        <p style="font-weight: bold; color: #007bff;">➡️ Q8: How can I suggest new features or report bugs?</p>
        <p style="margin-left: 20px;">Please reach out to the creator via their LinkedIn profile mentioned in the 'About Creator' section.</p>
    </div>
    """, unsafe_allow_html=True)
