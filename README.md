# spacer_minery.py, v_0.5 
# Detection of CRISPR spacers and protospacers

The aim of this project is to provide an all-in-one pipeline for the detection and extraction of spacer sequences by detecting the CRISPR repeat unit, then automatically search for cognate protospacers using BLAST to provide information on targets, sequence requirements, and potential PAM sequences.

# About

Author: Jake Bourgeois

Affiliation: Tufts Graduate School of Biomedical Sciences, Camilli Laboratory

---------------------------------

# Requirements

- Python 3.7+ (May be compatible with older versions, but this was coded in 3.7)
- BioPython
- BLAST command line tools, requires blastN (https://www.ncbi.nlm.nih.gov/books/NBK279671/)
- A working internet connection

---------------------------------

# Installation

Simply download the repository and extract to a directory.

The folder should contain:
- spacer_miner.py
- spacer_miner_creds.txt

Before running the script, put your email and api_key for [Entrez NCBI](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) in the spacer_miner_creds.txt file. If you do not, you may encounter errors, slowdowns, and even IP bans for repeated use of NCBI E-utilities.

---------------------------------

# Running the Script

spacer_miner.py uses a BLAST XML of the result from BLASTing a desired repeat sequence against a database. After arriving at the result page, select the sequence alignments you want to process, then export to XML. Currently, the script only extracts spacers in between perfect repeats to ensure spacer integrity, so you do not need to manually screen for repeat polymorphisms or other errors in this result.

After obtaining the BLAST XML, save the file in the same directory as spacer_miner.py.

The command is as follows:

`python3 spacer_miner.py BLAST_XML REPEAT_SEQ -p pam_length -c query_coverage_threshold -s query_similarity_theshold' -o output_prefix`

where the BLAST_XML is the path to the BLAST report described above, and REPEAT_SEQ is the repeat unit used in the BLAST search and downstream analysis.

Optional arguments:

| Argument | Description |
| --- | --- |
| -p, --pam_length | The number of nucleotides directly adjacent to the detected protospacer to extract. Takes from both 5' and 3' direction |
| -c, --query_coverage_threshold | The minimum coverage a protospacer must cover over the repeat sequenece. For example, 96% coverage over a 33bp spacer sequence allows for a single nucleotide mismatch on the terminal ends |
| -s, --query_similarity_threshold | The minimum similarity a protospacer must be identical to the repeat sequence. |
| -o, --output_prefix | Output prefix for analysis files |


Note: The largest bottleneck in the script is BLASTing over the web to detect putative protospacers in non-redundant. Please be patient while the script finishes; it could take several days for large amounts of spacers.

---------------------------------

# Output

spacer_miner.py outputs the following files in tab-delimited format:

**crispr_hit_overview.tsv**: Overview of accessions harboring the repeat sequence. Number of validated repeats refers to the number of repeats matching perfectly to the canonical 28 bp repeat sequence.

**spacer_overview.tsv**: Overview of mined spacers. Spacer number is oriented in the order of acquisition, ie. the most recently acquired spacer is the largest number.

**crispr_repeat_mined_data.tsv**: Overview of mined protospacers detected by BLASTing mined spacers to the NCBI non-redundant database.

**spacer_histogram.tsv**: Frequency of unique mined spacers across the samples

**potential_new_crispr_arrays.tsv**: A quick list of novel accessions detected during protospacer detection that may contain even more novel spacers. A future version may loop the script to recursively search these novel sources to even further expand the database.

**perfect_crispr_repeat_mined_data.tsv**: A subset of crispr_repeat_mined_data.tsv that only contains hits of 100% coverage and identity.

---------------------------------

Please direct questions and comments to my email, jacob.bourgeois@tufts.edu.

If you use this script and find it useful in your analysis, please cite us at: 
Bourgeois J, Lazinski DW, Camilli A. 2020. Identification of spacer and protospacer sequence requirements in the Vibrio cholerae type I-E CRISPR/Cas system. mSphere 5:e00813-20. https://doi.org/10.1128/mSphere.00813-20.

Thank you and happy hunting!

