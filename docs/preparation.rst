II. Preparation
===============

Input Data
----------

1. To analyze your spatial proteomics datasets, you need the proteomics report file(s) derived from your spectral search software, such as MaxQuant, Spectronaut, DIANN, or others. Your data must be reported as a **pivot report table**, meaning that your table includes one column per sample, as well as additional columns for further information. The necessary columns are:

   a. One column per sample (fraction).
   b. One column containing a **unique identifier** (e.g., protein groups, protein ID, etc.).
   c. One column containing key names that match the key names in your marker list (usually gene names). Ensure these keys are compatible, including case sensitivity.

2. Furthermore, you need a file containing your marker proteins. C-COMPASS provides prepared marker lists from previous publications, or you can use a custom export from a database relevant to your project. This file must include at least two columns:

   a. A column containing key names matching those in your dataset (usually gene names, see II.1.c).
   b. A column containing **class annotations** (for spatial proteomics experiments, this should represent the compartments where the marker proteins are located).

3. An additional dataset containing the total proteomes of the fractionation samples (proteomes derived from whole cell/tissue lysate) can be provided for **class-centric analysis** of compartments. This file should contain:

   a. One column per total proteome sample.
   b. One column containing the **same unique identifier** as used in the fractionation samples (see II.1.b).


Additional Notes
----------------

A) All input files must be **tab-delimited** (.tsv or .txt).  
B) If using an export file from **Perseus**, ensure that the file does not contain a second-layer header.  
C) Input datasets (for both fractionation and total proteome) can be stored in the same file or split across different files. If they are split, ensure that the **identifiers** are consistent.
