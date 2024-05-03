# Proteomics accelerator
A simple step-by-step guide to accelerate your proteomics study, even with no replicates. It is suitable for exploratory searches for potential proteins of interest as preliminary data.

## Introduction
This workflow is designed for proteomic analyses and includes only key steps.
Workflow Components
1. Data Imputation
2. Normalization
3. Enrichment Analysis
4. Visualization

## Package Citations
Please make sure to properly cite the packages and tools used in this workflow.

## Disclaimer
Proteomic analysis involves complexities that may require specialized knowledge.
This workflow is intended as a general guide, and users are strongly encouraged to consult domain experts for state-of-the-art approaches and detailed insights.

## Info
Created March 3rd, 2024 by Thomas Huang at home.
Contact: **r12633002@ntu.edu.tw**

## Prepare the data
Your data should contain at least:
1. Protein abundances (not normalised)
2. Unique peptide (Used for filtering out potential technical errors)
3. Protein accession (recommended: NCBI refseq, ENSEMBL, SYMBOL, ENTREZID, UNIPROTID; but any format is okay as long as you know what it is)

## Prepare R environment
Make sure you execute *0_setup.R* and library all required packages.
