# Structural-Analysis-FES-Kinase
R pipeline for the structural characterization of FES kinase using AlphaFold models, pLDDT quality assessment, and molecular descriptor calculation.

# Automated Structural Characterization of FES Kinase

## Project Overview
This repository contains a specialized bioinformatics pipeline developed for the structural analysis of the **FES (Feline Sarcoma oncogene) Tyrosine Kinase**. [cite_start]Developed as a technical demonstration for the **ISPA-FINBA Computational Biochemistry Unit (Ref: 05 2026 PTAI-TECTITN2)**, this project bridges transcriptomics-based gene identification with 3D structural proteomics[cite: 14, 16].

[cite_start]The pipeline automates the assessment of protein stability and structural quality, utilizing high-resolution models to evaluate potential therapeutic targets[cite: 16].

## Key Features & Methodology
[cite_start]In alignment with the requirements for the **INVESTIGO program**, this project implements the following[cite: 16, 81]:

1.  **Structural Modeling**: Integration of high-resolution structural predictions from the **AlphaFold Protein Structure Database**.
2.  **Quality Assessment (pLDDT)**: Systematic analysis of **Predicted Local Distance Difference Test (pLDDT)** scores to identify functional domains and assess model reliability.
3.  **Molecular Descriptors**: 
    * **Radius of Gyration (Rg)**: Manual implementation in R for robust assessment of protein compactness and folding state.
    * **Secondary Structure Analysis**: Identification of catalytic regions within the kinase domain.
4.  [cite_start]**Pipeline Automation**: A complete workflow that processes raw PDB data into refined structural outputs and sequence alignments (FASTA)[cite: 16].
5.  [cite_start]**HPC & Linux Readiness**: Scripting architecture designed for non-interactive execution in **Linux/UNIX** environments and **High-Performance Computing (HPC)** clusters[cite: 44, 56].



## Technical Stack
* [cite_start]**Programming**: R (utilizing the `bio3d` library)[cite: 46, 58].
* **Structural Input**: AlphaFold2 PDB Models & AlphaMissense Pathogenicity Scores.
* [cite_start]**Automation**: Bash scripting for pipeline deployment in scientific computing environments[cite: 44].
* **Visualization**: High-resolution rendering of protein domains and residue flexibility.

## Repository Contents
* [cite_start]`FES_structure_pipeline.R`: Core automated analysis script[cite: 16].
* `FES_Structural_Confidence.png`: Visual profile of structural confidence (pLDDT) across the protein backbone.
* `FES_sequence.fasta`: Extracted amino acid sequence for downstream evolutionary analysis.
* [cite_start]`FES_CAlpha_Model.pdb`: Processed structural model (C-Alpha trace) optimized for molecular dynamics and structural alignment[cite: 16, 48].

## Execution (Linux/HPC Terminal)
To run the automated characterization:

```bash
# Clone repository
git clone [https://github.com/](https://github.com/)[Your-Anonymous-ID]/Structural-Analysis-FES.git

# Execute R-based pipeline
Rscript FES_structure_pipeline.R
