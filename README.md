# Structural Analysis of FES Kinase

This project provides an automated R script to extract and analyze structural features from **AlphaFold2** protein models. It is designed to evaluate model confidence, protein compactness, and residue-level flexibility.

## Features

The pipeline performs the following tasks:

1.  **Sequence Extraction**: Exports the protein sequence to FASTA format.
2.  **Quality Assessment**: Calculates the average **pLDDT** score to measure the overall confidence of the predicted model.
3.  **Protein Compactness**: Implements a manual calculation of the **Radius of Gyration (Rg)** to evaluate the protein's folding state.
4.  **Flexibility Analysis**: Generates a structural confidence plot (pLDDT profile) using C-Alpha atoms.
5.  **Model Simplification**: Exports a refined PDB file containing only C-Alpha atoms for downstream analysis.



## Requirements

* **R**
* **bio3d** library (`install.packages("bio3d")`)

## How to Use

1.  Ensure your AlphaFold2 PDB file is in the project directory.
2.  Run the script in R:
    ```R
    source("FES_analysis.R")
    ```



## Output Files

* `FES_sequence_native.fasta`: Primary protein sequence.
* `FES_Structural_Confidence.png`: Plot showing the confidence score per residue.
* `FES_CAlpha_Model.pdb`: Simplified structural model (C-Alpha only).
