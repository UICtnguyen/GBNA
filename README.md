# Geodesic Distance Analysis for Functional Connectomes and Diffusion Weighted Imaging Data

GBNA is short for "Geodesic Brain Network Analysis" 

This pipeline is still under development 

## Summary 

This project performs the analysis of functional connectivity (FC) and diffusion-weighted imaging (DWI) data by computing geodesic distances between brain regions. The main tasks of the pipeline include:
- Preprocessing: Filter and prepare data by ensuring consistent subject lists across FC and DWI datasets.
- Geodesic Distance Calculation: Compute pairwise geodesic distances for both FC and DWI data after computing eigenvalues with respective methods for each.
- Visualization: Generate heatmaps for the geodesic distance matrices.
- Export: Save the subject names and geodesic distance matrices to text and CSV files for further analysis.

## Package folder includes

The project directory includes the following files:

- example.py ## Example script for importing and running the analysis pipeline

- GBNA.py ## Main analysis code (includes class definitions and methods)

- README.md             ## This README file

- After running the example.py file, three more files will populate in the folder corresponding to 1) participants in order of analysis 2) a .txt file of the similarity matrix for functional connectivity ready to be uploaded and inputted into the SNF package 3) a .txt file of the similarity matrix for structural connectivity ready to be uploaded and inputted into the SNF package


## Usage

1.	Prepare Data: Ensure your input data (both FC and DWI) are ready and organized in your working directory. The dataset should be formatted in a way that each subject’s data is accessible through file paths and all subjects are in one main folder for either FC or DWI.
2.	Configure the Script:
- In GBNA.py, adjust the file paths and ensure the correct format for your input data. 
- The compute_geodesic_distances method computes pairwise geodesic distances for the provided datasets.
3.	Run the Analysis: To execute the analysis, simply run the example.py script. 


# Funded/Supported by the CoNeCT Lab at the University of Illinois, Chicago
Developed by Theresa Nguyen, with gratitude to Paul Thomas for guidance
