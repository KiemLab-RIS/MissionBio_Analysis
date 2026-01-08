# MissioBio_Analysis
This repository contains the scripts to analyze the MissionBio DNA + Protein sequencing data to quantify gene-editing and assign cell types based on protein antibody.

The analysis involves 2 steps:
1. Analysing the fastq files. The C++ script mbFilter and all the associated .cpp, .hpp, .o files are used for this step.
The C++ code mbFilter reads the MissionBio DNA and Protein fastq files, identifies the MissionBio fixed sequences in the fastq files and also uses
a primer file (primers3.tsv) for identifying start primer, stop primer and whole amplicon.
Special code is in place for looking at gene editing as well as EGFP expression barcodes. Protein reads are used to identify cell type.
Output is created in a excel-like tsv file and a long csv file.

A batch script(.sh) is created for each sample to run the mbFilter code. A meta_file.tsv is required for every sample which lists the DNA and Protein fastq files.
mbFilter primers3.tsv output_file_name meta_file.tsv whitelist.csv

2. Generating a report for each sample which assigns celltype and calculates the editing for each cell. The report.py file is used for this step.
   python report.py out_table.tsv(generated in the previous step.)
