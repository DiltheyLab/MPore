# Tool1
Mpore is a user-friendly method for  database-driven identification of active methyltransferases in prokaryotic genomes from Oxford Nanopore R10 sequencing data
Detailed information can be found in [Publication reference]

## Features 
- Basecalling supports POD5 for high performance with [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- Identifies candidate methyltransferases by homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Genome wide methylation signals is inferred form Nanopore signal data 
- Activity assesment of candidate methyltransferases by statistical modelling 
- MPores creates visulizations that show identified candidate enzymes, their recognition motifs, and the methylation status at the relevant genomic positions

## Instalation 
1. **Create a new Conda environment**
   ```bash
   conda create -n Bacterial_context1
   conda actiavte Bacterial_context1
2. **Download the Snakemake_Tool with** 
   ```bash
   git clone https://github.com/AzlanNI/Tool1.git
   cd Tool1

## Initialization
1. **Setup CSV File**
   
After Instalation move into the Ordner1 Folder and create a CSV file containing the columns File_name, Reference_Path and Pod5_path. The first column is used as name in the downstream analysis for a isolate. And example for a CSV file would be 

```csv
File_name,Reference_path,pod5_path
12256U,/home/azlan/Myco_Data/ref/12256U.fasta,/home/azlan/Myco_Data/pod5s/12256U
8958VA,/home/azlan/M_hominis/ref/8958VA.fasta,/home/azlan/Myco_Data/pod5s/8958VA
```
Here File_names are the Isolates 12256U and 8958VA with their responding Reference and Pod5 paths. 





   


   







