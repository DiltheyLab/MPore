# Tool1
Mpore is a user-friendly method for  database-driven identification of active methyltransferases in prokaryotic genomes from Oxford Nanopore R10 sequencing data
Detailed information can be found in [Publication reference]

## Features 
- Basecalling supports POD5 for high performance with [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- Identifies candidate methyltransferases by homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Genome wide methylation signals is generated form Nanopore signal data 
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
   
After Instalation move into the Ordner1 Folder and create a CSV file containing the columns File_name, Reference_Path and Pod5_path. The first column is used as name in the downstream analysis for a isolate. An example for a CSV file would be 

```csv
File_name,Reference_path,pod5_path
12256U,/home/azlan/Myco_Data/ref/12256U.fasta,/home/azlan/Myco_Data/pod5s/12256U
8958VA,/home/azlan/M_hominis/ref/8958VA.fasta,/home/azlan/Myco_Data/pod5s/8958VA
```
Here File_names are the Isolates 12256U and 8958VA with their responding Reference and Pod5 paths. 

2. **Setup Environment variables**

Now setup variables for the command with ur paths and directories. Before Setting the variables up download [dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md) and check its directory.

```bash
export INPUT_CSV=Data_Test_Myco
export OUTPUT_DIR=/home/azlan/Myco_Data/Output
export DORADO_PATH=/home/azlan/Tools/dorado-0.8.0-linux-x64
export USER_MOTIF_LIST=Motifs_6mA.txt
export INCLUDE_REBASE_MOTIFS=true
export TSV_data=TSV_Enzyme.csv
export SPLIT=true
export Log_analysis=True
export REBASE_Motifs=TSV_REBASE_data.tsv
```
- INPUT_CSV is the csv file created in step 1
- OUTPUT_DIR is the path where the results should be saved to
- DORADO_PATH is the directory in which dorado is saved
- USER_MOTIF_LIST list of Motifs of interest
- INCLUDE_REBASE includes motifs and findings from REBASE next t the Motifs of interest
- TSV_DATA includes methylatransferases their recognition sites and catalyzed methylation from [REBASE](http://rebase.neb.com/rebase/rebase.html)
- SPLIT toggles on a memory efficient workflow at the cost of runtime
- LOG_ANALYSIS toggles on MPores statistical modelling

It is recommended to toggle on LOG_ANALYSIS to activate MPores activity assesment for candidate methyltransferses. SPLIT should also be toggled on if the user is unsure about RAM capacity. 

3. **Run Command**
```bash
Snakemake Command
```
















   


   







