# MPore
Mpore is a user-friendly method for  database-driven identification of active methyltransferases in prokaryotic genomes from Oxford Nanopore R10 sequencing data
Detailed information can be found in [Publication reference]

## Features 
- Basecalling supports POD5 for high performance modification calls with [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- Identifies candidate methyltransferases by homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Genome wide methylation signals is generated form Nanopore data 
- Activity assesment of candidate methyltransferases by L1 regularized logistic regression
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
For visualization purposes it is recommended to not use to long File_names, futhermore it is recommended to use  string elements without any whitespaces. 

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
export TSV_REBASE_data=TSV_REBASE_data.tsv
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
- TSV_REBASE_data includes methyltransferases their recognition sites and catalyzed methylation from [REBASE](http://rebase.neb.com/rebase/rebase.html) in an concatenated format 
- SPLIT toggles on a memory efficient workflow at the cost of runtime
- LOG_ANALYSIS toggles on MPores statistical modelling

It is recommended to toggle on LOG_ANALYSIS to activate MPores activity assesment for candidate methyltransferses. SPLIT should  be toggled on if the user is unsure about RAM capacity. 

3. **Run Command**
```bash
/usr/bin/time -v snakemake -s Snkemake_entire_thing \
  --conda-frontend conda --cores all --use-conda --resources dorado=1 \
  --rerun-incomplete \
  > snakemake_pipeline_output.txt 2> snakemake_resource_usage_erros.txt
```
- snakemake_pipeline_putput includes the standard output of the pipeline (echo, cat, executed snakemake rules)
- snakemake_resource_usage_errors includes the errors and resource usage for the current usage
It is also possible to use --Config as parameter to set the environment variables if one does not want to export them at the start

4. **Output**

The Output of MPore includes following files: 
- Bam files resulting from dorado basecalling and their corresponding bed files including site specific information (including number of modified reads, coverage etc.)
- All found CDS from PROKKA are included in the file_name directories 
- BlASTP results txt files including alignment results for the identified CDS against the REBASE data
- All_Isolates_gene_loci.csv presenting the Enzymes with an e_value < e-25 including methylases used for the downstream analysis 
- Beta_coef_p_values_{methyltype} where methyltype can be 4mC, 5mC or 6mA. Includes enzymes their beta coefficient estimate from the L1-regularized logistic regression and the origin gene loci for the enzymes
- Context_influence_{methyltype}.xlsx showing the influence of the flanking genomic context on the average genome methylation
- MTase_presence_e_25_values.csv showing identified enzymes and their corresponding e-value in isolates they were found in
- Sample_DF_{file_name}_{methyltype}.csv showing all analyzed motifs (User defined + REBASE motifs of identified enzymes) and their average methylation score in the Nanopore data 
- Sample_DF_detailed_{file_name}_{methyltype}.csv showing all individual sites for the analyzed motifs (User defined + REBASE motifs of identified enzymes) and the methylation score at each individual position
- plots_{methyltype} this directory includes all Boxplots comparing the nanopore methylation score across isolates with a data under the boxplots if enzymes have been identified which would be influcening the positions showed in the boxplots
- multipanel_plot_{file_name} the multipanel Plot for each isoalte shown in 5.
- heatmap_methylation_Score_{context}.png showing the overall methylation signal for motifs of identified methylases in a heatmap

5. **Visualisation**















   


   







