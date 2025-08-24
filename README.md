# MPore
MPore is a user-friendly method for database-driven identification of active methyltransferases in prokaryotic genomes using Oxford Nanopore R10 sequencing data.
Detailed information can be found in [Publication reference].


## Features 
- Basecalling supports POD5 for high-performance modification calls with [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- Identifies candidate methyltransferases by homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Generates genome-wide methylation signals from Nanopore data
- Assesses the activity of candidate methyltransferases using L1-regularized logistic regression
- MPore creates visualizations showing identified candidate enzymes, their recognition motifs, and the methylation status at the relevant genomic positions.

## Instalation 
1. **Create a new Conda environment**
   ```bash
   conda create -n Bacterial_context1
   conda actiavte Bacterial_context1
2. **Download the Snakemake_Tool** 
   ```bash
   git clone https://github.com/AzlanNI/Tool1.git
   cd Tool1

## Initialization
1. **Prepare the CSV file and user motif list**
   
After installation, move into the `MPore` Folder and create a CSV file containing the columns `File_name`, `Reference_Path` and `Pod5_path`. The first column is used as name in the downstream analysis for a isolate. An example for a CSV file would be
For visualization purposes, it is recommended not to use overly long `File_names` entries. Also, avoid whitespace characters and instead use continuous strings

```csv
File_name,Reference_path,pod5_path
12256U,/home/azlan/Myco_Data/ref/12256U.fasta,/home/azlan/Myco_Data/pod5s/12256U
8958VA,/home/azlan/M_hominis/ref/8958VA.fasta,/home/azlan/Myco_Data/pod5s/8958VA
```
Here, `File_name` corresponds to the isolates 12256U and 8958VA with their respective reference and POD5 paths
In addition to this CSV file, the user should also provide a text file containing motifs of interest.
If no motifs are of interest, simply provide a text file with the following format:

```motif-list
GATC
```
In this case, `GATC` will be used as a dummy motif.
By setting `INCLUDE_REBASE_MOTIFS=true`, all motifs of the candidate methyltransferases will be considered in addition to `GATC`.
*An empty motif file should not be used as input.*

2. **Setup environment variables**

Now, set up the required variables with your paths and directories.
Before doing so, download [dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md) and verify its installation path.

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
export MODE=2
```
- `INPUT_CSV`: the CSV file created in step 1
- `OUTPUT_DIR`: directory where the results will be saved
- `DORADO_PATH`: path to the Dorado installation
- `USER_MOTIF_LIST`: list of motifs of interest
- `INCLUDE_REBASE`: includes motifs and results from REBASE in addition to the motifs of interest
- `TSV_DATA`: [REBASE](http://rebase.neb.com/rebase/rebase.html)-derived file listing methyltransferases, their recognition sites, and associated methylation types 
- `TSV_REBASE_data`: concatenated [REBASE](http://rebase.neb.com/rebase/rebase.html) file with methyltransferases, recognition sites, and methylation types
- `SPLIT`: enables a memory-efficient workflow at the cost of longer runtime, while creating splitted results
- `LOG_ANALYSIS`: enables MPore’s statistical modeling
- `MODE=2`: initiates isolate-specific analysis, where a regularized logistic regression is fitted for each isolate (default mode without this variable is cross-isolate analysis)

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

The output of **MPore** includes the following files:  

- **BAM and BED files**  
  Generated from Dorado basecalling.  
  The BED files contain site-specific methylation information, including coverage and the number of modified reads.  

- **PROKKA annotations**  
  All predicted CDS (coding sequences) for each isolate are stored in the corresponding `file_name` directories.  

- **BLASTP result files**  
  Text files showing alignment results of the identified CDS against the REBASE database.  

- **All_Isolates_gene_loci.csv**  
  Contains all enzymes with an e-value < 1e-25, including methyltransferases used for downstream analyses and their gene loci of origin.  

- **Beta_coef_p_values_{methyltype}.csv**  
  Lists enzymes with their beta coefficient estimates from L1-regularized logistic regression.  
  `{methyltype}` can be `4mC`, `5mC`, or `6mA`.  
  The file also includes the origin gene loci for each enzyme.  

- **Context_influence_{methyltype}.xlsx**  
  Shows the influence of the flanking genomic context on the average genome-wide methylation.  

- **MTase_presence_e_25_values.csv**  
  Summarizes identified methyltransferases (MTases) and their corresponding e-values across isolates.  

- **Sample_DF_{file_name}_{methyltype}.csv**  
  Summarizes all analyzed motifs (both user-defined and REBASE-derived) with their average methylation scores per motif.  

- **Sample_DF_detailed_{file_name}_{methyltype}.csv**  
  Provides per-site methylation scores for each analyzed motif.  

- **plots_{methyltype}/ directory**  
  Contains boxplots comparing methylation scores across isolates.  
  Associated enzyme data is shown beneath each plot when applicable.  

- **multipanel_plot_{file_name}.png**  
  A combined visualization showing the relevant plots (see section 5) for each isolate.  

- **heatmap_methylation_Score_{context}.png**  
  A heatmap summarizing the global methylation signal across motifs of identified methyltransferases for a given genomic context. 

5. **Workflow and Multipanel** 

<img width="9629" height="4605" alt="MPore_workflow_bigger3_ad_2" src="https://github.com/user-attachments/assets/87c21dc3-539c-45f4-b28e-530aac24d909" />


In this image we show the general workflow of MPore and the datastructure used for L1-regularized logistic regression. 
The Barplot shows methylase findings for our benchmark dataset from [Link]. And a multipanel plot example for M.hominis created by MPore.  
For detailed information please look up our application note: [Link]
Contact: Azlan@uni-duesseldorf.de














   


   







