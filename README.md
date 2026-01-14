# MPore
MPore is a pipeline for database-driven identification and activity assessment of methyltransferases in prokaryotic genomes using Oxford Nanopore R10 sequencing data. 

MPore is designed to be user-friendly: it includes automated basecalling, motif extraction, statistical modeling, and visualization of methylation patterns. 

For detailed methodology and benchmarking, see [Publication reference, DOI].


## Features
- Performs basecalling from POD5 files using [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md) for high-accuracy modification calls
- Identifies candidate methyltransferases via homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Generates genome-wide methylation signals from Nanopore sequencing data
- Evaluates MTase activity using L1-regularized logistic regression
- Produces visualizations of enzymes, their recognition motifs, and site-specific methylation patterns


## External tools
These tools should be already downloaded to the system or should be downloaded before usage of MPore 
- [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- [PROKKA](https://github.com/tseemann/prokka)
- [BLASTP](https://github.com/blast-io/blast)


## Installation 
1. **Create a new Conda environment**
   ```bash
   conda create -c conda-forge -c bioconda -n Bacterial_context1 snakemake=8.20.5
   conda activate Bacterial_context1
   snakemake -help
2. **Clone the MPore respository** 
   ```bash
   git clone [https://github.com/DiltheyLab/MPore.git]
   cd MPore

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
An example to verify ur installation path could be the following: 
```bash
ls -l $(command -v dorado)
```
If u download dorado as self-contained release from [dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md) u should also verify that all methylation detection models are downloaded in the `dorado_directory/lib`. The models `dna_r10.4.1_e8.2_400bps_hac@v5.0.0_4mC_5mC@v1`,`dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA@v1` and `dna_r10.4.1_e8.2_400bps_hac@v5.0.0`should be included in ur models.
```bash
find $HOME -type d -name "dorado-*"
dorado download --model all 
```
Instead of downloading all models u can also download specific models 
```bash
dorado download --model <dna_r10.4.1_e8.2_400bps_hac@v5.0.0_4mC_5mC@v1> <dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA@v1> <dna_r10.4.1_e8.2_400bps_hac@v5.0.0>
```
To setup the variables move into the config.yaml file in the MPore directory. Inside the MPore directory u can open and adjust the config file for example with

```bash
cd /MPore
nano config.yaml
```
This config file includes all variables needed for MPore. And can look like this:
```
input_csv: Data_Test.csv
output_dir: /mnt/azlan/Nanomotif_data/Outpu
dorado_path: /home/azlan/Tools/dorado-0.8.0-linux-x64
user_motif_list: Motifs.txt
tsv_data: TSV_Enzyme.csv
tsv_rebase_data: TSV_REBASE_data.tsv
split: true
log_analysis: true
mode: 1
heatmap: true
```
Here is a outline and explanation for each variable being set for MPore:
- `INPUT_CSV`: the CSV file created in step 1
- `OUTPUT_DIR`: directory where the results will be saved
- `DORADO_PATH`: path to the Dorado installation
- `USER_MOTIF_LIST`: list of motifs of interest
- `TSV_DATA`: [REBASE](http://rebase.neb.com/rebase/rebase.html)-derived file listing methyltransferases, their recognition sites, and associated methylation types (should not be changed)
- `TSV_REBASE_data`: concatenated [REBASE](http://rebase.neb.com/rebase/rebase.html) file with methyltransferases, recognition sites, and methylation types (should not be changed)
- `SPLIT`: enables a memory-efficient workflow at the cost of longer runtime, while creating splitted results
- `LOG_ANALYSIS`: enables MPore’s statistical modeling
- `MODE=2`: initiates isolate-specific analysis, where a regularized logistic regression is fitted for each isolate (default mode without this variable is cross-isolate analysis)
- `heatmap`: enable heatmap plotting for figures like in step 5

It is recommended to toggle on LOG_ANALYSIS to activate MPores activity assesment for candidate methyltransferses. SPLIT should  be toggled on if the user is unsure about RAM capacity. 

3. **Run command**
The `/usr/bin/time -v -o snakemake_resource_usage.txt` is used to see how long snakemake runs and how much resources it used. The user can remove this informational line and look into the logs file found in `MPore/logs`. These logs file are helpful to see when each individual step started and ended. 
Also for debugging the log files should be advised since they give insight to detaled information.

```bash
/usr/bin/time -v -o snakemake_resource_usage.txt \
snakemake -s Snkemake_entire_thing \
  --configfile config.yaml \
  --conda-frontend conda \
  --cores all \
  --use-conda \
  --resources dorado=1 \
  > snakemake_pipeline_out.txt \
  2> snakemake_pipeline_err.txt
```

- `configfile` allows our pipeline to include all environmental variables set in step 2
- `conda-fronted conda` use conda (not mamba) to create and manage environments found in `/MPre/Enviornments`
- `cores all` allow snakemake to use all available CPU cores
- `resources dorado =1` limit dorado jobs runs to 1 at a time (useful fpr GPU and memory-intensive basecalling)
- `snakemake_pipeline_out` includes the standard output of the pipeline (echo, cat, executed snakemake rules)
- `snakemake_pipeline_err.txt` includes the errors for the current run


4. **Output**

The output of **MPore** includes the following files:  

- **BAM and BED files**  
  Generated from Dorado basecalling.  
  The BED files contain site-specific methylation information, including coverage and the number of modified reads.  

- **PROKKA annotations**  
  All predicted CDS (coding sequences) for each isolate are stored in the corresponding `file_name` directories.  

- **BLASTP result files**  
  Text files showing alignment results of the identified CDS against the REBASE database.  

- **All_isolates_gene_loci.csv**  
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

- **Plots_{methyltype}/ directory**  
  Contains boxplots comparing methylation scores across isolates.  
  Associated enzyme data is shown beneath each plot when applicable.  

- **Multipanel_plot_{file_name}.png**  
  A combined visualization showing the relevant plots (see section 5) for each isolate.  

- **Heatmap_methylation_score_{context}.png**  
  A heatmap summarizing the global methylation signal across motifs of identified methyltransferases for a given genomic context. 

5. **Workflow and multipanel** 

<img width="9629" height="4605" alt="MPore_workflow_bigger3_ad_2" src="https://github.com/user-attachments/assets/87c21dc3-539c-45f4-b28e-530aac24d909" />


In this image we show the general workflow of MPore and the datastructure used for L1-regularized logistic regression. 
The Barplot shows methylase findings for our benchmark dataset from [Link]. And a multipanel plot example for M.hominis created by MPore.  
For detailed information please look up our application note: [Link]
Contact: Azlan@uni-duesseldorf.de














   


   







