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


## External Tools

MPore relies on external tools that must be installed **before** running the pipeline.  
Please ensure that the following software is available on your system.

### Dorado

[Dorado](https://github.com/nanoporetech/dorado) is required for basecalling and methylation detection.

#### Verify installation

Make sure [Dorado](https://github.com/nanoporetech/dorado) is available in your `PATH`:

```bash
ls -l $(command -v dorado)
```
If this command does not return a valid path, Dorado is not correctly installed.
#### Model availability
If Dorado was installed using the self-contained release, ensure that the required methylation detection models are present in:
```bash
<dorado_directory>/lib
```
The following models must be available:
```bash
dna_r10.4.1_e8.2_400bps_hac@v5.0.0
dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA@v1
dna_r10.4.1_e8.2_400bps_hac@v5.0.0_4mC_5mC@v1
```
If u installed dorado as container you can locate your Dorado installation with:
```bash
find $HOME -type d -name "dorado-*"
```
#### Download models
To download all available models:
```bash
dorado download --model all
```
Alternatively, you can download only the required models:
```bash
dorado download --model \
dna_r10.4.1_e8.2_400bps_hac@v5.0.0 \
dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA@v1 \
dna_r10.4.1_e8.2_400bps_hac@v5.0.0_4mC_5mC@v1
```
For detailed installation instructions, please refer to the official Dorado documentation.

## Installation 
1. **Create a Conda environment**  
MPore is designed to run within a dedicated Conda environment.
   ```bash
   conda create -c conda-forge -c bioconda -n Bacterial_context1 snakemake=8.20.5
   conda activate Bacterial_context1
   snakemake -help
2. **Clone the MPore repository** 
   ```bash
   git clone [https://github.com/DiltheyLab/MPore.git]
   cd MPore

## Initialization
1. **Prepare the CSV file and user motif list**  
After installation, navigate to the `MPore` directory and create a CSV file containing the following columns:
- `File_name`
- `Reference_path`
- `pod5_path`
The `File_name` column is used as the isolate identifier throughout the downstream analysis.

For visualization and compatibility purposes, it is recommended to:
- avoid overly long `File_name` entries
- avoid whitespace characters
- use continuous, descriptive strings instead

An example input CSV file is shown below:

```csv
File_name,Reference_path,pod5_path
12256U,/home/azlan/Myco_Data/ref/12256U.fasta,/home/azlan/Myco_Data/pod5s/12256U
8958VA,/home/azlan/M_hominis/ref/8958VA.fasta,/home/azlan/Myco_Data/pod5s/8958VA
```
In this example, File_name corresponds to the isolates 12256U and 8958VA, each linked to their respective reference genome and POD5 directory.
#### Reference contig name formatting  
It is strongly recommended to verify and, if necessary, adjust contig names in the reference FASTA files.
Whitespace characters or additional annotations in contig headers can cause issues with several tools used in the pipeline.
A simple check of contig headers can be performed using:
```bash
grep ">" 12256U.fasta
```
Example output:
```bash
>12256U Mycoplasma hominis, complete genome
```
A properly formatted contig header should contain only a single identifier.  
An example of how to reformat the FASTA file is shown below:
```bash
zcat 12256U.fa.gz \ | awk '/^>/{print $1; next} {print}' \ > 12256U_formatted.fasta
```
Verification of the reformatted file:
```bash
grep ">" 12256U_formatted.fasta
```
Expected output:
```bash
>12256U 
```
#### Motif list
In addition to the input CSV file, the user must provide a text file containing DNA motifs of interest, with one motif per line.
If no specific motifs are of interest, a dummy motif file can be provided, for example:
```motif-list
GATC
```
In this case, `GATC` will be used as a placeholder motif.
all motifs associated with candidate methyltransferases will be included in the analysis in addition to the user-defined motifs.
#### Note: An empty motif file must not be used as input.

2. **Set up environment variables**

Next, configure the required paths and runtime options via the `config.yaml` file located in the MPore directory.
Navigate to the MPore folder and open the configuration file, for example using `nano`:
```bash
cd /MPore
nano config.yaml
```
The config.yaml file contains all variables required to run MPore and may look like the following:
```yaml
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
#### Configuration parameters
Below is an overview and explanation of each configuration parameter:

- `INPUT_CSV`  
  Path to the input CSV file created in Step 1
- `OUTPUT_DIR`  
  Directory where all MPore results will be written.
- `DORADO_PATH`  
  Path to the Dorado installation directory.
- `USER_MOTIF_LIST`  
  Text file containing user-defined DNA motifs of interest
- `TSV_DATA`  
[REBASE](http://rebase.neb.com/rebase/rebase.html)-derived file listing methyltransferases, their recognition sites, and associated methylation types (should not be modified).
- `TSV_REBASE_data`  
  Concatenated REBASE file containing methyltransferases, recognition sites, and methylation types (should not be modified).
- `SPLIT`  
  Enables a memory-efficient workflow at the cost of increased runtime. When enabled, intermediate and result files are split into smaller chunks.
- `LOG_ANALYSIS`  
  Enables MPore’s statistical modeling and activity assessment for candidate methyltransferases
- `MODE`  
  Analysis mode selection:
   - `mode: 1` — cross-isolate analysis (default), fitting a regularized logistic regression model across isolates
   - `mode: 2` — isolate-specific analysis, fitting a regularized logistic regression model per isolate  
- `heatmap`  
  Enables heatmap generation for summary figures (see Step 5)

#### Recommended settings
It is recommended to enable:  
`log_analysis: true`  
to activate MPore’s activity assessment for candidate methyltransferases.  
`split: true`  
if available RAM is limited or unknown, to reduce memory usage at the expense of longer runtime.

## Run the pipeline


MPore is executed using Snakemake.  
The optional prefix `/usr/bin/time -v -o snakemake_resource_usage.txt` records the total runtime and resource usage of the workflow.

This prefix can be removed if runtime and memory profiling is not required.  
For routine inspection and debugging, the log files located in `MPore/logs` are recommended, as they provide detailed information on when each individual pipeline step started and finished.

```bash
/usr/bin/time -v -o snakemake_resource_usage.txt \
snakemake -s Snakemake_entire_thing \
  --configfile config.yaml \
  --conda-frontend conda \
  --cores all \
  --use-conda \
  --resources dorado=1 \
  > snakemake_pipeline_out.txt \
  2> snakemake_pipeline_err.txt
```
#### Command options
- `-s Snakemake_entire_thing`  
  Specifies the main Snakemake workflow file
- `--configfile config.yaml`  
  Loads all configuration variables defined in Step 2
- `conda-fronted conda`  
  Uses Conda (instead of Mamba) to create and manage environments located in `MPore/Environments`
- `cores all`  
  Allows Snakemake to use all available CPU cores
- `use-conda`  
  Enables automatic creation and activation of rule-specific Conda environments
- `resources dorado =1`  
  Limits Dorado execution to a single job at a time. This is recommended for GPU- and memory-intensive basecalling steps
- `> snakemake_pipeline_out.txt`  
  Redirects standard output (e.g. executed rules, echo statements) to a file
- `2> snakemake_pipeline_err.txt`  
  Redirects standard error messages to a separate file for debugging

#### Logging and debugging
- Snakemake log files are written to the `MPore/logs` directory
- These logs provide detailed, rule-specific execution information and are highly recommended for debugging and performance assessment

## Output

The output generated by **MPore** consists of the following files and directories:

- **BAM and BED files**  
  Generated during Dorado basecalling  
  BED files contain site-specific methylation information, including coverage and the number of modified reads  

- **PROKKA annotations**  
  Predicted coding sequences (CDS) for each isolate  
  Annotation files are stored within the corresponding `file_name` output directories
  
- **BLASTP result files**  
  Text files containing BLASTP alignment results of predicted CDS against the REBASE database

- **All_isolates_gene_loci.csv**  
  Lists all enzymes with an e-value < `1e-25`, including candidate methyltransferases used for downstream analyses and their gene loci

- **Beta_coef_p_values_{methyltype}.csv**  
  Contains enzymes and their beta coefficient estimates derived from L1-regularized logistic regression `{methyltype}` can be one of `4mC`, `5mC`, or `6mA`    
  The file also includes the originating gene loci for each enzyme
  
- **Context_influence_{methyltype}.xlsx**  
  Summarizes the influence of flanking genomic context on average genome-wide methylation levels

- **MTase_presence_e_25_values.csv**  
  Provides an overview of identified methyltransferases (MTases) and their corresponding e-values across all analyzed isolates

- **Sample_DF_{file_name}_{methyltype}.csv**  
  Summarizes all analyzed motifs (both user-defined and REBASE-derived), including average methylation scores per motif

- **Sample_DF_detailed_{file_name}_{methyltype}.csv**  
  Contains per-site methylation scores for each analyzed motif

- **Plots_{methyltype}/**  
  Directory containing boxplots comparing methylation scores across isolates    
  Associated enzyme information is displayed beneath each plot when applicable

- **Multipanel_plot_{file_name}.png**  
  Combined visualization showing all relevant plots for a given isolate (see Section 5)

- **Heatmap_methylation_score_{context}.png**  
  Heatmap summarizing the global methylation signal across motifs of identified methyltransferases for a given genomic context

5. **Workflow and multipanel** 

<img width="9629" height="4605" alt="MPore_workflow_bigger3_ad_2" src="https://github.com/user-attachments/assets/87c21dc3-539c-45f4-b28e-530aac24d909" />


This figure illustrates the overall workflow of MPore and the data structure used for L1-regularized logistic regression.
The bar plot summarizes methyltransferase detection results for the benchmark dataset ([Link]).
In addition, a representative multipanel plot generated by MPore for Mycoplasma hominis is shown.
For a detailed description of the methodology and analysis workflow, please refer to the MPore application note ([Link]).
Contact: Azlan@uni-duesseldorf.de














   


   







