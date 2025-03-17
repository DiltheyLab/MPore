# Tool1
Mpore is a user-friendly method for  database-driven identification of active methyltransferases in prokaryotic genomes from Oxford Nanopore R10 sequencing data
Detailed information can be found in [Publication reference]

## Features 
- Basecalling supports POD5 for high performance with [Dorado](https://github.com/nanoporetech/dorado/blob/release-v0.9/README.md)
- Identifies candidate methyltransferases by homology search against [REBASE](http://rebase.neb.com/rebase/rebase.html) using [PROKKA](https://github.com/tseemann/prokka) and [BLASTP](https://github.com/blast-io/blast)
- Genome wide methylation signals is inferred form Nanopore signal data 
- Activity assesment of candidate methyltransferases by statistical modelling 
- MPores creates visulizations that show identified candidate enzymes, their recognition motifs, and the methylation status at the relevant genomic positions







