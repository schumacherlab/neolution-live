# Neolution-live
## Pipeline for neo-antigen prediction

### Introduction

This pipeline performs live predictions using netMHCpan and netChop to predict peptide affinity and proteasomal processing. RNA expression data and 'similarity-to-self' filters can used to further increase the precision of the predictions. Recently, a random forest model trained on mass spectrometry data was added, which integrates predicted affinity rank, proteasomal processing scores and RNA expression values to a combined model score. This score (0-1) gives the probability of peptide presentation and has been shown to substantially increase the precision and sensitivity over the conventional use of binary cutoff values for the various predicted parameters. This is now the preferred method. Suggested model score cutoff values are:

* 0.01 for TIL screens (more inclusive; we picked up low magnitude TIL hits at low prob scores)
* 0.02 for PBMC screens (more stringent; unlikely to pick up low magnitude responses with low prob scores in peripheral blood)

As input, the pipeline expects a tsv file with affected germline and tumor transcripts. Additional variant/transcript information can be provided and will be transferred into the output. See **Input file format** paragraph for more information regarding input file generation.

### Minimal usage example

`Rscript performPredictions.R -f /path/to/variants.tsv -m A0201 -d 0.02 -l 9`

The call should be run from the neolution script directory from the Terminal and will start neo-antigen predictions for **variants.tsv**, __HLA-A*02:01__, **9-mer** peptides, applying a model prediction cutoff score of 0.02. 

**By default, netMHCpan v4.0 will be used and 1/6th of the available cores are used per run for parallel computations.  
Make sure you *nice* your runs and don't exceed max. HPC load (max. load = # cores)!!**

**IMPORTANT:** For additional information regarding the required commandline arguments, read segment below.

### Commandline arguments

#### Required

1. full input file path
2. hla/mhc type (e.g. A0201)
3. model prediction or affinity(rank) cutoff value
4. peptide length (e.g. 9) 

#### Optional

1. netChop processing cutoff
2. RNA expression cutoff
3. single sequence input (fasta input: not paired tumor-normal, no rna expression)
4. structural variant predictions
5. simple self-similarity check (9-, 10-, 11-mers)
6. extended self-similarity check (9-mers only)
7. use self-epitope list
8. use database for peptide affinity lookups (9-mers, netMHCpan-2.4 only)
9. netMHCpan version

**NOTE: self-similarity checking requires predicted self-epitope lists of matching HLA & peptide length**

### Input file format

Required input format is a wide, tab-separated table with the following columns. [Neolution-prep](https://gitlab.nki.nl/l.fanchi/neolution-prep) generates all required data and provides required input file.

#### Variant

| variant\_id | chromosome | start\_position | end\_position | variant\_strand | ref\_allele | alt\_allele |
|-------------|------------|-----------------|---------------|-----------------|-------------|-------------|

#### Gene/Transcript

| gene\_id | transcript\_id | transcript\_strand | hugo\_symbol | variant\_classification | transcript\_remark | transcript\_extension | nmd\_status | nmd\_remark |
|----------|----------|---------|----------|-----------|---------|------------|-------------|-------------|

#### DNA/RNA

| dna\_ref\_read\_count | dna\_alt\_read\_count | dna\_total\_read\_count | dna\_vaf | rna\_ref\_read\_count | rna\_alt\_read\_count | rna\_total\_read\_count | rna\_vaf | rna\_alt\_expression | rna\_expression |
|----------|-----------|------------|------------|----------|-----------|-----------|------------|--------|--------|

#### Effect

| codon\_germline | codon\_tumor | aa\_germline | aa\_tumor | aa\_pos\_germline | aa\_pos\_tumor\_start | aa\_pos\_tumor\_stop | peptidecontextnormal | peptidecontexttumor |
|----------|-----------|------------|-------------|-----------|------------|-----------|----------|--------|

---

`Rscript performPredictions.R --help`  

**Usage: performPredictions.R [OPTIONS]**

**OPTIONS**  
`-f FILE, --file=FILE`  
*Full path to file containing variant calls (required)*

`-m MHC, --mhc=MHC`  
*MHC/HLA allele, formatted as follows: A0201 (required)*

`-l LENGTH, --length=LENGTH`  
*Peptide length (required)*

`-d MODEL, --model=MODEL`  
*Random forest model score cutoff (optional, suggested cutoff 0.02 for PBMC, 0.01 for TIL screens)*

`-r RANK, --rank=RANK`  
*netMHCpan rank cutoff (optional, suggested is 2.8 ~500nM for A\*02:01)*

`-a AFFINITY, --affinity=AFFINITY`  
*netMHCpan affinity cutoff (optional)*

`-p PROCESSING, --processing=PROCESSING`  
*netChop processing score cutoff (optional, default: >= 0.5)*

`-e EXPRESSION, --expression=EXPRESSION`  
*RNA expression cutoff (optional, default: > 0; use -1 for no filtering)*

`--single`  
*Single sequence predictions (not paired normal-tumor) (optional, default: FALSE)*

`--structural`  
*Structural variant predictions (optional, default: FALSE)*

`--selfsim`  
*Perform simple self-similarity check; compatible with 9-, 10-, 11-mers (optional, default: FALSE)*

`--extselfsim`  
*Perform extended self-similarity check; only compatible with 9-mers (optional, default: FALSE)*

`--selflist`  
*Add predicted self-epitopes to self-similarity check, requires length- & HLA-matched selflist (optional, default: FALSE)*

`--fasdb`  
*Look up peptide affinity in FASdb, predict if not found; only compatible with 9-mers & netMHCpan-2.4 (optional, default: FALSE)*

`--panversion`  
*Use different version of netMHCpan; must be installed in path specified in runConfig.R (optional, default: 4.0)*

`--verbose`  
*Be chatty (optional, default: FALSE)*

`-h, --help`  
*Show this help message and exit*
