# syntenyPlotteR
R package to draw synteny plots in three different styles

## To install:
`devtools::install_github("marta-fb/syntenyPlotteR")`


## Input file

Please provide a file with the synteny blocks information, following this format, separated by tabs.  

  DO NOT include a header line
  
  ref - name of reference species  
  chr - chromosome of reference species -- it ONLY accepts numbers  
  start - start of block in reference species  
  end - end of block in reference species  
  tarChr - chromosome of target species  
  tarSt - start of block in target species  
  orient - orientation of the block (1 or -1)  
  tar - name of target species  

## Evolution Highway style

### Usage

`draw.eh("input file","output file name",chromosomeRange)`


outputName - name of output file   
chromosomeRange - range of chromosome numbers in reference species.


### Example:  
`library(syntenyPlotteR)`

`draw.eh("testData.txt","testOut",1:5)`

## inferCARs style

### Usage

`library(syntenyPlotteR)`
`draw.ideogram(<synteny_file>, <target_chr_size>, <reference_chr_size>)`

  target_chr_size Lenghts of target chromosomes, scaffolds or contigs. Format: ID Length  
  reference_chr_size Lenghts of reference chromosomes, scaffolds or contigs. Format: ID Length  

### Example output:

![alt text](https://github.com/marta-fb/syntenyPlotteR/images/exampleInferCARs.png?raw=true)

## To come soon:

1. Update vignette
2. Alluvial style synteny plots
3. inferCARs style synteny plots
4. Preprint in Biorxiv

### Citation:
While I'm preparing the publication, please cite:  
Farré M, Kim J, et al. Evolution of gene regulation in ruminants differs between evolutionary breakpoint regions and homologous synteny blocks. Genome Research 2019 Apr;29(4):576-589
