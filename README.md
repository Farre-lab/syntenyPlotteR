# syntenyPlotteR
R package to draw synteny plots in three different styles

This package has been developed by Joana Damas (joanadamas@gmail.com) and Marta Farré (m.farre-belmonte@kent.ac.uk)

## To install:
`devtools::install_github("marta-fb/syntenyPlotteR")`


## Input file

Please provide a file with the synteny blocks information, following this format, separated by tabs.  

  DO NOT include a header line
  
* chr - chromosome of reference species -- it ONLY accepts numbers  
* start - start of block in reference species  
* end - end of block in reference species  
* tarChr - chromosome of target species  
* tarSt - start of block in target species 
* tarEnd - end of block in target species 
* orient - orientation of the block (1 or -1) 
* tarSpecies - name of the target species   

## Evolution Highway style

### Usage

`draw.eh("input file","output file name",chromosomeRange)`


* outputName - name of output file   
* chromosomeRange - range of chromosome numbers in reference species.


### Example:  
`library(syntenyPlotteR)`

`draw.eh("syntenicBlocks_oneChromosome_severalSps.txt","testOut",1:1)`

![alt text](https://github.com/marta-fb/syntenyPlotteR/blob/master/images/exampleEH.png?raw=true)  

## inferCARs style


### Usage

`library(syntenyPlotteR)`  
`draw.ideogram("syntenicBlocks.txt", "reference_chr_size", "target_chr_size")`

* reference_chr_size Lenghts of reference chromosomes, scaffolds or contigs. Format: ID Length  
* target_chr_size Lenghts of target chromosomes, scaffolds or contigs. Format: ID Length  


### Example output:

![alt text](https://github.com/marta-fb/syntenyPlotteR/blob/master/images/exampleInferCARs.png?raw=true)


## Linear pairwise style

### Usage

`library(syntenyPlotteR)`  
`draw.pairwise("syntenicBlocks.txt", "outputFileName", "reference_chr_size", "target_chr_size", "refName", "tarName")`

* outputFilename: Name of the output file
* target_chr_size: Lenghts of target chromosomes, scaffolds or contigs. Format: ID Length  
* reference_chr_size: Lenghts of reference chromosomes, scaffolds or contigs. Format: ID Length  
* refName: name of the reference genome
* tarName: name of the target genome

### Example output:

![alt text](https://github.com/marta-fb/syntenyPlotteR/blob/master/images/exampleLinearPair.jpg?raw=true)


## To come soon:

1. Update vignette
2. Preprint in Biorxiv

### Citation:
While we're preparing the publication, please cite:  
Farré M, Kim J, et al. Evolution of gene regulation in ruminants differs between evolutionary breakpoint regions and homologous synteny blocks. Genome Research 2019 Apr;29(4):576-589
