# syntenyPlotteR
R package to draw synteny plots in three different styles

This package has been developed by Joana Damas (joanadamas@gmail.com), Sarah Quigley (slq4@kent.ac.uk), Denis Larkin (dlarkin@rvc.ac.uk) and Marta Farré (m.farre-belmonte@kent.ac.uk)

It has been tested in Windows 10 R v4.0.3, Windows 11 with R v4.1.1 and v4.2.2 and MAC OS Ventura 13.1

## Table of contents

1. [To install](https://github.com/Farre-lab/syntenyPlotteR#to-install)
2. [Input files](https://github.com/Farre-lab/syntenyPlotteR#Input-files)
3. [Reformatting alignment data](https://github.com/Farre-lab/syntenyPlotteR#reformatting-alignment-data)
5. [Evolution Highway style](https://github.com/Farre-lab/syntenyPlotteR#evolution-highway-style)
6. [Chromosome painting style](https://github.com/Farre-lab/syntenyPlotteR#chromosome-painting-style)
7. [Linear style](https://github.com/Farre-lab/syntenyPlotteR#linear-style)

## To install:
`devtools::install_github("Farre-lab/syntenyPlotteR")`

### Package Requirements
This package requires the package ggplot2 to be able to run this must be installed and loaded in the R environment

This can be done by:

`install.packages("ggplot2")`
`library(ggplot2)`

alternatively ggplot2 is available in the tidyverse packages here:

`install.packages("tidyverse")`
`library(tidyverse)`


The reformat.syntenyData function also requires the stringR package included in the tidyverse collection of packages shown above or alternatively install using:

`install.packages("stringr")`
`library(stringr)`

---
## Input files 

### Input Alignment file

We provide a funtion to modify the output of DESCHRAMBLER or any other synteny tool that creates this type of file:

**Alignment output file format for reformat.syntenyData function**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example.deshrambler.output.png?raw=true" width="300" />

**Otherwise, you can provide your input Alignment file separated by tabs**

  DO NOT include a header line

* Reference chromosome ID 
* Reference start position
* Reference end position
* Target chromosome
* Target start position
* Target end position
* Orientation
* Reference species ID
* Target species ID

**Example input alignment file format for the three functions**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example.alignment.input.png?raw=true" width="400"/>


### Chromosome Length file for draw.ideogram and draw.linear functions

Please provide a file containing all aligned species in order from newest species – top of file – to ancestor – end of file, following this format, separated by tabs

  DO NOT include a header line
  
* Chromosome ID 
* Chromosome length
* Species ID

**Example file format**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example.lengths.input.png?raw=true" width="300"/>

---

## Reformatting alignment data

The syntenyPlotteR package includes a function to reformat alignmenty synteny data such as from DESCHRAMBLER or inferCARs - this does not curate files only reformats it

The function outputs a text file to the given file path or working directory containing the reformatted alignment data

### Usage 

`library(syntenyPlotteR)`

`reformat.syntenyData("file_data","filename",reference.species = "reference.sps",target.species = "target.sps")`

* file_data - alignment output file such as from DESCHRAMBLER/infeCARs
* filename - output file name for reformatted data in a character string


There are optional parameters for some customization of this function:

* reference.species - reference species ID as character string for final output table 
* target.species - target species ID as character string for final output table 

**Example code using data files in data folder**

`reformat.syntenyData("example_map_1.map","reformatted.data")`

---
## Evolution Highway style


### Usage
`library(syntenyPlotteR)`

`draw.eh("output",chrRange,...,fileformat = "png",colour = "lightblue",inverted.colour = "lightpink", sex.chromosome ="X",w=5.5,h=10,ps=10)`


* output - string assigned to name of output file name
* chrRange - the numeric range of chromosomes of the reference species this is entered as either a single number i.e. 1 or a range of numbers i.e. 1:22 
* ... - list of files containing the syntenic blocks (one file per alignment, in order from most recent species alignment file to ancestor alignment file)

There are optional parameters for some customization of this function:

* fileformat - format for saving the image i.e. png or pdf, parameter use: fileformat = "pdf" (the default value is "png")
* colour - colour of the syntenic blocks (not inverted blocks), parameter use: colour = "red" (the default value is "lightblue", see Rcolour pallette for colour options)
* colour - colour of the inverted syntenic blocks, parameter use: inverted.colour = "blue" (the default value is "lightpink", see Rcolour pallette for colour options)
* sex.chromosomes - The numeric range cannot accept characters, if sex chromosomes are required they can be added using: sex.chromosome ="X" or  sex.chromosome =c("X","Y") (the default value is "X")
* w -  The width of the image created can be changed by using: w = 5.5
* h -  The height of the image created can be changed by using: h = 10
* ps - The point size of the image created can be changed by using: ps = 10

 
**Example code using data files in data folder**

`draw.eh("example.eh",1:22, "example_alignment_1.txt","example_alignment_2.txt","example_alignment_3.txt")`

**Example output**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example.eh.png?raw=true" width="300"/>

---

## Chromosome painting style


### Usage

`library(syntenyPlotteR)`  

`draw.ideogram("file_data","sizefile","output",fileformat = "png",colours = colours.default,w=8.5,h=10,ps=5)`

* file_data - file containing the syntenic blocks
* sizefile - tab separated file of all chromosome, scaffold, or contig lengths and the species identifier (in order from newest species – top of file – to ancestor – end of file)
* output  -  string assigned to the output file name 

There are optional parameters for some customization of this function:

* fileformat - format for saving the image i.e. png or pdf, parameter use: fileformat = "pdf" (the default value is "png")
* colours - colours to assign to the ideograms in a concatenated string of chromosome IDs with assigned colour values which can be found with R colour Pallette, paramter use: colours = c("1" = "red", "2" = "blue", "3" = "green","4" = "orange", "5" = "purple","X" = "grey") if no colours are assigned default values will be used but colours MUST be assigned to all chromosomes
* w -  The width of the image created can be changed by using: w = 5.5
* h -  The height of the image created can be changed by using: h = 10
* ps - The point size of the image created can be changed by using: ps = 10

Target is the species which chromosomes will be painted. Reference will be used for painting and diagonals.
Chromosomes will be in the same order as in the target sizes file. 


**Example code using data files in data folder**

`draw.ideogram("example_alignment_1.txt","example_lengths.txt","example.ideogram")`


**Example output**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example.ideogram.png?raw=true" width="300"/>

---

## Linear style

### Usage

`library(syntenyPlotteR)`

`draw.linear(output,sizefile,...,fileformat = "png",colours = colours.default,w=13,h=5,opacity = .5)`

* output - string assigned to the output file name 
* sizefile - tab separated file of all chromosome, scaffold, or contig lengths and the species identifier (in order from newest species – top of file – to ancestor – end of file)
* ... - list of files containing the syntenic blocks (one file per alignment, in order from most recent species alignment file to ancestor alignment file)

Please ensure any species identifiers used between length and alignment files are matching (same identifiers and caseing)


There are optional parameters for some customization of this function:

* fileformat - format for saving the image i.e. png or pdf, parameter use: fileformat = "pdf" (the default value is "png")
* colours - colours to assign to the bands between ideograms in a concatenated string of chromosome IDs with assigned colour values which can be found with R colour Pallette, paramter use: colours = c("1" = "red", "2" = "blue", "3" = "green","4" = "orange", "5" = "purple","X" = "grey") if no colours are assigned default values will be used but colours MUST be assigned to all chromosomes
* w -  The width of the image created can be changed by using: w = 5.5
* h -  The height of the image created can be changed by using: h = 10
* opacity - the opacity of the ribbons can be changes using inputting: opacity = .5


**Example code using data files in data folder**

`draw.linear("example_linear","example_lengths.txt","example_alignment_1.txt","example_alignment_2.txt","example_alignment_3.txt")`


**Example output**

<img src="https://github.com/Farre-lab/syntenyPlotteR/blob/master/vignettes/images/example_linear.png?raw=true" width="600"/>

---



### Citation:
While we're preparing the publication, please cite:  
Farré M, Kim J, et al. Evolution of gene regulation in ruminants differs between evolutionary breakpoint regions and homologous synteny blocks. Genome Research 2019 Apr;29(4):576-589


