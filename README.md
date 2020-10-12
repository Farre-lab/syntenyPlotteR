# chromDraw
R package to draw Evolution Highway style synteny plots

It creates a PDF file with one page per reference chromosome.

## To install:
`devtools::install_github("marta-fb/chromDraw")`

## To use:

`draw.eh(input,outputName,chromosomeRange)`

input file - synteny blocks following this format (please include this line as a header as well):  
  ref,chr,start,end,tarChr,tarSt,tarEnd,orient,tar
  
  ref - name of reference species  
  chr - chromosome of reference species   
  start - start of block in reference species  
  end - end of block in reference species  
  tarChr - chromosome of target species  
  tarSt - start of block in target species  
  orient - orientation of the block (1 or -1)  
  tar - name of target species  

outputName - name of output file   
chromosomeRange - range of chromosome numbers in reference species for autosomes

NOTE: ChrX is hardcoded, if you don't have an X please ignore the last error message.

### Example:  
`library(chromDraw)`

`draw.eh(testData.csv,testOut.pdf,1:5)`

## To come soon:

1. Update vignette
2. Alluvial style synteny plots
3. inferCARs style synteny plots
4. Preprint in Biorxiv
