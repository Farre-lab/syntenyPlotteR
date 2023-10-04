#' Reformat synteny data
#'
#' This function takes output from alignment softwares such as deschrambler and inferCARs and re-formats it for syntenyPlotteR - this does not curate files only re-formats it
#'
#' It requires as input:
#'
#' 1. The data output from deschrambler or inferCARs
#'
#' 2. the desired output file name
#'
#' There are optional parameters for some customization of this function:
#'
#' 1. reference.species allows you to set the reference species identifier that will be set in the final output table i.e. `reference.species = "ref"`
#'
#' 2. target.species allows you to set the target species identifier that will be set in the final output table i.e. `target.species = "tar"`
#'
#' 3. The directory where the text file should be saved, as default the file is saved to temporary directory, change by inputting: `directory = "path/to/directory"`
#'
#' Example: `reformat.syntenyData("deschrambler.output", "reformatted.data", directory = "path/to/directory", reference.species = "ref", target.species = "tar" )`
#'
#' @title Reformat synteny data
#' @param file_data input file name for descrambler .map data
#' @param filename output file name for reformatted data
#' @param directory string containing file path to chosen directory to save text file
#' @param reference.species reference species identifier as a character string
#' @param target.species target species identifier as a character string
#' @return A text file with the reformatted data
#' @examples
#'
#' # Create object containing file path to external dataset
#' # (see vignette to follow examples with personal data)
#'
#' file <- system.file("extdata", "example_map_1.map", package = "syntenyPlotteR")
#'
#' # -----------------------------------------------------------------------------------
#'
#' # Run reformat.syntenyData function
#' # To run example and save file to working directory
#' # add directory parameter and set working directory
#' # To run example with personal data see vignette
#'
#' reformat.syntenyData(file, "outputName")
#' @export
#'

reformat.syntenyData <- function(file_data, filename, directory = NULL, reference.species = reference.sps, target.species = target.sps) {

  if (is.null(directory)) {
    directory <- tempdir()
  }

  desch <- utils::read.delim(file_data, header = FALSE) # input .map data from deschrambler
  line1 <- desch[seq(2, nrow(desch), 3), ]
  line2 <- desch[seq(3, nrow(desch), 3), ]
  combined <- as.data.frame(paste(line1, line2))
  names(combined) <- "combined"

  split <- as.data.frame(stringr::str_split_fixed(combined$combined, " ", 4))
  split$V1 <- gsub("[[:punct:]]", " ", split$V1)
  split$V3 <- gsub("[[:punct:]]", " ", split$V3)

  reference <- as.data.frame(stringr::str_split_fixed(split$V1, " ", 4))
  names(reference) <- c("reference.species", "reference.chromosome", "reference.start", "reference.stop")
  reference.no.chr <- gsub("chr", "", reference$reference.chromosome)
  reference$reference.chromosome <- reference.no.chr

  target <- as.data.frame(stringr::str_split_fixed(split$V3, " ", 4))
  names(target) <- c("target.species", "target.chromosome", "target.start", "target.stop")
  target.no.chr <- gsub("chr", "", target$target.chromosome)
  target$target.chromosome <- target.no.chr

  dataframe <- as.data.frame(cbind(
    reference$reference.chromosome, reference$reference.start, reference$reference.stop,
    target$target.chromosome, target$target.start, target$target.stop, split$V4,
    reference$reference.species, target$target.species
  ))

  reference.sps <- unique(dataframe$V8)
  target.sps <- unique(dataframe$V9)
  dataframe$V8 <- reference.species
  dataframe$V9 <- target.species
  message(paste0("Writing reformatted text file to ", directory))
  utils::write.table(dataframe, file = paste0(directory, "/", filename, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
