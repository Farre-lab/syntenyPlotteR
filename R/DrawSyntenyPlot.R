#' Draw Evolution Highway Plots
#'
#' This function draws Evolution Highway style plots.
#'
#' It requires as input:
#'
#' 1. Desired output file name
#'
#' 2. The range of chromosomes of the reference species this is entered as either a single number i.e. 1 or a range of numbers i.e. 1:22.
#' *Note: if you are inputting sex chromosomes or chromosomes with characters in the name input a concatenated string i.e. c(1,2,3,"X")*
#'
#' 3. File containing the syntenic blocks of alignments of one or more target species aligned to a single reference; following this format:
#' reference chromosome, reference start position, reference end position, target chromosome,
#' target start position, target end position, orient, reference species identifier, target species identifier
#'
#' There are optional parameters for some customization of this function:
#'
#' 1. The format for saving the image i.e. png or pdf can be altered by inputting: `fileformat = "pdf"` (the default value is "png")
#'
#' 2. The colour of the syntenic blocks (not inverted blocks) can be changed by inputting: `colour = "red"` (the default value is "lightblue", see Rcolour pallette for colour options)
#'
#' 3. The colour of the inverted syntenic blocks can be changed by inputting: `inverted.colour = "blue"` (the default value is "lightpink", see Rcolour pallette for colour options)
#'
#' 5. The width of the image created can be changed by inputting: `w = 5.5`
#'
#' 6. The height of the image created can be changed by inputting: `h = 10`
#'
#' 7. The point size of the image created can be changed by inputting: `ps = 10`
#'
#' 8. The directory where the image file should be saved, as default the image is saved to temporary directory, change by inputting: `directory = "path/to/directory"`
#'
#' The function works creating a graph for each reference chromosome using their start and end positions to create a block for the reference
#' and the target chromosome positions are used to colour the region where synteny was identified
#'
#' Example: `draw.eh("outputName",c(17,"X"), "example_eh_alignments_2.txt", directory = "path/to/directory", fileformat = "pdf")`
#'
#'
#' @title Evolution Highway style plot
#' @param output output file name
#' @param chrRange range of chromosome numbers in the reference as numbers i.e. 1:29
#' @param data_file file containing the syntentic blocks from the alignments
#' @param directory string containing file path to chosen directory to save image file
#' @param fileformat output file format desired using the format `fileformat = "png"` (default is "png")
#' @param colour set colour for non-inverted syntenic blocks using the format `colour = "red"` (default is "lightblue")
#' @param inverted.colour set colour for inverted syntenic blocks using the format `inverted.colour = "blue"` (default is "lightpink")
#' @param w width of output image using the format `w = 5.5` (default)
#' @param h height of output image using the format `h = 10` (default)
#' @param ps point size of output image using the format `ps = 10` (default)
#' @return An image file showing the comparative drawings
#' @examples
#'
#' # Create object containing file path to external dataset
#' # (see vignette to follow examples with personal data)
#'
#' file1 <- system.file("extdata", "example_eh_alignments_2.txt", package = "syntenyPlotteR")
#'
#' # -----------------------------------------------------------------------------------
#'
#' # Run draw.eh function
#' # To run example and save file to working directory
#' # add directory parameter and set working directory
#' # To run example with personal data see vignette
#'
#' draw.eh("outputName", c(17, "X"), file1, fileformat = "pdf")
#' @export
#'

draw.eh <- function(output,
                    chrRange,
                    data_file,
                    directory = NULL,
                    fileformat = "png",
                    colour = "lightblue",
                    inverted.colour = "lightpink",
                    w = 5.5,
                    h = 10,
                    ps = 10) {

  if (is.null(directory)) {
    directory <- tempdir()
  }

  colours <- c("1" = colour, "-1" = inverted.colour)

  outfile <- chr <- start <- end <- tarChr <- tarSt <- tarEnd <- orient <- tar <- text_size2 <- NULL
  alignments <- utils::read.table(data_file, header = FALSE)
  colnames(alignments) <- c("chr", "start", "end", "tarChr", "tarSt", "tarEnd", "orient", "ref", "tar")
  alignments$tar <- as.factor(alignments$tar)
  alignments$ref <- as.factor(alignments$ref)
  alignments$start <- as.numeric(gsub(",", "", alignments$start))
  alignments$end <- as.numeric(gsub(",", "", alignments$end))
  alignments$tarSt <- as.numeric(gsub(",", "", alignments$tarSt))
  alignments$tarEnd <- as.numeric(gsub(",", "", alignments$tarEnd))
  alignments$orient[alignments$orient == "+"] <- "1"
  alignments$orient[alignments$orient == "-"] <- "-1"
  alignments$orient <- factor(alignments$orient, levels = c("1", "-1"))
  alignments$text_size2 <- 80 * ((alignments$end - alignments$start) / 100000)


  plots <- ggplot2::ggplot()
  for (ID in c(chrRange)) {
    message(paste0("Saving eh image for chromosome ", ID, " to ", directory))
    subsetChr1 <- subset(alignments, chr == ID, select = c(chr, start, end, tarChr, tarSt, tarEnd, orient, tar, text_size2))
    min <- min(subsetChr1$start)
    max <- max(subsetChr1$end)
    plots <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = subsetChr1, mapping = ggplot2::aes(
          xmin = start, xmax = end, ymin = 0,
          ymax = 0.5, fill = orient, group = tar
        ),
        color = "white", alpha = 1, linewidth = 0.1
      ) +
      ggplot2::geom_rect(data = subsetChr1, ggplot2::aes(xmin = min, xmax = max, ymin = 0, ymax = 0.5), linewidth = 0.3, color = "black", fill = "NA") +
      ggplot2::facet_grid(~tar) +
      ggplot2::coord_flip() +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(
        panel.spacing.y = ggplot2::unit(c(-0.5, -0.5), "lines"),
        panel.spacing.x = ggplot2::unit(0, "lines"),
        panel.background = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::geom_text(data = subsetChr1, ggplot2::aes(x = start + (end - start) / 2, y = 0 + (0.5 - 0) / 2, label = tarChr, size = text_size2))

    print(plots)

    ggplot2::ggsave(paste0(directory,"/",output, ".", ID, ".", fileformat), plots, device = fileformat, width = w, height = h, pointsize = ps)
  }
}

#' Draw Linear Synteny Plots
#'
#' This function draws linear synteny plots.
#'
#' It requires:
#'
#' 1. The desired output file name;
#'
#' 2. Tab separated file of all chromosome, scaffold, or contig lengths and the species identifier,
#' in order from first target species in the alignment files followed by the first reference species in the alignment files
#' -- top of file -- to the last target species and reference species in the alignment files -- end of file.
#' in this format:
#' chromosome ID, chromosome length, species identifier
#'
#' 3. files containing the syntenic blocks - one file per alignment, in order from first target/reference
#' (most recent species pairwise alignment in ancestral reconstruction data) alignment file
#' to last target/reference (ancestor pairwise alignment in ancestral reconstruction data) alignment file
#' following this format:
#' reference chromosome, reference start position, reference end position, target chromosome,
#' target start position, target end position, orient, reference species identifier, target species identifier
#'
#' Please separate files by tab and ensure any species identifiers used between length and alignment files are matching (same identifiers and caseing)
#'
#'
#' There are optional parameters for some customization of this function:
#'
#' 1. The format for saving the image i.e. png or pdf can be altered by inputting: `fileformat = "pdf"` (the default value is "png")
#'
#' 2. The colour of the synteny bands can be altered by inputting a concatenated string of chromosome IDs with assigned colour values which can be found with R colour Pallette
#' e.g. `colours = c("1" = "red", "2" = "blue", "3" = "green","4" = "orange", "5" = "purple","X" = "grey")` if no colours are assigned default values will be used but colours MUST be assigned to all chromosomes
#'
#' 3. The width of the image created can be changed by inputting: `w = 13` (default)
#'
#' 4. The height of the image created can be changed by inputting: `h = 5` (default)
#'
#' 5. The opacity of the ribbons can be changed by inputting: `opacity = .5` (default)
#'
#' 6. The directory where the image file should be saved, as default the image is saved to temporary directory, change by inputting: `directory = "path/to/directory"`
#'
#' The function works using the chromosome length file to order the Y axis and provide chromosome lengths to draw chromosome ideograms and the alignment files provides coordinates to draw the alignment bands between ideograms
#'
#' Example: `draw.linear("outputname", "example_lengths.txt", "example_alignment_1.txt", "example_alignment_2.txt", "example_alignment_3.txt", directory = "path/to/directory", fileformat = "pdf")`
#'
#' @title Linear synteny plot
#' @param output output file name
#' @param sizefile Chromosome Size file
#' @param ... synteny files (any number of alignment files can be entered)
#' @param directory string containing file path to chosen directory to save image file
#' @param fileformat output file format specified using the format `fileformat = "pdf"` (the default is "png")
#' @param colours concatenated string of chromosome IDs and assigned colours if desired using the format `colours = c("1" = "red", "2" = "blue", "3" = "green", "X" = "grey")` if the no colours are assigned default values will be used
#' @param w width of output image using the format `w = 13` (default)
#' @param h height of output image using the format `h = 5` (default)
#' @param opacity opacity of syntenic bands using the format `opacity = .5` (default)
#' @return An image file showing the linear comparison drawings
#' @examples
#'
#' # Create objects containing file paths to external dataset
#' # (see vignette to follow examples with personal data)
#'
#' length.file <- system.file("extdata", "example_lengths.txt", package = "syntenyPlotteR")
#' file1 <- system.file("extdata", "example_alignment_1.txt", package = "syntenyPlotteR")
#' file2 <- system.file("extdata", "example_alignment_2.txt", package = "syntenyPlotteR")
#' file3 <- system.file("extdata", "example_alignment_3.txt", package = "syntenyPlotteR")
#'
#' # -----------------------------------------------------------------------------------
#'
#' # Run draw.linear function
#' # To run example and save file to working directory
#' # add directory parameter and set working directory
#' # To run example with personal data see vignette
#'
#' draw.linear("outputName", length.file, file1, file2, file3, fileformat = "pdf")
#' @export
#'
draw.linear <- function(output, sizefile, ..., directory = NULL, fileformat = "png", colours = colours.default, w = 13, h = 5, opacity = .5) {

  if (is.null(directory)) {
    directory <- tempdir()
  }

  synteny.data.reframing <- function(data, tar.y, ref.y, compiled.size) {
    synteny <- data.frame()
    for (i in c(1:nrow(data))) {
      reference <- data[i, "ref.species"]
      target <- data[i, "tar.species"]
      tar_chr <- data[i, "tarchr"]
      ref_chr <- data[i, "refchr"]
      dir <- data[i, "dir"]
      tar_sizes <- compiled.size[compiled.size$species == target, ]
      names(tar_sizes) <- c("tarchr", "size", "species", "xstart", "xend")
      ref_sizes <- compiled.size[compiled.size$species == reference, ]
      names(ref_sizes) <- c("refchr", "size", "species", "xstart", "xend")
      tar_add <- tar_sizes[as.character(tar_sizes$tarchr) == as.character(tar_chr), ]$xstart
      ref_add <- ref_sizes[as.character(ref_sizes$refchr) == as.character(ref_chr), ]$xstart
      tar_y <- tar.y
      ref_y <- ref.y
      tar_xstart <- data[i, "tarstart"] + tar_add
      tar_xend <- data[i, "tarend"] + tar_add
      ref_xstart <- data[i, "refstart"] + ref_add
      ref_xend <- data[i, "refend"] + ref_add

      inverted <- grepl("-", dir, fixed = TRUE)
      if (inverted == TRUE) {
        df <- data.frame(
          x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      } else {
        df <- data.frame(
          x = c(tar_xstart, ref_xstart, ref_xend, tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      }
      synteny <- rbind(synteny, df)
    }
    return(synteny)
  }

  colours.default <- c(
    "1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
    "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
    "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white", "LG2" = "black", "LG5" = "#CC99CC"
  )

  xstart <- xend <- refchr <- tarchr <- x <- y <- group <- fill <- chromosome <- species <- NULL
  sizes <- utils::read.delim(sizefile, header = FALSE) # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))

  count <- 0
  compiled.size <- data.frame()
  for (i in unique(sizes$species)) {
    size.intermediate <- sizes[sizes$species == i, ]
    for (x in c(1:nrow(size.intermediate))) {
      if (x == 1) {
        total_start <- 1
        total_end <- size.intermediate[x, "size"]
      } else {
        total_start <- total_end + 6000000
        total_end <- total_start + size.intermediate[x, "size"]
      }
      size.intermediate[x, "xstart"] <- total_start
      size.intermediate[x, "xend"] <- total_end
    }
    compiled.size <- rbind(compiled.size, size.intermediate)
  }

  for (z in unique(compiled.size$species)) {
    compiled.size$y[compiled.size$species == z] <- count
    count <- count + 2
  }

  list.of.files <- list()
  for (i in list(...)) {
    list.of.files[[i]] <- i
  }

  listsynt <- list()
  for (i in 1:length(list.of.files)) {
    num <- i
    file <- list.of.files[[num]]
    dataTMP <- utils::read.delim(file, header = FALSE)
    data2 <- dataTMP[, c(4, 5, 6, 1, 2, 3, 7, 8, 9)]
    colnames(data2) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir", "ref.species", "tar.species")
    data2$tarstart <- as.numeric(gsub(",", "", data2$tarstart))
    data2$tarend <- as.numeric(gsub(",", "", data2$tarend))
    data2$refstart <- as.numeric(gsub(",", "", data2$refstart))
    data2$refend <- as.numeric(gsub(",", "", data2$refend))
    reference <- data2[1, "ref.species"]
    target <- data2[1, "tar.species"]
    ref_y <- compiled.size[compiled.size$species == reference, "y"]
    tar_y <- compiled.size[compiled.size$species == target, "y"]
    if (tar_y[1] > ref_y[1]){
      ref_y <- ref_y[1] + 0.1
      tar_y <- tar_y[1]
    } else{
      ref_y <- ref_y[1]
      tar_y <- tar_y[1] + 0.1
    }
    x <- synteny.data.reframing(data2, tar_y, ref_y, compiled.size)
    x$fill <- as.factor(x$fill)
    listsynt[[i]] <- x
  }

  compiled.size$chromosome <- as.factor(compiled.size$chromosome)

  p <- ggplot2::ggplot()

  for (i in 1:length(listsynt)) {
    data <- listsynt[[i]]
    reference <- data[1, "ref"]
    target <- data[1, "tar"]
    ref_sizes <- compiled.size[compiled.size$species == reference, ]
    tar_sizes <- compiled.size[compiled.size$species == target, ]
    p <- p + ggplot2::geom_rect(
      data = ref_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10, fill = chromosome),
      color = "black", alpha = 0.85, linewidth = 0.2
    ) +
      ggplot2::geom_text(data = ref_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      ggplot2::geom_text(data = ref_sizes, mapping = ggplot2::aes(x = 2, y = y, label = species), size = 3, hjust = 1) +
      ggplot2::geom_rect(
        data = tar_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10), fill = "grey85",
        color = "black", alpha = 0.85, linewidth = 0.2
      ) +
      ggplot2::geom_text(data = tar_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      ggplot2::geom_text(data = tar_sizes, mapping = ggplot2::aes(x = 2, y = y, label = species), size = 3, hjust = 1) +
      ggplot2::geom_polygon(data = data, alpha = opacity, ggplot2::aes(x = x, y = y, group = group, fill = fill))
  }

  p <- p + ggplot2::scale_fill_manual(values = colours) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  message(paste0("Saving linear image to ", directory))
  print(p)
  ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), p, device = fileformat, width = w, height = h)
}

#' Draw synteny ideograms in Chromosome painting style
#'
#' This function draws pairwise synteny plots in chromosome painting style.
#'
#' Inputs are tab separated files;
#'
#' It requires as input:
#'
#' 1. File containing the syntenic blocks following this format:
#' reference chromosome, reference start position, reference end position, target chromosome,
#' target start position, target end position, orient, reference species identifier, target species identifier
#'
#' 2. Tab separated file of all chromosome, scaffold, or contig lengths and the species identifier,
#' in order from first target species in the alignment files followed by the first reference species in the alignment files
#' -- top of file -- to the last target species and reference species in the alignment files -- end of file.
#' in this format:
#' chromosome ID, chromosome length, species identifier
#'
#' 3. The desired output file name
#'
#' Please separate files by tab and ensure any species identifiers used between length and alignment files are matching (same identifiers and caseing)
#'
#' There are optional parameters for some customization of this function:
#'
#' 1. The format for saving the image i.e. png or pdf can be altered by inputting: `fileformat = "pdf"` (the default value is "png")
#'
#' 2. The colour of the ideograms can be altered by inputting a concatenated string of chromosome IDs with assigned colour values which can be found with R colour Pallette
#' e.g. `colours = c("1" = "red", "2" = "blue", "3" = "green","4" = "orange", "5" = "purple","X" = "grey")` if no colours are assigned default values will be used but colours MUST be assigned to all chromosomes
#'
#' 3. The width of the image created can be changed by inputting: `w = 8.5` (default)
#'
#' 4. The height of the image created can be changed by inputting: `h = 10` (default)
#'
#' 5. The point size of the image created can be changed by inputting: `ps = 5` (default)
#'
#' 6. The directory where the image file should be saved, as default the image is saved to temporary directory, change by inputting: `directory = "path/to/directory"`
#'
#' Target is the species which chromosomes will be painted. Reference will be used for painting and diagonals.
#' Chromosomes will be in the same order as in the target chromosomes in the chromosome length file
#'
#' Example: `draw.ideogram("example_alignment_1.txt", "example_lengths.txt", "outputname", directory = "path/to/directory", fileformat = "pdf")`
#'
#'
#' @title Draw ideograms in chromosome painting style
#' @param file_data Path to the syntenic blocks file
#' @param sizefile Chromosome size file
#' @param output output file name
#' @param directory string containing file path to chosen directory to save image file
#' @param fileformat output file format specified using `fileformat = "pdf"` (the default is "png")
#' @param colours concatenated string of chromosome IDs and assigned colours if desired using the format `colours = c("1" = "red", "2" = "blue", "3" = "green","X" = "grey")` if the no colours are assigned default values will be used
#' @param w width of output image using `w = 8.5` (default)
#' @param h height of output image using `h = 10` (default)
#' @param ps point size of output image using `ps = 5` (default)
#' @return An image file showing the ideogram
#' @examples
#'
#' # Create objects containing file paths to external dataset
#' # (see vignette to follow examples with personal data)
#'
#' length.file <- system.file("extdata", "example_lengths.txt", package = "syntenyPlotteR")
#' file1 <- system.file("extdata", "example_alignment_1.txt", package = "syntenyPlotteR")
#'
#' # -----------------------------------------------------------------------------------
#'
#' # Run draw.ideogram function
#' # To run example and save file to working directory
#' # add directory parameter and set working directory
#' # To run example with personal data see vignette
#'
#' draw.ideogram(file1, length.file, "outputName", fileformat = "pdf")
#' @export

draw.ideogram <- function(file_data, sizefile, output, directory = NULL, fileformat = "png", colours = colours.default, w = 8.5, h = 10, ps = 5) {

  if (is.null(directory)) {
    directory <- tempdir()
  }

  colours.default <- c(
    "1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
    "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
    "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white"
  )


  size <- tarstart <- tarend <- refchr <- ystart <- yend <- NULL
  data <- utils::read.delim(file_data, header = FALSE)

  colnames(data) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "orien", "tar", "ref")
  data$tarstart <- as.numeric(gsub(",", "", data$tarstart))
  data$tarend <- as.numeric(gsub(",", "", data$tarend))
  data$refstart <- as.numeric(gsub(",", "", data$refstart))
  data$refend <- as.numeric(gsub(",", "", data$refend))

  sizes <- utils::read.delim(sizefile, header = FALSE) # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))

  ref <- unique(data$ref)
  tar <- unique(data$tar)

  tar_sizes <- sizes[sizes$species == tar, ]
  ref_sizes <- sizes[sizes$species == ref, ]

  colnames(tar_sizes) <- c("tarchr", "size")
  colnames(ref_sizes) <- c("refchr", "size")

  data$tarchr <- factor(data$tarchr, levels = tar_sizes$tarchr)
  data$refchr <- factor(data$refchr, levels = ref_sizes$refchr)


  for (i in c(1:nrow(data))) {
    dir <- data[i, "orien"]
    chr <- data[i, "refchr"]
    full_len <- ref_sizes[ref_sizes$refchr == chr, 2]
    y1 <- round(data[i, "refstart"] / full_len, digits = 4)
    y2 <- round(data[i, "refend"] / full_len, digits = 4)

    inverted <- grepl("-", dir, fixed = TRUE)
    if (inverted == TRUE) {
      data[i, "ystart"] <- y2
      data[i, "yend"] <- y1
    } else {
      data[i, "ystart"] <- y1
      data[i, "yend"] <- y2
    }
  }



  plots <- ggplot2::ggplot(size = 0.2, font = 10, data = data) +
    ggplot2::geom_rect(
      data = tar_sizes, mapping = ggplot2::aes(xmin = 1, xmax = size, ymin = -0.1, ymax = 1.1),
      fill = "white", color = "black", alpha = 0.85, linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = data, mapping = ggplot2::aes(xmin = tarstart, xmax = tarend, ymin = -0.1, ymax = 1.1, fill = refchr),
      color = "black", alpha = 0.85, linewidth = 0.2
    ) +
    ggplot2::geom_segment(data = data, mapping = ggplot2::aes(x = tarstart, y = ystart, xend = tarend, yend = yend), linewidth = 0.2) +
    ggplot2::facet_grid(as.factor(tarchr) ~ .) +
    ggplot2::labs(fill = "Reference", x = "Chomosome length (Mb)", size = 10) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(size = 0.2),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.title.align = 0.5
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_x_continuous(
      breaks = c(
        0, 2.5e+07, 5e+07, 7.5e+07, 1e+08, 1.25e+08, 1.5e+08, 1.75e+08,
        2e+08, 2.25e+08, 2.5e+08, 2.75e+08, 3e+08, 3.25e+08, 3.5e+08
      ),
      labels = c(
        "0", "25", "50", "75", "100", "125", "150", "175",
        "200", "225", "250", "275", "300", "325", "350"
      )
    )
  message(paste0("Saving ideogram image to ", directory))
  print(plots)
  ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), plots, device = fileformat, width = w, height = h, pointsize = ps)
}
