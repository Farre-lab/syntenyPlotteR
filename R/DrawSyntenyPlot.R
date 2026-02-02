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
#' 
#' @param adj_file Character string or NULL. Optional file with adjacency scores
#'   (e.g. converted from DESCHRAMBLER). If NULL, no adjacency panel is drawn.
#' @param adj_point_size Numeric. Size of points in the adjacency score panel.
#' @param adj_line Logical. Whether to draw a connecting line between adjacency score points.
#' @param adj_colour_low Character string. Low colour for adjacency score gradient.
#' @param adj_colour_high Character string. High colour for adjacency score gradient.
#' @param adj_panel_frac Numeric. Fraction of total plot width allocated to the
#'   adjacency score panel (between 0 and 1).
#' @param adj_score_transform Character. Optional transformation applied to
#'   adjacency scores ("identity", "sqrt", or "log10").
#'
#' @param strip_angle Numeric. Angle of species (facet) labels.
#' @param strip_text_size Numeric. Text size of species (facet) labels.
#'
#' @param shorten_block_labels Logical. Whether to shorten chromosome labels
#'   inside EH blocks.
#' @param shorten_n Integer. Number of characters kept from the right when
#'   shortening block labels.
#' @param label_fixed_size Numeric. Base font size for labels inside EH blocks.
#' @param label_scale_by_width Logical. Whether to scale block label size by the
#'   number of species panels.
#' @param label_width_exponent Numeric. Exponent controlling how strongly label
#'   size decreases with increasing number of species panels.
#' @param label_min Numeric. Minimum allowed block label font size.
#' @param label_max Numeric. Maximum allowed block label font size.
#'
#' @return Invisibly returns NULL. The function is called for its side effect
#'   of saving one image per chromosome.
#'
#' @details
#' Adjacency scores are expected to be numeric values between 0 and 1 and are
#' plotted as a heatmap (orange to red) along genomic coordinates. One-to-many
#' mappings are retained in the adjacency panel.
#'
#' @seealso \code{\link{draw.linear}}
#'
#' @export

draw.eh <- function(output,
                    chrRange,
                    data_file,
                    directory = NULL,
                    fileformat = "png",
                    colour = "lightblue",
                    inverted.colour = "lightpink",
                    w = 5.5,
                    h = 10,
                    ps = 10,

                    # ---- adjacency panel options ----
                    adj_file = NULL,                 # AdjS: ancestor, chr, pos, score
                    adj_point_size = 1.1,
                    adj_line = TRUE,
                    adj_colour_low = "#FDB863",      # orange
                    adj_colour_high = "#B2182B",     # red
                    adj_panel_frac = 0.02,          # narrow panel - % of total panel
                    adj_score_transform = c("identity", "sqrt", "log10"),

                    # ---- facet strip (species) labels ----
                    strip_angle = 90,
                    strip_text_size = 10,

                    # ---- in-block label controls ----
                    shorten_block_labels = TRUE,
                    shorten_n = 3,                   # keep last N characters
                    label_fixed_size = 4,          # constant base size
                    label_scale_by_width = TRUE,     # shrink if many species columns
                    label_width_exponent = 0.5,      # 0.5 = 1/sqrt(n_tar)
                    label_min = 2,
                    label_max = 10) {

  if (is.null(directory)) directory <- tempdir()
  adj_score_transform <- match.arg(adj_score_transform)

  colours <- c("1" = colour, "-1" = inverted.colour)

  # ---- Read EH alignment input ----
  alignments <- utils::read.table(data_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(alignments) <- c("chr", "start", "end", "tarChr",
                            "tarSt", "tarEnd", "orient", "ref", "tar")

  alignments$tar <- as.factor(alignments$tar)
  alignments$ref <- as.factor(alignments$ref)

  alignments$start  <- as.numeric(gsub(",", "", alignments$start))
  alignments$end    <- as.numeric(gsub(",", "", alignments$end))
  alignments$tarSt  <- as.numeric(gsub(",", "", alignments$tarSt))
  alignments$tarEnd <- as.numeric(gsub(",", "", alignments$tarEnd))

  alignments$orient[alignments$orient == "+"] <- "1"
  alignments$orient[alignments$orient == "-"] <- "-1"
  alignments$orient <- factor(alignments$orient, levels = c("1", "-1"))

  # ---- Build shortened in-block labels (guaranteed "last N from the right") ----
  alignments$tarChr_short <- trimws(as.character(alignments$tarChr))

  if (isTRUE(shorten_block_labels)) {
    n <- as.integer(shorten_n)
    if (!is.finite(n) || n < 1) n <- 3

    alignments$tarChr_short <- vapply(
      alignments$tarChr_short,
      FUN = function(s) {
        if (is.na(s)) return(NA_character_)
        s <- trimws(s)
        if (nchar(s) <= n) return(s)
        substr(s, nchar(s) - n + 1, nchar(s))
      },
      FUN.VALUE = character(1)
    )
  }

  # ---- Optional: read adjacency scores ----
  adj_all <- NULL
  if (!is.null(adj_file)) {
    adj_all <- utils::read.table(adj_file, header = FALSE, sep = "\t",
                                 stringsAsFactors = FALSE)
    if (ncol(adj_all) < 4)
      stop("adj_file must have >=4 columns: ancestor, chr, pos, score")

    colnames(adj_all)[1:4] <- c("ancestor", "chr", "pos", "score")
    adj_all$chr   <- as.character(adj_all$chr)
    adj_all$pos   <- as.numeric(gsub(",", "", adj_all$pos))
    adj_all$score <- as.numeric(gsub(",", "", adj_all$score))

    if (adj_score_transform == "sqrt") {
      adj_all$score_t <- sqrt(pmax(adj_all$score, 0))
    } else if (adj_score_transform == "log10") {
      adj_all$score_t <- log10(pmax(adj_all$score, 1e-12))
    } else {
      adj_all$score_t <- adj_all$score
    }
  }

  # ---- Plot per chromosome ----
  for (ID in c(chrRange)) {

    message(paste0("Saving eh image for chromosome ", ID, " to ", directory))

    subsetChr1 <- subset(
      alignments, chr == ID,
      select = c(chr, start, end, tarChr_short, orient, tar)
    )

    if (nrow(subsetChr1) == 0) next

    ymin_gen <- min(subsetChr1$start, na.rm = TRUE)
    ymax_gen <- max(subsetChr1$end, na.rm = TRUE)

    # ---- Constant label sizes (optionally width-scaled only) ----
    n_tar <- length(levels(subsetChr1$tar))
    if (!is.finite(n_tar) || n_tar < 1) n_tar <- 1

    size <- rep(label_fixed_size, nrow(subsetChr1))
    if (isTRUE(label_scale_by_width)) {
      size <- size * (1 / (n_tar ^ label_width_exponent))
    }
    size <- pmin(pmax(size, label_min), label_max)
    subsetChr1$text_size2 <- size

    # ---- MAIN EH PANEL ----
    p_main <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = subsetChr1,
        ggplot2::aes(
          xmin = 0, xmax = 0.5,
          ymin = start, ymax = end,
          fill = orient,
          group = tar
        ),
        color = "white",
        linewidth = 0.1
      ) +
      ggplot2::geom_rect(
        data = subsetChr1,
        ggplot2::aes(xmin = 0, xmax = 0.5,
                     ymin = ymin_gen, ymax = ymax_gen),
        linewidth = 0.3,
        color = "black",
        fill = NA
      ) +
      ggplot2::geom_text(
        data = subsetChr1,
        ggplot2::aes(
          x = 0.25,
          y = start + (end - start) / 2,
          label = tarChr_short,
          size = text_size2
        )
      ) +
      ggplot2::scale_size_identity(guide = "none") +   # <-- KEY FIX (no rescaling)
      ggplot2::facet_grid(~tar) +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
      ggplot2::scale_y_reverse(limits = c(ymax_gen, ymin_gen),
                               expand = c(0, 0)) +
      ggplot2::theme(
        panel.spacing.y = grid::unit(c(-0.5, -0.5), "lines"),
        panel.spacing.x = grid::unit(0, "lines"),
        panel.background = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(
          angle = strip_angle,
          size = strip_text_size,
          margin = ggplot2::margin(b = 4),
          vjust = 0.5,
          hjust = 0
        ),
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )

    # ---- ADJACENCY PANEL ----
    p_adj <- NULL
    if (!is.null(adj_all)) {
      adj_chr <- adj_all[
        adj_all$chr == as.character(ID) &
          is.finite(adj_all$pos) &
          is.finite(adj_all$score_t),
        , drop = FALSE
      ]

      if (nrow(adj_chr) > 0) {

        adj_chr$score_c <- pmin(pmax(adj_chr$score_t, 0), 1)

        p_adj <- ggplot2::ggplot(adj_chr,
                                 ggplot2::aes(x = score_t, y = pos)) +
          ggplot2::geom_vline(
            xintercept = c(0, 0.5, 1),
            colour = "grey70",
            linewidth = 0.3,
            linetype = "dashed"
          ) +
          { if (adj_line)
              ggplot2::geom_path(
                colour = "grey50",
                linewidth = 0.4,
                alpha = 0.8
              )
          } +
          ggplot2::geom_point(
            ggplot2::aes(colour = score_c),
            size = adj_point_size
          ) +
          ggplot2::scale_colour_gradient(
            low = adj_colour_low,
            high = adj_colour_high,
            limits = c(0, 1),
            name = NULL
          ) +
          ggplot2::scale_y_reverse(limits = c(ymax_gen, ymin_gen),
                                   expand = c(0, 0)) +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
            plot.margin = ggplot2::margin(t = 5, r = 4, b = 5, l = 0)
          )
      }
    }


# ---- COMBINE PANELS ----
    final_plot <- p_main

    if (!is.null(p_adj)) {

    # clamp fraction safely
      f <- adj_panel_frac
      if (!is.finite(f)) f <- 0.10
      f <- max(min(f, 0.9), 0.02)   # keep it sane

      main_w <- 1 - f
      adj_w  <- f

      if (requireNamespace("patchwork", quietly = TRUE)) {
        final_plot <- p_main + p_adj +
        patchwork::plot_layout(widths = c(main_w, adj_w))
      } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        final_plot <- gridExtra::arrangeGrob(
          p_main, p_adj,
          nrow = 1,
          widths = grid::unit(c(main_w, adj_w), "null")
        )
      }
    }


    ggplot2::ggsave(
      filename = paste0(directory, "/", output, ".", ID, ".", fileformat),
      plot = final_plot,
      device = fileformat,
      width = w,
      height = h,
      pointsize = ps
    )
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
