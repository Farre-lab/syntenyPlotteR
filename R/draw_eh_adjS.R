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
