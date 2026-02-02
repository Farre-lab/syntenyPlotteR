#' Convert DESCHRAMBLER adjacency scores to adjS format
#'
#' @param apcf_file Path to Ancestor.APCF
#' @param adjs_file Path to Ancestor.ADJS
#' @param block_list_file Path to SFs/block_list.txt
#' @param out_file Optional output file
#' @param ancestor_name Optional ancestor name override
#' @param include_ends Include 0-end adjacencies
#' @param chr_prefix Prefix added to APCF IDs
#' @return data.frame with columns: ancestor, chr, pos, score
#' @export


deschrambler_to_adjS <- function(apcf_file,
                                 adjs_file,
                                 block_list_file,
                                 out_file = NULL,
                                 ancestor_name = NULL,   # <-- NEW
                                 include_ends = FALSE,
                                 chr_prefix = "") {

  # ---- helper: robust SF token parser ----
  parse_sf_token <- function(tok) {
    if (length(tok) == 0 || is.na(tok)) return(NA_integer_)
    tok <- trimws(tok)
    if (tok == "") return(NA_integer_)

    if (grepl("^[+-]?\\d+$", tok)) {
      return(as.integer(tok))
    }

    if (grepl("^\\d+[+-]$", tok)) {
      s <- substr(tok, nchar(tok), nchar(tok))
      n <- as.integer(substr(tok, 1, nchar(tok) - 1))
      if (s == "-") n <- -n
      return(n)
    }

    cleaned <- gsub("[^0-9+-]", "", tok)
    if (grepl("^[+-]?\\d+$", cleaned)) {
      return(as.integer(cleaned))
    }

    NA_integer_
  }

  # ---- read block_list.txt ----
  bl <- utils::read.table(block_list_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(bl) < 5) {
    stop("block_list_file must have >=5 columns (DESCHRAMBLER block_list.txt format).")
  }

  colnames(bl)[1:5] <- c("ref_chr", "start0", "end1", "ori", "sfid")

  bl$sfid   <- as.integer(bl$sfid)
  bl$start0 <- as.numeric(bl$start0)
  bl$end1   <- as.numeric(bl$end1)

  bl$len <- bl$end1 - bl$start0
  if (any(!is.finite(bl$len) | bl$len <= 0)) {
    stop("Non-positive or NA SF lengths detected in block_list.txt.")
  }

  sf_len <- bl$len
  names(sf_len) <- as.character(bl$sfid)

  # ---- read Ancestor.ADJS ----
  adjs <- utils::read.table(adjs_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(adjs) < 3) {
    stop("adjs_file must have >=3 columns: SF1 SF2 score.")
  }

  colnames(adjs)[1:3] <- c("sf1", "sf2", "score")
  adjs$sf1   <- as.integer(adjs$sf1)
  adjs$sf2   <- as.integer(adjs$sf2)
  adjs$score <- as.numeric(adjs$score)

  key12 <- paste(adjs$sf1, adjs$sf2, sep = "\t")
  key21 <- paste(adjs$sf2, adjs$sf1, sep = "\t")
  score_map <- setNames(adjs$score, key12)
  score_map[key21] <- adjs$score

  # ---- parse Ancestor.APCF ----
  ap_lines <- readLines(apcf_file, warn = FALSE)
  if (length(ap_lines) == 0) stop("apcf_file is empty.")

  if (is.null(ancestor_name)) {
    if (!startsWith(ap_lines[1], ">")) {
      stop("Ancestor.APCF must start with '>' or provide ancestor_name explicitly.")
    }
    ancestor_name <- strsplit(sub("^>\\s*", "", ap_lines[1]), "\\s+")[[1]][1]
  }

  current_apcf <- NULL
  apcf_tokens <- list()

  for (ln in ap_lines[-1]) {
    if (startsWith(ln, "#")) {
      m <- regmatches(ln, regexec("^#\\s*APCF\\s+(\\S+)", ln))[[1]]
      if (length(m) >= 2) {
        current_apcf <- m[2]
        apcf_tokens[[current_apcf]] <- integer(0)
      } else {
        current_apcf <- NULL
      }
      next
    }

    if (is.null(current_apcf)) next
    if (trimws(ln) == "") next

    toks <- strsplit(trimws(ln), "\\s+")[[1]]
    vals <- vapply(toks, parse_sf_token, integer(1))
    vals <- vals[!is.na(vals)]

    if (length(vals) > 0) {
      apcf_tokens[[current_apcf]] <- c(apcf_tokens[[current_apcf]], vals)
    }
  }

  if (length(apcf_tokens) == 0) {
    stop("No APCF sections parsed from apcf_file.")
  }

  # ---- build adjS rows ----
  out_rows <- list()

  for (apcf_id in names(apcf_tokens)) {

    sfs_signed <- apcf_tokens[[apcf_id]]
    if (length(sfs_signed) < 2) next

    sfs_abs <- as.character(abs(sfs_signed))
    missing <- sfs_abs[!sfs_abs %in% names(sf_len)]
    if (length(missing) > 0) {
      stop(sprintf(
        "APCF %s references SF IDs not found in block_list.txt: %s",
        apcf_id, paste(unique(missing), collapse = ", ")
      ))
    }

    lengths <- sf_len[sfs_abs]
    cum_end <- cumsum(lengths)

    for (i in seq_len(length(sfs_signed) - 1)) {
      a <- sfs_signed[i]
      b <- sfs_signed[i + 1]

      scr <- score_map[[paste(a, b, sep = "\t")]]
      if (is.null(scr) || is.na(scr)) {
        scr <- score_map[[paste(abs(a), abs(b), sep = "\t")]]
      }

      out_rows[[length(out_rows) + 1]] <- data.frame(
        ancestor = ancestor_name,
        chr      = paste0(chr_prefix, apcf_id),
        pos      = as.numeric(cum_end[i]),
        score    = as.numeric(scr),
        stringsAsFactors = FALSE
      )
    }
  }

  out_df <- do.call(rbind, out_rows)
  if (is.null(out_df) || nrow(out_df) == 0) {
    out_df <- data.frame(
      ancestor = character(),
      chr      = character(),
      pos      = numeric(),
      score    = numeric(),
      stringsAsFactors = FALSE
    )
  }

  out_df$score <- pmin(pmax(out_df$score, 0), 1)

  if (!is.null(out_file)) {
    utils::write.table(
      out_df,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  return(out_df)
}
