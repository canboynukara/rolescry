.to_char_df <- function(x) {
  as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE, check.names = FALSE)
}

.read_raw_excel <- function(path, sheet) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    m <- tryCatch(
      as.data.frame(
        readxl::read_excel(path, sheet = sheet %||% 1, col_names = FALSE, .name_repair = "minimal"),
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(m)) return(.to_char_df(m))
  }
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    raw <- tryCatch(
      openxlsx::read.xlsx(
        path, sheet = sheet %||% 1, colNames = FALSE,
        skipEmptyRows = FALSE, skipEmptyCols = FALSE, na.strings = character(0)
      ),
      error = function(e) NULL
    )
    if (!is.null(raw)) return(.to_char_df(raw))
  }
  stop("read_data(): reading Excel files needs the 'readxl' or 'openxlsx' package. ",
       "Install one with install.packages(\"readxl\").")
}

#' Read a data file with automatic header detection
#'
#' Reads a tabular file into a \code{data.frame}, detecting the header row with
#' [detect_header()] (the data is first read with no header so the header row is
#' visible as data). Delimited text (.csv/.tsv) is read with base R and always
#' works; spreadsheet and statistical formats use optional packages and degrade
#' gracefully with an actionable error if the package is absent.
#'
#' Supported: \code{.csv}, \code{.tsv}/\code{.tab} (base);
#' \code{.xlsx}/\code{.xls}/\code{.xlsm} (Suggests: readxl or openxlsx);
#' \code{.sav}/\code{.sas7bdat}/\code{.dta} (Suggests: haven, header intrinsic);
#' \code{.rds} (base, returned as stored).
#'
#' @param path Path to the file.
#' @param header Optional integer giving the 1-based header row to use directly,
#'   bypassing detection. \code{NULL} (default) auto-detects.
#' @param sheet Optional sheet name/index for Excel files.
#' @param na_strings Character vector of tokens mapped to \code{NA} before type
#'   conversion.
#' @param verbose Logical; emit the detected header row via \code{message()}.
#'
#' @return A \code{data.frame} with detected column names and per-column types
#'   inferred via \code{\link[utils]{type.convert}}.
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' writeLines(c("age,sex,score", "34,M,8.1", "51,F,7.4"), tmp)
#' read_data(tmp)
#' file.remove(tmp)
#'
#' @seealso [detect_header()], [detect_roles()]
#' @export
read_data <- function(path, header = NULL, sheet = NULL,
                      na_strings = c("", "NA", "N/A", "n/a", "na", "NULL", "null", "."),
                      verbose = FALSE) {
  if (!file.exists(path)) stop(sprintf("read_data(): file not found: %s", path))
  ext <- tolower(sub("^.*\\.", "", basename(path)))

  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext %in% c("sav", "zsav", "sas7bdat", "dta")) {
    if (!requireNamespace("haven", quietly = TRUE)) {
      stop("read_data(): reading .", ext, " files needs the 'haven' package.")
    }
    rd <- switch(ext,
      sav = , zsav = haven::read_sav(path),
      sas7bdat = haven::read_sas(path),
      dta = haven::read_dta(path)
    )
    return(as.data.frame(rd, stringsAsFactors = FALSE))
  }

  raw <- if (ext %in% c("xlsx", "xls", "xlsm")) {
    .read_raw_excel(path, sheet)
  } else {
    sep <- if (ext %in% c("tsv", "tab")) "\t" else ","
    utils::read.csv(
      path, header = FALSE, sep = sep, colClasses = "character",
      na.strings = character(0), check.names = FALSE, stringsAsFactors = FALSE,
      strip.white = TRUE
    )
  }
  if (is.null(raw) || nrow(raw) == 0L) stop("read_data(): no rows read from ", path)

  hr <- if (is.null(header)) {
    detect_header(raw, verbose = verbose)
  } else {
    hh <- as.integer(header)
    if (is.na(hh) || hh < 1L || hh > nrow(raw)) stop("read_data(): `header` out of range.")
    list(header_row = hh, names = .extract_header(raw, hh))
  }

  if (hr$header_row >= nrow(raw)) {
    body <- raw[0L, , drop = FALSE]
  } else {
    body <- raw[seq.int(hr$header_row + 1L, nrow(raw)), , drop = FALSE]
  }
  names(body) <- hr$names
  out <- as.data.frame(
    lapply(body, function(col) {
      col[col %in% na_strings] <- NA
      utils::type.convert(col, as.is = TRUE)
    }),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(out) <- NULL
  out
}
