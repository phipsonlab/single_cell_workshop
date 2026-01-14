#' Download Workshop Data from Zenodo
#'
#' Downloads the workshop dataset files from Zenodo to the local data directory.
#' The dataset contains single-nucleus RNA-seq data from human heart tissue
#' across three developmental stages (foetal, young, adult).
#'
#' @param dest_dir Character. Destination directory for data files.
#'   Default: "data"
#' @param overwrite Logical. Overwrite existing files? Default: FALSE
#'
#' @return Invisible NULL. Messages indicate download progress.
#'
#' @details
#' This function downloads the following files:
#' \itemize{
#'   \item \code{heart-counts.Rds}: Sparse count matrix (33,939 genes x 54,140 cells)
#'   \item \code{cellinfo_updated.Rds}: Cell metadata with sample, group, and sex info
#' }
#'
#' Data source: Sim et al. (2021) "Sex-Specific Control of Human Heart
#' Maturation by the Progesterone Receptor", Circulation.
#' DOI: 10.1161/CIRCULATIONAHA.120.051921
#'
#' @examples
#' \dontrun{
#' # Download to default 'data/' directory
#' download_data()
#'
#' # Download to custom directory
#' download_data(dest_dir = "my_data")
#'
#' # Force re-download
#' download_data(overwrite = TRUE)
#' }
#'
#' @export
download_data <- function(dest_dir = "data", overwrite = FALSE) {

    # Zenodo record ID (update when published)
    zenodo_record <- "18237749"
    base_url <- paste0("https://zenodo.org/records/", zenodo_record, "/files/")

    # Files to download
    files <- c(
        "heart-counts.Rds",
        "cellinfo_updated.Rds"
    )

    # Create destination directory
    if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
        message("Created directory: ", dest_dir)
    }

    message("Downloading workshop data from Zenodo (record: ", zenodo_record, ")")
    message("Destination: ", normalizePath(dest_dir, mustWork = FALSE))
    message("")

    # Download each file
    for (f in files) {
        dest_file <- file.path(dest_dir, f)

        if (file.exists(dest_file) && !overwrite) {
            message("- ", f, ": already exists, skipping (use overwrite=TRUE to re-download)")
            next
        }

        url <- paste0(base_url, f, "?download=1")
        message("- Downloading: ", f, " ...")

        tryCatch({
            utils::download.file(url, dest_file, mode = "wb", quiet = FALSE)
            message("  Saved to: ", dest_file)
        }, error = function(e) {
            warning("Failed to download ", f, ": ", e$message,
                    "\n  Check your internet connection or try again later.")
        })
    }

    message("")
    message("Download complete!")
    message("You can now proceed with the workshop modules.")

    invisible(NULL)
}
