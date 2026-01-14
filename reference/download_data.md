# Download Workshop Data from Zenodo

Downloads the workshop dataset files from Zenodo to the local data
directory. The dataset contains single-nucleus RNA-seq data from human
heart tissue across three developmental stages (foetal, young, adult).

## Usage

``` r
download_data(dest_dir = "data", overwrite = FALSE)
```

## Arguments

- dest_dir:

  Character. Destination directory for data files. Default: "data"

- overwrite:

  Logical. Overwrite existing files? Default: FALSE

## Value

Invisible NULL. Messages indicate download progress.

## Details

This function downloads the following files:

- `heart-counts.Rds`: Sparse count matrix (33,939 genes x 54,140 cells)

- `cellinfo_updated.Rds`: Cell metadata with sample, group, and sex info

Data source: Sim et al. (2021) "Sex-Specific Control of Human Heart
Maturation by the Progesterone Receptor", Circulation. DOI:
10.1161/CIRCULATIONAHA.120.051921

## Examples

``` r
if (FALSE) { # \dontrun{
# Download to default 'data/' directory
download_data()

# Download to custom directory
download_data(dest_dir = "my_data")

# Force re-download
download_data(overwrite = TRUE)
} # }
```
