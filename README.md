# RUN THIS FIRST!!

## Large Data Files

Two large data files are not included in this repository and must be downloaded separately from Zenodo. **You only need to do this once.**

Run the following code in R from the root of this project:

```r
# Download and extract large data files from Zenodo
# Run this once after cloning the repository

zenodo_url <- "https://zenodo.org/records/17423417/files/First_Report_supplementary_code.zip?download=1"
zip_path <- tempfile(fileext = ".zip")

message("Downloading data from Zenodo...")
download.file(zenodo_url, destfile = zip_path, mode = "wb")

message("Extracting files...")
files_to_extract <- c(
  "First_Report_supplementary_code/data/NormBurnRatioPlusDiff.tif",
  "First_Report_supplementary_code/data/canopy_metrics_bdf_30m.tif"
)

unzip(zip_path, files = files_to_extract, exdir = tempdir())

# Move files to the local data/ folder
dir.create("data", showWarnings = FALSE)
file.copy(
  from = file.path(tempdir(), files_to_extract),
  to   = "data/",
  overwrite = TRUE
)

# Clean up
unlink(zip_path)
message("Done! Files are in the data/ folder.")
```

The full dataset and project archive are available on Zenodo at:  
[https://zenodo.org/records/17423417](https://zenodo.org/records/17423417)
