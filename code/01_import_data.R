if (grepl("\\.gz$", file_path)) {
  counts <- read.csv(gzfile(file_path),
                     header = TRUE,
                     check.names = FALSE)
} else {
  counts <- read.csv(file_path,
                     header = TRUE,
                     check.names = FALSE)
}

dim(counts)
head(counts)
