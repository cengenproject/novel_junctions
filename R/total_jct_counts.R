# Count the total number of SJ, for reference (note it may take some memory)


# inits ----


message("Starting", date())


input_dir <- "/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_junctions"
output_dir <- "data/intermediates/filtered_by_sample"


outliers_to_ignore <- readLines("data/outliers_to_ignore.txt")

read_sj_file <- function(cur_path){
  readr::read_tsv(cur_path,
                  col_names = c("chr", "start", "end", "strand",
                                "motif","annotated","nb_uniq", "nb_multi", "overhang"),
                  show_col_types = FALSE)
}

# fn must be a function operating on the rows of a matrix (e.g. rowSums, rowMaxs, etc...)
combine_sj <- function(sj_file, fn){
  dplyr::bind_rows(sj_file) |>
    dplyr::mutate(nb_uniq = sum(nb_uniq),
                  nb_multi = sum(nb_multi),
                  .by = c("chr", "start", "end", "strand", "motif", "annotated"))
}



message("read and process samples")

all_files <- dplyr::tibble(path = list.files(input_dir, full.names = TRUE),
                           replicate = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                           sample_id = stringr::str_split_fixed(replicate, "t", 2)[,1]) |>
  dplyr::filter(! sample_id %in% outliers_to_ignore) |>
  dplyr::mutate(sj_file = purrr::map(path, read_sj_file, .progress = TRUE)) |>
  dplyr::summarize(sj_file_combined = list(combine_sj(sj_file, rowSums)),
                   .by = sample_id)


per_sample <- map2_dfr(all_files$sj_file_combined, all_files$sample_id,
                       ~ {.x |>
          add_column(sample_id = .y)}) |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{1,3}$")[,2]) |>
  distinct()




per_sample |>
  select(chr,start,end,strand,motif,annotated) |>
  distinct()
#> # A tibble: 1,140,108 × 4
#>    chr   start   end strand
#>    <chr> <dbl> <dbl>  <dbl>
#>  1 I      2713  3202      2
#>  2 I      4359  5194      2
#>  3 I      5297  6036      2
# >>same with or without $motif and $annotated

per_sample |>
  select(chr,start,end,strand,motif,annotated) |>
  distinct() |>
  pull(annotated) |>
  table()
#>       0       1
#> 1026619  113489
#> 




