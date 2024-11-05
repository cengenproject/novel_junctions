
library(tidyverse)

# Import results

candidate_sjs <- readxl::read_xlsx("data/export/240418_novel_junctions.xlsx")



# number of neighbor genes per junction
candidate_sjs |>
  pull(neighbor_genes) |>
  str_split(",") |>
  sapply(length) |>
  table()

# number of genes that can appear as neighbors
candidate_sjs |>
  pull(neighbor_genes) |>
  str_split(",") |>
  unlist() |>
  unique() |>
  length()

# keeping up to 1 gene per junction
nb_genes <- replicate(1000,{
  candidate_sjs |>
    pull(neighbor_genes) |>
    str_split(",") |>
    sapply(\(x) sample(x, 1)) |>
    unique() |>
    length()
})
hist(nb_genes, breaks = 50)
table(nb_genes)


neighbors_per_sj <- candidate_sjs |>
  pull(neighbor_genes) |>
  str_split(",")

# upper bound: number of genes that appear at least once
neighbors_per_sj |>
  unlist() |>
  unique() |>
  length()
#> 1430

# lower bound: number of genes that are alone in their neighborhood
neighbors_per_sj[ sapply(neighbors_per_sj, length) == 1] |>
  unlist() |>
  unique() |>
  length()
#> 1289










