# using results of `filter_novel_junctions` (that filtered the individual samples)

library(GenomicFeatures)
library(tidyverse)

library(wbData)




fl <- list.files("data/intermediates/filtered_by_sample/")

fl <- fl[! startsWith(fl, "Refr")]


per_sample <- map_dfr(fl,
                      ~read_tsv(file.path("data/intermediates/filtered_by_sample/", .x),
                                col_types = cols(
                                  chr = col_character(),
                                  start = col_integer(),
                                  end = col_integer(),
                                  strand = col_factor(levels = c("1","2")),
                                  nb_uniq = col_integer()
                                )) |>
                        add_column(sample_id = .x)) |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{1,3}$")[,2]) |>
  # note: mistake in previous script: when t1 and t2, appears twice (but count correct)
  distinct()

stopifnot(all( ! is.na(per_sample$neuron_id) ))


per_neur <- per_sample |>
  mutate(nb_samples_neur = length(unique(sample_id)), .by = "neuron_id") |>
  summarize(nb_samples_detected = n(),
            .by = c(chr,start,end,strand,neuron_id, nb_samples_neur)) |>
  mutate(prop_samples_detected = nb_samples_detected/nb_samples_neur) |>
  mutate(coordinates = paste0(chr,":",start,"-",end))


## explore
# barplot(table(per_neur$nb_samples_detected))
# hist(per_neur$prop_samples_detected)
# 
# per_neur |> filter(nb_samples_detected > 5)
# per_neur |> filter(nb_samples_detected > 2)


# filter junctions ----
candidate_sjs <- per_neur |>
  filter(prop_samples_detected > .5,
         nb_samples_detected > 1) |>
  mutate(strand = recode(strand,
                         `1` = "+",
                         `2` = "-",
                         `0` = "*")) |>
  select(chr, start, end, strand, neuron_id, coordinates) |>
  summarize(number_neurons = n(),
            neurons = paste(neuron_id, collapse = ", "),
            .by = c(chr, start, end, strand, coordinates))



# find corresponding gene ----

txdb <- wb_load_TxDb(289)
gtable <- wb_load_gene_ids(289)

spl_genes_list <- gtable$gene_id[(gtable$biotype=="protein_coding_gene" |
                                    gtable$biotype=="lincRNA_gene" |
                                    gtable$biotype=="pseudogene")]
spliceable_genes_gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = spl_genes_list))

get_matching_genes <- function(i){
  matching_rows <- from(which_spl_genes_matching) == i
  matching_gene_rows <- to(which_spl_genes_matching)[matching_rows]
  matching_gene_ids <- spliceable_genes_gr$gene_id[matching_gene_rows]
  matching_gene_ids
}



gr_candidate_sjs <- GRanges(
  seqnames = candidate_sjs$chr,
  ranges = IRanges(start = candidate_sjs$start,
                   end = candidate_sjs$end),
  strand = candidate_sjs$strand
)

which_spl_genes_matching <- findOverlaps(gr_candidate_sjs,
                                         spliceable_genes_gr,
                                         maxgap = 60,
                                         ignore.strand = TRUE)

candidate_sjs <- candidate_sjs |>
  mutate(neighbor_genes = map_chr(row_number(),
                                  ~ get_matching_genes(.x) |>
                                    i2s(gtable,
                                        warn_missing = TRUE) |>
                                    paste(collapse = ", ")),
         .after = "coordinates")


# Export results ----

arrange(candidate_sjs, chr, start, end) |>
  writexl::write_xlsx("data/export/240418_novel_junctions.xlsx")















