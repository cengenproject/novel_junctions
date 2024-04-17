
# inits ----
library(GenomicFeatures)
library(wbData)


message("Starting", date())


input_dir <- "/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_junctions"
output_dir <- "data/intermediates/filtered_by_sample"


outliers_to_ignore <- readLines("data/outliers_to_ignore.txt")



# gene coordinates for matching
txdb <- wb_load_TxDb(289)
gtable <- wb_load_gene_ids(289)

spl_genes_list <- gtable$gene_id[(gtable$biotype=="protein_coding_gene" |
                                    gtable$biotype=="lincRNA_gene" |
                                    gtable$biotype=="pseudogene")]
spliceable_genes_gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = spl_genes_list))

rrna_gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = gtable$gene_id[gtable$biotype=="rRNA_gene"]))



# functions ----


filter_one_sample <- function(sample_sj){
  
  candidate_sjs <- sample_sj |>
    dplyr::filter(annotated == 0, nb_uniq > 1) |>
    dplyr::filter(end - start +1 < 1000) |>
    dplyr::filter(motif != 0 & strand != 0) |>
    dplyr::mutate(coord = paste0(chr,":",start,"..",end))
  
  
  
  
  
  gr_candidates_sj <- GRanges(seqnames = candidate_sjs$chr,
                              ranges = IRanges(start = candidate_sjs$start,
                                               end = candidate_sjs$end),
                              strand = dplyr::recode(candidate_sjs$strand,
                                                     `1` = "+",
                                                     `2` = "-",
                                                     `0` = "*"))
  
  # match genes
  nb_spl_genes_matching <- countOverlaps(gr_candidates_sj,
                                         spliceable_genes_gr,
                                         maxgap = 60,
                                         ignore.strand = TRUE)
  
  nb_rrna_matching <- countOverlaps(gr_candidates_sj,
                                    rrna_gr,
                                    maxgap = 60,
                                    ignore.strand = TRUE)
  
  
  
  filtered_candidate_sjs <- candidate_sjs |>
    dplyr::mutate(nb_spl_genes_matching = nb_spl_genes_matching,
                  nb_rrna_matching = nb_rrna_matching) |>
    dplyr::filter(nb_rrna_matching == 0,
                  nb_spl_genes_matching > 0) |>
    dplyr::filter(nb_uniq > 2)
  
  
  
  # Local thresholding
  
  
  # max sj count per gene
  gr_sample_sj <- GRanges(seqnames = sample_sj$chr,
                          ranges = IRanges(start = sample_sj$start,
                                           end = sample_sj$end),
                          strand = dplyr::recode(sample_sj$strand,
                                                 `1` = "+",
                                                 `2` = "-",
                                                 `0` = "*"))
  
  max_sj_cnt_per_gene <- lapply(as(spliceable_genes_gr, "GRangesList"),
                                \(gr) findOverlaps(gr,
                                                   gr_sample_sj)) |>
    lapply(to) |>
    lapply(\(rows) sample_sj$nb_uniq[rows]) |>
    sapply(max)
  
  
  
  
  # per candidate SJ: max of neighbor genes' SJ
  gr_filtered_candidate_sjs <- GRanges(
    seqnames = filtered_candidate_sjs$chr,
    ranges = IRanges(start = filtered_candidate_sjs$start,
                     end = filtered_candidate_sjs$end),
    strand = dplyr::recode(filtered_candidate_sjs$strand,
                           `1` = "+",
                           `2` = "-",
                           `0` = "*")
  )
  
  which_spl_genes_matching <- findOverlaps(gr_filtered_candidate_sjs,
                                           spliceable_genes_gr,
                                           maxgap = 60,
                                           ignore.strand = TRUE)
  
  
  get_matching_genes <- function(i){
    matching_rows <- from(which_spl_genes_matching) == i
    matching_gene_rows <- to(which_spl_genes_matching)[matching_rows]
    matching_gene_ids <- spliceable_genes_gr$gene_id[matching_gene_rows]
    matching_gene_ids
  }
  
  
  sjs_to_keep <- filtered_candidate_sjs |>
    dplyr::mutate(neighbor_genes = lapply(dplyr::row_number(), get_matching_genes),
                  max_neighbor_counts = sapply(neighbor_genes,
                                               \(gene_ids) max(max_sj_cnt_per_gene[gene_ids]))) |>
    dplyr::mutate(perc_of_neighborhood = nb_uniq / max_neighbor_counts) |>
    dplyr::filter(perc_of_neighborhood > .2)
  
  
  sjs_to_keep |>
    dplyr::select(chr,start,end, strand, nb_uniq)
}



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


 # |>
 #  dplyr::mutate(sample_id = cur_sample) |>
 #  readr::write_tsv(file.path(output_dir, cur_sample) |> paste0(".tsv"))



message("read and process samples")

all_files <- dplyr::tibble(path = list.files(input_dir, full.names = TRUE),
                    replicate = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                    sample_id = stringr::str_split_fixed(replicate, "t", 2)[,1]) |>
  dplyr::filter(! sample_id %in% outliers_to_ignore) |>
  dplyr::mutate(sj_file = lapply(path, read_sj_file)) |>
  dplyr::summarize(sj_file_combined = list(combine_sj(sj_file, rowSums)),
                   .by = sample_id) |>
  dplyr::mutate(sj_file_combined = purrr::map(sj_file_combined, filter_one_sample, .progress = TRUE)) |>
  dplyr::mutate(out_path = paste0(output_dir, "/", sample_id, ".tsv"))

message("save files")

purrr::walk2(all_files$out_path, all_files$sj_file_combined,
             ~ readr::write_tsv(.y, .x))
  



message( "done")

