
# Principle


STAR junction file gives a junction count and whether it's an annotated or novel junction. We have one file per technical replicate (so 1 or 2 files per sample)., we group them if necessary to have 1 file per biological replicate.


In `filter_novel_junctions.R`, read each the STAR junction files, for each biological replicate filter the novel junctions and save them.

Then in `detect_consistent_junctions.R`, load all the junctions filtered from individual samples, and filter them based on appearing in enough samples. Export resulting list.


## Details

In the following, we define "is in the neighborhood of this gene" as being within 60 bp of that gene, ignoring the strand.

For sample filtering, we keep splice junctions that fulfill 4 criteria:
* is flanked by canonical splice site motifs
* no longer than 1 kb
* at least 2 supporting reads in that sample
* is not in the neighborhood of an rRNA gene
* is in the neighborhood of a protein-coding gene, long-non-coding RNA gene, or pseudogene
* has at least 20% as many reads as the most highly detected splice junction from the neighbor genes


For the second filtering, we keep splice junctions that:
* are detected in at least two samples of the same neuron type
* are detected in at least half the samples of the same neuron type


These thresholds were chosen by examining candidate splice junctions in the genome browser.


## Estimate number of genes with novel junctions

In `R/nb_genes_with_novel_sj.qmd` (and the resulting html file), we estimate bounds on the number of genes that can contain novel junctions.

