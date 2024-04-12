
# Principle


STAR junction file gives a junction count and wether it's an annotated or novel junction. We have one file per technical replicate (so 1 or 2 files per sample). 


In `filter_novel_junctions.R`, read all the STAR junction files, for each one filter the novel junctions and save them.

Then load all the junctions filtered from individual samples, and filter them based on appearing in enough samples.



