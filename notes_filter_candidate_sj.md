# Some notes about filtering candidate junctions


Take sj files from STAR:
  
  
```
chr     start   end strand   motif annotated nb_uniq nb_multi overhang
I       4359    5194     2       2         1       3        0       33
I       5297    6036     2       2         1       6        0       49
I       5474    6019     1       1         0       1        0       33
I       6328    9726     2       2         1       3        0       44
I       9847    10094    2       2         1       2        0       18
```


## Choice of criteria

unannotated, enough overhang, enough reads in two samples from the same neuron.

Also need to eliminate too long.

Note: in unannotated junctions, the lowest overhang is 12 (default filter `outSJfilterOverhangMin`).

Examples in AFDr39:
* I:811,083..819,364 (4 reads, 44 overhang) likely real (has many reads in other samples)
* I:7428781..7522115 (48 reads, 12 overhang), but 93 kb long: likely not real
* V:9,706,254..9,706,346 (4 reads, 15 overhang) likely not real, in a gene with extremely high expression levels
* II:5206985..5207035 (9 reads, 18 overhang) likely real but not expressed in neurons (highly expressed gene)
* 

Test criteria: at least 5 reads, not longer than 1 kb.
`filter(xx, annotated == 0, nb_uniq > 5) |> filter(end - start +1 < 1000) |> mutate(coord = paste0(chr,":",start,"..",end)) |> as.data.frame()|> head()`


```
 chr   start     end strand motif annotated nb_uniq nb_multi overhang
1   I  881202  881278      2     2         0       8        0       48
2   I 1567020 1567600      1     1         0      11        1       50
3   I 2062064 2062117      2     2         0      10        0       49
4   I 5306448 5306517      2     2         0      29        0       50
5   I 5364121 5364670      1     1         0      10        0       45
6   I 5705453 5705855      1     1         0      15        0       49
               coord
1   I:881202..881278
2 I:1567020..1567600
3 I:2062064..2062117
4 I:5306448..5306517
5 I:5364121..5364670
6 I:5705453..5705855
```

* 1: not real, in a highly expressed gene
* 2: real
* 3: not real, in a highly expressed gene
* 4: real
* 5: real
* 6: real

So, it would be useful to filter out those in highly expressed genes.


# Threshold on nb reads

```
candidate_jcts <- filter(xx, annotated == 0) |>
filter(end - start +1 < 1000) |>
mutate(coord = paste0(chr,":",start,"..",end))


candidate_jcts |> filter(nb_uniq < 5) |> slice_sample(n = 4) |> as.data.frame()
  chr    start      end strand motif annotated nb_uniq nb_multi overhang
1  II 12727313 12728170      2     2         0       0        1       29
2 III   913446   913595      2     2         0       1        0       23
3   X 13911549 13912279      2     2         0       1        0       27
4 III  6069971  6070026      1     1         0       1        0       41
                  coord
1 II:12727313..12728170
2    III:913446..913595
3  X:13911549..13912279
4  III:6069971..6070026
```


* 1: not real? Unclear why it's in the list, no evidence in other samples either
* 2: count as real (weak evidence)
* 3: not real
* 4: not real

so 0 or 1 read not enough to support it.

```
  chr    start      end strand motif annotated nb_uniq nb_multi overhang
1   X 15264112 15264314      1     1         0       4        0       48
2  II  5974867  5974911      2     2         0       2        0       34
3  IV 13882283 13882437      1     1         0       2        0       37
4   V 15534740 15535182      1     1         0       2        0       43
5   X 11755621 11756072      1     1         0       2        0       25
6   V 20783611 20783737      1     1         0       3        0       37
                  coord
1  X:15264112..15264314
2   II:5974867..5974911
3 IV:13882283..13882437
4  V:15534740..15535182
5  X:11755621..11756072
6  V:20783611..20783737
```

* 1: real
* 2: unclear (weak evidence), in a highly expressed gene
* 3: in a highly expressed gene
* 4: not real 
* 5: in a highly expressed gene
* 6: in a highly expressed gene


So the fact that we have only 2 reads may not be a problem, the high expression of the containing gene is much worse



Also note that the motif and strand tend to be undefined together:
```
table(candidate_jcts$motif, candidate_jcts$strand)
   
       0    1    2
  0 1495    0    0
  1    0 3097    0
  2    0    0 2826
  3    0  139    0
  4    0    0  150
  5    0   36    0
  6    0    0   30
```

and looking at a few, they don't look real:
```

candidate_jcts |>
  dplyr::filter(motif == 0 | strand == 0) |>
  dplyr::mutate(coord = paste0(chr,":",start,"..",end)) |>
  dplyr::slice_sample(n = 3) |> 
  as.data.frame()
  chr    start      end strand motif annotated nb_uniq nb_multi overhang
1   V 10792159 10792552      0     0         0       6        0       48
2  IV  4701296  4701853      0     0         0       9        0       59
3  IV  4805577  4805626      0     0         0       4        0       30
                 coord
1 V:10792159..10792552
2  IV:4701296..4701853
3  IV:4805577..4805626
```



# Looking at overlaps (ignoring strand and allowing 60 bp gap)
```
table(nb_genes_matching)
nb_genes_matching
   0    1    2    3    4    5    6    7 
 398 5272  530   54   15    5    3    1
```

Check the ones that don't overlap any gene

```
> candidate_jcts |>
+   dplyr::mutate(nb_genes_matching = nb_genes_matching) |>
+   dplyr::filter(nb_genes_matching == 0) |>
+   dplyr::slice_sample(n = 6) |>
+   as.data.frame()
  chr    start      end strand motif annotated nb_uniq nb_multi overhang                coord nb_genes_matching
1  IV  1455920  1456286      1     1         0       5        0       36  IV:1455920..1456286                 0
2  II   612508   612557      2     2         0       3        0       33    II:612508..612557                 0
3  IV  1446429  1446727      2     2         0       4        0       70  IV:1446429..1446727                 0
4   I  3156024  3156078      2     2         0       2        0       37   I:3156024..3156078                 0
5   X 11228730 11228768      2     2         0      14        0       19 X:11228730..11228768                 0
6   V 17309454 17310248      1     1         0       2       54       70 V:17309454..17310248                 0
```


* 1: intergenic, too low
* 2: could be real, at 40 bp (old version, only allowing 20 bp gap) of clec-118; low expression though, weak evidence
* 3: not real
* 4: not real
* 5: not real
* 6: real junction, probably not real gene

Looking specifically at some with high junction counts, a few may be real (including real new genes), some correspond either to markers (rol-6) or unclear.




## matching no gene of any kind but high count
```
> candidate_jcts |>
+   dplyr::mutate(nb_protcod_genes_matching = nb_genes_matching,
+                 nb_tot_genes_matching = all_nb_genes_matching,
+                 match_rrna = match_rrna) |>
+   dplyr::filter(match_rrna == 0,
+                 nb_genes_matching == 0,
+                 nb_tot_genes_matching == 0) |>
+   dplyr::filter(nb_uniq > 10) |>
+   dplyr::slice_sample(n = 6) |>
+   as.data.frame()
  chr    start      end strand motif annotated nb_uniq
1   X 15338691 15338743      1     1         0      63
2 III  2583617  2583678      2     2         0      20
3 III  7445263  7445541      1     1         0      12
4 III 13295594 13295635      2     2         0      22
5   I  5980864  5980911      2     2         0      16
6 III 13295432 13295483      2     2         0      22
  nb_multi overhang                  coord
1        0       71   X:15338691..15338743
2        0       57   III:2583617..2583678
3        0       69   III:7445263..7445541
4        5       58 III:13295594..13295635
5        0       71     I:5980864..5980911
6        0       72 III:13295432..13295483
  nb_protcod_genes_matching nb_tot_genes_matching
1                         0                     0
2                         0                     0
3                         0                     0
4                         0                     0
5                         0                     0
6                         0                     0
  match_rrna
1          0
2          0
3          0
4          0
5          0
```

* 1: real, too far from annotated gene, would need to reconstruct tx
* 2: too far from annotated gene, not real anyway
* 3: weird thing, not real
* 4: weird thing, not real (or novel gene)
* 5: real, too far from annotated gene, would need to reconstruct tx (exple: simr-1)
* 6: same as 4


So ignoring loses some real signal, but would need to reconstruct transcripts


# match ncRNA but not spliceable RNA

```
candidate_jcts |>
  dplyr::mutate(spl_genes_matching = spl_genes_matching,
                nb_tot_genes_matching = all_nb_genes_matching,
                match_rrna = match_rrna) |>
  dplyr::filter(match_rrna == 0,
                spl_genes_matching == 0,
                nb_tot_genes_matching > 0) |>
  dplyr::filter(nb_uniq > 10) |>
  dplyr::slice_sample(n = 6) |>
  as.data.frame()
  chr    start      end strand motif annotated nb_uniq nb_multi overhang
1   X 11992009 11992053      2     2         0      91        0       69
2   X 16937710 16938464      2     2         0      11        0       67
3   V 10332208 10332276      1     1         0      11        0       25
4   V 17316010 17316280      1     1         0      28        0       63
5 III  8211650  8211695      2     2         0      15        0       14
6   V 11427769 11427853      1     5         0      13        0       51
                 coord spl_genes_matching nb_tot_genes_matching match_rrna
1 X:11992009..11992053                  0                     1          0
2 X:16937710..16938464                  0                     1          0
3 V:10332208..10332276                  0                     2          0
4 V:17316010..17316280                  0                     1          0
5 III:8211650..8211695                  0                     1          0
6 V:11427769..11427853                  0                     1          0
```

* 1: likely linc-15 longer (real, too far), there is a ncRNA closer
* 2: weird, in ncRNA with lots of SJ
* 3: noise
* 4: too far from gene (probably real), there is a ncRNA closer
* 5: not real (noise)
* 6: not real (noise)


We can ignore those for now


## several matching genes

```
candidate_jcts |>
  dplyr::mutate(spl_genes_matching = spl_genes_matching,
                nb_tot_genes_matching = all_nb_genes_matching,
                match_rrna = match_rrna) |>
  dplyr::filter(match_rrna == 0,
                spl_genes_matching > 1) |>
  dplyr::filter(nb_uniq > 10) |>
  dplyr::slice_sample(n = 6) |>
  as.data.frame()
  chr    start      end strand motif annotated nb_uniq nb_multi overhang                 coord
1  IV 13371022 13371249      1     1         0      15        0       51 IV:13371022..13371249
2  IV  1549053  1549405      2     2         0      11        0       53   IV:1549053..1549405
3   I  6253315  6253700      1     1         0      18        2       61    I:6253315..6253700
4  IV  9403748  9403797      2     2         0     174        9       75   IV:9403748..9403797
5   I  4174775  4175452      2     2         0      15        0       43    I:4174775..4175452
6  II  3215128  3216068      1     1         0      24        0       61   II:3215128..3216068
  spl_genes_matching nb_tot_genes_matching match_rrna
1                  2                     2          0
2                  2                     2          0
3                  2                     2          0
4                  2                     2          0
5                  2                     2          0
6                  2                     2          0
```


* 1: two genes overlap (but high expr genes)
* 2: gene in intron (but high expr gene)
* 3: link to neighbor gene (high expr)
* 4: gene in intron, real (annotated in pseudogene)
* 5: link to neighbor gene (high expr)
* 6: link to neighbor gene (could be real)


So we want to keep those and submit them to the local thresholding



## More spliceable genes matching

previous was mostly for 2 matches. There are 8 additional SJ that match 3 or more spliceable genes.
```
candidate_jcts |>
  dplyr::mutate(spl_genes_matching = spl_genes_matching,
                nb_tot_genes_matching = all_nb_genes_matching,
                match_rrna = match_rrna) |>
  dplyr::filter(match_rrna == 0,
                spl_genes_matching > 2) |>
  dplyr::filter(nb_uniq > 10) |>
  as.data.frame()
```

we don't exclude they could be real (though not convincing).

## Summary of selection

We remove anything in rRNA. We remove SJ that don't match genes (loosing a few interesting cases). We keep those with spliceable genes matching (no matter how many).

We can also do a first basic filtering and remove if we don't have at least 3 reads.





# Filtering on local counts

```

xx <- filtered_candidate_sjs |>
  dplyr::mutate(neighbor_genes = lapply(dplyr::row_number(), get_matching_genes),
                max_neighbor_counts = sapply(neighbor_genes,
                                             \(gene_ids) max(max_sj_cnt_per_gene[gene_ids]))) 




xx |>
  dplyr::select(coord, nb_uniq, max_neighbor_counts) |>
  dplyr::mutate(perc = nb_uniq / max_neighbor_counts) |>
  dplyr::filter(perc > .2,
                perc < .3) |>
  dplyr::slice_sample(n = 3)
```

Under 10% we should discard; over 20% it usually looks real (though in many cases gene lowly expressed), in 10-20% not totally clear, usually looks not real: discard.




