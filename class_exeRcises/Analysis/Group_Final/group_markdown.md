Investigation of Promoter Window Size on Genomic DBP Activity
================
Breakout 2

We were interested in looking at how the window size around a promoter
affects the structure of the bimodal density plot observed in class RMD
20. Specifically we broke down our results by both biotype (mRNA,
lncRNA) and by counting method (“count”, “occurence”). First we imported
the libraries and local functions that we needed to complete this task.

Next we imported the Gencode and consensus peak data that we generated
in the class 11 RMD.

``` r
# Code from class 20 RMD
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2022/data/genomes/gencode.v32.annotation.gtf")
# Loading in files via listing and rtracklayer import
consensus_fl <- list.files("/scratch/Shares/rinnclass/CLASS_2022/EM/CLASS_2022/class_exeRcises/Analysis/11_r_functions/consensus_peaks", full.names = T)
# importing (takes ~5min)
consensus_peaks <- lapply(consensus_fl, rtracklayer::import)
# cleaning up file names
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2022/EM/CLASS_2022/class_exeRcises/Analysis/11_r_functions/consensus_peaks/|_consensus_peaks.bed","", consensus_fl)
# Filtering consensus peaks to those DBPs with at least 250 peaks
num_peaks_threshold <- 250
num_peaks <- sapply(consensus_peaks, length)
filtered_consensus_peaks <- consensus_peaks[num_peaks > num_peaks_threshold]
# Result: these were the DBPs that were filtered out.
filtered_dbps <- consensus_peaks[num_peaks < num_peaks_threshold]
names(filtered_dbps)
```

    ##  [1] "CEBPZ"    "GPBP1L1"  "H3K27me3" "HMGA1"    "IRF3"     "MLLT10"  
    ##  [7] "MYBL2"    "NCOA5"    "RNF219"   "RORA"     "ZBTB3"    "ZFP36"   
    ## [13] "ZFP62"    "ZMAT5"    "ZNF10"    "ZNF17"    "ZNF260"   "ZNF382"  
    ## [19] "ZNF48"    "ZNF484"   "ZNF577"   "ZNF597"   "ZNF7"

``` r
# We have this many remaining DBPs
length(filtered_consensus_peaks)
```

    ## [1] 460

We checked the names of dbps and the number of consensus peaks to ensure
that the data had been imported correctly.

First, we looked at both biotypes (mRNA, lncRNA) using the “counts”
method of count\_peaks\_per\_feature. This method counts each and every
binding event to a promoter window such that if the same protein binds
multiple times, each binding event is added to the total. This method
does not allow for the differentiation of number of unique DBPs that
bound inside the promoter window.

``` r
numbers = seq(0,5000, by = 250)

# Create matrixes to store the data
mrna_ngenes <- 19965
mrna_mat = matrix(numeric(), ncol = 1000, nrow = length(numbers))

lncrna_ngenes <- 16849
lncrna_mat = matrix(numeric(), ncol = 1000, nrow = length(numbers))

# subset promoters to save loop time (size 1000 random sample)
gene_samples_lncRNA <- gencode_gr[gencode_gr$type == "gene"]
gene_samples_lncRNA <- gene_samples_lncRNA[gene_samples_lncRNA$gene_type %in% "lncRNA"]
gene_samples_lncRNA <- gene_samples_lncRNA[sample(1:length(gene_samples_lncRNA), 1000)]

gene_samples_mRNA <- gencode_gr[gencode_gr$type == "gene"]
gene_samples_mRNA <- gene_samples_mRNA[gene_samples_mRNA$gene_type %in% "protein_coding"]
gene_samples_mRNA <- gene_samples_mRNA[sample(1:length(gene_samples_mRNA), 1000)]

for (i in 1:length(numbers)){ # loop over all window sizes 
  
  # get the regions for each biotype.
  lncrna_promoters <- get_promoter_regions(gene_samples_lncRNA, biotype = "lncRNA",     
                                           upstream = numbers[i]/2, downstream = numbers[i]/2)
  mrna_promoters <- get_promoter_regions(gene_samples_mRNA, biotype = "protein_coding",     
                                           upstream = numbers[i]/2, downstream = numbers[i]/2)
  
  # get gene names and count the number of peaks
  names(mrna_promoters) <- mrna_promoters$gene_id
  mrna_promoter_peak_counts <- count_peaks_per_feature(mrna_promoters, filtered_consensus_peaks, type = "counts")
  
  names(lncrna_promoters) <- lncrna_promoters$gene_id
  lncrna_promoter_peak_counts <- count_peaks_per_feature(lncrna_promoters, filtered_consensus_peaks, type = "counts")
  
  # save to matrix
  mrna_mat[i,] <- colSums(mrna_promoter_peak_counts)
  lncrna_mat[i,] <- colSums(lncrna_promoter_peak_counts)
}
```

We saved the results of count\_peaks\_per\_feature to a large matrix and
named the rows and columns for readability. In this matrix, rows are
different window sizes, and columns are different genes. Each entry in
this matrix represents the number of DBPs that bound with that window
size to that gene.

``` r
# assign row names as window sizes
rownames(mrna_mat) <- as.character(numbers)
rownames(lncrna_mat) <- as.character(numbers)

# assign column names as genes
colnames(mrna_mat) <- mrna_promoters$gene_id
colnames(lncrna_mat) <- lncrna_promoters$gene_id
```

Next, we employed pivot\_longer in order to make our data easily
plottable. This transformed our matrix into a three-column data frame
where the first column was biotype, the second column was window size,
and the third column was DBP count.

``` r
# transpose and cast to data frame
t_mrna_mat <- as.data.frame(t(mrna_mat))

# pivot_longer for plotting ease
mrna_mat_plotting <- t_mrna_mat %>%                          
  pivot_longer(colnames(t_mrna_mat)) %>% 
  as.data.frame()

# makes the legend ordered based on window size
mrna_mat_plotting$name <- as.numeric(mrna_mat_plotting$name)
```

``` r
# same operations with lncRNA data
t_lncrna_mat <- as.data.frame(t(lncrna_mat))

lncrna_mat_plotting <- t_lncrna_mat %>%                          
  pivot_longer(colnames(t_lncrna_mat)) %>% 
  as.data.frame()

lncrna_mat_plotting$name <- as.numeric(lncrna_mat_plotting$name)
```

Then we bound both of the dataframes together and proceeded to plot all
data on the same figure.

``` r
# Add dummy variable to differentiate biotype
mrna_mat_dummy <- data.frame(mrna_mat_plotting, 'dummy' = rep('mRNA'))

lncrna_mat_dummy <- data.frame(lncrna_mat_plotting, 'dummy' = rep('lncRNA'))

# Mash the data together
merged_mrna_lncrna <- rbind(mrna_mat_dummy, lncrna_mat_dummy)

# Plot it pretty
ggplot(data = merged_mrna_lncrna) +
  geom_density(aes(x = value, color = as.factor(name))) + 
  facet_wrap(~dummy) + 
  labs(color = 'Window Length') +
  labs(title = "Count of DBP Peaks Across Biotype and Promoter Window Size") + 
  theme(legend.position = 'bottom') +
  xlab("Number of DBP Peaks") + 
  scale_color_manual(values = get_palette(palette = 'RdYlBu', k=length(numbers)))+
  xlim(1,700)
```

![](group_markdown_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave('figures/peak_counts.png')
```

This figure shows the distinct bimodal feature of both biotypes of the
“counts” method. As window size increases, the shape of each of the
curves changes. With lncRNA, we see that the second peak of the bimodal
plot almost disappears as window size increases. The tail of the density
is stretched and made more heavy as window size increases. This is in
contrast to mRNA, where both peaks are clearly visible at all window
sizes. The tail is heavier with a larger window size for this density as
well. This was to be expected, as with a larger window size there is
more area on the genome for DBPs to bind and therefore one would expect
more binding events to be observed. This is consistent with what these
plots show for both mRNA and lncRNA.

Next, we repeated the same figure using the “occurence” method of
count\_peaks\_per\_feature. This method only counts unique DBPs that
bound inside the promoter window. Therefore, even if a DBP bound
multiple times inside the window, it would only be counted once.

``` r
# same procedure with "occurence"
numbers = seq(0,5000, by = 250)

# Create matrixes to store the data
mrna_mat = matrix(numeric(), ncol = 1000, nrow = length(numbers))
lncrna_mat = matrix(numeric(), ncol = 1000, nrow = length(numbers))

for (i in 1:length(numbers)){ # loop over all window sizes 
  
   # get the regions for each biotype.
  lncrna_promoters <- get_promoter_regions(gene_samples_lncRNA, biotype = "lncRNA",     
                                           upstream = numbers[i]/2, downstream = numbers[i]/2)
  mrna_promoters <- get_promoter_regions(gene_samples_mRNA, biotype = "protein_coding",     
                                           upstream = numbers[i]/2, downstream = numbers[i]/2)
  
  # get gene names and count the number of peaks
  names(mrna_promoters) <- mrna_promoters$gene_id
  mrna_promoter_peak_counts <- count_peaks_per_feature(mrna_promoters, filtered_consensus_peaks, type = "occurrence")
  
  names(lncrna_promoters) <- lncrna_promoters$gene_id
  lncrna_promoter_peak_counts <- count_peaks_per_feature(lncrna_promoters, filtered_consensus_peaks, type = "occurrence")
  
  # save to matrix
  mrna_mat[i,] <- colSums(mrna_promoter_peak_counts)
  lncrna_mat[i,] <- colSums(lncrna_promoter_peak_counts)
}
```

We performed the same data transformation on the “occurence” data as
with the “counts” data.

``` r
# assign row names as window sizes
rownames(mrna_mat) <- as.character(numbers)
rownames(lncrna_mat) <- as.character(numbers)

# assign column names as genes
colnames(mrna_mat) <- mrna_promoters$gene_id
colnames(lncrna_mat) <- lncrna_promoters$gene_id
```

``` r
# same operations with lncRNA data
t_mrna_mat <- as.data.frame(t(mrna_mat))

mrna_mat_plotting <- t_mrna_mat %>%                         
  pivot_longer(colnames(t_mrna_mat)) %>% 
  as.data.frame()

mrna_mat_plotting$name <- as.numeric(mrna_mat_plotting$name)
```

``` r
# same operations with lncRNA data
t_lncrna_mat <- as.data.frame(t(lncrna_mat))

lncrna_mat_plotting <- t_lncrna_mat %>%                          
  pivot_longer(colnames(t_lncrna_mat)) %>% 
  as.data.frame()

lncrna_mat_plotting$name <- as.numeric(lncrna_mat_plotting$name)
```

``` r
# Add dummy variable to differentiate biotype
mrna_mat_dummy <- data.frame(mrna_mat_plotting, 'dummy' = rep('mRNA'))

lncrna_mat_dummy <- data.frame(lncrna_mat_plotting, 'dummy' = rep('lncRNA'))

# Mash the data together
merged_mrna_lncrna <- rbind(mrna_mat_dummy, lncrna_mat_dummy)

# Plot it pretty
ggplot(data = merged_mrna_lncrna) +
  geom_density(aes(x = value, color = as.factor(name))) + 
  facet_wrap(~dummy) + 
  labs(color = 'Window Length') +
  labs(title = "Unique DBP Across Biotype and Promoter Window Size") + 
  theme(legend.position = 'bottom') +
  xlab("Number of DBP Proteins") + 
  scale_color_manual(values = get_palette(palette = 'RdYlBu', k=length(numbers)))+
  xlim(1,500)
```

![](group_markdown_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave('figures/peak_occurences.png')
```

We see that the bimodal distribution conserves its character from the
“counts” plot to the “occurrence” plot. This could indicate that
proteins that have high binding have also a diverse number of proteins
binding. We see the same pattern for increasing window size, where
increasing window size makes the tail of the distribution heavier.

Our data indicates that the window size does not affect the character of
the bimodal distribution of DNA binding proteins. This suggests that
there are two classes of types of promoters, one where very few proteins
bind and one where many proteins bind. Between biotypes, we see that
there are many more promoters with high number of DNA binding proteins
in the mRNA promoters than the lncRNA promoters. With the lncRNA
promoters, there are many promoters with very few DNA binding proteins
and few with lots of DNA binding proteins.

It is important to note that all the ChIP data used in this analysis are
from HEPG2, an immortalized liver cancer cell line. This is atypical
cell line with mutations and abnormal chromosome number. It would be
interesting to see if these result held across immortalized cell lines
and healthy cell lines.

In general, it would be interesting to explore what characterizes
promoters that constitute these two distinct peaks. A future direction
of this research would be to look at promoter motifs and to determine if
there are conserved nucleotides corresponding to the promoters that are
represented by the two different peaks. Looking at this pattern across
biotypes may reveal a pattern where conserved nucleotides that
correspond to the second peak appear in both biotypes. Another possible
direction this research could take would be to look at the types of DBPs
that are binding at each of these peaks and see if there are proteins
that are characteristic of one of these two promoter categories.

Looking at RNAseq transcription levels at each of these promoters may
give additional insight into what differentiates these two different
classes of promoters. This could evidence the idea that having more DBPs
bind to a promoter allows for more (or less) transcription to take
place. A final interesting direction to take this research would be to
perform this same analysis with cells that had been grown in different
conditions in order to determine whether the set of genes that have high
expression are conserved.

Determining why there is such a distinctive peak of 200-300 proteins of
which many promoters would be essential to understanding how DBPs are
promoting or inhibiting expression of certain genes. Comparing similar
expression levels across the two peaks would provide insight into what
kinds of DBPs and promoter motifs preferentially increase expression.

We found that there are two distinct categories of promoters in terms of
DBP affinity. One set of promoters has low affinity, and another has
high affinity. The distinctiveness of this bimodal distribution is
conserved across promoter window size and biotype. In conclusion, there
are many directions that this research could be taken in. This includes
characterizing the two distinct categories of promoters and looking at
how expression changes between these two different categories.
