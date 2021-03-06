# 16S-Cas9-Seq

Julian Liber

Analysis of 16S-Cas9-Seq approach.

Three methodologies were tested:
1) 515F-806R no digestion
 - Samples were amplified using 16S V4 primers 515F-806R
 - Gel purification performed

2) 515F-806R with Cas9 digestion
- Samples were amplified using 16S V4 primers 515F-806R
- Gel purification performed
- Cas9 nuclease and 2 guide RNAs, targeting the host's chloroplast and mitochondria, were used to digest host amplicons
- Products were gel purified, reamplified with the same primers, and gel purified again.

3) 799F-1193R
- Samples were amplified using 16S V5-6 primers 799F-1193R
- The lower gel band was excised after gel electrophoresis to exclude host chloroplast and mitochondrial amplicons

All were products were normalized to 2 ng/uL and submitted for library preparation and MiSeq 2x250 sequencing using v2 500 chemistry.

Libraries were returned as demultiplexed and split paired-end FASTQ files. These were [filtered, trimmed, and denoised](https://github.com/liberjul/16S-Cas9-Seq/blob/main/scripts/01a_filter_trim_denoise.sh) in QIIME2, [exported](https://github.com/liberjul/16S-Cas9-Seq/blob/main/scripts/02a_export.sh), then CONSTAX2 was used to [classify](https://github.com/liberjul/16S-Cas9-Seq/blob/main/scripts/03a_classify_constax.sh) ASVs with SILVA138 SSURef database. Analysis and plotting was [performed in R](https://github.com/liberjul/16S-Cas9-Seq/blob/main/scripts/04a_ecological_analysis.R).

Some brief results:

### Read counts were notably different between methods
<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/read_count_boxplot.png" alt="Boxplot showing that Cas9 sequencing produced the highest reads, followed by no digest 515F-806r and yet lower 799F-1193R" width="60%"/>

### Cas9 did not effectively remove host chloroplast and mitochondria
<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/high_level_barplot_with_host.png" alt="Stacked barchart showing that Cas9 sequencing and no-digest 515F-806R libraries were dominated by chloroplast reads, while 799F-1193R showed decently diverse libraries." width="60%"/>

### Cas9 appears to have reduced the amount of mitochondrial reads contaminating the samples
<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/high_level_boxplot_with_host.png" alt="Boxplot showing that Cas9 sequencing and no-digest 515F-806R libraries were dominated by chloroplast reads, while 799F-1193R showed decently diverse libraries. Notably, mitochondria (but also Proteobacteria) are less abundant in the Cas9 digested libraries" width="60%"/>

### After removing host reads, both no-digest and Cas9 515F-806R libraries had extremely low microbial diversity
#### Phylum-level
<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/high_level_barplot_without_host.png" alt="Stacked barchart showing that Cas9 sequencing and no-digest 515F-806R libraries were dominated by Proteobacteria reads, while 799F-1193R showed decently diverse libraries." width="60%"/>

#### Genus-level

Genus was inferred using CONSTAX2 classifications, which allows for comparison across libraries which were independently processed and thus ASVs are not equivalent between libraries.

<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/genus_barchart.png" alt="Stacked barchart showing that Cas9 sequencing and no-digest 515F-806R libraries were dominated by Pseudomonas reads, while 799F-1193R showed decently diverse libraries." width="80%"/>

<img src="https://github.com/liberjul/16S-Cas9-Seq/blob/main/plots/genus_level_diversity.png" alt="Boxplot showing that Cas9 sequencing and no-digest 515F-806R libraries were extremely low in diversity across multiple metrics, while 799F-1193R showed decently diverse libraries." width="60%"/>
