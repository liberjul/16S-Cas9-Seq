library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)
library(rbiom)
library(vegan)
library(RVAideMemoire)
library(indicspecies)
library(ggpubr)
library(MASS)
library(bbmle)
library(lmtest)

setwd("C:/Users/julia/Dropbox (Duke Bio_Ea)/He Lab/Julian_Liber/16S-Cas9-Seq")

otu_dat_nd58 <- read.biom("./exported-feature-table/table_1_nd58/feature-table.biom")$counts
otu_dat_nd58_mat <- as.matrix(otu_dat_nd58)
otu_dat_nd58_df <- as.data.frame(otu_dat_nd58_mat)%>%
  rownames_to_column(var = "OTU_ID") %>%
  mutate(OTU_name = paste("ASV", row_number(), sep="_"))

tax_dat_nd58 <- read.delim("./exported-feature-table/rep-seqs_1_nd58/constax_taxonomy.txt", stringsAsFactors = F)
comb_nd58 <- left_join(otu_dat_nd58_df, tax_dat_nd58, by = "OTU_ID")

otu_dat_cas9d58 <- read.biom("./exported-feature-table/table_1_cas9d58/feature-table.biom")$counts
otu_dat_cas9d58_mat <- as.matrix(otu_dat_cas9d58)
otu_dat_cas9d58_df <- as.data.frame(otu_dat_cas9d58_mat)%>%
  rownames_to_column(var = "OTU_ID") %>%
  mutate(OTU_name = paste("ASV", row_number() + dim(otu_dat_nd58_df)[1], sep="_"))

tax_dat_cas9d58 <- read.delim("./exported-feature-table/rep-seqs_1_cas9d58/constax_taxonomy.txt", stringsAsFactors = F)
comb_cas9d58 <- left_join(otu_dat_cas9d58_df, tax_dat_cas9d58, by = "OTU_ID")

otu_dat_gp711 <- read.biom("./exported-feature-table/table_1_gp711/feature-table.biom")$counts
otu_dat_gp711_mat <- as.matrix(otu_dat_gp711)
otu_dat_gp711_df <- as.data.frame(otu_dat_gp711_mat) %>%
  rownames_to_column(var = "OTU_ID") %>%
  mutate(OTU_name = paste("ASV", row_number() + dim(otu_dat_nd58_df)[1] + dim(otu_dat_cas9d58_df)[1], sep="_"))

tax_dat_gp711 <- read.delim("./exported-feature-table/rep-seqs_1_gp711/constax_taxonomy.txt", stringsAsFactors = F)
comb_gp711 <- left_join(otu_dat_gp711_df, tax_dat_gp711, by = "OTU_ID")

comb_tax <- rbind(tax_dat_nd58, tax_dat_cas9d58, tax_dat_gp711)
comb_otus <- left_join(otu_dat_nd58_df, otu_dat_cas9d58_df, by = "OTU_ID") %>%
  left_join(otu_dat_gp711_df, by="OTU_ID")

comb_all <- full_join(comb_nd58, comb_cas9d58,
                      by=c("OTU_name", "OTU_ID", colnames(comb_nd58)[9:20])) %>%
  full_join(comb_gp711,
            by=c("OTU_name", "OTU_ID", colnames(comb_nd58)[9:20]))
comb_all
method_df <- data.frame(sample = 1:19, method = factor(c(rep("No digest 515F-806R", 6), 
                                                         rep("Cas9 digest 515F-806R", 6),
                                                         rep("Gel purify 799F-1193R", 7)),
                                                       levels = c("No digest 515F-806R",
                                                                  "Cas9 digest 515F-806R",
                                                                  "Gel purify 799F-1193R")),
                        sample_name = c(rep(c("Endo-61", "Endo-62", "Endo-63", "Total-161", "Total-162", "Total-163"), 3),
                                        "Neg Con")) %>%
  mutate(library_name = paste("JL", str_pad(sample, 2, pad="0"), sep = "_")) %>%
  dplyr::select(library_name, sample_name, method)

hl_table <- comb_all %>%
  select(c(High_level_taxonomy, starts_with("JL_"))) %>%
  mutate(High_level_taxonomy = str_remove(High_level_taxonomy, "Bacteria_")) %>%
  group_by(High_level_taxonomy) %>%
  summarise(across(starts_with("JL_"), ~ sum(.x, na.rm = T)))
hl_table

hl_table_prop <- hl_table %>%
  mutate(across(starts_with("JL_"), ~ .x / sum(.x)))


## Read count analysis
comb_all %>%
  dplyr::select(starts_with("JL_")) %>%
  mutate(across(starts_with("JL_"), ~ replace_na(.x, 0))) %>%
  colSums() %>%
  as.data.frame() %>%
  rename("read_count" = ".") %>%
  rownames_to_column(var = "library_name") %>%
  left_join(method_df, by = "library_name") -> read_counts
var(read_counts$read_count)
mean(read_counts$read_count)
# negative binomial due to overdispersion (var >>> mean)
m1 <- glm.nb(read_count ~ method, data = read_counts)
summary(m1)

read_counts %>%
  ggplot(aes(x = method, y = read_count, color = method)) +
  geom_boxplot(alpha=0) +
  geom_jitter() +
  scale_y_continuous(limits = c(0, 1.4e5)) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Read count") -> read_count_boxplot
read_count_boxplot
ggsave("./plots/read_count_boxplot.png", read_count_boxplot)

hl_table_prop %>%
  pivot_longer(cols = starts_with("JL_"),
               names_to = "library_name") %>%
  left_join(method_df, by="library_name") %>%
  ggplot(aes(fill=High_level_taxonomy, y = value, x = sample_name)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(~method,
             switch = "x",
             scales = "free_x",
             space = "free_x") +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0),
        legend.position = "right") +
  labs(x = "Sample",
       y = "Proportion of reads",
       fill = "High-Level Classification") -> hl_bar_all
hl_bar_all
ggsave("./plots/high_level_barplot_with_host.png", hl_bar_all, width=8, height=4)

hl_table_prop %>%
  filter(! High_level_taxonomy %in% c("Chloroplast", "Mitochondria")) %>%
  pivot_longer(cols = starts_with("JL_"),
               names_to = "library_name") %>%
  left_join(method_df, by="library_name") %>%
  ggplot(aes(fill=High_level_taxonomy, y = value, x = sample_name)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(~method,
             switch = "x",
             scales = "free_x",
             space = "free_x") +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0),
        legend.position = "right") +
  labs(x = "Sample",
       y = "Proportion of reads",
       fill = "Phylum") -> hl_bar_nonhost
hl_bar_nonhost
ggsave("./plots/high_level_barplot_without_host.png", hl_bar_nonhost, width=8, height=4)

hl_table_prop %>%
  pivot_longer(cols = starts_with("JL_"),
               names_to = "library_name") %>%
  left_join(method_df, by="library_name") %>%
  ggplot(aes(fill=High_level_taxonomy, color = High_level_taxonomy, y = value, x = method)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Method",
       y = "Proportion of reads",
       fill = "High-Level Classification",
       color = "High-Level Classification") -> hl_boxplot_all
hl_boxplot_all
ggsave("./plots/high_level_boxplot_with_host.png", hl_boxplot_all)

## Analysis by Genus
gen_table <- comb_all %>%
  filter(! High_level_taxonomy %in% c("Chloroplast", "Mitochondria")) %>% # remove host
  select(c(Rank_6, starts_with("JL_"))) %>% # genus and counts
  mutate(Genus = str_remove(Rank_6, "_1"), # remove suffix
         Genus = str_replace(Genus, "Escherichia-shigella", "Escherichia-Shigella")) %>% # fix name
  group_by(Genus) %>% # group by genus
  summarise(across(starts_with("JL_"), ~ sum(.x, na.rm = T))) # sum by genus
gen_table$Genus[gen_table$Genus == ""] <- "Other" # No genus assigned = other
gen_table %>% # weighted genus counts (by reads per sample)
  mutate(across(starts_with("JL_"), ~ .x / sum(.x))) %>%
  dplyr::select(starts_with("JL_")) %>%
  rowSums() -> gen_counts
gen_table$Total <- gen_counts # add column to allow for sorting
gen_table %>%
  arrange(-Total) -> gen_sorted
# take top 10 named genera
if ("Other" %in% gen_sorted$Genus[1:10]){
  gen_sorted %>%
    filter(Genus != "Other") %>%
    .[1:10,] %>%
    dplyr::select(-Total) -> gen_top_10
  gen_sorted %>%
    .[12:dim(.)[1],] -> gen_else
  gen_sorted %>%
    filter(Genus == "Other") %>%
    rbind(., gen_else) -> gen_else
} else {
  gen_sorted %>%
    .[1:10,] %>%
    dplyr::select(-Total) -> gen_top_10
  gen_sorted %>%
    .[11:dim(.)[1],] -> gen_else
}
# combine others into single row
gen_else %>%
  mutate(Genus = "Other") %>%
  group_by(Genus) %>%
  summarize(across(starts_with("JL_"), ~ sum(.x, na.rm = T))) -> gen_else
gen_top <- rbind(gen_top_10, gen_else)

gen_top %>%
  mutate(across(starts_with("JL_"), ~ .x / sum(.x))) %>%
  pivot_longer(cols = starts_with("JL_"),
               names_to = "library_name") %>%
  left_join(method_df, by="library_name") %>%
  ggplot(aes(fill=Genus, color = Genus, y = value, x = factor(sample_name, levels = c("Endo-61", "Endo-62", "Endo-63", "Total-161", "Total-162", "Total-163", "Neg Con")))) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~method,
             switch = "x",
             scales = "free_x",
             space = "free_x") +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0),
        legend.position = "right") +
  labs(x = "Sample",
       y = "Proportion of reads",
       fill = "Genus",
       color = "Genus") -> genus_bar_top
genus_bar_top
ggsave("./plots/genus_barchart.png", genus_bar_top, width = 10, height = 4)

## now for diversity
# make a matrix with genera names as row names
gen_table %>%
  dplyr::select(starts_with("JL_")) %>%
  as.matrix() -> gen_mat
rownames(gen_mat) <- gen_table$Genus
# Create plotable df for diversity metrics
div_df <- data.frame(Shannon = diversity(t(gen_mat), index = "shannon"),
                     Simpson = diversity(t(gen_mat), index = "simpson"),
                     Richness = colSums(gen_mat > 0)) %>%
  rownames_to_column(var = "library_name") %>%
  left_join(method_df, by="library_name")
div_df %>%
  pivot_longer(cols=c(Shannon, Simpson, Richness), names_to = "Metric") %>%
  ggplot(aes(x = method, y = value, color = method)) +
  facet_wrap(. ~ Metric, scales = "free")+
  geom_boxplot() +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = -60, hjust = 0),
        legend.position = "none", plot.margin = margin(1, 50, 1, 1)) +
  labs(x = "Method",
       y = "Diversity at Genus Rank") -> div_boxplot
div_boxplot
ggsave("./plots/genus_level_diversity.png", div_boxplot)
