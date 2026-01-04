library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)

### Aedes aegypti siRNA (<22nt) and piRNA (>22nt) mappers to viral genomes  # note that we did not analyse viral siRNAs in this study.
table=read.table("Aedes.aegypti.sRNAseq.virus-mappers.counts", header=F)
table_d=dcast(table,V1~V2,value.var="V3")
table_d[is.na(table_d)] <- 0

# Subset data for Liao ning virus segment 1 (NC_007736.1) to segment 12 (NC_007747.1) 
subset_data <- table_d %>%
  filter(grepl("^NC_00773[6-9]\\.1@p@|^NC_00774[0-7]\\.1@p@", V1)) %>%
  select(V1, `Aedes_aegypti_ovaries_sRNA_rep1`, `Aedes_aegypti_eggs_sRNA_rep1`, `Aedes_aegypti_ovaries_sRNA_rep2`, `Aedes_aegypti_eggs_sRNA_rep1`) %>%
  pivot_longer(cols = -V1, names_to = "sample", values_to = "value") %>%
  mutate(
    label = sub(".*@p@", "", V1),           # extract p+ / p-
    log_value = log10(value + 1)            # optional log transform
  )

# Sort rows: @p@- first, then @p@+
subset_data <- subset_data %>%
  mutate(
    V1 = factor(V1, levels = unique(V1[order(label)])),
    sample = factor(sample, levels = c(
      "Aedes_aegypti_ovaries_sRNA_rep1",
      "Aedes_aegypti_eggs_sRNA_rep1",
      "Aedes_aegypti_ovaries_sRNA_rep2",
      "Aedes_aegypti_eggs_sRNA_rep1"
    ))
  )

s<-ggplot(subset_data, aes(x = sample, y = V1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("navy", "skyblue", "yellow", "red"),
    values = scales::rescale(c(0, 2, 5, 100)),
    name = "CPM"
  ) +
  labs(
    x = NULL,
    y = "Viral ID",
    title = "Expression heatmap (custom gradient)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file="Aaeg_ovary-vs-egg_virus-mappers_LNV_heatmap.pdf",width=6,height=6)
s
dev.off()
