library(tidyverse)
library(ggplot2)

#edited_all.txt is numbers summarized extracted from reads_collapsed.fa file output by from miRDeep collapse reads script
pre_mirdeep <- read.delim("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/edited_all_reads.txt", header = TRUE, sep=" ") %>%
  filter(index.length <= 45) %>%
  mutate(prop2 = prop*100)

ggplot(pre_mirdeep, aes(x=index.length, y=prop2, fill=start_nt)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "black")) +
  xlab("Read Length") +
  ylab("Percent of Total Reads (%)") +
  labs(fill = "Nucleotide") +
  theme(legend.position="top")


  

