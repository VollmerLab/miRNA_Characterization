library(tidyverse)
library(broom)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(BipartiteModularityMaximization)
library(bipartite)
library(rgexf)

#read in bonafides file
bonafides_acerv <- read_csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/bonafides_master.csv")
colnames(bonafides_acerv)[1] <- "microRNA"

#read in PITA results
results_file_3 <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/pita_results.tab", header=TRUE)

acerv_gff <- read.gff("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/k2_structuralAnnotations_wUTRs.gff3")

#remove 5' UTRs
UTR_keep <- acerv_gff %>%
  filter(type == "three_prime_UTR" | type == "CDS") %>%
  mutate(number = paste(start, end, sep = "-")) %>%
  mutate(UTR = paste(seqid, number, sep = ":")) %>%
  select(UTR, type)

#subsetting results file according to pita instructions (7-8 pair, 1 G:U wobbles, no mismatch)
results_file_sub_3 <- results_file_3 %>%
  filter(Seed == "7:0:0" | Seed == "7:0:1" |Seed == "8:0:0" | Seed == "8:0:1") %>%
  filter(ddG <= -10) %>%
  filter(microRNA %in% bonafides_acerv$microRNA) %>%
  filter(UTR %in% UTR_keep$UTR) %>%
  left_join(UTR_keep, by="UTR")

#process to get extended target matches
target_sequence_3 <- results_file_sub_3 %>%
  select(UTR, Start, End) %>%
  mutate(UTR2 = UTR) %>%
  separate(UTR2, into=c("seqid", "numbers"), sep = ":") %>%
  separate(numbers, into=c("start","end"), sep = "-") %>%
  mutate(seq_start = as.numeric(start)) %>%
  mutate(target_end = as.numeric(Start)) %>%
  mutate(target_start = as.numeric(End) - 18) %>%
  mutate(tool_start = seq_start + target_start) %>%
  mutate(tool_end = seq_start + target_end) %>%
  select(seqid, tool_start, tool_end)

write.table(target_sequence_3, "~/Desktop/miRNA_Characterization/miRNA_Characterization/Intermediate_Files/targets_3.bed", row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

#in bash
#bedtools getfasta  -fi k2_genome.fasta -bed targets_3.bed -fo targets_ext.fasta

fasta_to_table_pita <- function(input_file) {
  input_file %>%
    as_tibble() %>%
    dplyr::rename(variable = V1) %>%
    group_by(grp = str_c('variable', rep(1:2, length.out = n()))) %>%
    mutate(rn = row_number()) %>%
    ungroup %>%
    pivot_wider(names_from = grp, values_from = variable) %>%
    select(-rn) %>%
    dplyr::rename(utr_id = variable1, utr_target_sequence = variable2) %>%
    separate(col=utr_id, into=c("carrot", "utr_id"), sep=">") %>%
    select(-carrot)
}

write.table(utr[,1:2],"~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/targets_ext.txt")

utrs_3 <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/targets_ext.txt")
utr <- cbind(utrs_3, results_file_sub_3)


####
mirna_sequence <- bonafides_acerv %>%
  select(microRNA, miRNA_common_name, miRNA_category, consensus_mature_sequence) 

# file of miRNAs converted to T instead of U
bonafides_t <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/bonafides_t.txt") 
colnames(bonafides_t) <- "microRNA_sequence"

mirna_sequence2 <- mirna_sequence %>%
  select(microRNA, miRNA_common_name, miRNA_category) %>%
  cbind(bonafides_t)
###

library(stringi)
library(stringdist)

results_mirna <- results_file_sub_3 %>%
  select(UTR, microRNA, Start)

utr2 <- utr %>%
  select(UTR, utr_target_sequence, Start) %>%
  left_join(results_mirna, by=c("UTR", "Start")) %>%
  unique() %>%
  left_join(mirna_sequence2, by="microRNA") %>%
  mutate(utr_reverse = stri_reverse(utr_target_sequence)) %>%
  mutate(mirna_complement = chartr("ACGT", "TGCA", microRNA_sequence)) %>%
  mutate(mirna_to_match = substr(mirna_complement, 2, (nchar(mirna_complement) - 1))) %>%
  mutate(seed11_sequence = substr(mirna_to_match, 1, 13)) %>%
  mutate(leftover_seq = str_remove(mirna_to_match, seed11_sequence)) %>%
  mutate(mirna_len = str_length(mirna_to_match)) %>%
  mutate(target_seq = substr(utr_reverse, 2, mirna_len)) %>%
  mutate(utr11_sequence = substr(target_seq, 1, 13)) %>%
  mutate(utrleftover_seq = str_remove(target_seq, utr11_sequence)) %>%
  mutate(lev_good = stringdist(seed11_sequence, utr11_sequence, method = "hamming") <= 2) %>%
  mutate(cleave_good10 = stringdist(substr(seed11_sequence, 9, 9), substr(utr11_sequence, 9, 9), method = "hamming") == 0) %>%
  mutate(cleave_good11 = stringdist(substr(seed11_sequence, 10, 10), substr(utr11_sequence, 10, 10), method = "hamming") == 0) %>%
  filter(lev_good == TRUE) %>%
  filter(cleave_good10 != FALSE | cleave_good11 != FALSE)

good_targets <- utr2 %>%
  select(UTR, microRNA) %>%
  unique() %>%
  filter(microRNA %in% bonafides_acerv$microRNA)

kegg <- read_csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/kegg_annotations.csv.gz")

mRNA_annotations <- acerv_gff %>%
  filter(grepl("product=", attributes)) %>%
  separate(col=attributes, into=c("name_id","name_parent", "gene", "other"), sep=";") %>%
  select(-other) %>%
  separate(col=gene, into=c("remove", "gene_name"), sep="=") %>%
  select(-remove)

kegg_merge <- as_tibble(acerv_gff) %>%
  mutate(gene_id = str_extract(attributes, 'Acer_[0-9]+')) %>%
  filter(type == 'mRNA' | type == 'three_prime_UTR' | type == 'CDS') 

gene_table <- as.data.frame(table(kegg_merge$gene_id)) %>%
  filter(Freq != 1)

kegg_merge_3UTR <- kegg_merge %>%
  filter(gene_id %in% gene_table$Var1) %>%
  filter(type == 'mRNA') %>%
  select(seqid, gene_id, start, end) %>%
  left_join(kegg, by = 'gene_id')

results_file_annotations <- good_targets %>%
  separate(col=UTR, into=c("seqid", "range"), sep=":") %>%
  separate(col=range, into=c("start", "end"), sep="-") %>%
  select('seqid', 'start', 'end', 'microRNA') %>%
  distinct() %>%
  mutate_at(c("start", "end"), as.integer) %>%
  left_join(kegg_merge_3UTR, by = "seqid") %>%
  filter(start.x >= start.y & end.x <= end.y)

#Creating table of target hits
mirna_names <- bonafides_acerv %>%
  select(microRNA, miRNA_common_name, miRNA_category)

results_each_miRNA <- results_file_annotations %>%
  nest(.by = microRNA) %>%
  mutate(num_targets = map(data, ~ nrow(.))) %>%
  mutate(num_swiss_annotations = map(data, ~ length(na.omit(.$uniprot_id)))) %>%
  mutate(num_kegg_annotations = map(data, ~ length(na.omit(.$kegg_orthology)))) %>%
  left_join(mirna_names, by="microRNA") %>%
  filter(miRNA_common_name != "NA") 

target_numbers <- results_each_miRNA %>%
  dplyr::select(miRNA_common_name, num_targets, num_swiss_annotations, num_kegg_annotations, miRNA_category) %>%
  arrange(miRNA_category) %>%
  unnest(cols = c(num_targets, num_swiss_annotations, num_kegg_annotations))

mean(target_numbers$num_targets)
sd(target_numbers$num_targets)/sqrt(length(target_numbers$num_targets))

mean(target_numbers$num_swiss_annotations)
sd(target_numbers$num_swiss_annotations)/sqrt(length(target_numbers$num_swiss_annotations))

mean(target_numbers$num_kegg_annotations)
sd(target_numbers$num_kegg_annotations)/sqrt(length(target_numbers$num_kegg_annotations))

ac <- target_numbers %>% filter(miRNA_category == "acroporids") 
mean(ac$num_targets)
sd(ac$num_targets)/sqrt(length(ac$num_targets))

mean(ac$num_swiss_annotations)
sd(ac$num_swiss_annotations)/sqrt(length(ac$num_swiss_annotations))

mean(ac$num_kegg_annotations)
sd(ac$num_kegg_annotations)/sqrt(length(ac$num_kegg_annotations))

cni <- target_numbers %>% filter(miRNA_category == "cnidarian") 
mean(cni$num_targets)
sd(cni$num_targets)/sqrt(length(cni$num_targets))

mean(cni$num_targets)
sd(cni$num_targets)/sqrt(length(cni$num_targets))

mean(cni$num_swiss_annotations)
sd(cni$num_swiss_annotations)/sqrt(length(cni$num_swiss_annotations))

mean(cni$num_kegg_annotations)
sd(cni$num_kegg_annotations)/sqrt(length(cni$num_kegg_annotations))

uni <- target_numbers %>% filter(miRNA_category == "unique") 
mean(uni$num_targets)
sd(uni$num_targets)/sqrt(length(uni$num_targets))

mean(uni$num_targets)
sd(uni$num_targets)/sqrt(length(uni$num_targets))

mean(uni$num_swiss_annotations)
sd(uni$num_swiss_annotations)/sqrt(length(uni$num_swiss_annotations))

mean(uni$num_kegg_annotations)
sd(uni$num_kegg_annotations)/sqrt(length(uni$num_kegg_annotations))

#violin plot of target numbers
target_numbers$miRNA_category <- factor(target_numbers$miRNA_category, levels = c("unique", "acroporids", "cnidarian"))

total <- ggplot(target_numbers, aes(x=miRNA_category, y=num_targets, fill=miRNA_category)) + 
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Number Targets", y = "Counts of Targets") +
  ylim(0,350) +
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  labs(y="Counts of Targets") +
  ggtitle("Total Number of Targets") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE) 

swiss <- ggplot(target_numbers, aes(x=miRNA_category, y=num_swiss_annotations, fill=miRNA_category)) + 
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Swiss-Prot Annotations", y = "") +
  ylim(0,350) +
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  ggtitle("Number of Targets with Swiss-Prot Annotations") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE) 


kegg <- ggplot(target_numbers, aes(x=miRNA_category, y=num_kegg_annotations, fill=miRNA_category)) + 
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="KO Annotations", y = "") +
  ylim(0,350) +
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  ggtitle("Number of Targets with KO Annotations") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE) 


ggarrange(total, swiss, kegg, ncol = 3)


#read in target table with match number and region type
results_file_annotations_wmatch <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/results_file_annotations_wmatch.txt")
#merging amount of seed match

#comparing CDS to UTR targets
CDS_info <- results_file_annotations_wmatch %>%
  select(-kegg_path_id, -kegg_path_name, -minor_category, -major_category) %>%
  distinct() %>%
  #filter(type == "CDS") %>%
  select(seqid, start.x, end.x, type) %>%
  distinct()

CDSi <- CDS_info %>%
  filter(type == "CDS")

UTRi <- CDS_info %>%
  filter(type == "three_prime_UTR")

merges <- merge(CDSi, UTRi, by=c("seqid", "start.x", "end.x"))


results_file_annotations_wcat <- results_file_annotations %>%
  select(-kegg_pathway, -kegg_brite, -kegg_module) %>% 
  left_join(select(kegg_paths, gene_id, kegg_gene, kegg_orthology,
                   kegg_path_id, name, contains('category')) %>%
              rename(kegg_path_name = name), 
            by=c("gene_id", 'kegg_gene', 'kegg_orthology')) %>%
  left_join(mirna_names, by="microRNA")

results_file_annotations_wtype <- results_file_annotations_wcat %>%
  select(-microRNA, -miRNA_category) %>%
  mutate(UTR_sub = paste(seqid, start.x, sep=":")) %>%
  mutate(UTR = paste(UTR_sub, end.x, sep = "-")) %>%
  left_join(results_file_good, by=c("UTR", "miRNA_common_name"))



#do three-prime UTR sites have more miRNAs vs CDS
target_times <- results_file_annotations_wtype %>%
  select(-kegg_path_id, -kegg_path_name, -major_category, -minor_category) %>%
  unique() %>%
  group_by(UTR,type) %>%
  summarise(count = n()) 

library(car)
hits_glm2 <- glm(count == 1 ~ type, family = "binomial", data=target_times)
summary(hits_glm2)
car::Anova(hits_glm2)

#threeprimeUTR = 1, CDS = 0
#how to the count change for UTR relative to CDS target

#log odds count = 1.799 + -0.6754type
exp(coefficients(hits_glm2))

#the odds of only 1 miRNA hitting a target of a UTR is half that of the odds of hitting a target in CDS region.
confint(hits_glm2)

#chi-square proporition test for CDS vs UTR numbers
prop.test(x = c(3355, 1097), n = c(28059, 11168))

#create data file for network 
pathways_allannot <- results_file_annotations_wcat %>%
  as_tibble() %>%
  filter(!grepl("Human Diseases", major_category)) %>%
  select(microRNA, uniprot_id, kegg_orthology, gene_id) %>%
  distinct() %>%
  replace_na(list(kegg_orthology = "unknown", uniprot_id = "unknown")) %>%
  mutate(new_cat = case_when(kegg_orthology != "unknown" ~ kegg_orthology,
                             kegg_orthology == "unknown" & uniprot_id != "unknown" ~ uniprot_id,
                             kegg_orthology == "unknown" & uniprot_id == "unknown" ~ gene_id)) %>%
  select(microRNA, new_cat) %>%
  left_join(mirna_names, by="microRNA") %>%
  select(-microRNA) %>%
  drop_na() %>%
  group_by(miRNA_common_name, new_cat) %>%
  summarise(cat_count = n()) %>%
  ungroup %>%
  pivot_wider(names_from=new_cat, values_from = cat_count) %>%
  replace(is.na(.), 0) %>%
  remove_rownames %>% 
  column_to_rownames(var="miRNA_common_name") %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) %>%
  t() %>%
  as.data.frame() %>%
  #each pathway has to have at least 2 hits
  filter(.[,ncol(.)] >= 2) %>%
  select(-(last_col(offset = 0):last_col())) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) %>%
  t() %>%
  as.data.frame() %>%
  #each miRNA has hit at least 1 category
  filter(.[,ncol(.)] >= 2) %>%
  select(-(last_col(offset = 0):last_col())) %>%
  as.matrix() %>%
  t()

#network measurements with bipartite
d <- dfun(t(pathways_allannot))
specializes <- data.frame(miRNA_common_name = attr(d$dprime, which="names"), dprime = d$dprime) %>%
  left_join(bonafides_acerv, by="miRNA_common_name") %>%
  select(miRNA_common_name, miRNA_category, dprime) %>%
  mutate(dprime_round = round(dprime, 2)) 

specializes$miRNA_category <- factor(specializes$miRNA_category, levels = c("unique", "acroporids", "cnidarian"))

between_degree <- specieslevel(pathways_allannot, level="higher", index=c("betweenness", "degree")) %>%
  rownames_to_column(var ="miRNA_common_name") %>%
  left_join(bonafides_acerv, by="miRNA_common_name") %>%
  select(miRNA_common_name, miRNA_category, betweenness, degree) 

between_degree$miRNA_category <- factor(between_degree$miRNA_category, levels = c("unique", "acroporids", "cnidarian"))

spec_v <- specializes %>%  
  ggplot(aes(x=fct_infreq(miRNA_category), y=dprime,fill=miRNA_category)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  labs(y="Degree of Specialization (d')") +
  ggtitle("Specialization") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE)

bet_v <- between_degree %>%  
  ggplot(aes(x=fct_infreq(miRNA_category), y=betweenness,fill=miRNA_category)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  labs(y="Betweenness") +
  ggtitle("Betweenness") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE)


deg_v <- between_degree %>%  
  ggplot(aes(x=fct_infreq(miRNA_category), y=degree,fill=miRNA_category)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, fill="white")+
  ylim(0,140)+
  scale_fill_manual(values=cols) + 
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  labs(y="Degree") +
  ggtitle("Degree") +
  scale_x_discrete(labels=c("unique" = "Unique", "acroporids" = "Acroporid",
                            "cnidarian" = "Cnidarian")) +
  guides(fill=FALSE)

ggarrange(deg_v, bet_v, spec_v, ncol = 3)

#anovas on centrality measures
between_degree$bet_0 <- between_degree$betweenness + 0.000001
between_degree$log_bet <- log(between_degree$bet_0)
between_degree$log_degree <- log(between_degree$degree)
specializes$log_spec <- log(specializes$dprime)

hist(specializes$log_spec)
hist(between_degree$log_bet)
hist(between_degree$log_degree)

spec_aov <- aov(log_spec ~ miRNA_category,
             data = specializes
)

bet_aov <- aov(log_bet ~ miRNA_category,
             data = between_degree
)

deg_aov <- aov(log_degree ~ miRNA_category,
             data = between_degree
)

#create modules
set.seed(1234)
bip2 <- bipmod(pathways_allannot, ITER = 10)
bip2$MODULARITY

targets2 <- bip2$ASSIGN[1:743]
mis2 <- bip2$ASSIGN[744:810]

targets_dat <- data.frame(module = targets2, target_id = rownames(pathways_allannot)) %>%
  arrange(module)

#creating modules from network
module1_ko <- targets_dat %>%
  filter(module == 1) %>%
  select(target_id)

module2_ko <- targets_dat %>%
  filter(module == 2) %>%
  select(target_id)

module3_ko <- targets_dat %>%
  filter(module == 3) %>%
  select(target_id)

module4_ko <- targets_dat %>%
  filter(module == 4) %>%
  select(target_id)

module5_ko <- targets_dat %>%
  filter(module == 5) %>%
  select(target_id)

module6_ko <- targets_dat %>%
  filter(module == 6) %>%
  select(target_id)

module7_ko <- targets_dat %>%
  filter(module == 7) %>%
  select(target_id)

module8_ko <- targets_dat %>%
  filter(module == 8) %>%
  select(target_id)

module9_ko <- targets_dat %>%
  filter(module == 9) %>%
  select(target_id)

module10_ko <- targets_dat %>%
  filter(module == 10) %>%
  select(target_id)

module11_ko <- targets_dat %>%
  filter(module == 11) %>%
  select(target_id)

module12_ko <- targets_dat %>%
  filter(module == 12) %>%
  select(target_id)

module13_ko <- targets_dat %>%
  filter(module == 13) %>%
  select(target_id)

module14_ko <- targets_dat %>%
  filter(module == 14) %>%
  select(target_id)

module15_ko <- targets_dat %>%
  filter(module == 15) %>%
  select(target_id)

module16_ko <- targets_dat %>%
  filter(module == 16) %>%
  select(target_id)

mis_dat <- data.frame(module = mis2, mi_id = colnames(pathways_allannot)) %>%
  arrange(module)

module1_mi <- mis_dat %>%
  filter(module == 1) %>%
  select(mi_id)

module2_mi <- mis_dat %>%
  filter(module == 2) %>%
  select(mi_id)

module3_mi <- mis_dat %>%
  filter(module == 3) %>%
  select(mi_id)

module4_mi <- mis_dat %>%
  filter(module == 4) %>%
  select(mi_id)

module5_mi <- mis_dat %>%
  filter(module == 5) %>%
  select(mi_id)

module6_mi <- mis_dat %>%
  filter(module == 6) %>%
  select(mi_id)

module7_mi <- mis_dat %>%
  filter(module == 7) %>%
  select(mi_id)

module8_mi <- mis_dat %>%
  filter(module == 8) %>%
  select(mi_id)

module9_mi <- mis_dat %>%
  filter(module == 9) %>%
  select(mi_id)

module10_mi <- mis_dat %>%
  filter(module == 10) %>%
  select(mi_id)

module11_mi <- mis_dat %>%
  filter(module == 11) %>%
  select(mi_id)

module12_mi <- mis_dat %>%
  filter(module == 12) %>%
  select(mi_id)

module13_mi <- mis_dat %>%
  filter(module == 13) %>%
  select(mi_id)

module14_mi <- mis_dat %>%
  filter(module == 14) %>%
  select(mi_id)

module15_mi <- mis_dat %>%
  filter(module == 15) %>%
  select(mi_id)

module16_mi <- mis_dat %>%
  filter(module == 16) %>%
  select(mi_id)

#over-representation analysis
ora_test <- function(x, k, m, N, direction = 'two.sided'){
  ##https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
  data.frame(significant = c(x, m - x),
             not_significant = c(k - x, N - m - (k - x))) %>%
    fisher.test(alternative = direction, simulate.p.value=T) %>%
    tidy
}

#number for each annotation
#x = number in hub miRNA module
#m = total number in all modules
#N = total for entire dataset all anotations 
#k = total for entire dataset all anotations in hub module

results_by_annot <-  results_file_annotations_wcat %>%
  as_tibble() %>%
  filter(!grepl("Human Diseases", major_category)) %>%
  select(microRNA, uniprot_id, kegg_orthology, gene_id, kegg_path_name) %>%
  distinct() %>%
  replace_na(list(kegg_orthology = "unknown", uniprot_id = "unknown")) %>%
  mutate(new_cat = case_when(kegg_orthology != "unknown" ~ kegg_orthology,
                             kegg_orthology == "unknown" & uniprot_id != "unknown" ~ uniprot_id,
                             kegg_orthology == "unknown" & uniprot_id == "unknown" ~ gene_id)) %>%
  mutate(module = case_when(new_cat %in% module1_ko$target_id ~ "module1",
                            new_cat %in% module2_ko$target_id ~ "module2",
                            new_cat %in% module3_ko$target_id ~ "module3",
                            new_cat %in% module4_ko$target_id ~ "module4",
                            new_cat %in% module5_ko$target_id ~ "module5",
                            new_cat %in% module6_ko$target_id ~ "module6",
                            new_cat %in% module7_ko$target_id ~ "module7",
                            new_cat %in% module8_ko$target_id ~ "module8",
                            new_cat %in% module9_ko$target_id ~ "module9",
                            new_cat %in% module10_ko$target_id ~ "module10",
                            new_cat %in% module11_ko$target_id ~ "module11",
                            new_cat %in% module12_ko$target_id ~ "module12",
                            new_cat %in% module13_ko$target_id ~ "module13",
                            new_cat %in% module14_ko$target_id ~ "module14",
                            new_cat %in% module15_ko$target_id ~ "module15", 
                            new_cat %in% module16_ko$target_id ~ "module16")) %>%
  filter(module != "NA") %>%
  filter(kegg_path_name != "NA") %>%
  mutate(kegg_path_new = gsub("- fly", "", kegg_path_name)) %>%
  mutate(kegg_path_new2 = gsub("- animal", "", kegg_path_new)) %>%
  mutate(kegg_path_new3 = gsub("- other", "", kegg_path_new2)) %>%
  mutate(kegg_path_new4 = gsub("- yeast", "", kegg_path_new3)) %>%
  mutate(kegg_path_new5 = gsub("- worm", "", kegg_path_new4)) %>%
  nest(.by=kegg_path_new5) %>%
  mutate(num_genes = map(data, ~ nrow(.))) %>%
  #change to module you are interested in
  mutate(num_mod6 = map(data, ~ nrow(filter(., module == "module2")))) #%>%

total_mod6 <- rep(sum(as.numeric(results_by_annot$num_mod6)), nrow(results_by_annot))
total_gene <- rep(sum(as.numeric(results_by_annot$num_genes)), nrow(results_by_annot))

results_by_annot$total_mod6 <- total_mod6
results_by_annot$total_gene <- total_gene

ora_mod <- results_by_annot %>% 
  rowwise(kegg_path_new5) %>% 
  summarise(ora_test(x = num_mod6, m = num_genes, N = total_gene, k = total_mod6, 'greater'),
            .groups = 'drop') %>%
  mutate(p_adj = p.adjust(p.value, 'fdr'))

ora_mod %>% filter(p_adj <= 0.05)

ora_keggs <- results_file_annotations_wmatch %>%
  filter(miRNA_common_name %in% module1_mi$mi_id) %>%
  filter(kegg_orthology %in% module1_ko$target_id)

#output file to visualize in gephi
path_geph_wtype <- results_file_annotations_wtype %>%
  as_tibble() %>%
  filter(!grepl("Human Diseases", major_category)) %>%
  select(microRNA, uniprot_id, kegg_orthology, gene_id, type, minor_category) %>%
  distinct() %>%
  replace_na(list(kegg_orthology = "unknown", uniprot_id = "unknown")) %>%
  mutate(new_cat = case_when(kegg_orthology != "unknown" ~ kegg_orthology,
                             kegg_orthology == "unknown" & uniprot_id != "unknown" ~ uniprot_id,
                             kegg_orthology == "unknown" & uniprot_id == "unknown" ~ gene_id)) %>%
  select(microRNA, new_cat, type, minor_category) %>%
  left_join(bonafide_names, by="microRNA") %>%
  select(-microRNA) %>%
  replace_na(list(minor_category = "unknown")) %>%
  mutate(new_cat2 = case_when(grepl("Acer", new_cat) ~ "unknown",
                              grepl("up:", new_cat) ~ "swissprot",
                              grepl("ko:", new_cat) & minor_category == "unknown" ~ "KEGG", 
                              grepl("ko:", new_cat) & minor_category != "Immune system" ~ "KEGG",
                              grepl("ko:", new_cat) & minor_category == "Immune system" ~ "Immune KEGG")) %>%
  select(miRNA_common_name, new_cat, miRNA_category, new_cat2) %>%
  distinct()

path_input <- path_geph_wtype %>%
  select(miRNA_common_name, new_cat)

nodes <- data.frame(Id = path_geph_wtype$miRNA_common_name, Category = path_geph_wtype$new_cat2) %>%
  distinct() %>%
  mutate(kind = rep("miRNA", nrow(.))) 

edges <- data.frame(Id = path_geph_wtype$new_cat, Category = path_geph_wtype$new_cat2) %>%
  distinct() %>%
  mutate(kind = rep("target", nrow(.)))

edges$Category <- as.character(edges$Category)

new_data <- data.frame(Id = c(nodes$Id, edges$Id), category = c(nodes$Category, edges$Category), kind = c(nodes$kind, edges$kind))

new_data2 <- data.frame(Id = new_data$Id, label = new_data$Id )

category <- as.data.frame(new_data$category)
colnames(category) <- "category"

kind <- as.data.frame(new_data$kind)
colnames(kind) <- "kind"

attributes <- cbind(category, kind)

write.gexf(new_data2, path_input, nodesAtt=attributes, output="~/Desktop/miRNA_Characterization/miRNA_Characterization/Intermediate_Files/gephi_file_wtype3.gexf")
  
        