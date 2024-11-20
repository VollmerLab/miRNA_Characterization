library(tidyverse)
library(dplyr)
library(seqinr)
library(stringdist)
library(purrr)
library(ape)
library(msa)
library(ggplot2)

#Getting conserved miRNAs from mature sequences using criteria (Wheeler et al., 2009)
#1. nt 2-7 are exact match -> seed
#2. length is +/- 2
#3. remainder of alignment only can contain up to 3 mismatches

fasta_to_table <- function(input_file) {
  input_file %>%
    as_tibble() %>%
    dplyr::rename(variable = V1) %>%
    group_by(grp = str_c('variable', rep(1:2, length.out = n()))) %>%
    mutate(rn = row_number()) %>%
    ungroup %>%
    pivot_wider(names_from = grp, values_from = variable) %>%
    select(-rn) %>%
    dplyr::rename(mirna_id = variable1, mature_sequence = variable2) %>%
    separate(col=mirna_id, into=c("carrot", "mirna_id"), sep=">") %>%
    select(-carrot)
}

mirna_seq_files <- list.files(path = '~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/raw_mirna_seq_cnidarians/', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(species = str_remove(file, dirname(file)) %>% 
           str_remove('//') %>% str_remove('_mirnas_t.fa') %>% str_remove('_t.fa') %>% str_remove('praher_')) %>%
  rowwise() %>%
  summarise(fasta_to_table(read.table(file)), .groups = 'drop') %>%
  mutate(seed_sequence = substr(mature_sequence, 2, 7)) %>%
  mutate(leftover_seq = str_remove(mature_sequence, seed_sequence)) %>%
  mutate(mirna_len = str_length(mature_sequence)) %>%
  mutate(mirna_id2 = mirna_id) %>%
  separate(mirna_id2, into=c("species", "other_id"), sep = "-") %>%
  select(-other_id)

get_conserved <- function(specie, file) {
  target_species <- file %>% filter(species == specie)
  
  other_mirnas <- file %>% filter(species != specie)
  
  seed_seq <- target_species$seed_sequence
  len_seq <- target_species$mirna_len
  leftover_seq <- target_species$leftover_seq
  
  #nt 2-7 of seed match with 1 wobble to account for first nt switch
  find_seed <- function(seed){match_data <- filter(other_mirnas, stringdist(seed_sequence, seed, method = "lv") <= 1)}
  
  #length is within +/- 2
  find_len <- function(mir_length){match_data <- filter(other_mirnas, mirna_len <= mir_length + 2 & mirna_len >= mir_length - 2)}
  
  #no more than 3 mismatches in sequence without seed
  find_lev <- function(leftover){match_data <- filter(other_mirnas, stringdist(leftover, leftover_seq, method = "lv") <= 3)}
  
  match_list_seed <- lapply(seed_seq, find_seed)
  match_list_len <- lapply(len_seq, find_len)
  match_list_left <- lapply(leftover_seq, find_lev)
  
  target_species$seed_matches <- match_list_seed
  target_species$len_matches <- match_list_len
  target_species$lev_good <- match_list_left
  
  targets_species_final <- target_species %>%
    mutate(conserved = map2(seed_matches, lev_good, inner_join)) #%>%
}

stylophora <- get_conserved("spi", mirna_seq_files)
cervicornis <- get_conserved("ace", mirna_seq_files)
millepora <- get_conserved("ami", mirna_seq_files)
digitifera <- get_conserved("adi", mirna_seq_files)
nematostella <- get_conserved("nve", mirna_seq_files)
exaiptsia <- get_conserved("epa", mirna_seq_files)
carnea <- get_conserved("eca", mirna_seq_files)
callimorphus <- get_conserved("sca", mirna_seq_files)
senile <- get_conserved("mse", mirna_seq_files)
viridis <- get_conserved("avi", mirna_seq_files)
jardinei <- get_conserved("cja", mirna_seq_files)
coerulea <- get_conserved("hco", mirna_seq_files)
aurita <- get_conserved("aau", mirna_seq_files)
malayensis <- get_conserved("sma", mirna_seq_files)
esculentum <- get_conserved("res", mirna_seq_files)
hydra <- get_conserved("hma", mirna_seq_files)

all_species <- rbind(stylophora, cervicornis, millepora ,digitifera, nematostella , exaiptsia, carnea, callimorphus, senile, viridis, jardinei, coerulea, aurita, malayensis, esculentum, hydra) %>%
  filter(map_lgl(conserved, ~ nrow(.) != 0))

all_species_shared <- all_species %>%
  select(-seed_matches, -lev_good) %>%
  unnest(cols = conserved, names_sep = "_") %>%
  mutate(base_call = paste(species, mirna_id, sep = "_")) %>%
  mutate(conserved_call = paste(conserved_species, conserved_mirna_id, sep = "_")) %>%
  select(base_call, conserved_call) %>%
  pivot_wider(names_from = base_call, values_from = conserved_call, values_fill = NA, values_fn = list) %>%
  pivot_longer(everything()) %>% 
  mutate(value = map(value, `length<-`, max(lengths(value)))) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  unnest(everything())

write.csv(all_species_shared, "~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/all_species_shared.csv")

#Base tree that is then edited in FigTree
myTree <- ape::read.tree(text = "((Bilateria), (((Hma),((Sma, Aau), Res)), ((Hco), (((((Cja), (Spi), ((Ace), (Adi, Ami))))), (((Avi),(Epa, Mse)), (Eca, Sca, Nve))))));")
plot(myTree)

ape::write.nexus(myTree, file='~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/base_cnidarian_tree.nex')

#miRNA 100 alignment
mirna_100 <- readRNAStringSet("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/mir100u.fa", format="fasta")
mir100_aln <- msa(mirna_100)

msaPrettyPrint(mir100_aln, output="pdf", showNames="right", showNumbering = "none", showConsensus = "none", showLogo = "none", furtherCode=c("\\defconsensus{.}{lower}{upper}"), showLegend = FALSE, askForOverwrite = FALSE)

#acerv mature miRNA graph
acerv_mirnas <- all_mirdeep_files %>%
  filter(species == "adi") %>%
  mutate(mirna_length = str_length(consensus_mature_sequence)) %>%
  rowwise() %>%
  mutate(start_nt = str_sub(consensus_mature_sequence, 1, 1)) %>%
  select(mirna_length, start_nt) %>%
  table() %>%
  as.data.frame() %>%
  mutate(prop = (Freq/sum(Freq))*100)

ggplot(acerv_mirnas, aes(x=mirna_length, y=prop, fill=start_nt)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "black")) +
  xlab("Read Length") +
  ylab("Percent of Total Reads (%)") +
  labs(fill = "Nucleotide") +
  theme(legend.position="right")