pita_analysis
================
Brecia Douglas
2024-02-09

Load libraries

``` r
library(tidyverse)
library(ape)
library(tidyr)
library(dplyr)
library(kableExtra)
```

``` r
results_file <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/tst_A_cerv_pita_results.tab", header=TRUE)

#targets file takes results file and collapses all miRNAs sites that have a target of a single UTR and sums the energy score for that miRNA
#targets_file <- read.table("~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/tst_A_cerv_pita_results_targets.tab", header=TRUE)

#read in bonafides file
bonafides_acerv <- read_csv("~/Desktop/mirna_1Nov/bonafides_master.csv")

#subsetting results file according to pita instructions
results_file_sub <- results_file %>%
  filter(Seed == "7:0:0" | Seed == "7:0:1" |Seed == "8:0:0" | Seed == "8:0:1") %>%
  filter(ddG <= -10) %>%
  filter(microRNA %in% bonafides_acerv$provisional_id) 
```

Getting categories of pathways

``` r
kegg_annotations <- read_csv('~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/kegg_annotations.csv.gz', show_col_types = FALSE)


#### KEGG Pathway Identification ####
if(file.exists('~/Desktop/RNA_Analysis_Lab/RNA_Analysis/PITA_files/kegg_orthogroup_pathways.csv.gz')){
  kegg_paths <- read_csv('~/Desktop/RNA_Analysis_Lab/RNA_Analysis/PITA_files/kegg_orthogroup_pathways.csv.gz', show_col_types = FALSE)
} else {
  kegg_paths <- kegg_annotations %>%
    filter(!is.na(kegg_orthology)) %>%
    pivot_longer(cols = c(starts_with('kegg'), -kegg_orthology, -kegg_gene),
                 names_prefix = 'kegg_',
                 names_to = 'kegg_type',
                 values_to = 'kegg_path_id') %>%
    filter(!is.na(kegg_path_id)) %>%
    mutate(kegg_path_id = str_split(kegg_path_id, ';;')) %>%
    unnest(kegg_path_id) %>%
    filter(kegg_type == 'pathway') %>%
    
    nest_by(kegg_type, kegg_path_id) %>%
    mutate(kegg_api = possibly(keggGet, otherwise = list(NULL), quiet = FALSE)(kegg_path_id)) %>%
    
    unnest_wider(kegg_api) %>%
    select(-contains(c('ENTRY', 'REFERENCE', 'DBLINKS', 'DISEASE', 'MODULE', 'DRUG', 'KO_PATHWAY'))) %>%
    mutate(name_length = lengths(NAME)) %>%
    rowwise %>%
    mutate(DESCRIPTION = str_c(DESCRIPTION, collapse = ';;;;'),
           DESCRIPTION = if_else(DESCRIPTION == "", NA_character_, DESCRIPTION),
           
           #May not want to do this - related pathways - network??
           REL_PATHWAY = str_c(REL_PATHWAY, collapse = ';;;;'),
           REL_PATHWAY = if_else(REL_PATHWAY == "", NA_character_, REL_PATHWAY)) %>%
    mutate(NAME = case_when(name_length == 2 ~ list(str_c(unlist(NAME), collapse = ': ')), 
                            name_length == 0 ~ list(NA_character_),
                            TRUE ~ list(NAME))) %>%
    ungroup %>%
    select(-name_length) %>%
    unnest(NAME) %>%
    janitor::clean_names() %>%
    separate(class, sep = '; ', into = c('major_category', 'minor_category')) %>%
    unnest(data)
  
  write_csv(kegg_paths, '~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/kegg_orthogroup_pathways.csv.gz')
}
```

``` r
acerv_gff <- read.gff("~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/k2_structuralAnnotations_wUTRs.gff3")

 
table(acerv_gff$type)
```

    ## 
    ##             CDS            exon  five_prime_UTR            gene            mRNA 
    ##          170894          182660           11532           33794           28059 
    ##     start_codon      stop_codon three_prime_UTR            tRNA 
    ##           26951           26364           11168            5735

``` r
#33,794 total 28,059 validated genes
#11,168 three prime UTRs

kegg <- read_csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/kegg_annotations.csv.gz")

mRNA_annotations <- acerv_gff %>%
  filter(type == "mRNA") %>%
  filter(grepl("product=", attributes)) %>%
  separate(col=attributes, into=c("name_id","name_parent", "gene", "other"), sep=";") %>%
  select(-other) %>%
  separate(col=gene, into=c("remove", "gene_name"), sep="=") %>%
  select(-remove)

kegg_merge <- as_tibble(acerv_gff) %>%
  mutate(gene_id = str_extract(attributes, 'Acer_[0-9]+')) %>%
  filter(type == 'mRNA') %>%
  select(seqid, gene_id, start, end) %>%
  left_join(kegg, by = 'gene_id')

results_file_annotations <- results_file_sub %>%
  separate(col=UTR, into=c("seqid", "range"), sep=":") %>%
  separate(col=range, into=c("start", "end"), sep="-") %>%
  select('seqid', 'start', 'end', 'microRNA') %>%
  distinct() %>%
  mutate_at(c("start", "end"), as.integer) %>%
  left_join(kegg_merge, by = "seqid") %>%
  filter(start.x >= start.y & end.x <= end.y)

uniprot_id <- results_file_annotations %>%
  select(uniprot_id) %>%
  unique() %>%
  separate(uniprot_id, into = c("prefix", "identifier"), sep=":") %>%
  na.omit() %>%
  mutate(accession_id = paste0("accession:", identifier))

ids <- uniprot_id$identifier

write(ids, file = "~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/ids.txt")

uniprot_maps <- read.csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/pita_files/acerv_idmappings.csv") %>%
  mutate(uniprot_id = paste0("up:", Entry))

results_file_annotations_wswiss <- results_file_annotations %>%
  left_join(uniprot_maps, by = "uniprot_id")

results_file_annotations_wcat <- results_file_annotations %>%
  select(-kegg_pathway, -kegg_brite, -kegg_module) %>% 
  left_join(select(kegg_paths, gene_id, kegg_gene, kegg_orthology,
                   kegg_path_id, name, contains('category')) %>%
              rename(kegg_path_name = name), 
            by=c("gene_id", 'kegg_gene', 'kegg_orthology')) 
```

Creating table of target hits

``` r
#table for each miRNA
#          category   number targets    annotated   unannotated   swissprot   kegg    uniq kegg
#mirna 1   cnidarian
#mirna 2   unique

mirna_names <- bonafides_acerv %>%
  select(provisional_id, miRNA_common_name, miRNA_category) %>%
  rename(microRNA = provisional_id)

results_each_miRNA <- results_file_annotations_wswiss %>%
  nest(.by = microRNA) %>%
  mutate(num_targets = map(data, ~ nrow(.))) %>%
  mutate(num_swiss_annotations = map(data, ~ length(na.omit(.$Protein.names)))) %>%
  mutate(num_kegg_annotations = map(data, ~ length(na.omit(.$kegg_gene)))) %>%
  left_join(mirna_names, by="microRNA")

target_numbers <- results_each_miRNA %>%
  dplyr::select(miRNA_common_name, num_targets, num_swiss_annotations, num_kegg_annotations, miRNA_category) %>%
  arrange(miRNA_category)

kable(target_numbers)
```

| miRNA_common_name        | num_targets | num_swiss_annotations | num_kegg_annotations | miRNA_category |
|:-------------------------|:------------|:----------------------|:---------------------|:---------------|
| miRNA_7                  | 102         | 69                    | 24                   | acroporids     |
| miRNA_33                 | 306         | 228                   | 62                   | acroporids     |
| miRNA_5_3p               | 845         | 600                   | 179                  | acroporids     |
| miRNA_2_3p               | 227         | 169                   | 64                   | acroporids     |
| miRNA_1                  | 85          | 62                    | 15                   | acroporids     |
| miRNA_13                 | 531         | 379                   | 127                  | acroporids     |
| miRNA_4_3p               | 220         | 166                   | 45                   | acroporids     |
| miRNA_5                  | 387         | 265                   | 79                   | acroporids     |
| miRNA_4                  | 134         | 104                   | 23                   | acroporids     |
| miRNA_8                  | 237         | 174                   | 56                   | acroporids     |
| miRNA_19                 | 233         | 170                   | 46                   | acroporids     |
| miRNA_19                 | 233         | 170                   | 46                   | acroporids     |
| miRNA_2                  | 170         | 126                   | 37                   | acroporids     |
| miRNA_17                 | 66          | 43                    | 11                   | acroporids     |
| miRNA_29                 | 57          | 42                    | 10                   | acroporids     |
| miRNA_29                 | 57          | 42                    | 10                   | acroporids     |
| miRNA_24                 | 47          | 36                    | 10                   | acroporids     |
| miRNA_10                 | 82          | 58                    | 18                   | acroporids     |
| miRNA_10                 | 82          | 58                    | 18                   | acroporids     |
| miRNA_2022               | 335         | 236                   | 83                   | cnidarian      |
| miRNA_2037               | 352         | 254                   | 93                   | cnidarian      |
| miRNA_2025               | 254         | 179                   | 55                   | cnidarian      |
| miRNA_100                | 71          | 48                    | 12                   | cnidarian      |
| miRNA_2023               | 356         | 227                   | 72                   | cnidarian      |
| miRNA_2030               | 34          | 22                    | 10                   | cnidarian      |
| miRNA_2036               | 62          | 44                    | 11                   | cnidarian      |
| miRNA_2050               | 524         | 364                   | 108                  | scleractinia   |
| miRNA_14                 | 211         | 141                   | 46                   | scleractinia   |
| miRNA_9425               | 107         | 78                    | 23                   | scleractinia   |
| Acerv_scaffold_112_40019 | 1072        | 783                   | 247                  | unique         |
| Acerv_scaffold_54_26073  | 159         | 106                   | 36                   | unique         |
| Acerv_scaffold_20_11433  | 236         | 161                   | 57                   | unique         |
| Acerv_scaffold_2_1277    | 294         | 219                   | 68                   | unique         |
| Acerv_scaffold_159_44966 | 130         | 107                   | 34                   | unique         |
| Acerv_scaffold_16_9855   | 122         | 92                    | 28                   | unique         |
| Acerv_scaffold_23_14072  | 159         | 112                   | 29                   | unique         |
| Acerv_scaffold_50_24648  | 427         | 296                   | 87                   | unique         |
| Acerv_scaffold_57_26966  | 227         | 154                   | 50                   | unique         |
| Acerv_scaffold_84_34634  | 223         | 152                   | 47                   | unique         |
| Acerv_scaffold_120_41029 | 459         | 338                   | 86                   | unique         |
| Acerv_scaffold_9_6107    | 73          | 52                    | 19                   | unique         |
| Acerv_scaffold_9_6292    | 73          | 52                    | 19                   | unique         |
| Acerv_scaffold_0_318     | 275         | 191                   | 51                   | unique         |
| Acerv_scaffold_25_14740  | 124         | 93                    | 33                   | unique         |
| Acerv_scaffold_33_18101  | 124         | 93                    | 33                   | unique         |
| Acerv_scaffold_14_8906   | 203         | 152                   | 53                   | unique         |
| Acerv_scaffold_16_9758   | 203         | 152                   | 53                   | unique         |
| Acerv_scaffold_227_47575 | 203         | 152                   | 53                   | unique         |
| Acerv_scaffold_73_32210  | 203         | 152                   | 53                   | unique         |
| Acerv_scaffold_57_27039  | 50          | 36                    | 7                    | unique         |
| Acerv_scaffold_8_5139    | 73          | 57                    | 16                   | unique         |
| Acerv_scaffold_23_14118  | 203         | 158                   | 51                   | unique         |
| Acerv_scaffold_116_40356 | 65          | 44                    | 21                   | unique         |
| Acerv_scaffold_34_18301  | 92          | 62                    | 19                   | unique         |
| Acerv_scaffold_14_8954   | 513         | 377                   | 128                  | unique         |
| Acerv_scaffold_3_2245    | 124         | 94                    | 27                   | unique         |
| Acerv_scaffold_151_44408 | 116         | 85                    | 25                   | unique         |
| Acerv_scaffold_36_19235  | 167         | 115                   | 41                   | unique         |
| Acerv_scaffold_79_33834  | 421         | 308                   | 113                  | unique         |
| Acerv_scaffold_51_24751  | 453         | 317                   | 105                  | unique         |
| Acerv_scaffold_15_9107   | 384         | 276                   | 76                   | unique         |
| Acerv_scaffold_28_16127  | 463         | 332                   | 104                  | unique         |
| Acerv_scaffold_28_15991  | 327         | 236                   | 76                   | unique         |
| Acerv_scaffold_47_23866  | 327         | 236                   | 76                   | unique         |
| Acerv_scaffold_3_1884    | 68          | 41                    | 12                   | unique         |
| Acerv_scaffold_28_16363  | 134         | 94                    | 27                   | unique         |
| Acerv_scaffold_1_839     | 317         | 234                   | 83                   | unique         |
| Acerv_scaffold_56_26942  | 88          | 55                    | 19                   | unique         |
| Acerv_scaffold_4_2537    | 95          | 73                    | 20                   | unique         |
| Acerv_scaffold_7_4123    | 98          | 69                    | 24                   | unique         |
| Acerv_scaffold_64_30361  | 43          | 31                    | 10                   | unique         |
| Acerv_scaffold_72_32054  | 140         | 104                   | 32                   | unique         |
| Acerv_scaffold_16_9509   | 89          | 58                    | 23                   | unique         |
| Acerv_scaffold_95_37471  | 89          | 58                    | 23                   | unique         |
| Acerv_scaffold_33_18099  | 57          | 39                    | 12                   | unique         |
| Acerv_scaffold_38_19930  | 143         | 109                   | 30                   | unique         |
| Acerv_scaffold_17_10230  | 53          | 36                    | 13                   | unique         |
| Acerv_scaffold_48_24090  | 27          | 20                    | 6                    | unique         |
| Acerv_scaffold_23_14008  | 112         | 73                    | 25                   | unique         |
| Acerv_scaffold_20_11778  | 111         | 86                    | 33                   | unique         |
| Acerv_scaffold_26_14789  | 17          | 12                    | 5                    | unique         |
| Acerv_scaffold_32_17956  | 29          | 20                    | 4                    | unique         |
