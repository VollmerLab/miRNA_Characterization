mirdeep_analysis
================
Brecia Douglas
2024-02-02

Loading in R packages

``` r
library(tidyverse)
library(dplyr)
library(kableExtra)
library(seqinr)
library(stringdist)
library(purrr)
```

Read in mirdeep files, bonafides and find shared miRNAs

``` r
#need to add in 3 jellyfish
#read in bonafides file
bonafides_acerv <- read_csv("~/Desktop/mirna_1Nov/bonafides_master.csv")

all_mirdeep_files <- list.files(path = '~/Desktop/miRNA_Characterization/miRNA_Characterization/mirdeep_files/', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(species = str_remove(file, dirname(file)) %>% str_remove('//') %>% str_remove('result_') %>% str_remove('.csv')) %>%
  rowwise(species) %>%
  summarise(read_csv(file), .groups = 'drop') %>%
  filter(provisional_id %in% bonafides_acerv$provisional_id)

shared_mirnas <- all_mirdeep_files %>%
  select("species", "provisional_id", "example_miRBase_miRNA_with_the_same_seed") %>%
  filter(example_miRBase_miRNA_with_the_same_seed != "-") %>%
  pivot_wider(names_from = species, values_from = example_miRBase_miRNA_with_the_same_seed) %>%
  filter(provisional_id %in% bonafides_acerv$provisional_id)

common_mirnas <- shared_mirnas %>%
  mutate_all(funs(replace(., is.na(.), "-"))) %>%
  mutate(n_missing = apply(. == "-", 1, sum)) %>%
  filter(n_missing <= 8) 

kable(common_mirnas)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
provisional_id
</th>
<th style="text-align:left;">
adi
</th>
<th style="text-align:left;">
ami
</th>
<th style="text-align:left;">
avi
</th>
<th style="text-align:left;">
cja
</th>
<th style="text-align:left;">
eca
</th>
<th style="text-align:left;">
epa
</th>
<th style="text-align:left;">
hco
</th>
<th style="text-align:left;">
hma
</th>
<th style="text-align:left;">
mse
</th>
<th style="text-align:left;">
nve
</th>
<th style="text-align:left;">
sca
</th>
<th style="text-align:left;">
spi
</th>
<th style="text-align:right;">
n_missing
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Acerv_scaffold_1_1115
</td>
<td style="text-align:left;">
adi-mir-100_5p
</td>
<td style="text-align:left;">
ami-nve-f-mir-100-5p
</td>
<td style="text-align:left;">
avi-nve-f-mir-100
</td>
<td style="text-align:left;">

- </td>
  <td style="text-align:left;">
  eca-nve-f-mir-100
  </td>
  <td style="text-align:left;">
  epa-mir-100
  </td>
  <td style="text-align:left;">

  - </td>
    <td style="text-align:left;">

    - </td>
      <td style="text-align:left;">
      mse-nve-f-mir-100
      </td>
      <td style="text-align:left;">

      - </td>
        <td style="text-align:left;">
        sca-nve-f-mir-100
        </td>
        <td style="text-align:left;">
        spi-mir-temp-1
        </td>
        <td style="text-align:right;">
        4
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        Acerv_scaffold_49_24339
        </td>
        <td style="text-align:left;">
        adi-mir-2023_3p
        </td>
        <td style="text-align:left;">
        ami-nve-f-mir-2023-3p
        </td>
        <td style="text-align:left;">
        avi-nve-mir-2023
        </td>
        <td style="text-align:left;">

        - </td>
          <td style="text-align:left;">
          eca-nve-f-mir-2023
          </td>
          <td style="text-align:left;">
          epa-mir-2023
          </td>
          <td style="text-align:left;">
          hco_scaf371_26412
          </td>
          <td style="text-align:left;">

          - </td>
            <td style="text-align:left;">
            mse-nve-f-mir-2023
            </td>
            <td style="text-align:left;">
            nve-mir-2023-3p
            </td>
            <td style="text-align:left;">
            sca-nve-mir-2023-3p
            </td>
            <td style="text-align:left;">
            spi-mir-temp-4
            </td>
            <td style="text-align:right;">
            2
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            Acerv_scaffold_40_20559
            </td>
            <td style="text-align:left;">
            adi-mir-p-novel-1-3p
            </td>
            <td style="text-align:left;">
            ami-mir-p-novel-1-3p
            </td>
            <td style="text-align:left;">

            - </td>
              <td style="text-align:left;">

              - </td>
                <td style="text-align:left;">

                - </td>
                  <td style="text-align:left;">
                  epa-mir-12434
                  </td>
                  <td style="text-align:left;">

                  - </td>
                    <td style="text-align:left;">

                    - </td>
                      <td style="text-align:left;">

                      - </td>
                        <td style="text-align:left;">
                        nve-mir-9437
                        </td>
                        <td style="text-align:left;">

                        - </td>
                          <td style="text-align:left;">

                          - </td>
                            <td style="text-align:right;">
                            8
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Acerv_scaffold_151_44415
                            </td>
                            <td style="text-align:left;">
                            adi-mir-2030_5p
                            </td>
                            <td style="text-align:left;">
                            ami-nve-f-mir-2030-5p
                            </td>
                            <td style="text-align:left;">
                            avi-nve-f-mir-2030
                            </td>
                            <td style="text-align:left;">
                            cja-scaffold_14_292623
                            </td>
                            <td style="text-align:left;">
                            eca-nve-f-mir-2030
                            </td>
                            <td style="text-align:left;">
                            epa-mir-2030
                            </td>
                            <td style="text-align:left;">
                            hco_scaf68_13789
                            </td>
                            <td style="text-align:left;">
                            hma-mir-2030
                            </td>
                            <td style="text-align:left;">
                            mse-nve-f-mir-2030
                            </td>
                            <td style="text-align:left;">
                            nve-mir-2030-5p
                            </td>
                            <td style="text-align:left;">
                            sca-nve-f-mir-2030
                            </td>
                            <td style="text-align:left;">
                            spi-mir-temp-40
                            </td>
                            <td style="text-align:right;">
                            0
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;">
                            Acerv_scaffold_98_38275
                            </td>
                            <td style="text-align:left;">
                            adi-mir-9425_5p
                            </td>
                            <td style="text-align:left;">
                            ami-nve-f-mir-9425-5p
                            </td>
                            <td style="text-align:left;">

                            - </td>
                              <td style="text-align:left;">
                              cja-scaffold_3_48436
                              </td>
                              <td style="text-align:left;">

                              - </td>
                                <td style="text-align:left;">

                                - </td>
                                  <td style="text-align:left;">

                                  - </td>
                                    <td style="text-align:left;">

                                    - </td>
                                      <td style="text-align:left;">

                                      - </td>
                                        <td style="text-align:left;">
                                        nve-mir-9425
                                        </td>
                                        <td style="text-align:left;">

                                        - </td>
                                          <td style="text-align:left;">

                                          - </td>
                                            <td style="text-align:right;">
                                            8
                                            </td>
                                            </tr>
                                            <tr>
                                            <td style="text-align:left;">
                                            Acerv_scaffold_9_6161
                                            </td>
                                            <td style="text-align:left;">
                                            adi-mir-2025_3p
                                            </td>
                                            <td style="text-align:left;">
                                            ami-nve-f-mir-2025-3p
                                            </td>
                                            <td style="text-align:left;">
                                            avi-nve-f-mir-2025
                                            </td>
                                            <td style="text-align:left;">

                                            - </td>
                                              <td style="text-align:left;">
                                              eca-nve-f-mir-2025
                                              </td>
                                              <td style="text-align:left;">
                                              epa-mir-2025
                                              </td>
                                              <td style="text-align:left;">

                                              - </td>
                                                <td style="text-align:left;">

                                                - </td>
                                                  <td style="text-align:left;">
                                                  mse-nve-f-mir-2025
                                                  </td>
                                                  <td style="text-align:left;">
                                                  nve-mir-2025-3p
                                                  </td>
                                                  <td style="text-align:left;">
                                                  sca-nve-f-mir-2025
                                                  </td>
                                                  <td style="text-align:left;">

                                                  - </td>
                                                    <td style="text-align:right;">
                                                    4
                                                    </td>
                                                    </tr>
                                                    <tr>
                                                    <td style="text-align:left;">
                                                    Acerv_scaffold_1_691
                                                    </td>
                                                    <td style="text-align:left;">
                                                    adi-mir-2036_3p
                                                    </td>
                                                    <td style="text-align:left;">
                                                    ami-nve-f-mir-2036-3p
                                                    </td>
                                                    <td style="text-align:left;">
                                                    avi-apa-b-mir-2036
                                                    </td>
                                                    <td style="text-align:left;">

                                                    - </td>
                                                      <td style="text-align:left;">
                                                      eca-nve-f-mir-2036
                                                      </td>
                                                      <td style="text-align:left;">
                                                      epa-mir-2036
                                                      </td>
                                                      <td style="text-align:left;">

                                                      - </td>
                                                        <td style="text-align:left;">

                                                        - </td>
                                                          <td style="text-align:left;">
                                                          mse-nve-f-mir-2036
                                                          </td>
                                                          <td style="text-align:left;">
                                                          nve-mir-2036-3p
                                                          </td>
                                                          <td style="text-align:left;">
                                                          sca-nve-f-mir-2036
                                                          </td>
                                                          <td style="text-align:left;">
                                                          spi-mir-temp-30
                                                          </td>
                                                          <td style="text-align:right;">
                                                          3
                                                          </td>
                                                          </tr>
                                                          <tr>
                                                          <td style="text-align:left;">
                                                          Acerv_scaffold_43_22758
                                                          </td>
                                                          <td style="text-align:left;">
                                                          adi-mir-2022_3p
                                                          </td>
                                                          <td style="text-align:left;">
                                                          ami-nve-f-mir-2022-3p
                                                          </td>
                                                          <td style="text-align:left;">
                                                          avi-adi-mir-p-novel-8
                                                          </td>
                                                          <td style="text-align:left;">

                                                          - </td>
                                                            <td style="text-align:left;">
                                                            eca-nve-f-mir-2022a
                                                            </td>
                                                            <td style="text-align:left;">
                                                            epa-mir-2022a
                                                            </td>
                                                            <td style="text-align:left;">

                                                            - </td>
                                                              <td style="text-align:left;">
                                                              hma-mir-2022
                                                              </td>
                                                              <td style="text-align:left;">
                                                              mse-nve-f-mir-2022-3p
                                                              </td>
                                                              <td style="text-align:left;">
                                                              nve-mir-2022-3p
                                                              </td>
                                                              <td style="text-align:left;">

                                                              - </td>
                                                                <td style="text-align:left;">
                                                                spi-mir-temp-25
                                                                </td>
                                                                <td style="text-align:right;">
                                                                3
                                                                </td>
                                                                </tr>
                                                                <tr>
                                                                <td style="text-align:left;">
                                                                Acerv_scaffold_33_18099
                                                                </td>
                                                                <td style="text-align:left;">
                                                                adi-mir-2025_3p
                                                                </td>
                                                                <td style="text-align:left;">
                                                                ami-nve-f-mir-2025-3p
                                                                </td>
                                                                <td style="text-align:left;">
                                                                avi-nve-f-mir-2025
                                                                </td>
                                                                <td style="text-align:left;">

                                                                - </td>
                                                                  <td style="text-align:left;">
                                                                  eca-nve-f-mir-2025
                                                                  </td>
                                                                  <td style="text-align:left;">
                                                                  epa-mir-2025
                                                                  </td>
                                                                  <td style="text-align:left;">

                                                                  - </td>
                                                                    <td style="text-align:left;">

                                                                    - </td>
                                                                      <td style="text-align:left;">
                                                                      mse-nve-f-mir-2025
                                                                      </td>
                                                                      <td style="text-align:left;">
                                                                      nve-mir-2025-3p
                                                                      </td>
                                                                      <td style="text-align:left;">
                                                                      sca-nve-f-mir-2025
                                                                      </td>
                                                                      <td style="text-align:left;">

                                                                      - </td>
                                                                        <td style="text-align:right;">
                                                                        4
                                                                        </td>
                                                                        </tr>
                                                                        <tr>
                                                                        <td style="text-align:left;">
                                                                        Acerv_scaffold_11_7131
                                                                        </td>
                                                                        <td style="text-align:left;">

                                                                        - </td>
                                                                          <td style="text-align:left;">
                                                                          ami-nve-f-mir-2037-3p
                                                                          </td>
                                                                          <td style="text-align:left;">
                                                                          avi-nve-f-mir-2037
                                                                          </td>
                                                                          <td style="text-align:left;">

                                                                          - </td>
                                                                            <td style="text-align:left;">
                                                                            eca-nve-f-mir-2037
                                                                            </td>
                                                                            <td style="text-align:left;">
                                                                            epa-mir-2037
                                                                            </td>
                                                                            <td style="text-align:left;">

                                                                            - </td>
                                                                              <td style="text-align:left;">

                                                                              - </td>
                                                                                <td style="text-align:left;">
                                                                                mse-nve-f-mir-2037-3p
                                                                                </td>
                                                                                <td style="text-align:left;">
                                                                                nve-mir-2037-3p
                                                                                </td>
                                                                                <td style="text-align:left;">
                                                                                sca-nve-f-mir-2037-3p
                                                                                </td>
                                                                                <td style="text-align:left;">
                                                                                spi-mir-temp-20
                                                                                </td>
                                                                                <td style="text-align:right;">
                                                                                4
                                                                                </td>
                                                                                </tr>
                                                                                </tbody>
                                                                                </table>

Getting conserved miRNAs from mature sequences using criteria 1. nt 2-7
are exact match -\> seed 2. length is +/- 2 3. remainder of alignment
only can contain up to 3 mismatches

``` r
fasta_to_table <- function(input_file) {
  input_file %>%
  as_tibble() %>%
  dplyr::rename(variable = V1) %>%
  group_by(grp = str_c('variable', rep(1:2, length.out = n()))) %>%
  mutate(rn = row_number()) %>%
  ungroup %>%
  pivot_wider(names_from = grp, values_from = variable) %>%
  select(-rn) %>%
  rename(mirna_id = variable1, mature_sequence = variable2) %>%
  separate(col=mirna_id, into=c("carrot", "mirna_id"), sep=">") %>%
  select(-carrot)
}

mirna_seq_files <- list.files(path = '~/Desktop/miRNA_Characterization/miRNA_Characterization/raw_mirna_seq_cnidarians/', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(species = str_remove(file, dirname(file)) %>% 
           str_remove('//') %>% str_remove('_mirnas_t.fa') %>% str_remove('_t.fa') %>% str_remove('praher_')) %>%
  rowwise(species) %>%
  summarise(fasta_to_table(read.table(file)), .groups = 'drop') %>%
  mutate(seed_sequence = substr(mature_sequence, 2, 7)) %>%
  mutate(leftover_seq = str_remove(mature_sequence, seed_sequence)) %>%
  mutate(mirna_len = str_length(mature_sequence))

get_conserved <- function(specie, file) {
  target_species <- file %>% filter(species == specie)

  other_mirnas <- file %>% filter(species != specie)
  
  seed_seq <- target_species$seed_sequence
  len_seq <- target_species$mirna_len
  leftover_seq <- target_species$leftover_seq
  
  find_seed <- function(seed){match_data <- filter(other_mirnas, seed_sequence == seed)}
  
  find_len <- function(mir_length){match_data <- filter(other_mirnas, mirna_len <= mir_length + 2 & mirna_len >= mir_length - 2)}
  
  find_lev <- function(leftover){match_data <- filter(other_mirnas, stringdist(leftover, leftover_seq, method = "lv") < 3)}

  match_list_seed <- lapply(seed_seq, find_seed)
  match_list_len <- lapply(len_seq, find_len)
  match_list_left <- lapply(leftover_seq, find_lev)
  
  target_species$seed_matches <- match_list_seed
  target_species$len_matches <- match_list_len
  target_species$lev_good <- match_list_left

  targets_species_final <- target_species %>%
    mutate(temp = map2(seed_matches, len_matches, inner_join)) %>%
    mutate(conserved = map2(temp, lev_good, inner_join)) %>%
    select(-temp)
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
```
