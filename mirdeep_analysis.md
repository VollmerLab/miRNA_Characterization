mirdeep_analysis
================
Brecia Douglas
2024-02-02

Loading in R packages

``` r
library(tidyverse)
library(dplyr)
library(kableExtra)
library(UpSetR)
library(ggplot2)
library(ggupset)
library(ComplexUpset)
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
