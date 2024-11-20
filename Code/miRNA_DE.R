library(tidyverse)
library(edgeR)
library(lme4)
library(multidplyr)
library(lmerTest)
library(afex)
library(emmeans)
library(BiocParallel)
library(variancePartition)
library(magrittr)
library(qvalue)
library(ggh4x)

#miRNA differential expression with susceptibility
bonafides_acerv <- read_csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/bonafides_master.csv")
colnames(bonafides_acerv)[1] <- "microRNA"

miRNA_counts <- read.csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/mirna_counts.csv") %>%
  filter(microRNA %in% bonafides_acerv$microRNA)


samples_code <- read.csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/samples.csv")


miRNA_counts_mod <- miRNA_counts %>%
  as_tibble() %>%
  left_join(.,samples_code, by="sample_id") %>%
  dplyr::select(-c("sample_id")) %>%
  pivot_wider(names_from = sample_name, values_from = n) %>%
  replace(is.na(.), 0) 

#format counts table
miRNA_counts_mod_de <- miRNA_counts_mod %>%
  remove_rownames %>%
  column_to_rownames(var = "microRNA") %>%
  as.matrix(.)


#normalize count table
#create DGElist
counts_miR_dge <- DGEList(miRNA_counts_mod_de)

#calculate TMM normalization factors
counts_miR_norm_fac <- calcNormFactors(counts_miR_dge, method="TMM")

samples_names <- as.data.frame(colnames(miRNA_counts_mod_de))
colnames(samples_names) <- "sample"

sample_info <- samples_names %>%
  separate(sample,c("time", "exposure", "tank", "genotype", "susceptibility"), "_") %>%
  cbind(samples_names, .) %>%
  mutate(treat_outcome = case_when(exposure == 'H' ~ 'H',
                                   TRUE ~ str_c(exposure, susceptibility, sep = '_'))) %>% 
  mutate(final_disease_state = case_when(treat_outcome == "H" ~ "H",
                                         treat_outcome == "D_R" ~ "H",
                                         treat_outcome == "D_S" ~ "D",
                                         TRUE ~ as.character(treat_outcome))) %>%
  mutate(fragment_id = str_c(exposure, tank, genotype, final_disease_state, sep="_")) %>%
  mutate(across(where(is.character), as.factor))

# estimate weights using linear mixed model of dream
param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
vobjDream <- variancePartition::voomWithDreamWeights(counts = counts_miR_norm_fac, 
                                  formula = ~ time * exposure * susceptibility + 
                                    (1 | genotype) +
                                    (1 | fragment_id) +
                                    (1 | tank),
                                  
                                  data = sample_info %>%
                                    arrange(sample) %>%
                                    column_to_rownames('sample'),
                                  BPPARAM = param) 

model_data <- counts_miR_norm_fac %>% 
  cpm(log = TRUE, prior.count = 0.5,
      normalized.lib.sizes = TRUE) %>%
  as.data.frame %>%
  t %>%
  as_tibble(rownames = "sample") %>%
  full_join(sample_info, by="sample") %>%
  pivot_longer(cols = -any_of(colnames(sample_info)), 
               names_to = "miRNA_names", values_to = "value") %>%
  left_join(vobjDream$weights %>% 
              set_colnames(colnames(vobjDream$E)) %>%
              set_rownames(rownames(vobjDream$E)) %>%
              as_tibble(rownames = 'miRNA_names') %>%
              pivot_longer(cols = -miRNA_names,
                           names_to = 'sample',
                           values_to = 'weight'), by=c("miRNA_names", "sample")) %>%
  mutate(across(c(exposure, final_disease_state, susceptibility), factor)) 

cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 'stringr', 'tidyr'))
cluster_copy(cluster, c('alpha'))
safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

full_gene_models <- model_data %>%
  nest(data = -c(miRNA_names)) %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(model = list(lmer(value ~ time * exposure * susceptibility + 
                             (1 | genotype) + (1 | fragment_id) + 
                             (1 | tank),
                           weights = weight,
                           data = data)),
         
         anova_table = list(anova(model, 
                                  ddf = 'Kenward-Roger')),
         
         as_tibble(anova_table, rownames = 'term') %>%
           select(term, `Pr(>F)`) %>%
           pivot_wider(names_from = term, 
                       values_from = `Pr(>F)`,
                       names_prefix = 'p_') %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  collect() %>%
  ungroup %>% 
  mutate(across(starts_with('p_'),
                ~p.adjust(., method = 'fdr'), 
                .names = 'fdr_{.col}')) %>%
  rename_with(~str_replace(., 'fdr_p_', 'fdr_')) %>%
  mutate(across(starts_with('p_'), safe_qvalue,
                .names = 'q_{.col}')) %>%
  rename_with(~str_replace_all(., 'q_p_', 'q_'))

#graphing 3 significant miRNAs
sig_mirnas <- full_gene_models %>%
  filter(fdr_susceptibility<= 0.05 | fdr_timeXsusceptibility <= 0.05)

miRNA_33 <- emmeans(sig_mirnas$model[[1]], ~ time * susceptibility)
miRNA_2_3p <- emmeans(sig_mirnas$model[[2]], ~ time * susceptibility)
miRNA_2022 <- emmeans(sig_mirnas$model[[3]], ~ time * susceptibility)

line_graph <- sig_mirnas %>%
  mutate(posthoc = list(miRNA_33, miRNA_2_3p, miRNA_2022)) %>%
  select(miRNA_names, posthoc) %>%
  rowwise(miRNA_names) %>%
  reframe(as_tibble(posthoc)) %>%
  rename(microRNA = "miRNA_names") %>%
  left_join(bonafides_acerv, by="microRNA") %>%
  mutate(miRNA_common_name = factor(miRNA_common_name, levels = c("miRNA_33","miRNA_2-3p" ,"miRNA_2022"))) 

levels(line_graph$miRNA_common_name) <- c("miRNA-33", "miRNA-2-3p", "miRNA-2022")

line_graph %>%  
  ggplot(aes(x=time, y=emmean, group = susceptibility)) +
  geom_line(aes(color = susceptibility), position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, color = susceptibility), width=.1, position = position_dodge(width = 0.2)) +
  geom_point(aes(color = susceptibility), size = 5, fill = "white", position = position_dodge(width = 0.2)) + 
  #scale_shape_manual(values = c("No" = 21, "Yes" = 19)) +
  scale_color_manual(values = c("S" = "#e30e0e", "R" = "#22A7B6"), labels = c("Resistant", "Susceptible")) +
  labs(x="Time", y=expression(Normalized~Log[2]~(cpm))) +
  labs(color="Resistance") +
  scale_x_discrete(labels = c("T3" = "Day 3", "T7" = "Day 7")) +
  facet_grid2(.~ miRNA_common_name, independent = TRUE, scales = "free", labeller = label_wrap_gen(width=1),
              strip = strip_themed(background_x = list(element_rect(fill = "#de6a01"),
                                                       element_rect(fill = "#de6a01"),
                                                       element_rect(fill = "#8e5eb8")),
                                   text_x = element_text(face = "bold", size = 12))) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = 3), order = 2),
         color = guide_legend(override.aes = list(size = 3), order = 1)) 



anovas_add <- full_gene_models %>%
  select(miRNA_names, starts_with("fdr_")) %>%
  rename_with(~str_replace(., 'X', ':')) %>%
  rename_with(~str_replace(., 'X', ':')) %>%
  rename_with(~str_replace(., 'fdr_', '')) %>%
  rename_with(~str_replace(., 'susceptibility', 'resistance')) %>%
  pivot_longer(!miRNA_names, names_to= "main_effect", values_to = "fdr_p_value") %>%
  rename(microRNA = miRNA_names)
  
anovas <- full_gene_models %>%
  select(miRNA_names, anova_table, starts_with("fdr_")) %>%
  rowwise(miRNA_names) %>%
  reframe(as_tibble(anova_table)) %>%
  rename(ndf = NumDF,
         ddf = DenDF,
         f.value = `F value`,
         p.value = `Pr(>F)`,
         microRNA = miRNA_names) %>%
  left_join(bonafides_acerv[,1:3], by="microRNA") %>%
  mutate(main_effect = rep(c("time", "exposure", "resistance", "time:exposure", "time:resistance", "exposure:resistance", "time:exposure:resistance"), 67)) %>%
  left_join(anovas_add, by=c("microRNA", "main_effect"))  %>%
  select(-microRNA, -`Sum Sq`, -`Mean Sq`) %>%
  select(miRNA_common_name, miRNA_category, main_effect, ndf, ddf, f.value, p.value, fdr_p_value)

colnames(anovas) <- c("miRNA", "Type", "Main Effect", "NumDF", "DenDF", "F-value", "p-value", "FDR p-value")

write.csv(anovas, "~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/miRNA_main_effects.csv")

#two-way anova of counts
seq <- read.csv("~/Desktop/miRNA_Characterization/miRNA_Characterization/Files/miRNA_Sequencing.csv") 

seq2 <- seq %>%
  separate(Sample,c("time", "exposure", "tank", "genotype", "susceptibility"), "_") 

counts_aov <- aov(`QC...Decontaminated` ~ exposure * susceptibility , data = seq2)
summary(counts_aov)
  
  



