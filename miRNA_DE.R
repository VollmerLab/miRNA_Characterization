#miRNA differential expression with susceptibility

miRNA_counts <- read.csv("~/Desktop/RNA_Analysis_Lab/RNA_Analysis/PITA_files/mirna_counts.csv") %>% select(-X) %>%
  rename(microRNA = mirna_name) %>%
  filter(microRNA %in% bonafide_names$microRNA)

samples_code <- read.table("~/Desktop/RNA_Analysis_Lab/RNA_Analysis/PITA_files/config.txt")
colnames(samples_code) <- c("sample_name", "sample_id")

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
#normalize count table
#create DGElist
counts_miR_dge <- DGEList(miRNA_counts_mod_de)
#calculate TMM normalization factors
counts_miR_norm_fac <- calcNormFactors(counts_miR_dge, method="TMM")
#get the normalized counts
counts_miR_norm <- cpm(counts_miR_norm_fac, log=TRUE)

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
  mutate(treatment = str_c(exposure, final_disease_state, sep = '_')) %>%
  mutate(treatment_full = str_c(exposure, final_disease_state, susceptibility, sep = '_')) %>%
  mutate(across(where(is.character), as.factor))

# estimate weights using linear mixed model of dream
param <- SnowParam(parallel::detectCores() - 1, "SOCK", progressbar = TRUE)
vobjDream <- voomWithDreamWeights(counts = counts_miR_norm, 
                                  formula = ~ time * treatment_full + 
                                    (1 | genotype) +
                                    (1 | fragment_id) +
                                    (1 | tank),
                                  
                                  data = sample_info %>%
                                    arrange(sample) %>%
                                    column_to_rownames('sample'),
                                  BPPARAM = param) 

model_data <- counts_miR_norm %>% 
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
  mutate(across(c(exposure, final_disease_state, susceptibility, treatment), factor)) 

cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 'stringr', 'tidyr'))
cluster_copy(cluster, c('alpha'))
safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)


full_gene_models <- model_data %>%
  nest(data = -c(miRNA_names)) %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(model = list(lmer(value ~ time * treatment_full + 
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

library(treedataverse)
library(tidyverse)
library(magrittr)
library(lmerTest)
library(emmeans)
library(multidplyr)
library(ComplexUpset)
library(qvalue)
library(patchwork)
library(chisq.posthoc.test)
library(ggnested)
library(ggtext)

refit_models <- FALSE
alpha <- 0.05

run_posthoc <- function(data, model, a_priori, sig_values, alpha){
  names(sig_values) <- c('time', 'treatment_full', 'timeXtreatment_full')
  
  if(sig_values[['timeXtreatment_full']] < alpha | 
     (sig_values[['time']] < alpha & sig_values[['treatment_full']] < alpha)){
    
    posthoc <- emmeans(model, ~time * treatment_full)
    
  } else if(sig_values[['treatment_full']] < alpha){
    
    posthoc <- emmeans(model, ~treatment_full)
    
  } else if(sig_values[['time']] < alpha){
    
    posthoc <- emmeans(model, ~time)
    
  } else {
    posthoc <- emmeans(model, ~1)
  }
  
  contrasts <- contrast(posthoc, method = a_priori, adjust = 'none')
  
  tibble(posthoc = list(posthoc), contrasts = list(contrasts))
}

#contrast order  (time)_(susceptibility)
# c('T3_R', 'T7_R', 'T3_S', 'T7_S')

emmeans(full_gene_models$model[[1]], ~ time * treatment_full)
#contrast order (time)_(exposure)_(outcome)_susceptibility
# c('T3_D_D_S', 'T7_D_D_S', 'T3_D_H_R', 'T7_D_H_R', 'T3_H_H_R', 'T7_H_H_R', 'T3_H_H_S', 'T7_H_H_S' )

brecia_contrasts <- list(time = c(-1/4, 1/4, -1/4, 1/4, -1/4, 1/4, -1/4, 1/4),
                         exposure = c(1/4, 1/4, 1/4, 1/4, -1/4, -1/4, -1/4, -1/4),
                         outcome = c(1/2, 1/2, -1/6, -1/6, -1/6, -1/6, -1/6, -1/6),
                         susceptibility = c(-1/4, -1/4, 1/4, 1/4, 1/4, 1/4, -1/4, -1/4),
                         
                         t3_exposure = c(1/2, 0, 1/2, 0, -1/2, 0, -1/2, 0),
                         t3_outcome = c(1, 0, -1/3, 0, -1/3, 0, -1/3, 0),
                         t3_susceptibility = c(1/2, 0, -1/2, 0, -1/2, 0, 1/2, 0),
                         t7_exposure = c(0, 1/2, 0, 1/2, 0, -1/2, 0, -1/2),
                         t7_outcome = c(0, 1, 0, -1/3, 0, -1/3, 0, -1/3),
                         t7_susceptibility = c(0, 1/2, 0, -1/2, 0, -1/2, 0, 1/2))

full_gene_post <- full_gene_models %>%
  rename_with(.cols = starts_with('fdr_'), ~str_replace(., 'fdr_', 'sig_')) %>%
  filter(sig_treatment_full < alpha | sig_timeXtreatment_full < alpha) %>%
  
  
  mutate(significance_group = case_when(sig_timeXtreatment_full < alpha ~ 'interaction',
                                        sig_treatment_full < alpha & sig_time < alpha ~ 'interaction',
                                        sig_treatment_full < alpha ~ 'treatment_full',
                                        sig_time < alpha ~ 'time',
                                        TRUE ~ 'none')) %>%
  mutate(planned_contrast = case_when(sig_timeXtreatment_full < alpha ~ list(brecia_contrasts),
                                      sig_treatment_full < alpha & sig_time < alpha ~ list(brecia_contrasts),
                                      sig_treatment_full < alpha ~ list('pairwise'),
                                      sig_time < alpha ~ list('poly'),
                                      TRUE ~ list('identity'))) %>%
  
  rowwise %>%
  # partition(cluster) %>%
  mutate(run_posthoc(data = data, model = model,
                     a_priori = planned_contrast, 
                     sig_values = c_across(starts_with('sig_')), 
                     alpha = alpha)) %>%
  
  mutate(as_tibble(contrasts) %>%
           select(contrast, estimate, p.value) %>%
           rename(p = p.value,
                  logFC = estimate) %>%
           pivot_wider(names_from = contrast,
                       values_from = c('logFC', 'p'),
                       names_sep = '.')) %>%
  ungroup %>%
  mutate(across(starts_with('p.'),
                ~p.adjust(., method = 'fdr'), 
                .names = 'fdr.{.col}')) %>%
  rename_with(~str_replace(., 'fdr.p.', 'fdr.')) %>%
  mutate(across(starts_with('p.'), safe_qvalue,
                .names = 'q.{.col}')) %>%
  rename_with(~str_replace_all(., 'q.p.', 'q.'))
```
Graphing signficant miRNAs differential expression

```{r}
library(ComplexUpset)
library(ggplot2)
library(ggupset)

full_gene_post %>%
  select(miRNA_names, starts_with('fdr')) %>% 
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  
  pivot_longer(cols = -miRNA_names,
               names_to = c('term'),
               values_to = 'significance',
               names_prefix = 'fdr_') %>%
  filter(significance) %>%
  group_by(miRNA_names) %>%
  summarise(terms = list(term),
            .groups = 'drop') %>%
  
  ggplot(aes(x = terms)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)

sig <- full_gene_post %>%
  select(miRNA_names, starts_with("fdr.")) %>%
  rename(microRNA = "miRNA_names") %>%
  mutate(across(starts_with('fdr'), ~. < 0.05)) %>%
  mutate(sig_post_hocs = pmap_chr(across(where(is.logical)), ~ toString(names(c(...))[which(c(...))]))) %>%
  filter(sig_post_hocs != "") %>%
  left_join(bonafide_names, by="microRNA") 

#line graphs
line_data <- full_gene_post %>%
  select(miRNA_names, posthoc) %>%
  rowwise(miRNA_names) %>%
  reframe(as_tibble(posthoc)) %>%
  mutate(treatment_full_2 = treatment_full) %>%
  separate(col = treatment_full_2, into=c("exposure", "outcome", "susceptibility"), sep = "_") %>%
  rename(microRNA = "miRNA_names") %>%
  left_join(sig, by="microRNA") %>%
  mutate(header = paste(miRNA_common_name, sig_post_hocs, sep = ": "))

library(showtext)

line_data %>%  
  ggplot(aes(x=time, y=emmean, group = treatment_full)) +
  geom_line(aes(color = outcome, linetype = exposure), position=position_dodge(width=0.1)) +
  geom_point(aes(shape = susceptibility, color = outcome), position=position_dodge(width=0.1), size = 2)  + 
  scale_color_manual(values = c("firebrick2", "dodgerblue2"), labels = c("Disease", "Healthy")) +
  labs(x="Time", y=expression(Normalized~Log[2]~(cpm))) +
  labs(shape="Susceptibility", color="Outcome", linetype = "Exposure") +
  scale_linetype_discrete(labels = c("Disease", "Healthy")) +
  scale_shape_discrete(labels = c("Resistant", "Susceptible")) +
  scale_x_discrete(labels = c("T3" = "Day 3", "T7" = "Day 7")) +
  facet_wrap(.~ header, scale="free_y", labeller = label_wrap_gen(width=1)) 




tst <- line_data %>%
  filter(miRNA_names == "Acerv_scaffold_1_1115")


ggplot(tst, aes(x=time, y=emmean, group = treatment_full)) +
  geom_line(aes(color = final_disease_state, linetype = exposure), position=position_dodge(width=0.05)) +
  geom_point(aes(shape = susceptibility, color = final_disease_state), position=position_dodge(width=0.05))  + 
  scale_color_manual(values = c("red", "darkgreen"))
labs(x="Time", y="Emmean") #+
# guides(color=guide_legend("Treatment")) 