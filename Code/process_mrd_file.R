if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  mrd_file <- args[1]
  output_file <- args[2]
} else {
  mrd_file <- '../output.mrd'
  output_file <- '../mirna_counts.csv'
}


suppressMessages(suppressWarnings(library(tidyverse)))


# mrd_string <- read_lines(mrd_file)

find_start_end_lines <- function(mrd_string){
  # mrd_string <- read_lines(mrd_file)
  
  #all miRNA starts with a '>' as first character in line
  start_mirna_line <- str_which(mrd_string, '^>')
  
  #ending then 3 blank lines then start next miRNA
  #shift to left to align ending with correct miRNA
  #add in final line for last miRNA (subtract 3 for the extra lines)
  end_mirna_line <- c((start_mirna_line - 3 - 1)[-1], length(mrd_string) - 3)
  
  mirna_name <- str_remove(mrd_string[start_mirna_line], '^>')
  tibble(mirna_name, start_mirna_line, end_mirna_line)
}

# find_start_end_lines(mrd_string) %>%
#   slice(909)

# get_overall_results <- function(mrd_string, start, end){
#   mirDeep2 already outputs this info as a tsv file - uneeded
#   full_mirna <- mrd_string[start:end]
#   last_overall_line <- str_which(full_mirna, 'star read count')
#   
#   overall_info <- full_mirna[1:last_overall_line]
# }

get_mirna_score <- function(mrd_string, start, end){
  full_mirna <- mrd_string[start:end]
  
  str_subset(full_mirna, 'score total') %>%
    str_extract('[\\-0-9\\.]+') %>%
    parse_double()
}

get_total_mature <- function(mrd_string, start, end){
  full_mirna <- mrd_string[start:end]
  
  str_subset(full_mirna, 'mature read count') %>%
    str_extract('[\\-0-9\\.]+') %>%
    parse_double()
}

# get_mirna_score(mrd_string, start = start_mirna_line, end = end_mirna_line)

get_mirna_structure <- function(mrd_string, start, end){
  # full_mirna <- read_lines(mrd_file, skip = (start - 1), n_max = (end))
  full_mirna <- mrd_string[start:end]
  start_structure <- str_which(full_mirna, '^exp')
  end_structure <- str_which(full_mirna, '^pri_struct')
  
  full_mirna[start_structure:end_structure]
}

# mirna_structure <- get_mirna_structure(mrd_string, 14623, 22553) 

find_mature_sequence <- function(mirna_structure){
  mirna_structure <- str_subset(mirna_structure, '^exp') %>%
    str_remove('exp[ ]*')
  
  out <- str_locate(mirna_structure, 'M+') %>%
    as.integer()
  names(out) <- c('mature_start', 'mature_end')
  out
}

# mature_location <- find_mature_sequence(mirna_structure)# %>% bind_rows()

# start <- 509692; end <- 509740

get_sample_results <- function(mrd_string, start, end){
  # full_mirna <- read_lines(mrd_file, skip = (start - 1), n_max = (end))
  full_mirna <- mrd_string[start:end]
  start_samples <- str_which(full_mirna, '^pri_struct') + 1
  
  full_mirna[start_samples:length(full_mirna)] %>%
    as_tibble() %>%
    mutate(value = str_replace_all(value, '\\.', 'Z')) %>%
    separate(value, into = c('sample_id', 'mirdeep_code', 'count', 'sequence', 'mismatch')) %>%
    mutate(sequence = str_replace_all(sequence, 'Z', '.'),
           count = str_remove(count, '^x') %>% as.integer,
           mismatch = as.integer(mismatch))
}

# tmp <- get_sample_results(mrd_string, 512023, 512119)

parse_mature_sequence <- function(sequence, mature_start, mature_end){
  str_sub(sequence, start = mature_start, end = mature_end)
}

# parse_mature_sequence(tmp$sequence, mature_location[1], mature_location[2])

parse_mrd <- function(mrd_file, min_score = 0, min_mature_length = 18, max_mature_length = 25){
  mrd_string <- read_lines(mrd_file)
  
  find_start_end_lines(mrd_string) %>%
    # sample_n(5) %>%
    # slice(909) %>%
    rowwise(mirna_name) %>%
    mutate(score = get_mirna_score(mrd_string, start = start_mirna_line, end = end_mirna_line)) %>%
    filter(score >= min_score) %>%
    mutate(get_mirna_structure(mrd_string, start = start_mirna_line, end = end_mirna_line) %>%
             find_mature_sequence() %>%
             bind_rows(),
           mature_length = length(mature_start:mature_end)) %>%
    filter(mature_length >= min_mature_length) %>%
    filter(mature_length <= max_mature_length) %>%
    mutate(sample_results = list(get_sample_results(mrd_string, start = start_mirna_line, end = end_mirna_line))) %>%
    reframe(mutate(sample_results, sequence = parse_mature_sequence(sequence, mature_start, mature_end)))
}


parse_mrd(mrd_file) %>%
  filter(str_detect(sequence, '[a-zA-Z]')) %>%
  group_by(mirna_name, sample_id) %>%
  summarise(n = sum(count),
            .groups = 'drop') %>%
  write_csv(output_file)


