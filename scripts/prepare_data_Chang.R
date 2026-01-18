rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(writexl)
library(here)

data.folder <- here::here("data", "psycholinguistic", "Chang_Naming")

f.cols <- c(
  # "Stop", "Aspirated", "Voiced", "Bilabial", "Alveolar", "Palatoalveolar", "Alveolopalatal", 
  "LogCF", # log transformed character frequency
  "NS",    # number of strokes
  # "REG",   # regular (1) or irregular (0)
  # "UNP",   # unpronounceable (1) or pronounceable (0)
  "CON",   # phonological consistency
  "PC",    # phonetic combinability, number of characters that can be created by a phonetic radical
  "SC",    # semantic combinability, number of characters that can be created by a semantic radical
  "SAR",   # semantic ambiguity rating, measuring the number of meanings of a character
  "IMG",   # imageability, measuring how easily a mental image could be aroused by a character
  "AoA"
)

## Chang et al. (2020), the trial-wise dataset =================================

df.20.v0 <- readRDS(file.path(data.folder, "AoA_character_naming.rds")) 

df.20.v1 <- df.20.v0 %>%
  dplyr::rename(
    Char = character, 
    RT = rt, 
    zRT = z_rt
  ) 

for ( z in c(1, 2) ) {
  out.name.1 <- paste0("all_subjs_", c("raw", "zvars")[z], "_20 (indv; with infos).xlsx")
  out.name.2 <- paste0("all_subjs_", c("raw", "zvars")[z], "_20 (indv).xlsx")
  
  rt.col <- c("RT", "zRT")[z]
  
  df.20.v2 <- df.20.v1 %>% 
    select(all_of(c("Char", "subject_id", rt.col, f.cols)))
  
  if ( z == 2 ) {
    df.20.v2 <- df.20.v2 %>% 
      # group_by(subject_id) %>% 
      dplyr::mutate(across(
        .cols = all_of(f.cols),
        .fns = ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)
      )) # %>% 
    # ungroup() 
  }
  
  writexl::write_xlsx(
    df.20.v2, file.path(data.folder, out.name.1)
  )
  
  df.20.v2 %>% 
    select(all_of(c(f.cols, rt.col))) %>% 
    writexl::write_xlsx(file.path(data.folder, out.name.2))
  
  for ( sid in df.20.v2$subject_id ) {
    out.name.3 <- paste0(c("", "zscored_")[z], "sub_", sid, ".xlsx")
    out.path.3 <- file.path(data.folder, out.name.3)
    
    if ( ! file.exists(out.path.3) ) {
      df.20.v2 %>% 
        subset(subject_id == sid) %>% 
        dplyr::select(all_of(c(f.cols, rt.col))) %>% 
        writexl::write_xlsx(out.path.3, format_headers = TRUE)
    }
  }
}

## Chang & Lee (2016), the group-mean dataset ==================================

df.16.v0 <- read.csv(file.path(data.folder, "Chang_Lee_2016.csv"))

df.16.v1 <- df.16.v0 %>% 
  dplyr::mutate(
    LogCF = log(Frequency)
  ) %>% 
  dplyr::rename(
    Char = Character, 
    # NS = Stroke,
    # CON = `Consistency..token.`, 
    # PC = `Phonetic.Combinability`, 
    # SC = `Semantic.Combinability`, 
    # SAR = Semantic.Ambiguity.Rating, 
    # mean_ACC = `Naming.Acc`, 
    mean_RT = `Naming.RT`
  ) %>% 
  select(all_of(c(
    "Char", "mean_RT" # "LogCF", "NS", "CON", "PC", "SC", "SAR"
  )))

for ( z in c(1, 2) ) {
  out.name.4 <- paste0("all_subjs_", c("raw", "zvars")[z], "_16 (mean; with infos).xlsx")
  out.name.5 <- paste0("all_subjs_", c("raw", "zvars")[z], "_16 (mean).xlsx")
  
  if ( z == 2 ) {
    df.16.v1 <- df.16.v1 %>% 
      dplyr::mutate(across(
        .cols = all_of("mean_RT"),
        .fns = ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)
      )) %>% 
      dplyr::rename(
        z_mean_RT = mean_RT
      )
  }
  
  df.20.fp <- paste0("all_subjs_", c("raw", "zvars")[z], "_20 (indv; with infos).xlsx")
  df.20.reload <- readxl::read_excel(file.path(data.folder, df.20.fp)) %>% 
    select(-all_of(c("subject_id", c("RT", "zRT")[z])))
  
  df.merged <- list(df.16.v1, df.20.reload) %>% 
    purrr::reduce(dplyr::left_join, by = "Char") %>% 
    na.omit() %>% 
    unique()
  
  writexl::write_xlsx(
    df.merged, file.path(data.folder, out.name.4)
  )
  
  df.merged %>% 
    select(all_of(c(f.cols, c("mean_RT", "z_mean_RT")[z]))) %>% 
    writexl::write_xlsx(file.path(data.folder, out.name.5))
}




