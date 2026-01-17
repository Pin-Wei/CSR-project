rm(list = ls())

library(dplyr)
library(tidyr)
library(pracma)
library(here)

path.root <- here::here("data", "SRT", "BinShan")
data.path.list <- list.files(
  path = file.path(path.root, "Exp1_Spatial", "Data&Analysis"), 
  pattern = ".csv$", 
  recursive = TRUE, 
  full.names = TRUE
)
data.path <- data.path.list[2]
data.name <- basename(data.path)
sid <- sub(".*_([0-9]{2})\\.csv$", "\\1", data.name)
suffix <- sub("(.*)_[A-Z]{2}_[0-9]{2}\\.csv$", "\\1", data.name)
out.name <- paste0(sid, "_", suffix, ".csv")
out.path <- file.path(path.root, "Exp1_Spatial", "PreprocessedData", out.name)
  
data <- read.csv(data.path)
data.wide <- data %>% 
  dplyr::filter(
    State == "Sequential"
  ) %>% 
  dplyr::mutate(
    Item = rep(1:12, length.out = length(RT)), 
    Trial = ceiling(seq_along(RT) / 12)
  ) %>%
  dplyr::select(all_of(c(
    "Block", "Item", "RT", "Trial"
  ))) %>%
  tidyr::pivot_wider(
    names_from = Item,
    values_from = RT,
    names_prefix = "Item_"
  ) %>%
  arrange(Trial)

RTs.prev  <- unname(as.matrix(data.wide[1:nrow(data.wide)-1, 5:14]))
RTs.curr  <- unname(as.matrix(data.wide[2:nrow(data.wide)  , 5:14]))
RTs.1back <- unname(as.matrix(data.wide[2:nrow(data.wide)  , 4:13]))
RTs.2back <- unname(as.matrix(data.wide[2:nrow(data.wide)  , 3:12]))
k <- (nrow(data.wide) - 1) * 10

data.out <- data.frame(
  "Block" = c(
    rep(1, each=10*7), rep(2:13, each=10*8)),
  
  "n+1" = pracma::Reshape(
    RTs.curr - RTs.1back, k, 1), 
    # RT reduction from the previous item
  
  "n+2" = pracma::Reshape(
    RTs.curr - RTs.2back, k, 1), 
    # RT reduction from the item two trials ago
  
  "seqlearn" = pracma::Reshape(
    RTs.curr - RTs.prev, k, 1)
    # RT reduction from the last repetition
)
write.csv(data.out, row.names = FALSE, file = out.path)