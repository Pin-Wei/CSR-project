rm(list = ls())

require(dplyr)
require(pracma)
# require(tidyr)

setwd('C:/Users/PinWei/my_CSR/Data_SRT_BinShan')

file_paths <- list.files(
  path = file.path("Exp1_Spatial", "Data&Analysis"), 
  pattern = ".csv$", 
  recursive = TRUE, 
  full.names = TRUE
)

data <- read.csv(file_paths[2])

data_wide <- data %>% 
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
  pivot_wider(
    names_from = Item,
    values_from = RT,
    names_prefix = "Item_"
  ) %>%
  arrange(Trial)

RTs.prev  <- unname(as.matrix(data_wide[1:nrow(data_wide)-1, 5:14]))
RTs.curr  <- unname(as.matrix(data_wide[2:nrow(data_wide)  , 5:14]))
RTs.1back <- unname(as.matrix(data_wide[2:nrow(data_wide)  , 4:13]))
RTs.2back <- unname(as.matrix(data_wide[2:nrow(data_wide)  , 3:12]))
k <- (nrow(data_wide) - 1) * 10

data_out <- data.frame(
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