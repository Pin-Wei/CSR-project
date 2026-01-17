rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(trimr)
library(here)

path.root <- here::here("data", "SRT", "Joanne")

learned.subjs <- c(
  2, 3, 4, 5, 7, 16, 17, 21, 22, 32, 33, 36, 
  9, 10, 11, 12, 14, 23, 25, 26, 29, 30, 39, 40
)

for ( inst.type in c("explicit", "implicit") ) {
  
  data.paths <- list.files(
    path = file.path(path.root, "RawData", inst.type), 
    pattern = ".csv$", 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  for ( fp in data.paths ) {
    
    DF <- read.csv(fp) %>% 
      dplyr::filter( block %in% c(
        "R", paste0("Block", seq(1:14))
      )) %>% 
      dplyr::rename(c(
        "SID"      = "participant", 
        "Block"    = "block", 
        "State"    = "rule", 
        "Target"   = "Stimuli", 
        "Response" = "Stimuli_key_resp.keys", 
        "RT"       = "Stimuli_key_resp.rt"
      )) %>% 
      dplyr::mutate(
        Index = rownames(.), 
        Target = purrr::map_chr(
          Target, ~ gsub(".jpg", "", .x)), 
        Block = purrr::map_chr(
          Block, ~ gsub("Block", "", .x))
      ) %>% 
      dplyr::mutate(
        Hit = dplyr::case_when(
          Target == '1' & Response == 5 ~ 1, 
          Target == '2' & Response == 6 ~ 1, 
          Target == '3' & Response == 7 ~ 1, 
          Target == '4' & Response == 8 ~ 1, 
          .default = 0
        )) %>% 
      dplyr::mutate(across(
        .cols = c(
          "Block", "State", "Target", "Response", "Hit"), 
        .fns = as.factor
      )) %>% 
      dplyr::select(all_of(c(
        "Index", "SID", "Block", "State", 
        "Target", "Response", "RT", "Hit"
      ))) 
    
    # ## convert unit from second to millisecond
    # DF$RT <- DF$RT * 1000
    
    sid <- unique(DF$SID)[[1]]
    
    if ( sid %in% learned.subjs ) {
      subj.group <- "learned"
    } else {
      subj.group <- "not_learned"
    }
    
    out.dir <- file.path(
      path.root, "IntermediateData", inst.type, subj.group
    )
    
    if ( ! file.exists(out.dir) ) { 
      dir.create(out.dir, recursive = TRUE) 
    }
    
    write.csv(
      DF, row.names = FALSE, 
      file = file.path(out.dir, sprintf('sub_%02d.csv', sid))
    )
    
    ## Response Time Trimming (trimr):
    ## https://cran.r-project.org/web/packages/trimr/vignettes/overview.html
    
    DF.trimmed <- trimr::modifiedRecursive(
      data = DF,
      minRT = 0, 
      pptVar = "SID", 
      condVar = "Block",
      accVar = "Hit",
      rtVar = "RT",
      returnType = "raw"
    )
    
    write.csv(
      DF.trimmed, row.names = FALSE, 
      file = file.path(out.dir, sprintf('sub_%02d_cleaned.csv', sid))
    )
  }
}

## discarded -------------------------------------------------------------------
# ## modified from trimr::modifiedRecursive()
# 
# for ( cond in unique(DF[["Block"]]) ) {
#   
#   sub.DF <- subset(DF, Block == cond)
#   
#   repeat {
#     
#     sample.size <- nrow(sub.DF)
#     
#     if( sample.size <= 2 ){ break }
#     if( sample.size > 100 ){ sample.size <- 100 }
#     
#     # number of trials have been removed on each recursion 
#     removed.nT <- 0
#     
#     # 1-1. look up the SD criterion to use for the current sample size
#     sd.crit <- trimr::linearInterpolation$modifiedRecursive[
#       which(trimr::linearInterpolation$sampleSize == sample.size)
#     ]
#     
#     # 1-2. temporarily remove the slowest RT value from the distribution and than calculate mean ans SD
#     RT.max.idx <- which.max(sub.DF[["RT"]])
#     temp.DF <- sub.DF %>% 
#       filter(row_number() != RT.max.idx)
#     temp.mean <- mean(temp.DF[["RT"]])
#     temp.sd <- sd(temp.DF[["RT"]])
#     
#     # 1-3. determine cut-off values
#     RT.max.crit <- temp.mean + (sd.crit * temp.sd)
#     RT.min.crit <- temp.mean - (sd.crit * temp.sd)
#     rm(list = c("temp.DF", "temp.mean", "temp.sd"))
#     
#     # 2-1. remove data points that are above (slower than) the cut-off value
#     RT.max <- sub.DF[RT.max.idx, ][["RT"]]
#     if( RT.max > RT.max.crit ) {
#       sub.DF <- sub.DF %>% 
#         filter(row_number() != RT.max.idx)
#       removed.nT <- 1
#     }
#     
#     # 2-1. remove data points that are below (faster than) the cut-off value
#     RT.min.idx <- which.min(sub.DF[["RT"]])
#     RT.min <- sub.DF[RT.min.idx, ][["RT"]]
#     if( RT.min < RT.min.crit ) {
#       sub.DF <- sub.DF %>% 
#         filter(row_number() != RT.min.idx)
#       removed.nT <- 1
#     }
#     
#     # repeated until no outliers remain
#     if( removed.nT == 0 ){ break }
#     # or until the sample size drops below four
#     if( nrow(sub.DF) < 5 ){ break }
#   }
# }
