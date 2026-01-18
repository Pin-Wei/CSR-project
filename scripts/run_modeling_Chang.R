rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(lmerTest) 
library(stargazer)
# library(jtools)
library(car)
library(stats)
library(MuMIn)
library(openxlsx)
library(here)

## class & functions ===========================================================

setClass(
  "Config", 
  slots=c(
    use_z="logical", 
    rm_ex="logical", 
    rt_mode="character", # "individual" or "group"
    mdl_typ="character"  # "GLM", "CSR", or "gLMEM"
  ), 
  prototype = list(
    use_z=TRUE, 
    rm_ex=TRUE, 
    rt_mode="individual", 
    mdl_typ="CSR"
  )
)

get_input_path <- function(cfg) {
  stopifnot(is(cfg, "Config"))
  dat.root <- here::here("data", "psycholinguistic", "Chang_Naming")
  input_map <- list(
    "raw|individual"=file.path(dat.root, "all_subjs_raw_20 (indv; with infos).xlsx"), 
    "z|individual"  =file.path(dat.root, "all_subjs_zvars_20 (indv; with infos).xlsx"), 
    "raw|group"     =file.path(dat.root, "all_subjs_raw_16 (mean; with infos).xlsx"), 
    "z|group"       =file.path(dat.root, "all_subjs_zvars_16 (mean; with infos).xlsx")
  )
  dv_map <- list(
    "raw|individual"="RT", 
    "z|individual"  ="zRT", 
    "raw|group"     ="mean_RT", 
    "z|group"       ="z_mean_RT"
  )
  z_tag <- if (isTRUE(cfg@use_z)) "z" else "raw"
  key <- paste(z_tag, cfg@rt_mode, sep = "|")
  if ( is.null(input_map[[key]]) ) {
    stop(sprintf("invalid key: %s", key))
  }
  return(list(
    fp=input_map[[key]], 
    dv=dv_map[[key]]
  ))
}

get_output_path <- function(cfg) {
  stopifnot(is(cfg, "Config"))
  tags <- paste0(
    if (cfg@use_z) " (z-scored)" else "", 
    if (cfg@rm_ex) " (out-rm)" else ""
  )
  out.root <- here::here("output", "psycholinguistic", "Chang_Naming")
  out.name <- sprintf("[R] %s_reg_results%s.xlsx", cfg@mdl_typ, tags)
  return(file.path(out.root, cfg@rt_mode, out.name))
}

run_mdl <- function(mdl.type, dat, dv, f.cols) {
  switch(mdl.type, 
    "CSR" = {
      formula <- as.formula(paste(
        dv, "~", paste(f.cols, collapse = " + "), # linear terms
        "+", paste0("I(", f.cols, "^2)", collapse = " + "), # quadratic terms
        "+", paste(combn(f.cols, 2, function(x) paste(x, collapse = " * ")), collapse = " + ") # interaction terms
      ))
      mdl <- lm(formula, data = dat)
    }, 
    "GLM" = {
      formula <- as.formula(
        paste(dv, "~", paste(f.cols, collapse = " + "))
      )
      mdl <- lm(formula, data = dat)
    }, 
    "gLMEM" = {
      formula <- as.formula(paste(
        dv, "~", paste(f.cols, collapse = " + "), 
        "+ (1 | Char)", 
        # "+ (1 +", paste(f.cols, collapse = " + "), "|| subject_id)"
        "+ (1 | subject_id)"
      ))
      mdl <- lmerTest::lmer(
        formula = formula, 
        data = dat, 
        REML = FALSE, 
        na.action = na.exclude, 
        verbose = 1
      )
    }
  )
  return(mdl)
}

gen_res_df <- function(mdl) {
  is_lm <- inherits(mdl, "lm")
  
  if ( is_lm ) {
    coefs <- stats::coef(mdl)
    summ <- summary(mdl)
    r.df <- data.frame(
      R2 = summ$r.squared, 
      R2_adj = summ$adj.r.squared
    )
  } else {
    coefs <- lme4::fixef(mdl)
    rsq <- MuMIn::r.squaredGLMM(mdl)
    r.df <- data.frame(
      R2m = rsq[1],
      R2c = rsq[2]
    )
  }

  res.df <- cbind(
    t(as.data.frame(coefs)), 
    data.frame(
      nT = nobs(mdl), 
      LogLik = stats::logLik(mdl), 
      AIC = stats::AIC(mdl), 
      AICc = MuMIn::AICc(mdl), 
      BIC = stats::BIC(mdl)
    ), 
    r.df
  )
  return(res.df)
}

save_res_wb <- function(mdl, fp) {
  is_lm <- inherits(mdl, "lm")
  
  if ( is_lm ) {
    summ <- stargazer(mdl, align = TRUE, type = "text")
  } else {
    summ <- capture.output(summary(mdl))
  }
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Model Summary")
  openxlsx::writeData(wb, sheet = "Model Summary", 
                      x = summ)
  openxlsx::addWorksheet(wb, "VIF")
  openxlsx::writeData(wb, sheet = "VIF", 
                      x = capture.output(car::vif(mdl)))
  openxlsx::saveWorkbook(wb, fp, overwrite = TRUE)
}

## main ========================================================================

cfg.list <- list(
  new("Config", rt_mode="individual", mdl_typ="CSR"), 
  new("Config", rt_mode="individual", mdl_typ="GLM"), 
  new("Config", rt_mode="individual", mdl_typ="gLMEM"),
  new("Config", rt_mode="group", mdl_typ="CSR"), 
  new("Config", rt_mode="group", mdl_typ="GLM")
)

f.cols <- c("LogCF", "NS", "CON", "PC", "SC", "SAR", "IMG", "AoA")

for ( i in seq_along(cfg.list) ) {
  
  cfg         <- cfg.list[[i]]
  input.infos <- get_input_path(cfg)
  input.path  <- input.infos$fp
  dv          <- input.infos$dv
  output.path <- get_output_path(cfg)
  out.folder  <- dirname(output.path)
  
  if ( ! file.exists(out.folder)) {
    dir.create(out.folder, recursive=TRUE) 
  }
  
  if ( ! file.exists(output.path) ) {
    
    dat <- readxl::read_excel(input.path)
    
    if ( cfg@rm_ex == TRUE ) {
      dat <- dat %>% 
        dplyr::filter(if_all(all_of(f.cols), ~ abs(.x) <= 3))
    }
    
    if ((cfg@rt_mode == "individual") & (cfg@mdl_typ != "gLMEM")) {
      
      dat %>%
        split(., dat$subject_id) %>% 
        lapply(., function(x) {
          mdl <- run_mdl(cfg@mdl_typ, x, dv, f.cols)
          res_df <- gen_res_df(mdl)
          return(res_df)
        }) %>% 
        dplyr::bind_rows(.id = "SID") %>% 
        as.data.frame() %>% 
        writexl::write_xlsx(output.path)
      
    } else {
      
      mdl <- run_mdl(cfg@mdl_typ, dat, dv, f.cols)
      res.df <- gen_res_df(mdl)
      writexl::write_xlsx(res.df, output.path)
      
      output.path.2 <- sub(".xlsx", " v2.xlsx", output.path)
      save_res_wb(mdl, output.path.2)
    }
  }
}



    


      

