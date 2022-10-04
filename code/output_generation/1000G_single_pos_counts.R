# Generate spreadsheets with single position per sheet, single document per subtype
library(tidyverse)
library(writexl)

input_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"
out_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/excel/single_pos_counts/"

read_sp <- function(pop, st, rp, data_dir){
  f_name <- paste0(data_dir, pop, "/", st, "_rp", rp, ".csv")
  df <- read_csv(f_name, show_col_types = FALSE) %>%
    select(Nuc, singletons, controls, n_gw)
  return(df)
}

read_sp_kl <- function(pop, st, rp, data_dir){
  df <- read_sp(pop, st, rp, data_dir) %>%
    mutate(p_s = singletons / sum(singletons),
           p_c = controls / sum(controls)) %>%
    mutate(kl = p_s * log(p_s / p_c))
  return(sum(df$kl))
}

read_sp_re <- function(pop, st, rp, data_dir){
  df <- read_sp(pop, st, rp, data_dir) %>%
    select(Nuc, singletons, controls) %>%
    pivot_longer(singletons:controls, names_to = "status", values_to = "n")
  mod_obj <- glm(n ~ Nuc + status, data = df, family = poisson)
  df$fitted <- predict(mod_obj, type = "response")
  N <- sum(df$n)
  df$re <- (df$n / N) * log(df$n / df$fitted)
  return(deviance(mod_obj) / (2*N))
}

read_sp_stat <- function(pop, st, data_dir, r_start = 1){
  df <- data.frame(rp = numeric(), re = numeric(), kl = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- bind_rows(df,
                    data.frame(rp = i,
                               re = read_sp_re(pop, st, i, data_dir),
                               kl = read_sp_kl(pop, st, i, data_dir)))
  }
  return(df)
}

read_all_sp <- function(pop, st, data_dir){
  df <- read_sp(pop, st, -10, data_dir)
  df$rp <- -10
  for(i in c(-9:-1, 1:10)){
    df2 <- read_sp(pop, st, i, data_dir)
    df2$rp <- i
    df <- bind_rows(df, df2)
  }
  return(df)
}

for(pop in c("ALL", "AFR", "AMR", "EUR", "EAS", "SAS")){
  print(pop)
  for(st in c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
    print(st)
    df <- read_all_sp(pop, st, input_dir)
    out_file <- paste0(out_dir, pop, "_", st, ".csv")
    write_csv(df, out_file)
  }
}

# Single excel spreadsheet for each pop-subtype pair

for(pop in c("ALL", "AFR", "AMR", "EUR", "EAS", "SAS")){
  print(pop)
  for(st in c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
    print(st)
    df <- list(`rp-10` = read_sp(pop, st, -10, input_dir),
               `rp-9` = read_sp(pop, st, -9, input_dir),
               `rp-8` = read_sp(pop, st, -8, input_dir),
               `rp-7` = read_sp(pop, st, -7, input_dir),
               `rp-6` = read_sp(pop, st, -6, input_dir),
               `rp-5` = read_sp(pop, st, -5, input_dir),
               `rp-4` = read_sp(pop, st, -4, input_dir),
               `rp-3` = read_sp(pop, st, -3, input_dir),
               `rp-2` = read_sp(pop, st, -2, input_dir),
               `rp-1` = read_sp(pop, st, -1, input_dir),
               `rp1` = read_sp(pop, st, 1, input_dir),
               `rp2` = read_sp(pop, st, 2, input_dir),
               `rp3` = read_sp(pop, st, 3, input_dir),
               `rp4` = read_sp(pop, st, 4, input_dir),
               `rp5` = read_sp(pop, st, 5, input_dir),
               `rp6` = read_sp(pop, st, 6, input_dir),
               `rp7` = read_sp(pop, st, 7, input_dir),
               `rp8` = read_sp(pop, st, 8, input_dir),
               `rp9` = read_sp(pop, st, 9, input_dir),
               `rp10` = read_sp(pop, st, 10, input_dir))
    out_file <- paste0(out_dir, pop, "_", st, ".xls")
    write_xlsx(df, path = out_file)
  }
}

# Statistics
out_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/excel/single_pos_stat/"

for(pop in c("ALL", "AFR", "AMR", "EUR", "EAS", "SAS")){
  print(pop)
  for(st in c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
    print(st)
    if(str_starts(st, "cpg")){
      r_start <- 2
    } else {
      r_start <- 1
    }
    df <- read_sp_stat(pop, st, input_dir, r_start)
    out_file <- paste0(out_dir, pop, "_", st, ".csv")
    write_csv(df, out_file)
  }
}


# Now do BRIDGES
read_sp_b <- function(st, rp){
  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/single_pos_df/"
  f_name <- paste0(data_dir, st, "_rp", rp, ".csv")
  df <- read_csv(f_name, show_col_types = FALSE) %>%
    select(Nuc, singletons, controls, n_gw)
  return(df)
}

read_all_sp_b <- function(st){
  df <- read_sp_b(st, -10)
  df$rp <- -10
  for(i in c(-9:-1, 1:10)){
    df2 <- read_sp_b(st, i)
    df2$rp <- i
    df <- bind_rows(df, df2)
  }
  return(df)
}

out_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/excel/single_pos_counts/"

for(st in c("AT_CG", "AT_GC", "AT_TA",
            "GC_AT", "GC_TA", "GC_CG",
            "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
  print(st)
  df <- read_all_sp_b(st)
  out_file <- paste0(out_dir, st, ".csv")
  write_csv(df, out_file)
}

for(st in c("AT_CG", "AT_GC", "AT_TA",
            "GC_AT", "GC_TA", "GC_CG",
            "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
  print(st)
  df <- list(`rp-10` = read_sp_b(st, -10),
             `rp-9` = read_sp_b(st, -9),
             `rp-8` = read_sp_b(st, -8),
             `rp-7` = read_sp_b(st, -7),
             `rp-6` = read_sp_b(st, -6),
             `rp-5` = read_sp_b(st, -5),
             `rp-4` = read_sp_b(st, -4),
             `rp-3` = read_sp_b(st, -3),
             `rp-2` = read_sp_b(st, -2),
             `rp-1` = read_sp_b(st, -1),
             `rp1` = read_sp_b(st, 1),
             `rp2` = read_sp_b(st, 2),
             `rp3` = read_sp_b(st, 3),
             `rp4` = read_sp_b(st, 4),
             `rp5` = read_sp_b(st, 5),
             `rp6` = read_sp_b(st, 6),
             `rp7` = read_sp_b(st, 7),
             `rp8` = read_sp_b(st, 8),
             `rp9` = read_sp_b(st, 9),
             `rp10` = read_sp_b(st, 10))
  out_file <- paste0(out_dir, st, ".xls")
  write_xlsx(df, path = out_file)
}


read_sp_kl_b <- function(st, rp){
  df <- read_sp_b(st,rp) %>%
    mutate(p_s = singletons / sum(singletons),
           p_c = controls / sum(controls)) %>%
    mutate(kl = p_s * log(p_s / p_c))
  return(sum(df$kl))
}

read_sp_re_b <- function(st, rp){
  df <- read_sp_b(st, rp) %>%
    select(Nuc, singletons, controls) %>%
    pivot_longer(singletons:controls, names_to = "status", values_to = "n")
  mod_obj <- glm(n ~ Nuc + status, data = df, family = poisson)
  df$fitted <- predict(mod_obj, type = "response")
  N <- sum(df$n)
  df$re <- (df$n / N) * log(df$n / df$fitted)
  return(deviance(mod_obj) / (2*N))
}

read_sp_stat_b <- function(st, r_start = 1){
  df <- data.frame(rp = numeric(), re = numeric(), kl = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- bind_rows(df,
                    data.frame(rp = i,
                               re = read_sp_re_b(st, i),
                               kl = read_sp_kl_b(st, i)))
  }
  return(df)
}

out_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/excel/single_pos_stat/"
for(st in c("AT_CG", "AT_GC", "AT_TA",
            "GC_AT", "GC_TA", "GC_CG",
            "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
  print(st)
  if(str_starts(st, "cpg")){
    r_start <- 2
  } else {
    r_start <- 1
  }
  df <- read_sp_stat_b(st, r_start)
  out_file <- paste0(out_dir, st, ".csv")
  write_csv(df, out_file)
}
