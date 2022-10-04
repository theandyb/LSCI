#' Load single position results for a single flanking position
#'
#' @param population 1KGP super population code (AMR, AFR, EUR, EAS, SAS, ALL)
#' @param subtype Mutation subtype
#' @param rp Relative position
#' @param data_dir Location of count files
#' @return Data frame with nucleotide-level statistics ("residuals")
load_sp_results <- function(population, subtype, rp,
                            data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  if(population == "BRIDGES"){
    f_name <- paste0(data_dir, subtype, "_rp", rp, ".csv")
  } else{
    f_name <- paste0(data_dir, population, "/", subtype, "_rp", rp, ".csv")
  }
  df <- read_csv(f_name, col_types = cols())
  df_tall <- df %>%
    select(Nuc, singletons, controls) %>%
    pivot_longer(-Nuc) %>%
    rename(status = name, n = value)

  mod_obj <- glm(n ~ Nuc + status, data = df_tall, family = poisson())
  n_s <- df_tall %>% filter(status == "singletons") %>% pull(n) %>% sum()
  n_c <- df_tall %>% filter(status == "controls") %>% pull(n) %>% sum()
  n_t <- n_s + n_c
  df_tall$resid <- resid(mod_obj)
  df_tall$fitted <- mod_obj$fitted.values
  df_tall$re <- 2 * df_tall$n * log(df_tall$n / mod_obj$fitted.values) / (2 * sum(df_tall$n))
  #df_tall$pred <- predict(mod_obj, type = "response")
  new_dat <- data.frame(Nuc = c("A", "C", "G", "T"), status = "singletons")
  new_dat$pred <- predict(mod_obj, type = "response", newdata = new_dat)
  new_dat$deviance <- {df_tall %>%
      select(Nuc, resid) %>%
      group_by(Nuc) %>%
      summarize(resid = sum(resid^2)) %>%
      pull(resid)}
  df <- df%>%
    inner_join({new_dat %>% select(-status)}, by = "Nuc") %>%
    mutate(pct_pred = pred / sum(pred),
           pct_s = singletons / sum(singletons)) %>%
    mutate(kl_r = pct_s * log(pct_s / pct_pred)) %>%
    select(-starts_with("pct"))
  df <- df %>%
    inner_join({df_tall %>%
        select(Nuc, re) %>%
        group_by(Nuc) %>%
        summarize(re = sum(re))},
        by = "Nuc") %>%
    mutate(s_re = sign(kl_r) * re)
  return(df)
}

#' Load all single position results for a subtype
load_all_sp_results <- function(population, subtype, r_start = 1,
                                data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  df <- load_sp_results(population, subtype, -10, data_dir)
  df$rp <- -10
  for(i in c(-9:-1, r_start:10)){
    df2 <- load_sp_results(population, subtype, i, data_dir)
    df2$rp <- i
    df <- bind_rows(df, df2)
  }
  return(df)
}

load_all_results_all_sp <- function(subtype, r_start = 1,
                                    data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- load_all_sp_results("AFR", subtype, r_start = r_start, data_dir = data_dir)
  final$pop <- "AFR"

  for(pop in c("AMR", "EAS" ,"EUR", "SAS")){
    df <- load_all_sp_results(pop, subtype, r_start = r_start, data_dir = data_dir)
    df$pop <- pop
    final <- bind_rows(final, df)
  }

  return(final)
}

#' Relative entropy by position
#' This is the version not based on the loglinear model
#'
kl_by_position <- function(population ,subtype, r_start = 1,
                           data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- load_sp_results(population, subtype, i, data_dir)
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = sum(df$kl_r)))
  }
  return(final)
}

#' Zhu et al re by position
#'
zhu_by_pos <- function(population, subtype, r_start = 1,
                       data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- load_sp_results(population, subtype, i, data_dir)
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = sum(df$re, na.rm = T)))
  }
  return(final)
}

#' Relative Entropy from loglinear model
#'
dev_re <- function(population, subtype, rp, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  df <- load_sp_results(population, subtype, rp, data_dir) %>%
    select(Nuc, singletons, controls) %>%
    pivot_longer(singletons:controls, names_to = "status", values_to = "n")
  mod_obj <- glm(n ~ status + Nuc, family = "poisson", data = df)
  df <- mod_obj$data
  df$res <- residuals(mod_obj)^2
  n_s <- df %>%
    filter(status == "singletons") %>%
    pull(n) %>%
    sum()
  df <- df %>%
    mutate(re = ifelse(status == "singletons",  res / (2*n_s), res / (10*n_s)))
  return(sum(df$re))
}

#' Loglinear model relative entropy across positions
#'
dev_re_by_position <- function(population, subtype, r_start = 1,
                               data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = dev_re(population, subtype, i, data_dir)))
  }
  return(final)
}

#' Goodness of fit chi sq by position - control model
gof_by_position <- function(population, subtype, r_start = 1,
                            data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- load_sp_results(population, subtype, i, data_dir)
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = sum(df$chi_sq_ct, na.rm = T)))
  }
  return(final)
}

#' Relative entropy by position - control model
re_by_position <- function(population, subtype, r_start = 1,
                           data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- load_sp_results(population, subtype, i, data_dir)
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = sum(df$kl_r, na.rm = T)))
  }
  return(final)
}

# Deviance by position
dev_by_position <- function(population, subtype, r_start = 1,
                            data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  final <- data.frame(pos = numeric(), statistic = numeric())
  for(i in c(-10:-1, r_start:10)){
    df <- load_sp_results(population, subtype, i, data_dir)
    final <- bind_rows(final,
                       data.frame(pos = i,
                                  statistic = sum(df$deviance, na.rm = T)))
  }
  return(final)
}

#' Plot position x statistic
plot_pos_stat <- function(population, subtype, stat_func,
                          title_text, ylab_text,
                          r_start = 1,
                          data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/" ){
  df <- stat_func(population, subtype, r_start, data_dir)
  if(!str_starts(subtype, "cpg")){
    subtype2 <- str_replace(subtype, "_", " → ")
  } else {
    subtype <- str_sub(subtype, 5)
    subtype2 <- paste0("(cpg) ", str_replace(subtype, "_", " → "))
  }
  p <-df %>%
    ggplot(aes(x = pos, y = statistic)) +
    geom_point() +
    geom_line() +
    ggtitle(paste0(title_text, subtype2),
            paste0("Population: ", population)) +
    xlab("Relative Position") +
    ylab(ylab_text)
  return(p)
}

plot_nuc_level_bar <- function(population, subtype, stat_col,
                               title_text, r_start=1,
                               data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  df <- load_all_sp_results(population, subtype, r_start, data_dir) %>%
    replace_na(list(stat_col = 0)) %>%
    group_by(rp) %>%
    mutate(pct = !!rlang::sym(stat_col) / sum(!!rlang::sym(stat_col)))

  if(!str_starts(subtype, "cpg")){
    subtype2 <- str_replace(subtype, "_", " → ")
  } else {
    subtype <- str_sub(subtype, 5)
    subtype2 <- paste0("(cpg) ", str_replace(subtype, "_", " → "))
  }

  p <- df %>%
    ggplot(aes(x = rp, y = pct, fill = Nuc)) +
    geom_col() +
    xlab("Relative Position") +
    ylab("Proportion Contributed") +
    labs(fill = "Nucleotide")
 if(population != ""){
   sub_text = paste0("Subtype: ", subtype2, "; Population: ", population)
   p <- p + ggtitle(title_text, sub_text)
 }
  else {
    p <- p + ggtitle(title_text)
  }
  return(p)
}

plot_signed_nuc_by_pos <- function(population, subtype, r_start = 1,
                                   data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  df <- load_all_sp_results(population, subtype, r_start, data_dir = data_dir)
  df <- df %>%
    replace_na(list("chi_sq_ct" = 0)) %>%
    mutate(s_val = sign(singletons - exp_ct)) %>%
    mutate(signed_chi_sq = s_val * chi_sq_ct)
  g_title <- paste0("Nucleotide Level Signed ChiSq Residual: ", str_replace_all(subtype, "_", "-"))

  p <- df %>%
    ggplot(aes(x = rp, y = signed_chi_sq, colour = Nuc)) +
    geom_point() +
    geom_line() +
    xlab("Relative Position") +
    ylab("Signed ChiSq Contribution") +
    labs(fill = "Nucleotide") +
    ggtitle(g_title, paste0("Population: ", population))
  return(p)
}

plot_signed_re_by_pos <- function(population, subtype, r_start = 1,
                                   data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"){
  df <- load_all_sp_results(population, subtype, r_start, data_dir = data_dir)

  if(!str_starts(subtype, "cpg") & !str_starts(subtype, "all")){
    subtype2 <- str_replace(subtype, "_", " → ")
  } else {
    subtype <- str_sub(subtype, 5)
    if(str_starts(subtype, "cpg")){
      subtype2 <- paste0("(cpg) ", str_replace(subtype, "_", " → "))
    } else{
      subtype2 <- paste0(str_replace(subtype, "_", " → "))
    }
  }

  g_title <- paste0("Nucleotide Level Relative Entropy: ", subtype2)

  p <- df %>%
    ggplot(aes(x = rp, y = s_re, colour = Nuc)) +
    geom_point() +
    geom_line() +
    xlab("Relative Position") +
    ylab("Relative Entropy") +
    labs(fill = "Nucleotide") +
    ggtitle(g_title, paste0("Population: ", population))
  return(p)
}
