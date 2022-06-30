#' Get contingency table for two positions
#'
#' @param population 1KGP superpopulation code
#' @param subtype Basic mutation subtype
#' @param p1 First relative position
#' @param p2 Second relative position
#' @param data_dir Location of count tables
#' @return Count table
singleton_count_2pos <- function(population, subtype, p1, p2,
                                 data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0)
  df_tab <- xtabs(singletons ~ p1 + p2, data = df)
  return(df_tab)
}

#' Perform statistical test for singleton counts stratified by nucleotides at two flanking positions
#' @param population 1KGP superpopulation code
#' @param subtype Simple mutation subtype
#' @param stat_func Function to get statistic for pair of positions
#' @return Data.frame with statistics at each pair of positions
stat_pair_all <- function(population, subtype, stat_func, r_start=1,
                          data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/", ...){
  final <- data.frame(p1 = numeric(),
                      p2 = numeric(),
                      statistic = numeric())
  for(i in c(-10:-1,r_start:9)){
    for(j in (i+1):10){
      if(j == 0) next
      if(j == 1 & r_start > 1) next
      stat_val <- stat_func(population, subtype, i, j, data_dir, ...)
      final <- bind_rows(final, data.frame(p1 = i,
                                           p2 = j,
                                           statistic = stat_val))
    }
  }
  return(final)
}

#' Heatmap plot for two-position statistics
#'
#' @param population 1KGP super-population code
#' @param subtype Simple mutation subtype
#' @param statistic_function Function which returns data frame with statistic at pairs of positions
#' @param plot_title Text to begin plot title
#' @param fill_label Legend title
#' @return Figure with statistic values at each pair of positions
heatmap_plot <- function(population, subtype,
                         statistic_function, plot_title, fill_label, r_start = 1,
                         data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/", ...){
  df <- stat_pair_all(population, subtype, statistic_function, r_start, data_dir, ...)

  if(!str_starts(subtype, "cpg")){
    subtype2 <- str_replace(subtype, "_", " → ")
  } else {
    subtype <- str_sub(subtype, 5)
    subtype2 <- paste0("(cpg) ", str_replace(subtype, "_", " → "))
  }

  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    ggtitle(paste0(plot_title, subtype2),
            paste0("Population: ", population))+
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill = fill_label) +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}

#' Perform chi-sq test of independence for singleton counts stratified by nucleotides at two flanking positions
#' @param population 1KGP superpopulation code
#' @param subtype Simple mutation subtype
#' @return Chi square statistic
chi_sq_ind <- function(population, subtype, p1, p2, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  df_tab <- singleton_count_2pos(population, subtype, p1, p2, data_dir = data_dir)
  stat_val <- (chisq.test(df_tab)$statistic %>% unname())
  return(stat_val)
}

chi_sq_ind_plot <- function(population, subtype, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, chi_sq_ind, "Chi Sq Test of Independence: ", "Chi Sq Statistic", data_dir = data_dir)
  return(p)
}

chi_sq_gof_pair <- function(population, subtype, p1, p2, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/") {
  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0)
  df$p_c <- df$controls / sum(df$controls)
  df$e_s <- sum(df$singletons) * df$p_c
  df$chi_res <- (df$e_s - df$singletons)^2 / df$e_s
  df$p_s <- df$singletons / sum(df$singletons)
  df$kl_res <- df$p_s * log(df$p_s / df$p_c)
  return(sum(df$chi_res))
}

chi_sq_re_pair <- function(population, subtype, p1, p2,  data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    mutate(p_c = controls / sum(controls),
           p_s = singletons / sum(singletons)) %>%
    mutate(re = p_s * log(p_s / p_c))
  return(sum(df$re))
}

chi_sq_gof_plot <- function(population, subtype, data_dir ="/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, chi_sq_gof_pair,
                    "Control Di-Nucleotide Rate: ", "Chi Sq Statistic", data_dir = data_dir)
  return(p)
}

chi_sq_re_plot <- function(population, subtype, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, chi_sq_re_pair,
                    "Control Di-Nucleotide Rate: ", "Relative Entropy", data_dir = data_dir)
  return(p)
}

get_mar <- function(freq_tab){
  n_t <- sum(freq_tab)
  col_mar <- colSums(freq_tab) / n_t
  row_mar <- rowSums(freq_tab) / n_t

  prod_mat <- row_mar %*% t(col_mar)
  return(prod_mat)
}

get_int <- function(freq_tab){
  n_t <- sum(freq_tab)
  col_mar <- colSums(freq_tab) / n_t
  row_mar <- rowSums(freq_tab) / n_t

  prod_mat <- row_mar %*% t(col_mar)
  int_mat <- (freq_tab / n_t) / prod_mat
  return(int_mat)
}

control_int_mult_pair <- function(population, subtype, p1, p2,
                                  data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/",
                                  type = "gof") {
  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0)
  control_tab <- df %>%
    select(p1, p2, controls) %>%
    pivot_wider(names_from = p2, values_from = controls) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'p1')
  singleton_tab <- df %>%
    select(p1, p2, singletons) %>%
    pivot_wider(names_from = p2, values_from = singletons) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'p1')

  exp_s <- get_mar(singleton_tab) * get_int(control_tab)
  exp_p_s <- exp_s / sum(exp_s) # normalize
  exp_s <- exp_p_s * sum(singleton_tab)
  chi_sq <- (exp_s - singleton_tab)^2 / exp_s

  if(type == "gof"){
    stat_val <- sum(chi_sq)
  } else{
    p_s <- singleton_tab / sum(singleton_tab)
    re <- p_s * log(p_s / exp_p_s)
    stat_val <- sum(re)
  }

  return(stat_val)
}

control_mult_int_plot <- function(population, subtype, r_start = 1,
                                  data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/",  type = "gof"){
  if(type == "gof") {
    l_title <- "Chi Sq. Statistic"
  } else {
    l_title <- "Relative Entropy"
  }
  p <- heatmap_plot(population, subtype, control_int_mult_pair,
                    "Control Estimate Interaction: ",
                    l_title, r_start = r_start, data_dir = data_dir , type)
  return(p)
}

deviance_pair <- function(population, subtype, p1, p2, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())

  return(mod_obj %>% deviance)
}

deviance_pair_re <- function(population, subtype, p1, p2, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  n_singletons <- df %>% filter(status=="singletons") %>% pull(n) %>% sum()
  n_controls <- df %>% filter(status=="controls") %>% pull(n) %>% sum()
  df <- df %>%
    rowwise() %>%
    mutate(re.res = ifelse(status == "singletons", (res / (2*n_singletons)), (res / (2*n_controls)) ))
  re <- sum(df$re.res)
  return(re)
}

zhu_pair_re <- function(population, subtype, p1, p2, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  n_obs <- sum(df$n)
  re <- sum(df$res)
  re <- re / (2*n_obs)
  return(re)
}

deviance_plot <- function(population, subtype, r_start = 1, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, deviance_pair,
                    "Loglinear Model Deviance: ", "Deviance", r_start = r_start, data_dir = data_dir)
  return(p)
}

deviance_re_plot <- function(population, subtype, r_start = 1, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, deviance_pair_re,
                    "Loglinear Model Relative Entropy: ", "Relative Entropy", r_start = r_start, data_dir = data_dir)
  return(p)
}

zhu_re_plot <- function(population, subtype, r_start = 1, data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  p <- heatmap_plot(population, subtype, zhu_pair_re,
                    "Loglinear Model Relative Entropy: ", "Relative Entropy", r_start = r_start, data_dir = data_dir)
  return(p)
}

get_residuals_re <- function(subtype, population, p1, p2,
                             data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  df$fit <- predict(mod_obj, type = "response")
  n_singletons <- df %>% filter(status=="singletons") %>% pull(n) %>% sum()
  n_controls <- df %>% filter(status=="controls") %>% pull(n) %>% sum()
  df <- df %>%
    rowwise() %>%
    mutate(re.res = ifelse(status == "singletons", (res / (2*n_singletons)), (res / (2*n_controls)) ))
  df <- df %>%
    filter(status == "singletons") %>%
    mutate(dir = sign(n - fit)) %>%
    mutate(s.res = re.res * dir) %>%
    ungroup() %>%
    mutate(prop.res = sqrt(res) / sum(sqrt(res))) %>%
    arrange(desc(res))
  return(df)
}

get_residuals_re2 <- function(subtype, population, p1, p2,
                              data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  df$fit <- predict(mod_obj, type = "response")
  n_singletons <- df %>% filter(status=="singletons") %>% pull(n) %>% sum()
  n_controls <- df %>% filter(status=="controls") %>% pull(n) %>% sum()
  df <- df %>%
    rowwise() %>%
    mutate(re.res = ifelse(status == "singletons", (res / (2*n_singletons)), (res / (2*n_controls)) ))
  df <- df %>%
    select(p1, p2, re.res) %>%
    group_by(p1, p2) %>%
    summarize(re.res = sum(re.res)) %>%
    ungroup() %>%
    mutate(prop.re = re.res / sum(re.res)) %>%
    arrange(desc(re.res))
  return(df)
}

get_res_zhu <- function(subtype, population, p1, p2,
                                data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  df$fit <- predict(mod_obj, type = "response")
  n_obs <- sum(df$n)
  df$re.res <- df$res / (2*n_obs)
  df <- df %>%
    select(p1, p2, re.res) %>%
    group_by(p1, p2) %>%
    summarize(re.res = sum(re.res)) %>%
    ungroup() %>%
    mutate(prop.re = re.res / sum(re.res)) %>%
    arrange(desc(re.res))
  return(df)
}

get_res_zhu_sign <- function(subtype, population, p1, p2,
                             data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  df$res <- residuals(mod_obj) ^ 2
  df$fit <- predict(mod_obj, type = "response")
  df$sign <- sign(df$n - df$fit)
  n_obs <- sum(df$n)
  df$re.res <- df$res / (2*n_obs)
  df$s_re <- df$re.res * df$sign
  return(df)
}

heatmap_zhu <- function(subtype, pop, p1, p2,
                        data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  df <- get_res_zhu_sign(subtype, pop, p1, p2, data_dir) %>%
    filter(status == "singletons")

  if(!str_starts(subtype, "cpg")){
    subtype2 <- str_replace(subtype, "_", " → ")
  } else {
    subtype <- str_sub(subtype, 5)
    subtype2 <- paste0("(cpg) ", str_replace(subtype, "_", " → "))
  }

  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = s_re)) +
    geom_tile() +
    ggtitle(paste0("Interaction Relative Entropy: ", subtype2),
            paste0("Population: ", pop))+
    xlab(paste0("Relative Position: ", p2)) +
    ylab(paste0("Relative Position: ", p1)) +
    labs(fill = "Signed RE Residual") +
    scale_fill_distiller(palette = "RdBu", direction = -1)
  return(p)
}

heatmap_re_res <- function(subtype, pop, p1, p2, singletons_only = TRUE,
                           data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  if(singletons_only){
    df <- get_residuals_re(subtype, pop, p1, p2, data_dir)
  } else {
    df <- get_residuals_re2(subtype, pop, p1, p2, data_dir)
  }
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = re.res)) +
    geom_tile() +
    ggtitle(paste0("Interaction RE: ", subtype),
            paste0("Population: ", pop))+
    xlab(paste0("Relative Position: ", p2)) +
    ylab(paste0("Relative Position: ", p1)) +
    labs(fill = "RE Residual") +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}

heatmap_signed_re_res <- function(subtype, pop, p1, p2,
                                  data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  df <- get_residuals_re(subtype, pop, p1, p2, data_dir)
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = s.res)) +
    geom_tile() +
    ggtitle(paste0("Interaction RE: ", subtype),
            paste0("Population: ", pop))+
    xlab(paste0("Relative Position: ", p2)) +
    ylab(paste0("Relative Position: ", p1)) +
    labs(fill = "RE Residual") +
    scale_fill_distiller(palette = "RdBu", direction = -1)
  return(p)

}

control_int_mult_res <- function(population, subtype, p1, p2,
                                 data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/") {

  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0)
  control_tab <- df %>%
    select(p1, p2, controls) %>%
    pivot_wider(names_from = p2, values_from = controls) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'p1')
  singleton_tab <- df %>%
    select(p1, p2, singletons) %>%
    pivot_wider(names_from = p2, values_from = singletons) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'p1')

  exp_s <- get_mar(singleton_tab) * get_int(control_tab)
  exp_p_s <- exp_s / sum(exp_s) # normalize
  exp_s <- exp_p_s * sum(singleton_tab)
  #chi_sq <- (exp_s - singleton_tab)^2 / exp_s

  p_s <- singleton_tab / sum(singleton_tab)
  re <- p_s * log(p_s / exp_p_s)
  final <- rownames_to_column(re, "p1") %>%
    pivot_longer(-p1, names_to = "p2", values_to = "statistic")

  return(final)
}

plot_control_mult_res <- function(population, subtype, p1, p2,
                                  data_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"){
  df <- control_int_mult_res(population, subtype, p1, p2, data_dir)
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu") +
    ggtitle("Singleton di-nucleotide Enrichment/Depletion") +
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill="RE")
  return(p)
}
