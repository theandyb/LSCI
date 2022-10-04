# Two Position Model Functions

## 1KGP
deviance_pair <- function(population, subtype, p1, p2){

  data_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"
  f_name <- paste0(data_dir, population, "/", subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())

  return(mod_obj %>% deviance)
}

nuc_re <- function(population, subtype, p1, p2){

  data_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/all_count_2_pos/"
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
  return(df)
}

deviance_pair_re <- function(population, subtype, p1, p2){
  df <- nuc_re(population, subtype, p1, p2)
  return(sum(df$re.res))
}

deviance_all <- function(population, subtype, r_start = 1){
  final <- data.frame(p1 = numeric(),
                      p2 = numeric(),
                      statistic = numeric())
  for(i in c(-10:-1,r_start:9)){
    for(j in (i+1):10){
      if(j == 0) next
      if(j == 1 & r_start > 1) next
      stat_val <- deviance_pair(population, subtype, i, j)
      final <- bind_rows(final, data.frame(p1 = i,
                                           p2 = j,
                                           statistic = stat_val))
    }
  }
  return(final)
}

re_all <- function(population, subtype, r_start = 1){
  final <- data.frame(p1 = numeric(),
                      p2 = numeric(),
                      statistic = numeric())
  for(i in c(-10:-1,r_start:9)){
    for(j in (i+1):10){
      if(j == 0) next
      if(j == 1 & r_start > 1) next
      stat_val <- deviance_pair_re(population, subtype, i, j)
      final <- bind_rows(final, data.frame(p1 = i,
                                           p2 = j,
                                           statistic = stat_val))
    }
  }
  return(final)
}

deviance_plot <- function(population, subtype, r_start = 1){
  df <- deviance_all(population, subtype, r_start) %>% drop_na()
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    ggtitle(paste0("Deviance Interaction Test: ", subtype),
            paste0("Population: ", population,"; Min: ", round(min(df$statistic), 2), "; Max: ", round(max(df$statistic), 2)))+
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill = "Deviance") +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}

deviance_re_plot <- function(population, subtype, r_start = 1){
  df <- re_all(population, subtype, r_start) %>% drop_na()
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    ggtitle(paste0("Interaction RE: ", subtype),
            paste0("Population: ", population,"; Min: ", signif(min(df$statistic), 2), "; Max: ", signif(max(df$statistic), 2)))+
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill = "RE") +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}

## BRIDGES
deviance_pair_BRIDGES <- function(subtype, p1, p2){

  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/all_count_2_pos/"
  f_name <- paste0(data_dir, subtype, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    filter(singletons > 0) %>%
    select(p1, p2, singletons, controls) %>%
    gather(status, n, singletons:controls)
  mod_obj <- glm(n ~ (p1 + p2 + status)^2, data = df, family = poisson())
  return(mod_obj %>% deviance)
}

deviance_all_BRIDGES <- function(subtype, r_start = 1){
  final <- data.frame(p1 = numeric(),
                      p2 = numeric(),
                      statistic = numeric())
  for(i in c(-10:-1,r_start:9)){
    for(j in (i+1):10){
      if(j == 0) next
      if(j == 1 & r_start > 1) next
      stat_val <- deviance_pair_BRIDGES(subtype, i, j)
      final <- bind_rows(final, data.frame(p1 = i,
                                           p2 = j,
                                           statistic = stat_val))
    }
  }
  return(final)
}

deviance_plot_BRIDGES <- function(subtype, r_start = 1){
  df <- deviance_al_BRIDGES(subtype, r_start) %>% drop_na()
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    ggtitle(paste0("Deviance Interaction Test: ", subtype),
            paste0("Min: ", round(min(df$statistic), 2), "; Max: ", round(max(df$statistic), 2)))+
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill = "Deviance") +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}

nuc_re_BRIDGES <- function(subtype, p1, p2){
  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/all_count_2_pos/"
  f_name <- paste0(data_dir, subtype, "_p", p1, "_q", p2, ".csv")
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
  return(df)
}

deviance_pair_re_BRIDGES <- function(subtype, p1, p2){
  df <- nuc_re_BRIDGES(subtype, p1, p2)
  return(sum(df$re.res))
}

re_all_BRIDGES <- function(subtype, r_start = 1){
  final <- data.frame(p1 = numeric(),
                      p2 = numeric(),
                      statistic = numeric())
  for(i in c(-10:-1,r_start:9)){
    for(j in (i+1):10){
      if(j == 0) next
      if(j == 1 & r_start > 1) next
      stat_val <- deviance_pair_re_BRIDGES(subtype, i, j)
      final <- bind_rows(final, data.frame(p1 = i,
                                           p2 = j,
                                           statistic = stat_val))
    }
  }
  return(final)
}

deviance_re_plot_BRIDGES <- function(population, subtype, r_start = 1){
  df <- re_all_BRIDGES(population, subtype, r_start) %>% drop_na()
  p <- df %>%
    ggplot(aes(x = p2, y = p1, fill = statistic)) +
    geom_tile() +
    ggtitle(paste0("Interaction RE: ", subtype),
            paste0("Population: ", population,"; Min: ", signif(min(df$statistic), 2), "; Max: ", signif(max(df$statistic), 2)))+
    xlab("Relative Position 2") +
    ylab("Relative Position 1") +
    labs(fill = "RE") +
    scale_fill_distiller(palette = "Reds", direction = 1)
  return(p)
}
