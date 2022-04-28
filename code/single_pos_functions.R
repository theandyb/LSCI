# 1000G Functions
load_results <- function(subtype, rp, sp){
  data_dir <- "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/single_pos_df/"
  f_name <- paste0(data_dir, sp, "/", subtype, "_rp", rp, ".csv")
  df <- read_csv(f_name, col_types = cols())
  return(df)
}

load_all_results <- function(subtype, sp, r_start = 1){
  df <- load_results(subtype, -10, sp)
  df$rp <- -10
  for(i in c(-9:-1, r_start:10)){
    df2 <- load_results(subtype, i, sp)
    df2$rp <- i
    df <- bind_rows(df, df2)
  }
  return(df)
}

load_all_results_all_sp <- function(subtype, r_start = 1){
  df <- load_all_results(subtype, "AFR") %>%
    select(Nuc, rp, chi_sq_gw, chi_sq_ct) %>%
    rename_with(.fn = ~paste0(., ".AFR"), .cols = starts_with("chi"))
  for(sp in c("AMR", "EAS", "EUR", "SAS")){
    df2 <- load_all_results(subtype, sp, r_start) %>%
      select(Nuc, rp, chi_sq_gw, chi_sq_ct) %>%
      rename_with(.fn = ~paste0(., ".", sp), .cols = starts_with("chi"))
    df <- full_join(df, df2, by = c("Nuc", "rp"))
  }
  return(df)
}

statistic_by_position <- function(subtype, sp, r_start = 1){
  final <- data.frame(pos = numeric(), chi_sq = numeric(), singletons = numeric(),type = character())
  for(i in c(-10:-1, r_start:10)){
    df <- load_results(subtype, i, sp)
    final <- bind_rows(final,
                       data.frame(pos = c(i,i),
                                  chi_sq = c(sum(df$chi_sq_gw, na.rm = T), sum(df$chi_sq_ct, na.rm = T) ),
                                  singletons = c(sum(df$singletons, na.rm = T), sum(df$singletons, na.rm = T) ),
                                  type = c("Genome-wide", "Control")))
  }
  final <- final %>%
    mutate(cohen = sqrt(chi_sq / singletons))
  return(final)
}

plot_pos_stat <- function(subtype, sp, r_start = 1){
  df <- statistic_by_position(subtype, sp, r_start)
  p <-df %>%
    ggplot(aes(x = pos, y = chi_sq, colour = type)) +
    geom_point() +
    geom_line() +
    ggtitle(paste0("1000G Single-Position Results: ", str_replace_all(subtype, "_", "-")),
            paste0("Population: ", sp)) +
    xlab("Relative Position") +
    ylab("Chi Square Goodness of Fit Statistic") +
    labs(colour = "Background Rate")
  return(p)
}

plot_nuc_cont_by_position <- function(subtype, sp, type = "gw", r_start = 1){
  if(type == "gw"){
    df <- load_all_results(subtype, sp, r_start) %>%
      replace_na(list("chi_sq_gw" = 0)) %>%
      group_by(rp) %>%
      mutate(p_chi = chi_sq_gw / sum(chi_sq_gw))
    g_title <- paste0("Nucleotide Contribution to Genome-wide Statistic: ", str_replace_all(subtype, "_", "-"))
  } else{
    df <- load_all_results(subtype, sp, r_start) %>%
      replace_na(list("chi_sq_ct" = 0)) %>%
      group_by(rp) %>%
      mutate(p_chi = chi_sq_ct / sum(chi_sq_ct))
    g_title <- paste0("Nucleotide Contribution to Control-rate Statistic: ", str_replace_all(subtype, "_", "-"))
  }

  p <- df %>%
    ggplot(aes(x = rp, y = p_chi, fill = Nuc)) +
    geom_col() +
    xlab("Relative Position") +
    ylab("Proportion Contributed") +
    labs(fill = "Nucleotide") +
    ggtitle(g_title, paste0("Population: ", sp))
  return(p)
}

plot_signed_nuc_by_pos <- function(subtype, sp, type = "gw", r_start = 1){
  df <- load_all_results(subtype, sp, r_start)
  if(type == "gw"){
    df <- df %>%
      replace_na(list("chi_sq_gw" = 0)) %>%
      mutate(s_val = sign(singletons - exp_gw)) %>%
      mutate(signed_chi_sq = s_val * chi_sq_gw)
    g_title <- paste0("Nucleotide Genome-wide Signed ChiSq Residual: ", str_replace_all(subtype, "_", "-"))
  } else {
    df <- df %>%
      replace_na(list("chi_sq_ct" = 0)) %>%
      mutate(s_val = sign(singletons - exp_ct)) %>%
      mutate(signed_chi_sq = s_val * chi_sq_ct)
    g_title <- paste0("Nucleotide Control-Rate Signed ChiSq Residual: ", str_replace_all(subtype, "_", "-"))
  }
  p <- df %>%
    ggplot(aes(x = rp, y = signed_chi_sq, colour = Nuc)) +
    geom_point() +
    geom_line() +
    xlab("Relative Position") +
    ylab("Signed ChiSq Contribution") +
    labs(fill = "Nucleotide") +
    ggtitle(g_title, paste0("Population: ", sp))
  return(p)
}

### Functions for dealing with BRIDGES data

load_results_bridges <- function(subtype, rp){
  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/single_pos_df/"
  f_name <- paste0(data_dir, subtype, "_rp", rp, ".csv")
  df <- read_csv(f_name, col_types = cols())
  return(df)
}

load_all_results_bridges <- function(subtype, r_start = 1){
  df <- load_results_bridges(subtype, -10)
  df$rp <- -10
  for(i in c(-9:-1, r_start:10)){
    df2 <- load_results_bridges(subtype, i)
    df2$rp <- i
    df <- bind_rows(df, df2)
  }
  return(df)
}

statistic_by_position_bridges <- function(subtype, r_start = 1){
  final <- data.frame(pos = numeric(), chi_sq = numeric(), singletons = numeric(),type = character())
  for(i in c(-10:-1, r_start:10)){
    df <- load_results_bridges(subtype, i) %>%
      filter(singletons > 0)
    final <- bind_rows(final,
                       data.frame(pos = c(i,i),
                                  chi_sq = c(sum(df$chi_sq_gw), sum(df$chi_sq_ct) ),
                                  singletons = c(sum(df$singletons), sum(df$singletons) ),
                                  type = c("Genome-wide", "Control")))
  }
  final <- final %>%
    mutate(cohen = sqrt(chi_sq / singletons))
  return(final)
}
