library(dplyr)
library(tidyverse)

rm(list = ls())
df_c <- read.csv("./data/Sudan/Africa_data_except_Gambia_sort.csv")
df_t <- read.csv("./data/Sudan/Sudan_data.csv")

df_c <- df_c %>% select(-1, -2)

df_t <- df_t %>% select(-1)

df <- bind_rows(df_c, df_t)

df <- df %>%
  mutate(
    country_id = as.numeric(factor(country))
  ) %>%
  mutate(
    country_id = ifelse(country == "SDN", 0, country_id)
  ) %>%
  mutate(
    country = ifelse(country_id == 0, "Sudan", country)
  ) %>%
  select(last_col(), everything())

df %>% write_csv("./data/Sudan/africa_panel_data.csv")

weight_mat <- read.csv("./data/Sudan/weight_mat_avg.csv")

country_list <- df %>% select(country, country_id) %>% unique()

weight_mat <- weight_mat %>%
  inner_join(
    country_list,
    by = 'country'
  ) %>%
  select(-1) %>%
  select(last_col(), everything())

weight_mat %>% write_csv("./data/Sudan/weight_mat.csv")

weight_vec <- read.csv("./data/Sudan/weight_vec_avg.csv")

weight_vec <- weight_vec %>%
  inner_join(
    country_list,
    by = 'country'
  ) %>%
  select(-1) %>%
  select(last_col(), everything())

weight_vec %>% write_csv("./data/Sudan/weight_vec.csv")


df_t <- read_csv("./data/tabacco/smoking_treat.csv") %>% mutate(state_id = 0)
df_c <- read_csv("./data/tabacco/smoking_control.csv") %>%
  mutate(state_id = as.numeric(factor(state)))
df_c %>% select(-1)
df_t %>% select(-1)

df <- bind_rows(df_t, df_c) %>%
  select(-1)

df %>% write_csv("./data/tabacco/us_smoking_panel.csv")

state_list <- df %>% select(state_id, state) %>% unique()
weight_mat <- read_csv("./data/tabacco/adj_mat_smoking.csv") %>%
  rename('state' = '...1') %>%
  inner_join(state_list, by = 'state') %>%
  select(-1) %>%
  select(last_col(), everything())

weight_vec <- read_csv("./data/tabacco/adj_vec_smoking.csv") %>%
  rename('state' = '...1', 'value' = '0') %>%
  inner_join(state_list, by = 'state') %>%
  select(-1) %>%
  select(last_col(), everything())

weight_mat %>% write_csv("./data/tabacco/weight_mat.csv")
weight_vec %>% write_csv("./data/tabacco/weight_vec.csv")


df <- read_csv("./data/tabacco/us_smoking_panel.csv")
weight_mat <- read_csv("./data/tabacco/weight_mat.csv")
weight_vec <- read_csv("./data/tabacco/weight_vec.csv")

california_smoking <- list(
  panel = df,
  W = weight_mat,
  w = weight_vec
)

library(usethis)
use_data(california_smoking, overwrite = TRUE)

df <- read_csv("./data/Sudan/africa_panel_data.csv")
weight_mat <- read_csv("./data/Sudan/weight_mat.csv")
weight_vec <- read_csv("./data/Sudan/weight_vec.csv")

sudan_secession <- list(
  panel = df,
  W = weight_mat,
  w = weight_vec
)
use_data(sudan_secession, overwrite = TRUE)
