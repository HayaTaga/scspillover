library(dplyr)
library(tidyverse)

df <- read_csv("./data/tobacco.csv")

state_list <- df %>% select(state_name) %>% unique()


df
w <- read_csv("./data/aws.csv")
W <- read_csv("./data/rew.csv")

california_smoking$w


california_smoking_latest <- list(
  panel_df = df,
  w = w,
  W = w
)

library(usethis)
use_data(california_smoking_latest, overwrite = TRUE)
