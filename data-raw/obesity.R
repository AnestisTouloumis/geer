library("tidyverse")
library("readxl")


obesity <- read_excel("./data-raw/pcbi.1008108.s001.xls")
obesity <-
  obesity |>
  select(AGE, ZYGOSITY, HOST_SUBJECT_ID, FAMILY, ANCESTRY, OBESITYCAT, Bacteroides) |>
  arrange(FAMILY, HOST_SUBJECT_ID, AGE) |>
  filter(ZYGOSITY != "NA") |>
  rename(age = AGE,
         zygosity = ZYGOSITY,
         sid = HOST_SUBJECT_ID,
         fid = FAMILY,
         obesity = OBESITYCAT,
         bacteroides = Bacteroides,
         ancestry = ANCESTRY) |>
  relocate(fid, sid, obesity, age, zygosity, ancestry, bacteroides) |>
  mutate_at(c("fid", "sid", "obesity", "zygosity", "ancestry"), as.factor) |>
  as.data.frame()
usethis::use_data(obesity, overwrite = TRUE)
