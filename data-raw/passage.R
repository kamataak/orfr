## code to prepare `passage` dataset goes here

library(tidyverse)

# load("C:/Users/jnese/Desktop/BRT/GRANT-Aki/project/FinalPackage/jftn/ex_data2r.Rdata")

# Data preparation notes: ----
## Changes from passage.dat6
### only Year 5 (2018-19)
### only 100 students per grade
### "occasion" replaces "season"
### "id.passage" replaces "passage_id"
### "sec" replaces "secs"
### "numwords.pass" replaces "nwords.p"
### "id.student" replaces "student_id"; "id.student" is transformed
### only seasons 1 (as fall), 3 (as winter), and 4 (spring)

# stu_season_id2 == paste(student_id, year, season, sep = "_")
# In CODE and VIGNETTE, create a "stu_season_id2" ID by paste(student_id, occasion, sep = "_")

#Grade associated with student
#"season" be a numeric or character (unordered/leveled)

set.seed(3000)
passage <- passage.dat6 %>%
  filter(year == 5) %>% # Year 5 (2018-19) only
  group_by(grade, student_id) %>%
  nest() %>%
  group_by(grade) %>%
  slice_sample(n = 100) %>% # only 100 students per grade
  unnest(cols = "data") %>%
  ungroup() %>%
  mutate(
    occasion = case_when( # "occasion" replaces "season"
      season == 1 ~ "fall",
      season == 2 ~ "drop",
      season == 3 ~ "winter",
      season == 4 ~ "spring"),
  ) %>%
  rename(
    id.passage = passage_id, # "id.passage" replaces "passage_id"
    sec = secs, # "sec" replaces "secs"
    numwords.pass = nwords.p # "numwords.pass" replaces "nwords.p"
  ) %>%
  mutate(id.student = str_pad((student_id * 2) + 12, width = 4, side = "left", pad = "0")) %>%  # "id.student" replaces "student_id"; "id.student" is transformed
  filter(
    occasion != "drop" # seasons 1, 3, 4 only
  ) %>%
  select(id.student, occasion, grade, id.passage, numwords.pass, wrc, sec)

#usethis::use_data(passage, passage2, overwrite = TRUE)
