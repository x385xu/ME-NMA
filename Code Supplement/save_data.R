# install.packages("remotes") 
# remotes::install_version("nmadb", version = "1.2.0")
library(nmadb)
# install.packages("netmeta") 
library(netmeta)
# install.packages("dplyr") 
library(dplyr)
# install.packages("purrr") 
library(purrr)
source("functions/classify_nmadb.R")
source("functions/fit_netmeta.R")

if (!dir.exists("data")) {
  dir.create("data")
}

dat_nmadb <- getNMADB()

dat_nmadb <- dat_nmadb %>%
  select(Record.ID, Title, First.Author, Year,
         Number.of.Studies.., Number.of.Treatments,
         Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
         Primary.Outcome, Description.of.the.outcome.s.,
         Harmful.Beneficial.Outcome, dataset) %>%
  rename(recid = Record.ID) %>%
  mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")),
         .keep = "unused")

saveRDS(dat_nmadb, "data/dat_nmadb.rds")

# NOTE: This step may take several minutes to run
idx_nmadb <- classify_nmadb(dat = dat_nmadb)
save(idx_nmadb, file = "data/idx_nmadb.RData")

# Keep only two-arm networks with OR, RR, or MD
idx_twoarm <- idx_nmadb %>%
  filter(
    twoarm_only,
    effect_measure %in% c(
      "risk ratio",
      "odds ratio",
      "mean difference"
    )
  )

# Read and save datasets
# NOTE: This step may take several minutes to run
twoarm_data_list <- map(idx_twoarm$recid, function(id) {
  readByID(id)
})

# Name the list by recid
names(twoarm_data_list) <- idx_twoarm$recid

# Save combined object as a single .rds file


saveRDS(
  twoarm_data_list,
  file = "data/nmadb_twoarm_data_all.rds"
)
