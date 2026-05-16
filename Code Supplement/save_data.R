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
  select(Record.ID, Effect.Measure) %>%
  rename(recid = Record.ID)

saveRDS(dat_nmadb, file.path("data/dat_nmadb.rds"))

# Note: Download failed for multiple nmadb datasets.
# These failures occurred inside readByID(), because some records
# could not be retrieved from the external database or the exported
# .xlsx files were empty/corrupted. Failed records were skipped.
# NOTE: This step takes about 30 minutes to run
idx_nmadb <- classify_nmadb(dat = dat_nmadb)
save(idx_nmadb, file = file.path("data/idx_nmadb.RData"))

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
# NOTE: This step takes about 5 minutes to run
twoarm_data_list <- map(idx_twoarm$recid, function(id) {
  readByID(id)
})

# Name the list by recid
names(twoarm_data_list) <- idx_twoarm$recid

# Save combined object as a single .rds file
saveRDS(
  twoarm_data_list,
  file = file.path("data/nmadb_twoarm_data_all.rds")
)
