# newest version of netmeta does not have the function pairwise(), which is required while using the nmadb package
# remove.packages("netmeta")
# install.packages("netmeta_2.9-0.tar.gz", repos = NULL, type = "source")
library(netmeta)
# nmadb is removed from cran, can find it in the archive: https://cran.r-project.org/src/contrib/Archive/nmadb/?C=D;O=A
# https://github.com/cran/nmadb
# install.packages(c("devtools", "RCurl", "readxl", "jsonlite"))
# install.packages("nmadb_1.2.0.tar.gz", repos = NULL, type = "source")
library(nmadb)
library(dplyr)
library(purrr)
source("classify_nmadb.R")

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

saveRDS(dat_nmadb, "dat_nmadb.rds")

idx_nmadb <- classify_nmadb(dat = dat_nmadb)
save(idx_nmadb, file = "idx_nmadb.RData")

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

# Create folder
dir.create("nmadb_twoarm_data", showWarnings = FALSE)

# Read and save datasets
twoarm_data_list <- map(idx_twoarm$recid, function(id) {
  
  cat("Reading dataset:", id, "\n")
  
  dat <- readByID(id)
  
  saveRDS(
    dat,
    file = file.path(
      "nmadb_twoarm_data",
      paste0("nmadb_", id, ".rds")
    )
  )
  
  dat
})

# Name the list by recid
names(twoarm_data_list) <- idx_twoarm$recid

# Save combined object
saveRDS(
  twoarm_data_list,
  file = "nmadb_twoarm_data_all.rds"
)
