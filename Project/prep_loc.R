library(data.table)
library(ggplot2)
library(haven)

rm(list = ls())
ihme_loc <- "ETH"

root <- if (grepl("apple", sessionInfo()$platform)) "/Volumes/snfs/" else "/home/j/"
geospatial_u5m_dir <- paste0(root, "WORK/11_geospatial/02_processed data/U5M/")
fertilitypop_asfr_dir <- paste0(root, "WORK/02_fertilitypop/fertility/gbd_2017/data/microdata/cbh/ETH/")

# read in under 5 mortality data in order to get lat long
u5m_data <- readRDS(paste0(geospatial_u5m_dir, "Global/Global_CBH_1988_2017.RDS"))
u5m_data <- unique(u5m_data[country == "ETH", list(nid, psu = cluster_number, latnum, longnum)])
u5m_data[, psu := as.integer(psu)]

# read in prepped fertility data
asfr_files <- list.files(fertilitypop_asfr_dir, full.names = T)
asfr_data <- lapply(asfr_files, function(file) {
  data <- read_dta(file)
  data <- data.table(data)
  return(data)
})
asfr_data <- rbindlist(asfr_data, use.names = T, fill = T)
asfr_data[, psu := as.integer(psu)]

# merge on latitude and longitude
# asfr_data <- merge(asfr_data, u5m_data, by = c("nid", "psu"), all.x = T)

asfr_data <- asfr_data[nid == 218568, list(nid, ihme_loc_id = iso3, year_start, year_end,
                                           psu, hh_id, mother_id, child_id,
                                           interview_date_cmc, child_dob_cmc, mother_dob_cmc,
                                           strata, pweight)]

births <- copy(asfr_data)
births[, period_of_birth := interview_date_cmc - child_dob_cmc]
births <- births[between(period_of_birth, 1, 36)]
births[, mother_age := (child_dob_cmc - mother_dob_cmc) / 12]
births[, mother_age := plyr::round_any(mother_age, 5, floor)]
births <- births[, list(births = .N),
                 by = c("nid", "ihme_loc_id", "year_start", "year_end", "psu", "mother_age")]

exposure <- copy(asfr_data)
exposure[, age_end_period := (interview_date_cmc - 1) - mother_dob_cmc]
exposure <- exposure[, list(mother_age = c(plyr::round_any(age_end_period / 12, 5, floor),
                                           plyr::round_any(age_end_period / 12, 5, floor) - 5),
                            age_group = c("higher", "lower")),
                     by = c("nid", "ihme_loc_id", "year_start", "year_end", "psu", "age_end_period", "mother_id")]
exposure[age_group == "higher", person_months := age_end_period - (12 * mother_age) + 1]
exposure[age_group == "higher" & person_months > 36, person_months := 36]
exposure[age_group == "lower", person_months := 36 - (age_end_period - (12 * (mother_age + 5)) + 1)]
exposure[age_group == "lower" & person_months < 0, person_months := 0]
exposure[, person_years := person_months / 12]

exposure <- exposure[, list(person_years = sum(person_years)),
                     by = c("nid", "ihme_loc_id", "year_start", "year_end", "psu", "mother_age")]

cluster_data <- merge(births, exposure, all = T, by = c("nid", "ihme_loc_id", "year_start", "year_end", "psu", "mother_age"))
cluster_data <- cluster_data[between(mother_age, 15, 45)]

# dcast out and then melt in order to get rows for all combinations of psu and maternal age group
cluster_data <- dcast(cluster_data, nid + ihme_loc_id + year_start + year_end + psu ~ mother_age,
                      value.var = c("births", "person_years"))
cluster_data <- melt(cluster_data, id.vars = c("nid", "ihme_loc_id", "year_start", "year_end", "psu"),
               measure.vars = patterns("births", "person_years"),
               variable.name = "mother_age", value.name = c("births", "person_years"))
cluster_data[, mother_age := factor(mother_age, labels = seq(15, 45, 5))]
cluster_data[, mother_age := as.integer(as.character(mother_age))]
cluster_data[is.na(births), births := 0]
cluster_data[, asfr := births / person_years]

total <- cluster_data[, list(births = sum(births, na.rm = T),
                             person_years = sum(person_years, na.rm = T)),
                      c("nid", "ihme_loc_id", "year_start", "year_end", "mother_age")]
total[, asfr := births / person_years]
total[, psu := 1]
ggplot(cluster_data[psu %in% 1:5], aes(x = mother_age, y = asfr, group = psu, colour = as.factor(psu))) +
  geom_point(data = total, aes(x = mother_age, y = asfr), size = 5.0, colour = "red") +
  geom_line(data = total, aes(x = mother_age, y = asfr), size = 2.0, colour = "red") +
  geom_point() + geom_line() +
  theme_bw() + ylim(0, 1.0)

loc_xy <- u5m_data[nid == 218568 & !is.na(latnum) & !is.na(longnum)]


gaul_list <- 79
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

# Make triangulated mesh
mesh = inla.mesh.create(loc_xy)
spde = inla.spde2.matern(mesh)


cluster_data <- merge(cluster_data, u5m_data, by = c("nid", "psu"), all.x = T)
