rm(list = ls())
source("FILEPATH")
source("FILEPATH")
library(ggplot2)
library(RColorBrewer)
library(plotly)

out_dir <- "FILEPATH"

df <- fread("FILEPATH")
locs <- get_location_metadata(
  location_set_id = 35,
  gbd_round_id = 6,
  decomp_step = 'iterative'
)
draws <- names(df)[grepl("draw", names(df))]

df$prop <- rowMeans(df[, ..draws])
df[, (draws) := NULL]

# Add age group name and super region name
ages <- get_age_metadata(gbd_round_id = 6)
df <-
  merge(df, ages[, age_group_name, age_group_id], by = "age_group_id",
        all.x = TRUE)
stopifnot(!any(is.na(df$age_group_name)))
df <-
  merge(df, locs[, location_id, super_region_id], by = "location_id",
        all.x = TRUE)
stopifnot(!any(is.na(df$super_region_id)))
df <-
  merge(df,
        locs[, location_id, location_name_short],
        by.x = "super_region_id",
        by.y = 'location_id',
        all.x = TRUE)
stopifnot(!any(is.na(df$location_name_short)))
setnames(df, "location_name_short", "super_region_name")

for (measure in c(1, 6)) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  set.seed(42)
  colors <- sample(getPalette(length(unique(df$pathogen))))
  g <- ggplot(df[measure_id == measure],
              aes_string(fill = "pathogen", y = "prop", x = "super_region_name")) +
    geom_bar(position = "stack",
             stat = "summary",
             fun = mean) +
    scale_fill_manual(values = colors) +
    facet_wrap( ~ age_group_name) +
    theme(axis.text.x = element_text(angle = 90))
  fig <- ggplotly(g)
  htmlwidgets::saveWidget(fig, "FILEPATH")
}