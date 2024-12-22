required_packages <- c("data.table", "sf", "geobr", "gstat", "RColorBrewer", "ggplot2", "viridis")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)


prediction_results <- fread("data/oob-results-original.csv")
nrow(prediction_results)

prediction_results <- unique(prediction_results, by = c("latitude", "longitude"))
nrow(prediction_results)
head(prediction_results)

brazil_states <- geobr::read_state(year = 2018)
brazil_states <- sf::st_transform(brazil_states, crs = 4326)

rondonia_cities <- geobr::read_municipality(code_muni = "RO", year = 2018)
rondonia_cities <- sf::st_transform(rondonia_cities, crs = 4326)

prediction_sf <- sf::st_as_sf(
  prediction_results,
  coords = c("latitude", "longitude"),
  crs = 4326
)

prediction_sf <- sf::st_intersection(prediction_sf, rondonia_cities)
print(prediction_sf)

portovelho_predictions <- prediction_sf[prediction_sf$name_muni == "Porto Velho", ]

error_breaks <- c(0, 25, 50, 100)
error_palette <- colorRampPalette(c("blue", "black", "red"))(length(error_breaks) - 1)
plot(
  portovelho_predictions["observed"],
  reset = FALSE,
  pal = error_palette, breaks = error_breaks,
  pch = 20, cex = 0, key.pos = 1, key.length = 1,
  graticule = TRUE, axes = TRUE
)
plot(
  rondonia_cities[rondonia_cities$name_muni == "Porto Velho", "abbrev_state"],  
  main = "Spatial Distribution of Prediction Errors",
  col = "transparent",
  border = "lightgray",
  reset = FALSE, 
  add = TRUE
)
plot(
  portovelho_predictions["observed"],
  add = TRUE,
  pal = error_palette, breaks = error_breaks,
  pch = 21,
  cex = 0.5
)

cat("Calculando estatísticas descritivas dos valores observados de COS...\n")
observed_stats <- portovelho_predictions %>%
  summarise(
    mean_cos = mean(observed),
    median_cos = median(observed),
    min_cos = min(observed),
    max_cos = max(observed),
    sd_cos = sd(observed)
  )
predictions_utm <- sf::st_transform(portovelho_predictions, crs = 31983)
print(predictions_utm)

semivariogram <- gstat::variogram(
  observed ~ 1,
  locations = predictions_utm,
  cutoff = 16000, 
  width = 1600 
)
print(semivariogram)
plot(semivariogram)

nugget_value <- 130  # Nugget effect for short-range variability
partial_sill <- 80   # Maximum semivariance (partial sill)
semivariogram_range <- 6000  # Distance where the semivariogram reaches the sill
vgm_model <- gstat::vgm(psill = partial_sill, model = "Sph", range = semivariogram_range, nugget = nugget_value)
fit_vgm_model <- gstat::fit.variogram(semivariogram, model = vgm_model)
print(fit_vgm_model)

semivariogram_line <- gstat::variogramLine(fit_vgm_model, maxdist = 16000)
plot(
  x = semivariogram[["dist"]], 
  y = semivariogram[["gamma"]], 
  main = "Semivariogram of Prediction Errors",
  xlab = "Distance (m)",
  ylab = "Semivariance",
  pch = 20, panel.first = grid(), ylim = c(0, max(semivariogram$gamma))
)
lines(semivariogram_line[["dist"]], semivariogram_line[["gamma"]], col = "blue")

rondonia_limit <- sf::read_sf('./data/vector/RO_Municipios_2023/RO_Municipios_2023.shp')
cidade_portovelho <- rondonia_limit %>% dplyr::filter(NM_MUN == "Porto Velho")
cidade_utm <- sf::st_transform(cidade_portovelho, crs = 31983)
grid_points <- sf::st_make_grid(x = cidade_utm, n = c(10, 10), crs = 31983)

kriging_prediction <- gstat::krige(
  observed ~ 1,
  locations = predictions_utm,
  newdata = grid_points,
  model = fit_vgm_model
)

ggplot() +
  geom_tile(data = kriging_prediction, aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Kriging Prediction of Prediction Errors",
       x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Prediction Error (t/ha)")

kriging_sf <- st_as_sf(kriging_prediction, coords = c("x", "y"), crs = st_crs(cidade_utm))

kriging_clipped <- st_intersection(kriging_sf, cidade_utm)

ggplot() +
  geom_tile(data = kriging_clipped, aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Kriging Prediction in Porto Velho Area",
       x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Prediction Error (t/ha)")
