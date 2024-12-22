if (!require(pacman)) install.packages("pacman")
pacman::p_load(data.table, sf, geobr, gstat, RColorBrewer, ggplot2, viridis, dplyr)

prediction_results <- fread("data/oob-results-original.csv")
cat("Número de registros iniciais: ", nrow(prediction_results), "\n")

unique_prediction_results <- unique(prediction_results, by = c("latitude", "longitude"))
cat("Número de registros após remoção de duplicados: ", nrow(unique_prediction_results), "\n")

head(unique_prediction_results)

brazil_states <- geobr::read_state(year = 2018)
brazil_states <- sf::st_transform(brazil_states, crs = 4326)

sf_prediction_results <- sf::st_as_sf(
  unique_prediction_results,
  coords = c("latitude", "longitude"),
  crs = 4326
)

sf_prediction_results <- sf::st_intersection(sf_prediction_results, brazil_states)
print(sf_prediction_results)

breaks <- c(0, 25, 50, 100)
color_palette <- colorRampPalette(c("blue", "black", "red"))(length(breaks) - 1)

plot(
  sf_prediction_results["observed"],
  reset = FALSE,
  pal = color_palette, breaks = breaks,
  pch = 20, cex = 0, key.pos = 1, key.length = 1,
  graticule = TRUE, axes = TRUE
)
plot(
  brazil_states[brazil_states$abbrev_state == "RO", "name_state"],
  col = "transparent",
  border = "lightgray",
  reset = FALSE, add = TRUE
)
plot(
  sf_prediction_results["observed"],
  add = TRUE,
  pal = color_palette, breaks = breaks,
  pch = 21,
  cex = 0.5
)

cat("Calculando estatísticas descritivas dos valores observados de COS...\n")
observed_stats <- sf_prediction_results %>%
  summarise(
    mean_cos = mean(observed),
    median_cos = median(observed),
    min_cos = min(observed),
    max_cos = max(observed),
    sd_cos = sd(observed)
  )

sf_prediction_results_utm <- sf::st_transform(sf_prediction_results, crs = 31983)
print(sf_prediction_results_utm)

colnames(sf_prediction_results_utm)

ro_points <- sf_prediction_results_utm[["abbrev_state"]] == "RO"

cutoff_distance <- 1000  # Testar distâncias de 1000 m até 20000 m
semivariogram <- gstat::variogram(
  observed ~ 1,
  locations = sf_prediction_results_utm[ro_points, ],
  cutoff = cutoff_distance, # Distância máxima para estimar o semivariograma (em metros)
  width = cutoff_distance / 10 # Largura das caixas de distância (em metros)
)
print(semivariogram)

plot(
  x = semivariogram[["dist"]], # Distância
  y = semivariogram[["gamma"]], # Semivariança
  xlab = "Distância (m)",
  ylab = "Semivariança",
  pch = 20, panel.first = grid(), ylim = c(0, max(semivariogram$gamma))
)
grid()

psill <- 200 # Sill parcial é a máxima semivariança
range_distance <- 500 # A distância em que o semivariograma atinge o sill
nugget <- 0 # Valor do nugget (variabilidade a distâncias muito curtas)

variogram_model <- gstat::vgm(psill = psill, model = "Sph", range = range_distance, nugget = nugget)
fitted_variogram <- gstat::fit.variogram(semivariogram, model = variogram_model)
print(fitted_variogram)

variogram_model_line <- gstat::variogramLine(variogram_model, maxdist = cutoff_distance)
fitted_variogram_line <- gstat::variogramLine(fitted_variogram, maxdist = cutoff_distance)

plot(
  x = semivariogram[["dist"]],
  y = semivariogram[["gamma"]],
  xlab = "Distância (m)",
  ylab = "Semivariança",
  pch = 20, panel.first = grid(), ylim = c(0, max(semivariogram$gamma))
)
lines(
  x = variogram_model_line[["dist"]],
  y = variogram_model_line[["gamma"]],
  col = "red"
)
lines(
  x = fitted_variogram_line[["dist"]],
  y = fitted_variogram_line[["gamma"]],
  col = "blue"
)
legend(
  "bottomright",
  legend = c("Modelo Inicial", "Modelo Ajustado"),
  col = c("red", "blue"),
  lty = 1
)

municipality_limits <- sf::read_sf('./data/vector/RO_Municipios_2023/RO_Municipios_2023.shp')
municipality_limits <- st_transform(municipality_limits, crs = 31983)
print(municipality_limits)

grid <- sf::st_make_grid(x = municipality_limits, 
                         n = c(30, 30),
                         crs = 31983)
print(grid)

grid_points_utm <- st_transform(grid, crs = 31983)
grid_points_utm <- st_transform(grid_points_utm, crs = st_crs(sf_prediction_results_utm))

prediction <- gstat::krige(
  observed ~ 1,
  locations = sf_prediction_results_utm[ro_points, ],
  newdata = grid_points_utm,
  model = fitted_variogram
)

ggplot() +
  geom_tile(data = prediction, aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Previsão de Erros de Previsão pelo Método Kriging",
       x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Erro de Previsão (t/ha)")

prediction_sf <- st_as_sf(prediction, coords = c("x", "y"), crs = st_crs(municipality_limits))
prediction_cropped <- st_intersection(prediction_sf, municipality_limits)

ggplot() +
  geom_sf(data = municipality_limits, fill = "gray", color = "black", alpha = 0.2) +
  geom_sf(data = prediction_cropped, aes(fill = var1.pred), color = NA, alpha = 0.8) +
  scale_fill_viridis(name = "COS Preditivo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right") +
  labs(title = "Kriging de Valores Observados de COS para Rondônia",
       x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "COS Preditivo (t/ha)")

ggplot() +
  geom_sf(data = municipality_limits, fill = "gray", color = "black", alpha = 0.2) +
  geom_sf(data = prediction_cropped, aes(fill = var1.pred), color = NA, alpha = 0.8) +
  geom_sf(data = sf_prediction_results_utm[ro_points, ], aes(color = "Observado"), shape = 21, size = 0.5, alpha = 0.6) +
  scale_fill_viridis(name = "Erro Predito") +
  scale_color_manual(values = c("Observado" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right") +
  labs(title = "Kriging de Valores Observados de COS para Rondônia",
       x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "COS Preditivo (t/ha)",
       color = "Ponto de Observação")
