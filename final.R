#==================================================
# WoSIS AFRICA SOIL SPATIAL ANALYSIS PIPELINE
# Complete, Structured & Thesis Ready
#==================================================


#==================================================
# SECTION 0: LIBRARIES installing & GLOBAL SETTINGS
#==================================================
install.packages(c("dplyr", "sf", "ggplot2", "moments", "corrplot", "cluster", "FactoMineR", "factoextra", "tmap",
                 "gstat", "sp", "spdep", "terra", "patchwork"))
library(dplyr)
library(sf)
library(ggplot2)
library(moments)
library(corrplot)
library(cluster)
library(FactoMineR)
library(factoextra)
library(tmap)
library(gstat)
library(sp)
library(spdep)
library(terra)
library(patchwork)
library(tidyr)     # for pivot_longer

set.seed(123)


#==================================================
# SECTION 1: DATA INPUT
#==================================================
soil_data <- read.csv("C:/Users/DELL/Documents/wosis_2023_full_combined_africa_only.csv")


#==================================================
# SECTION 2: DATA QUALITY CONTROL
#==================================================

#--------------------------------------------------
# 2.1 Missing Value Detection
#--------------------------------------------------
null_percent <- colMeans(is.na(soil_data)) * 100
print(null_percent)

#--------------------------------------------------
# 2.2 Remove Variables with ≥50% Missing
#--------------------------------------------------
vars_keep <- names(null_percent)[null_percent < 50]
soil_clean <- soil_data %>% select(all_of(vars_keep))
print(names(soil_clean))
str(soil_clean)

#--------------------------------------------------
# 2.3 Imputation (Median / Mean)
#--------------------------------------------------
median_vars <- c( "NITKJD", "ORGC")
for (col in median_vars) {
  if (col %in% names(soil_clean)) {
    soil_clean[[col]][is.na(soil_clean[[col]])] <- median(soil_clean[[col]], na.rm = TRUE)
  }
}

mean_vars <- c("CLAY", "SILT", "SAND")
for (col in mean_vars) {
  soil_clean[[col]][is.na(soil_clean[[col]])] <- mean(soil_clean[[col]], na.rm = TRUE)
}

ph_cec_vars <- c("PHAQ")
for (col in ph_cec_vars) {
  if (col %in% names(soil_clean)) {
    soil_clean[[col]][is.na(soil_clean[[col]])] <- median(soil_clean[[col]], na.rm = TRUE)
  }
}

#--------------------------------------------------
# 2.4 Remove Rows Missing Key Variables
#--------------------------------------------------
key_vars <- c("CLAY", "SILT", "SAND", "ORGC", "PHAQ", "NITKJD")
key_vars <- key_vars[key_vars %in% names(soil_clean)]

soil_clean <- soil_clean %>%
  filter(rowMeans(is.na(select(., all_of(key_vars)))) == 0)

cat("Rows after cleaning:", nrow(soil_clean), "\n")

#--------------------------------------------------
# 3.3 Spatial Duplicate Check(same cooredinates)
#--------------------------------------------------
dup_coords <- soil_clean%>%
  st_drop_geometry() %>%
  group_by(longitude, latitude) %>%
  filter(n() > 1)

cat("Total observations:", nrow(soil_clean), "\n")
cat("Duplicate records:", nrow(dup_coords), "\n")

#==================================================
# SECTION 3: SPATIAL PREPARATION - HANDLE DUPLICATES (Average All Depths)
#==================================================
soil_unique <- soil_clean %>%
  group_by(longitude, latitude) %>%
  summarise(across(c(CLAY, SILT, SAND, ORGC, NITKJD, PHAQ), mean, na.rm = TRUE),
            lower_depth = mean(lower_depth),
            country_name = first(country_name),
            .groups = 'drop')

# Convert to sf using unique sites
soil_sf <- st_as_sf(soil_unique, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
soil_sf_utm <- st_transform(soil_sf, crs = 32637)

cat("Unique sites:", nrow(soil_sf_utm), "\n")

#--------------------------------------------------
# 3.1 Reproject to Africa UTM Zone 37N (meters)
#--------------------------------------------------
soil_sf_utm <- st_transform(soil_sf_utm, crs = 32637)
st_crs(soil_sf_utm)

#-------------------------------------------------
# Visualizing distributions of all soil properties before transformation
#-------------------------------------------------
soil_unique %>%
  select(CLAY, NITKJD, ORGC, PHAQ, SAND, SILT) %>%
  pivot_longer(cols = everything(), names_to = "Property", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Property, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Soil Properties", x = "Value", y = "Frequency")

#==================================================
# SECTION 4: STATISTICAL PREPROCESSING
#==================================================

#--------------------------------------------------
# 4.1 Skewness & Transformation
#--------------------------------------------------
numeric_vars <- c("CLAY", "SILT", "SAND", "PHAQ", "ORGC", "NITKJD")
numeric_vars <- numeric_vars[numeric_vars %in% names(st_drop_geometry(soil_sf_utm))]

numeric_df <- st_drop_geometry(soil_sf_utm) %>% select(all_of(numeric_vars))

skew_values <- sapply(numeric_df, function(x) skewness(x, na.rm = TRUE))
print("Skewness:")
print(round(skew_values, 3))

for (v in numeric_vars) {
  sk <- skew_values[v]
  if (abs(sk) > 1) {
    soil_sf_utm[[paste0(v, "_trans")]] <- log(soil_sf_utm[[v]] + 1)
    cat(v, "skewed:", round(sk, 2), "- log+1 applied\n")
  } else {
    soil_sf_utm[[paste0(v, "_trans")]] <- soil_sf_utm[[v]]
    cat(v, "not skewed:", round(sk, 2), "- kept original\n")
  }
}
transformed_vars <- paste0(numeric_vars, "_trans")

#--------------------------------------------------
# 4.2 Outlier Capping (IQR)
#--------------------------------------------------
cap_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQRv <- Q3 - Q1
  x[x < (Q1 - 1.5 * IQRv)] <- Q1 - 1.5 * IQRv
  x[x > (Q3 + 1.5 * IQRv)] <- Q3 + 1.5 * IQRv
  x
}

for (v in transformed_vars) {
  soil_sf_utm[[v]] <- cap_outliers(soil_sf_utm[[v]])
}

#--------------------------------------------------
# 4.3 Distribution Plots after transformation 
#--------------------------------------------------
par(mfrow = c(2, 4))
for (v in transformed_vars) {
  hist(soil_sf_utm[[v]], main = v, col = "skyblue")
}


#==================================================
# SECTION 5: MULTIVARIATE PREPARATION
#==================================================

#--------------------------------------------------
# 5.1 Scaling
#--------------------------------------------------
scaled_vars <- scale(st_drop_geometry(soil_sf_utm)[, transformed_vars])

#--------------------------------------------------
# 5.2 Correlation Analysis (Numbers + Plot)
#--------------------------------------------------
corr_matrix <- cor(scaled_vars, method = "spearman", use = "complete.obs")
print(round(corr_matrix, 3))

corrplot(corr_matrix,
         method = "color",
         addCoef.col = "black",
         type = "upper")

#--------------------------------------------------
# 5.3 Depth Summary
#--------------------------------------------------
soil_sf_utm$depth_class <- cut(soil_sf_utm$lower_depth,
                               breaks = c(0, 30, 60, 100, 200),
                               labels = c("0–30", "30–60", "60–100", ">100"),
                               include.lowest = TRUE)

soil_sf_utm %>%
  st_drop_geometry() %>%
  group_by(depth_class) %>%
  summarise(mean_ORGC = mean(ORGC_trans, na.rm = TRUE),
            var_ORGC = var(ORGC_trans, na.rm = TRUE))


#==================================================
# SECTION 6: CLUSTERING using k-mean clustering
#==================================================
#  -> Groups similar soil profiles/sites based on properties (e.g., texture, pH, carbon).
#  -> To cluster: Run as-is for 4 clusters. Change centers for 
#        different number (use elbow plot to choose optimal
#-------------------------------------------------
kmeans_res <- kmeans(scaled_vars, centers = 4, nstart = 25)
soil_sf_utm$cluster <- factor(kmeans_res$cluster)
fviz_cluster(kmeans_res, data = scaled_vars)

#==================================================
# PHASE 1: EXPLORATORY SPATIAL DATA ANALYSIS (ESDA)
#==================================================

#-------------------------------
# 0. Ensure ORGC_trans exists
#-------------------------------
if (!"ORGC_trans" %in% names(soil_sf_utm)) {
  soil_sf_utm$ORGC_trans <- log(soil_sf_utm$ORGC + 1)
}

#-------------------------------
# 1. Coordinate Reference System, Scale & Resolution
#-------------------------------
st_crs(soil_sf_utm)            # Check CRS
st_bbox(soil_sf_utm)           # Bounding box of data
coords <- st_coordinates(soil_sf_utm)
mean(dist(coords[1:1000, ]))   # Approximate average distance (sample of 1000 points)

#-------------------------------
# 2. Stationarity Check (Mean & Variance)
#-------------------------------
# By depth class
soil_sf_utm %>%
  st_drop_geometry() %>%
  group_by(depth_class) %>%
  summarise(
    mean_ORGC = mean(ORGC_trans, na.rm = TRUE),
    var_ORGC = var(ORGC_trans, na.rm = TRUE)
  )

# By latitude band
soil_sf_utm$lat_band <- cut(soil_sf_utm$latitude, breaks = 5)
soil_sf_utm %>%
  st_drop_geometry() %>%
  group_by(lat_band) %>%
  summarise(
    mean_ORGC = mean(ORGC_trans, na.rm = TRUE),
    var_ORGC = var(ORGC_trans, na.rm = TRUE)
  )

#-------------------------------
# 4. Trend Visualization
#-------------------------------
# North-South trend
ns_plot <- ggplot(st_drop_geometry(soil_sf_utm), aes(x = ORGC_trans, y = latitude)) +
  geom_point(alpha = 0.4, size = 1.2, color = "darkblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) +
  labs(title = "North-South Trend", x = "log(ORGC + 1)", y = "Latitude") +
  theme_minimal(base_size = 12)

# East-West trend
ew_plot <- ggplot(st_drop_geometry(soil_sf_utm), aes(x = longitude, y = ORGC_trans)) +
  geom_point(alpha = 0.4, size = 1.2, color = "darkgreen") +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) +
  labs(title = "East-West Trend", x = "Longitude", y = "log(ORGC + 1)") +
  theme_minimal(base_size = 12)

# Combine plots
ns_plot + ew_plot

#===============================
# 5. Detrending (CRITICAL STEP)
#===============================

trend_model <- lm(
ORGC_trans ~ longitude + latitude,
data = st_drop_geometry(soil_sf_utm)
)

soil_sf_utm$ORGC_resid <- residuals(trend_model)
#-----------------
# Residual summary (should be mean ≈ 0)
#-----------------
summary(soil_sf_utm$ORGC_resid)

#-----------------------------------------------------
# 2. Distribution Check: Histogram of transformed ORGC
#-----------------------------------------------------
ggplot(st_drop_geometry(soil_sf_utm), aes(x = ORGC_trans)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Distribution of log(ORGC + 1)",
       x = "log(ORGC + 1) [g/kg]",
       y = "Frequency") +
  theme_minimal(base_size = 14)
soil_sp <- as(soil_sf_utm, "Spatial")

#------------------------------
# stationarity residual check
#-----------------------------
ggplot(st_drop_geometry(soil_sf_utm),
aes(longitude, ORGC_resid)) +
geom_point(alpha = 0.4) +
geom_smooth(method = "loess", se = FALSE, color = "blue") +
labs(title = "Residuals After Detrending (East–West)",
x = "Longitude",
y = "ORGC residuals") +
theme_minimal()

ggplot(st_drop_geometry(soil_sf_utm),
aes(latitude, ORGC_resid)) +
geom_point(alpha = 0.4) +
geom_smooth(method = "loess", se = FALSE, color = "blue") +
labs(title = "Residuals After Detrending (North–South)",
x = "Latitude",
y = "ORGC residuals") +
theme_minimal()
#-------------------------------
# 5. Distribution Check of ORGC Residuals
#-------------------------------
ggplot(st_drop_geometry(soil_sf_utm),
aes(x = ORGC_resid)) +
geom_histogram(bins = 50,
fill = "steelblue",
color = "black",
alpha = 0.8) +
labs(title = "Distribution of ORGC Residuals",
x = "Residual value",
y = "Frequency") +
theme_minimal()
#-------------------------------
# 6. Visual Inspection (Interactive Maps)
#-------------------------------
tmap_mode("view")
for (v in transformed_vars) {
  map <- tm_shape(soil_sf_utm) +
    tm_dots(col = v,
            palette = "viridis",
            style = "quantile",
            size = 0.08,
            alpha = 0.8,
            popup.vars = c("Value" = v, "Country" = "country_name", "Depth" = "lower_depth"),
            title = v) +
    tm_layout(title = paste("Africa:", v),
              legend.position = c("left", "bottom"),
              frame = FALSE)
  print(map)
}
#===================================================
#phase 2:Spatial Dependence Modeling 
#===================================================

# Convert to sp for gstat (variogram needs SpatialPointsDataFrame)
soil_sp <- as(soil_sf_utm, "Spatial")

# 1. Experimental Variogram for ORGC_trans
v_emp <- variogram(ORGC_trans ~ 1, data = soil_sp, cutoff = 500000, width = 10000)  # cutoff ~500km, bin 10km

plot(v_emp, main = "Experimental Variogram - ORGC_trans")

# 2. Fit models (Spherical, Exponential, Gaussian)
v_fit_sph <- fit.variogram(v_emp, vgm("Sph"))
v_fit_exp <- fit.variogram(v_emp, vgm("Exp"))
v_fit_gau <- fit.variogram(v_emp, vgm("Gau"))

cat("Spherical fit:\n"); print(v_fit_sph)
cat("Exponential fit:\n"); print(v_fit_exp)
cat("Gaussian fit:\n"); print(v_fit_gau)

# Plot fitted
plot(v_emp, v_fit_sph, main = "Fitted Spherical Variogram")
plot(v_emp, v_fit_exp, main = "Fitted Exponential Variogram")
plot(v_emp, v_fit_gau, main = "Fitted Gaussian Variogram")

#-----------------------------------------------------
# 3. Spatial Weights Matrix - K-Nearest Neighbors (k=8)
#-----------------------------------------------------
# Fixed: Spatial Weights Matrix - K-Nearest Neighbors (k=8)

coords <- st_coordinates(soil_sf_utm)
knn8 <- knearneigh(coords, k = 8)
nb_knn8 <- knn2nb(knn8)
lw_knn <- nb2listw(nb_knn8, style = "W", zero.policy = TRUE)

# Plot connections
plot(st_geometry(soil_sf_utm), pch = 20, cex = 0.3, col = "gray")
plot(nb_knn8, coords, add = TRUE, col = "red", lwd = 0.8,
     main = "KNN Neighborhood (k=8)")

# Summary
summary(nb_knn8, zero.policy = TRUE)
#---------------------
# Moran's I
#---------------------
for (var in transformed_vars) {
  cat("\nMoran's I -", var, ":\n")
  print(moran.test(soil_sp[[var]], lw_knn))
}
#------------
# Nugget-sill
#------------
n_bins <- length(unique(v_emp$dist))
cat("Number of lag bins:", n_bins, "\n")
nugget <- v_emp$gamma[1]
sill <- max(v_emp$gamma)
ns_ratio <- nugget / sill
cat("Nugget:", nugget, "Sill:", sill, "Nugget/Sill ratio:", ns_ratio, "\n")
if(ns_ratio > 0.75){
  cat("High nugget effect: strong microscale variability or measurement error\n")
} else if(ns_ratio < 0.25){
  cat("Low nugget effect: strong spatial autocorrelation\n")
} else{
  cat("Moderate spatial autocorrelation\n")
}
#==================================================
# Phase 3: Interpolation and Prediction
#==================================================
# ------------------------------------------------------
# 1. Variogram (already fitted - just confirm)
# ------------------------------------------------------
soil_sp <- as(soil_sf_utm, "Spatial")
vg_exp <- variogram(ORGC_trans ~ 1, soil_sp, cutoff = 500000, width = 10000)
best_vgm <- fit.variogram(vg_exp, vgm("Sph"))  # Change to "Exp" or "Gau" if better

print(best_vgm)
#--------------------------
# model      psill  range
#1  Nug  0.06617613      0
#2  Sph  0.14187805 683302
#---------------------------

# ------------------------------------------------------
# 2. Create 50 km prediction grid
# ------------------------------------------------------
grid_pts <- st_make_grid(soil_sf_utm, cellsize = 50000, what = "centers") %>% 
  st_sf()

cat("Grid points:", nrow(grid_pts), "\n")  # ~29,700

# ------------------------------------------------------
# 3. Perform local Ordinary Kriging
# ------------------------------------------------------
soil_sp <- as(soil_sf_utm, "Spatial")

krige_pred <- krige(ORGC_trans ~ 1,
                    locations = soil_sp,
                    newdata = as(grid_pts, "Spatial"),
                    model = best_vgm,
                    nmax = 40,
                    maxdist = 200000)

# ------------------------------------------------------
# 4. Convert results & back-transform
# ------------------------------------------------------
krige_sf <- st_as_sf(krige_pred) %>%
  rename(prediction = var1.pred,
         variance   = var1.var) %>%
  mutate(ORGC_pred = exp(prediction) - 1)   # back to g/kg

#----------------------------------
cv <- krige.cv(ORGC_trans ~ 1,
               locations = soil_sp,
               model = best_vgm,
               nmax = 40,
               nfold = 10,
               verbose = FALSE)


cv_clean <- cv[!is.na(cv$var1.pred), ]
RMSE <- sqrt(mean(cv_clean$residual^2))
ME <- mean(cv_clean$residual)
MAE <- mean(abs(cv_clean$residual))
cat("Valid CV points:", nrow(cv_clean), "\n")
cat("RMSE:", round(RMSE, 3), "\n")
cat("ME:", round(ME, 3), "\n")
cat("MAE:", round(MAE, 3), "\n")
print(summary(cv))

# ------------------------------------------------------
# 5. Rasterize for nice maps
# ------------------------------------------------------
rast_template <- rast(ext(soil_sf_utm), resolution = 50000, 
                      crs = st_crs(soil_sf_utm)$wkt)

pred_rast <- rasterize(krige_sf, rast_template, field = "ORGC_pred", fun = "mean")
var_rast  <- rasterize(krige_sf, rast_template, field = "variance",   fun = "mean")

names(pred_rast) <- "ORGC_pred"
names(var_rast)  <- "variance"

# ------------------------------------------------------
# 6. Phase 3: Final Maps
# ------------------------------------------------------
tmap_mode("plot")

pred_map <- tm_shape(pred_rast) +
  tm_raster(col = "ORGC_pred",
            style = "quantile", n = 7,
            palette = "viridis",
            title = "Predicted ORGC (g/kg)") +
  tm_shape(soil_sf_utm) +
  tm_dots(size = 0.0, col = "gray", alpha = 0.4) +
  tm_layout(main.title = "Ordinary Kriging - Soil Organic Carbon (Africa)",
            main.title.size = 1.2,
            legend.outside = TRUE,
            legend.outside.position = "right")
print(pred_map)

#----------------------------------------------------
unc_map <- tm_shape(var_rast) + tm_raster(col = "variance",
     style = "quantile", n = 7, palette = "Reds", title = "Kriging Variance") + 
     tm_layout(main.title = "Kriging Uncertainty", main.title.size = 1.2,
     legend.outside = TRUE, legend.outside.position = "right") 

# Identify high-uncertainty regions (top 10% variance)

high_uncert <- krige_sf %>% 
  filter(variance >= quantile(variance, 0.90, na.rm = TRUE))
#-------------------------------------------------
# Add high-uncertainty overlay to your existing map
#--------------------------------------------------
unc_map_high <- unc_map + 
  tm_shape(high_uncert) +
  tm_borders(col = "black", lwd = 2) +        # Outline high-uncertainty regions
  tm_fill(col = "green", alpha = 0.2)           # Optional: semi-transparent fill for emphasis

# Print the enhanced map
print(unc_map_high)

#===============================================
# Prediction Map + uncertainty (full extent, clear legend/title)
#------------------------------------------------
pred_map <- tm_shape(pred_rast, bbox = st_bbox(soil_sf_utm)) +
  tm_raster(col = "ORGC_pred",
            palette = "viridis",
            style = "pretty",
            n = 10,
            title = "Predicted Soil Organic Carbon (g/kg)",
            legend.show = TRUE) +
  tm_layout(main.title = "Predicted Soil Organic Carbon across Africa",
            main.title.size = 1.5,
            legend.outside = TRUE,
            legend.outside.position = "right",
            frame = TRUE)

# Uncertainty Map (full extent, clear legend/title)
unc_map <- tm_shape(var_rast, bbox = st_bbox(soil_sf_utm)) +
  tm_raster(col = "variance",
            palette = "Reds",
            style = "pretty",
            n = 7,
            title = "Kriging Variance (Uncertainty)",
            legend.show = TRUE) +
  tm_layout(main.title = "Prediction Uncertainty across Africa",
            main.title.size = 1.5,
            legend.outside = TRUE,
            legend.outside.position = "right",
            frame = TRUE)

tmap_arrange(pred_map, unc_map, nrow = 1)


#==================================================
# Phase 4: Interpretation & Application
#==================================================
hotspots <- krige_sf %>% filter(ORGC_pred >= quantile(ORGC_pred, 0.90, na.rm = TRUE))
coldspots <- krige_sf %>% filter(ORGC_pred <= quantile(ORGC_pred, 0.10, na.rm = TRUE))
high_uncert <- krige_sf %>% filter(variance >= quantile(variance, 0.90, na.rm = TRUE))

cat("\n=== Interpretation ===\n")
cat("Hotspots (top 10% ORGC):", nrow(hotspots), "grid cells\n")
cat("Coldspots (bottom 10% ORGC):", nrow(coldspots), "grid cells\n")
cat("High uncertainty areas:", nrow(high_uncert), "grid cells\n")

#-------------------------------
# Hotspot & Coldspot Map
#-------------------------------
hc_map <- tm_shape(pred_rast) +
  tm_raster(col = "ORGC_pred", palette = "viridis", title = "ORGC (g/kg)", style = "quantile") +
  tm_shape(hotspots) + tm_borders(col = "red", lwd = 3) +
  tm_shape(coldspots) + tm_borders(col = "blue", lwd = 3) +
  tm_layout(main.title = "Soil Organic Carbon:Hotspots (red) & Coldspots (blue) of Soil Organic Carbon",
            legend.outside = TRUE)
print(hc_map)

#-----------------------------
# Uncertainty / Sampling Priority Map
#-----------------------------
samp_map <- tm_shape(var_rast) +
  tm_raster(col = "variance", palette = "Reds", title = "Kriging Variance") +
  tm_shape(high_uncert) + tm_borders(col = "black", lwd = 3) +
  tm_layout(main.title = "Priority Areas for Future Soil Sampling (black outline)",
            legend.outside = TRUE)
print(samp_map)


#-----------------------------
# Policy / Management Implications
#-----------------------------
cat("\n=== Policy & Management Implications ===\n")
cat("1. High-ORGC hotspots (red borders) should be conserved to maintain soil carbon stocks.\n")
cat("   Actions: restrict deforestation, promote agroforestry, and soil conservation practices.\n")
cat("2. Low-ORGC coldspots (blue borders) should be targeted for soil carbon improvement.\n")
cat("   Actions: introduce cover crops, organic amendments, reduced tillage, and support local farmers.\n")
cat("3. High-uncertainty zones (black outlines) should be prioritized for additional soil sampling.\n")
cat("   Actions: collect new soil data to improve model accuracy and reduce prediction errors.\n")
cat("4. Land-use planning and agricultural management should consider these maps for sustainable soil carbon management.\n")
cat("5. Policymakers can use these insights to design targeted interventions and carbon sequestration programs.\n")
#======================================
# End the analysis
#======================================

