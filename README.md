🌍 Spatial Data Project: Africa Soil Organic Carbon (ORGC) Analysis
Overview

This project performs a comprehensive spatial analysis of African soils using the WoSIS 2023 soil dataset. It is structured, complete, and thesis-ready, covering data cleaning, preprocessing, spatial analysis, and geostatistical prediction.

The analysis focuses on Soil Organic Carbon (ORGC) and includes both exploratory and predictive workflows, providing maps, uncertainty assessments, and policy insights.

Features

Data Cleaning & Imputation: Handles missing values, removes unreliable variables, and imputes using median/mean.

Spatial Preparation: Handles duplicate locations, averages soil properties by site, and converts data to spatial (sf) objects.

Exploratory Analysis (ESDA):

Property distributions

Depth-based and latitudinal summaries

Clustering of soil sites (k-means)

Trend visualization & detrending

Spatial Dependence Modeling:

Experimental variogram

Variogram model fitting (Spherical, Exponential, Gaussian)

Moran’s I for spatial autocorrelation

Interpolation & Prediction:

Ordinary Kriging for ORGC prediction

50 km resolution prediction grid

Rasterization for visualization

Uncertainty & Hotspot Analysis:

Maps of prediction uncertainty

Identification of ORGC hotspots & coldspots

Prioritization for future soil sampling

Policy & Management Insights:

Guidance for soil carbon conservation, improvement, and sustainable land management.

Requirements

R (≥4.0)

R Packages:

dplyr, sf, ggplot2, moments, corrplot, cluster, FactoMineR, factoextra, tmap, gstat, sp, spdep, terra, patchwork, tidyr

Usage

Clone the repository:

git clone https://github.com/haile6/spatial-data-project.git


Open the final.R script in RStudio.

Set the path to the WoSIS dataset:

soil_data <- read.csv("wosis_2023_full_combined_africa_only.csv")


Run the script section by section to reproduce:

Data cleaning & preprocessing

Exploratory spatial data analysis

Variogram modeling & Kriging prediction

Visualization of maps and uncertainty

Outputs

Predicted Soil Organic Carbon Map

Kriging Uncertainty Map

Hotspot & Coldspot maps

Policy & management recommendations
