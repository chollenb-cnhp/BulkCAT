#' Run BulkCAT analysis
#'
#' This function completes conservation status assessments (rarity ranks) using NatureServe
#' methodology. BulkCAT requires a user-curated multi-species point observations dataset with lat/lon columns.
#'
#' @param input_df Input multispecies observation points dataframe
#' @param sname Species name column. Default = "acceptedName"
#' @param lat Latitude column name. Default = "decimalLatitude"
#' @param lon Longitude column name. Default = "decimalLongitude"
#' @param eo_separation Minimum separation distance (m) for unique EO clusters. Default = 1000m.
#' @param grid_size Side length (m) for AOO grid cells. Default = 2000m.
#'
#' @return Dataframe with calculated rarity metrics for each species/element.
#' @export
run_bulkCAT <- function(input_df,
                        sname = "acceptedName",
                        lat = "decimalLatitude",
                        lon = "decimalLongitude",
                        eo_separation = 1000,
                        grid_size = 2000) {

  # ----------------------------------------------------------------
  # ----- Input validation -----
  # ----------------------------------------------------------------
  if (!is.data.frame(input_df)) stop("input_df must be a data frame.")

  required_cols <- c(sname, lat, lon)
  missing_cols <- setdiff(required_cols, names(input_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(input_df[[lat]]) || !is.numeric(input_df[[lon]])) {
    stop("Latitude and longitude columns must be numeric.")
  }

  df <- input_df[(!is.na(input_df[[lat]]) &
                   !is.na(input_df[[lon]]) &
                   !is.na(input_df[[sname]])),]

  if (nrow(df) == 0) stop("No rows remaining after filtering missing lat/lon/species.")

  if (all(is.na(input_df[[sname]]))) {
    stop("Species column has only missing values.")
  }

  # ----------------------------------------------------------------
  # ----- Prepare coordinate reference systems -----
  # ----------------------------------------------------------------
  equal_area_crs <- 6933  # World Cylindrical Equal Area projection
  wgs_crs <- 4326         # WGS84

  # ----------------------------------------------------------------
  # ----- Convert to sf object and transform to equal area -----
  # ----------------------------------------------------------------
  gdf <- sf::st_as_sf(df, coords = c(lat, lon), crs = wgs_crs)
  gdf_proj <- sf::st_transform(gdf, crs = equal_area_crs)

  # ----------------------------------------------------------------
  # ----- Identify species list -----
  # ----------------------------------------------------------------
  species_list <- sort(unique(gdf_proj[[sname]]))

  # create list object to store results
  results <- list()

  # ----------------------------------------------------------------
  # ----- Species loop -----
  # ----------------------------------------------------------------
  for (species in species_list) {
    cat("Processing:", species, "\n")

    # Match all descendant varieties/subspecies if just genus + species
    if (length(strsplit(species, "\\s+")[[1]]) == 2) {
      species_subset <- gdf_proj[startsWith(gdf_proj[[sname]], species), ]
    }
    else {
      # if trinomial name, only rank using the infraspecific observations
      species_subset <- gdf_proj[gdf_proj[[sname]] == species, ]
    }

    num_obs <- nrow(species_subset)     # count number of observations
    if (num_obs == 0) {
      warning("No observations for species: ", species, "; skipping.")
      next
    }

    ###### EOO: Convex Hull Area ######
    hull <- sf::st_convex_hull(sf::st_union(species_subset))
    eoo_area_km2 <- as.numeric(sf::st_area(hull)) / 1e6

    ###### AOO: 2x2 km Grid ######
    coords <- sf::st_coordinates(species_subset)
    bottomleftpoints <- floor(coords / grid_size)
    uniquecells <- unique(bottomleftpoints)
    aoo_cells <- nrow(uniquecells)

    ###### EO Cluster Count (1 km buffer) ######
    clustering <- dbscan::dbscan(coords, eps = eo_separation, minPts = 1)
    num_clusters <- length(unique(clustering$cluster))

    ###### Append to results ######
    results[[length(results) + 1]] <- data.frame(
      species = species,
      num_obs = num_obs,
      eoo_area_km2 = round(eoo_area_km2, 2),
      aoo_num_cells = aoo_cells,
      num_EOs = num_clusters,
      stringsAsFactors = FALSE
    )
  }

  final_results <- do.call(rbind, results)

  # ----------------------------------------------------------------
  # ----- RANKING ROLL‑UP -----
  # ----------------------------------------------------------------
  rules_df <- data.frame(
    EOOVal = c(100, 250, 1000, 5000, 20000, 200000, 2500000, 25000000, NA),
    EOOScore = c(0, 0.79, 1.57, 2.36, 3.14, 3.93, 4.71, 5.5, NA),
    AOOVal = c(1, 2, 5, 20, 125, 500, 5000, 50000, 10000000),
    AOOScore = c(0, 0.69, 1.38, 2.06, 2.75, 3.44, 4.13, 4.81, 5.5),
    NumVal = c(5, 20, 80, 300, 1200, 1000000, NA, NA, NA),
    NumScore = c(0, 1.38, 2.75, 4.13, 5.5, 5.5, NA, NA, NA),
    RankVal = c(1.5, 2.5, 3.5, 4.5, 6, NA, NA, NA, NA),
    RankScore = c("S1", "S2", "S3", "S4", "S5", NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )

  assign_points <- function(value, rules, metric) {
    val_col <- paste0(metric, "Val")
    score_col <- paste0(metric, "Score")
    idx <- which(value <= rules[[val_col]])[1]
    if (is.na(idx)) return(0)
    rules[[score_col]][idx]
  }

  score_eoo <- function(x) vapply(x, assign_points, numeric(1), rules = rules_df, metric = "EOO")
  score_aoo <- function(x) vapply(x, assign_points, numeric(1), rules = rules_df, metric = "AOO")
  score_num <- function(x) vapply(x, assign_points, numeric(1), rules = rules_df, metric = "Num")
  score_rank <- function(x) vapply(x, assign_points, character(1), rules = rules_df, metric = "Rank")

  final_results$Points <- (score_eoo(final_results$eoo_area_km2) +
                             2 * score_aoo(final_results$aoo_num_cells) +
                             score_num(final_results$num_EOs)) / 4

  final_results$SRank <- score_rank(final_results$Points)

  return(final_results)
}
