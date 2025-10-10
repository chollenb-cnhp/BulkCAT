#' Run BulkCAT analysis
#'
#' This function completes conservation status assessments (rarity ranks) using NatureServe
#' methodology. BulkCAT requires a user-curated multi-species point observations dataset with lat/lon columns.
#'
#' @param input_df Input multispecies observation points dataframe
#' @param sname Species name column. Default = "scientificName"
#' @param lat Latitude column name. Default = "decimalLatitude"
#' @param lon Longitude column name. Default = "decimalLongitude"
#' @param eo_separation Minimum separation distance (m) for unique EO clusters. Default = 1000m.
#' @param grid_size Side length (m) for AOO grid cells. Default = 2000m.
#'
#' @return Dataframe with calculated rarity metrics for each species/element.
#' @export
#'
run_bulkCAT <- function(input_df,
                        sname = "scientificName",
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

  # remove genus-only records (keep names with ≥ 2 parts)
  df <- df[sapply(strsplit(df[[sname]], "\\s+"), length) >= 2, ]

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
  gdf <- sf::st_as_sf(df, coords = c(lon, lat), crs = wgs_crs)
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
    } else {
      # trinomial: check if this is the only infraspecific taxon for the binomial
      species_parts <- strsplit(species, "\\s+")[[1]]
      binomial <- paste(species_parts[1:2], collapse = " ")

      # find all taxa in list starting with the binomial
      related_taxa <- grep(paste0("^", binomial, "\\b"), species_list, value = TRUE)
      if (length(related_taxa) == 2) {
        # only trinomial and binomial exists -> copy binomial results
        print("Only one infraspecies found, matching binomial results...")
        binomial_result <- results[nrow(results), ]
        results <- rbind(results, data.frame(
          species = species,
          num_obs = binomial_result$num_obs,
          eoo_area_km2 = binomial_result$eoo_area_km2,
          aoo_num_cells = binomial_result$aoo_num_cells,
          num_EOs = binomial_result$num_EOs,
          stringsAsFactors = FALSE
        ))
          next  # skip recalculation for this trinomial
        }
      else {
        # multiple infraspecific taxa -> only use exact match
        species_subset <- gdf_proj[gdf_proj[[sname]] == species, ]
      }}


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
    results <- rbind(results, data.frame(
      species = species,
      num_obs = num_obs,
      eoo_area_km2 = round(eoo_area_km2, 2),
      aoo_num_cells = aoo_cells,
      num_EOs = num_clusters,
      stringsAsFactors = FALSE
    ))
  }

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

  results$Points <- (score_eoo(results$eoo_area_km2) +
                             2 * score_aoo(results$aoo_num_cells) +
                             score_num(results$num_EOs)) / 4

  results$SRank <- score_rank(results$Points)

  return(results)
}
#' Deduplicate Records by Specified Columns
#'
#' Removes duplicate records from a data frame based on specified columns,
#' giving priority to institutions with more records.
#'
#' @param input_df A data frame containing the records to deduplicate.
#' @param cols A character vector of column names to check for duplicates.
#'   Defaults to c("recordedBy", "recordNumber", "scientificName", "eventDate").
#' @param institution_col A string specifying the column name used to prioritize
#'   records by institution count. Defaults to "institutionCode".
#'
#' @return A data frame with duplicates removed. Prints the number of records removed.
#' @export
#'
deduplicate <- function(input_df, cols = c("recordedBy", "recordNumber", "scientificName", "eventDate"), institution_col = "institutionCode") {
  # Step 1: count records per institution
  inst_counts <- table(input_df[[institution_col]])

  # Step 2: add counts back as a helper column
  input_df$inst_count <- inst_counts[input_df[[institution_col]]]

  # Step 3: order rows by inst_count (descending), then by original order
  ord <- order(-input_df$inst_count, seq_len(nrow(input_df)))
  input_df <- input_df[ord, ]

  # Step 4: drop duplicates based on user-defined columns
  dedup_df <- input_df[!duplicated(input_df[cols]), ]

  # Remove helper column
  dedup_df$inst_count <- NULL
  duplicates <- nrow(input_df) - nrow(dedup_df)
  cat("Observations removed via deduplication:", duplicates, "\n")

  return(dedup_df)
}
