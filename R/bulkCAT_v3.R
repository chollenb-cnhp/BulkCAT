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
#' @param threats Boolean if you would like threat options in the output. calc_threats() can also be run on the output later.
#' @param community Boolean if the input contains plant association (or non-scientific) names. Default = FALSE.
#' @param poly_layer Shapefile path or sf object for polygon layer to be used for calculating AOO. Default = NULL.
#' @return Dataframe with calculated rarity metrics for each species/element.
#' @export
#'
run_bulkCAT <- function(input_df,
                        sname = "scientificName",
                        lat = "decimalLatitude",
                        lon = "decimalLongitude",
                        eo_separation = 1000,
                        grid_size = 2000,
                        threats = FALSE,
                        community = FALSE,
                        poly_layer = NULL) {

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
  if (!is.null(poly_layer)) {
    if (inherits(poly_layer, "character")) {
      poly_sf <- st_read(poly_layer, quiet = TRUE)
    } else if (inherits(poly_layer, "sf")) {
      poly_sf <- poly_layer
    } else {
      stop("poly_layer must be either a path to a shapefile or an sf object.")
    }

    # ensure it's in World Equal Area
    poly_sf <- st_transform(poly_sf, crs = 6933)}
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
    if (community == FALSE){
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
      }}
    # if processing a plant community dataset
    else
    {
      species_subset <- gdf_proj[gdf_proj[[sname]] == species, ]
      num_obs <- nrow(species_subset)     # count number of observations
      if (num_obs == 0) {
        warning("No observations for species: ", species, "; skipping.")
        next
      }
    }

    ###### EOO: Convex Hull Area ######
    hull <- sf::st_convex_hull(sf::st_union(species_subset))
    eoo_area_km2 <- as.numeric(sf::st_area(hull)) / 1e6

    ###### AOO: 2x2 km Grid ######
    coords <- sf::st_coordinates(species_subset)
    if (!is.null(poly_layer)) {
      poly_subset <- poly_sf[poly_sf[[sname]]==species, ]

      # Fixed global origin (IUCN-style)
      origin <- c(0, 0)

      # Get bounding box
      bb <- sf::st_bbox(poly_subset)

      # Snap bbox to grid
      xmin <- floor((as.numeric(bb["xmin"]) - origin[1]) / grid_size) * grid_size + origin[1]
      ymin <- floor((as.numeric(bb["ymin"]) - origin[2]) / grid_size) * grid_size + origin[2]
      xmax <- ceiling((as.numeric(bb["xmax"]) - origin[1]) / grid_size) * grid_size + origin[1]
      ymax <- ceiling((as.numeric(bb["ymax"]) - origin[2]) / grid_size) * grid_size + origin[2]

      # Create snapped bbox polygon
      snapped_bbox <- sf::st_as_sfc(sf::st_bbox(
        c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
        crs = sf::st_crs(poly_subset)
      ))
      # Build grid
      grid <- sf::st_make_grid(
        snapped_bbox,
        cellsize = grid_size,
        square = TRUE
      )
      grid_sf <- sf::st_sf(geometry = grid)

      # Intersect grid with species geometry
      hits <- sf::st_intersects(grid_sf, poly_subset, sparse = FALSE)
      # Count occupied cells
      aoo_cells <- sum(apply(hits, 1, any))
    } else{

      bottomleftpoints <- floor(coords / grid_size)
      uniquecells <- unique(bottomleftpoints)
      aoo_cells <- nrow(uniquecells)
    }

    ###### EO Cluster Count (1 km buffer) ######
    clustering <- dbscan::dbscan(coords, eps = eo_separation, minPts = 1)
    num_clusters <- max(clustering$cluster)

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
    EOOVal = c(100, 250, 1000, 5000, 20000, 200000, 2500000, 1000000000000000, NA),
    EOOScore = c(0, 0.79, 1.57, 2.36, 3.14, 3.93, 4.71, 5.5, NA),
    AOOVal = c(1, 2, 5, 25, 125, 500, 2500, 12500, 1000000000000000),
    AOOScore = c(0, 0.69, 1.38, 2.06, 2.75, 3.44, 4.13, 4.81, 5.5),
    NumVal = c(5, 20, 80, 300, 1200, 1000000000000000, NA, NA, NA),
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
  if (threats){
    results <- calc_threats(results)
  }

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

#' Show the effects of different potential threat options on calculated ranks.
#'
#' Includes low, medium, and high threat options, calculating SRank_lowT, SRank_medT, SRank_highT
#' to assist with review of rarity-based ranks.
#'
#' @param input_df A data frame containing the records to deduplicate {usually an output of runBulkCAT()}
#' @param points_col A string specifying the column name which contains rarity-based points.
#'   Defaults to "Points".
#'
#' @return A data frame calculated SRanks and Points for different threat options.
#' @export
#'
calc_threats <- function(input_df, points_col = "Points") {
  # create rules_df based on NS methodology (same as in run_BulkCAT())
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
    if (length(idx)==0) return(0)
    rules[[score_col]][idx]
  }

  score_rank <- function(x) vapply(x, assign_points, character(1), rules = rules_df, metric = "Rank")


  # threat scenarios
  threats <- data.frame(
    suffix = c("_HighT", "_MedT", "_LowT"),
    value  = c(1.83, 3.67, 5.5),
    stringsAsFactors = FALSE
  )

  # apply threats
  for (i in seq_len(nrow(threats))) {

    threat <- threats$suffix[i]
    value  <- threats$value[i]

    points_col_new <- paste0(points_col, threat)
    srank_col_new  <- paste0("SRank", threat)

    input_df[[points_col_new]] <- input_df[[points_col]] * 0.7 + value * 0.3
    input_df[[srank_col_new]]  <- score_rank(input_df[[points_col_new]])
  }

  return(input_df)
}
