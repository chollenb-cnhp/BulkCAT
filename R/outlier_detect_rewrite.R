#' BulkCAT Outlier Detection for Multiple Species (sf-based)
#'
#' Identifies potentially erroneous or suspect herbarium / observation records
#' for rare species while preserving multiple valid, disjunct populations.
#' Each species is processed independently.
#'
#' @param data Data frame containing latitude and longitude columns.
#' @param sname Column name for species name (default: "scientificName").
#' @param lat Column name for latitude in decimal degrees (default: "decimalLatitude").
#' @param lon Column name for longitude in decimal degrees (default: "decimalLongitude").
#' @param elev_dem Optional file path to raster DEM. If NULL, elevation checks are skipped.
#' @param eco_shp Optional file path to ecoregion shapefile. If NULL, ecoregion checks are skipped.
#' @param suspicion_thresh Threshold for flagging a point as suspect (default: 0.8).
#' @param weights Named vector of weights for combining indicators
#' (default: c(isolation = 0.5, density = 0.3, ecoregion = 0.1, elevation = 0.2)).
#'
#' @return Data frame with original columns plus:
#' \item{nn_dist}{Nearest neighbor distance.}
#' \item{density_count}{Number of points found in specified density radius.}
#' \item{density_score}{Suspicion probability based on density count (0–1).}
#' \item{D2}{Distance to second nearest neighbor.}
#' \item{isolation_score}{Suspicion probability based on D2 (0-1).}
#' \item{eco_score}{Suspicion probability based on ecoregion (0–1).}
#' \item{combined_score}{Composite suspicion probability (0-1).}
#' \item{LEVEL3}{Level 3 EPA ecoregion ID.}
#' \item{elev}{Elevation in meters.}
#' \item{elev_score}{Suspicion probability based on elevation (0-1).}
#'
#' @import dbscan sf terra scales
#' @export
outlier_detect <- function(data,
                           sname = "scientificName",
                           lat = "decimalLatitude",
                           lon = "decimalLongitude",
                           elev_dem = NULL,      # raster path
                           eco_shp = NULL,       # shapefile path or sf object
                           suspicion_thresh = 0.8,
                           weights = c(isolation = 0.5, density = 0.3, ecoregion = 0.1, elevation = 0.2)) {

  # required packages
  if (!requireNamespace("sf", quietly = TRUE)) stop("Please install sf")
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Please install dbscan")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Please install scales")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Please install terra")

  library(terra)
  library(sf)
  library(dbscan)
  library(scales)

  # -------------------------
  # Input checks and setup
  # -------------------------
  if (!all(c(sname, lat, lon) %in% names(data))) {
    stop("Data must include columns: species name, latitude, and longitude.")
  }

  valid_idx <- which(!is.na(data[[lat]]) & !is.na(data[[lon]]))
  if (length(valid_idx) == 0) stop("No rows with valid coordinates found.")

  sf_all <- st_as_sf(data, coords = c(lon, lat), crs = 4326, remove = FALSE)
  sf_proj <- st_transform(sf_all, crs = 6933)  # Equal Earth (meters)
  coords_proj_all <- st_coordinates(sf_proj)

  # prepare output columns
  data$nn_dist <- NA_real_
  data$density_count <- NA_integer_
  data$density_score <- NA_real_
  data$D2 <- NA_real_
  data$isolation_score <- NA_real_
  data$eco_score <- NA_real_
  data$combined_score <- NA_real_
  data$suspicious <- NA_integer_

  species_list <- unique(data[[sname]])

  # -------------------------
  # Optional ecoregion join
  # -------------------------
  if (!is.null(eco_shp)) {
    if (inherits(eco_shp, "character")) {
      eco_sf <- st_read(eco_shp, quiet = TRUE)
    } else if (inherits(eco_shp, "sf")) {
      eco_sf <- eco_shp
    } else {
      stop("eco_shp must be either a path to a shapefile or an sf object.")
    }

    # ensure both have same CRS
    eco_sf <- st_transform(eco_sf, st_crs(sf_proj))
    sf_all <- suppressWarnings(st_join(sf_proj, eco_sf["LEVEL3"]))
    data$LEVEL3 <- sf_all$LEVEL3
  }

  # -------------------------
  # Optional elevation join
  # -----------------------------------------------

  if (!is.null(elev_dem)) {
    # read raster once outside loop ideally, but keep simple here
    if (!inherits(elev_dem, "SpatRaster")) {
      elev_ras <- terra::rast(elev_dem)
      print("loaded")
    } else {
      elev_ras <- elev_dem
      print("noLoad")
    }
    # extract elevations for projected points (transform to DEM CRS)
    pts_ras <- st_transform(sf_all, crs(elev_ras))
    elev_vals <- terra::extract(elev_ras, vect(pts_ras))[, 2]
    print(elev_vals)
    data$elev <- elev_vals

  }


  # -------------------------
  # loop per species
  # -------------------------
  for (sp in species_list) {
    idx_sp <- which(data[[sname]] == sp)
    coords_sp <- coords_proj_all[idx_sp, , drop = FALSE]
    n_sp <- nrow(coords_sp)

    # only analyze species with more than 2 points
    if (n_sp < 3) next

    dist_mat <- as.matrix(dist(coords_sp))
    diag(dist_mat) <- Inf
    nn_dist <- apply(dist_mat, 1, min, na.rm = TRUE)

    q1_nn <- quantile(nn_dist, 0.25, na.rm = TRUE)
    q3_nn <- quantile(nn_dist, 0.75, na.rm = TRUE)
    iqr_nn <- q3_nn - q1_nn
    max_nn <- max(nn_dist, na.rm = TRUE)

    # calculate clustering radius for density based on weighted combination of maxNN distance and IQR outliers
    r <- min((q3_nn + 1.5 * iqr_nn) * 0.3 + max_nn * 0.7, max_nn)
    if (!is.finite(r) || r <= 0) {
      r <- median(nn_dist, na.rm = TRUE)
      if (!is.finite(r) || r <= 0) r <- 1000
    }

    # calculate number of points within density radius
    density_counts <- apply(dist_mat, 1, function(row) sum(row <= r, na.rm = TRUE))
    max_density <- max(density_counts)
    density_score_sp <- 1 - density_counts / max_density

    # --- isolation (D2)
    D2 <- apply(dist_mat, 1, function(row) {
      vals <- sort(row)
      fin <- vals[is.finite(vals)]
      if (length(fin) >= 2) fin[2] else fin[1]
    })
    Q1_D2 <- quantile(D2, 0.25, na.rm = TRUE)
    Q3_D2 <- quantile(D2, 0.75, na.rm = TRUE)
    IQR_D2 <- Q3_D2 - Q1_D2
    denom_D2 <- Q3_D2 * 0.3 + max_nn * 0.7 - Q1_D2
    if (!is.finite(denom_D2) || denom_D2 <= 0) denom_D2 <- max(Q3_D2, 1e-9)
    isolation_score_sp <- (D2 - Q1_D2) / denom_D2
    isolation_score_sp[is.na(isolation_score_sp)] <- 0
    isolation_score_sp[isolation_score_sp < 0] <- 0
    isolation_score_sp[isolation_score_sp > 1] <- 1

    # --- optional ecoregion suspicion
    eco_score_sp <- rep(NA_real_, n_sp)

    if (!is.null(eco_shp) && "LEVEL3" %in% names(data)) {
      eco_subset <- data$LEVEL3[idx_sp]
      valid <- !is.na(eco_subset)

      if (any(valid)) {
        eco_tab <- table(eco_subset[valid])
        eco_names <- names(eco_tab)
        total_sp <- sum(eco_tab)

        counts <- eco_tab[match(eco_subset[valid], eco_names)]
        freq <- counts / total_sp

        # log-based scaling
        # transform freq (0 < freq <= 1) into 0..1 score
        # score is scaled such that 50-50 ecoregion split is 0 (no suspicion)
        score <- 1 - 2*freq
        score <- pmin(pmax(score, 0), 1)  # constrain to 0..1

        eco_score_sp[valid] <- score
      } else {
        eco_score_sp[] <- 0
      }
    } else {
      eco_score_sp[] <- 0
    }

    # --- optional elevation suspicion ---
    elev_score_sp <- rep(0, n_sp)
    if (!is.null(elev_dem) && "elev" %in% names(data)) {
      elev_subset <- data$elev[idx_sp]
      # valid elevation values
      valid <- !is.na(elev_subset)
    if (any(valid)) {
      elev_q1 <- quantile(elev_subset[valid], 0.25, na.rm = TRUE)
      elev_q3 <- quantile(elev_subset[valid], 0.75, na.rm = TRUE)
      elev_iqr <- elev_q3 - elev_q1
      lower <- elev_q1 - 1.5 * elev_iqr
      upper <- elev_q3 + 1.5 * elev_iqr

      # distance from range, scaled 0–1
      elev_score_sp[valid] <- 0
      elev_score_sp[elev_subset < lower] <- pmin(1, (lower - elev_subset[elev_subset < lower]) / elev_iqr)
      elev_score_sp[elev_subset > upper] <- pmin(1, (elev_subset[elev_subset > upper] - upper) / elev_iqr)
    }
    }

    # --- combined weighted score
    combined_sp <- (
      weights["isolation"] * isolation_score_sp +
        weights["density"]   * density_score_sp +
        weights["ecoregion"] * eco_score_sp +
        weights["elevation"] * elev_score_sp
    )
    suspicious_sp <- as.integer(combined_sp >= suspicion_thresh)

    # --- write back
    data$nn_dist[idx_sp] <- nn_dist
    data$density_count[idx_sp] <- density_counts
    data$density_score[idx_sp] <- density_score_sp
    data$D2[idx_sp] <- D2
    data$isolation_score[idx_sp] <- isolation_score_sp
    data$eco_score[idx_sp] <- eco_score_sp
    data$elev_score[idx_sp] <- elev_score_sp
    data$combined_score[idx_sp] <- combined_sp
    data$suspicious[idx_sp] <- suspicious_sp
  }

  return(data)
}
