#' Translate plant names to accepted names using multiple synonym sources
#'
#' This function standardizes scientific names by translating synonyms to
#' accepted names using a combination of:
#' \itemize{
#'   \item User-supplied lookup tables (LUTs)
#'   \item Built-in synonym tables (NatureServe, SEINet, USDA, WCVP)
#'   \item Optional fuzzy and exact matching via the World Checklist of Vascular Plants (rWCVP)
#' }
#'
#' The function proceeds in priority order:
#' \enumerate{
#'   \item Direct matches to the checklist
#'   \item User-supplied LUTs (in order)
#'   \item Built-in LUTs (in order given by \code{synonym_sources})
#'   \item Optional WCVP exact matches
#'   \item Optional WCVP fuzzy matches
#'   \item Optional binomial-only pass (if \code{ssp_mods = TRUE})
#' }
#'
#' A new column \code{acceptedName} is added to the data, along with
#' \code{translationSource} indicating how the name was translated.
#'
#' @param input_df A data frame containing taxon names to be standardized.
#' @param name_col Name of the column in \code{input_df} containing scientific names.
#' @param checklist Optional data frame of accepted names. If \code{NULL} or missing,
#'   the built-in NatureServe Network tracheophyte checklist is used. This may contain names which are not correct in your region.
#' @param checklist_name_col Column in \code{checklist} containing accepted names.
#' @param synonym_LUTs Optional list of user-supplied data frames with columns
#'   \code{inputName} and \code{outputName}.
#' @param synonym_sources Character vector giving the order of built-in synonym sources
#'   to apply. Any of \code{"NatureServe"}, \code{"SEINet"}, \code{"USDA"}, \code{"WCVP"}.
#' @param synonym_sources_rerun Logical. If TRUE, reload synonym sources (ex: with a new target checklist).
#' @param fuzzy Logical. If TRUE, attempt exact and fuzzy matching using \pkg{rWCVP}.
#' @param wcvp_rerun Logical. If TRUE, ignored cached \pkg{rWCVP} results and rerun. Note that you can pass saved LUTs from previous runs if desired in synonym_LUTs.
#' @param ssp_mods Logical. If TRUE, perform a second pass using only binomial names
#'   (genus + species) for unresolved infraspecific names.
#'
#' @details
#' The function never overwrites an already-resolved \code{acceptedName}.
#' Earlier sources in the priority chain take precedence over later ones.
#'
#' Built-in synonym tables are shipped with the package and loaded from
#' \code{inst/extdata}.
#'
#' @return The input data frame with two additional columns:
#' \itemize{
#'   \item \code{acceptedName} – the resolved accepted name (or NA if unresolved)
#'   \item \code{translationSource} – source used for the translation
#' }
#'
#' @seealso \code{\link[rWCVP]{wcvp_match_names}}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(scientificName = c("Pinus murrayana", "Quercus alba"))
#' synonymize(df)
#' }
#'
#' @export


synonymize <- function(input_df,
                       name_col = "scientificName",
                       checklist = NA,
                       checklist_name_col = "SNAME",
                       synonym_LUTs = list(),
                       synonym_sources = c("NatureServe", "SEINet", "USDA", "WCVP"),
                       synonym_sources_rerun = FALSE,
                       fuzzy = TRUE,
                       wcvp_rerun = FALSE,
                       ssp_mods = FALSE) {

  # -----------------------------------------
  # Load built in synonym sources
  # -----------------------------------------
  get_builtin_LUTs <- function(sources, synonym_sources_rerun = FALSE) {

    # Make sure cache exists
    if (!exists(".bulkcat_cache", envir = .GlobalEnv)) {
      print("Creating .bulkcat_cache")
      assign(".bulkcat_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
    }

    # Only load LUTs if not already in cache
    if (!exists("builtin_LUTs", envir = .bulkcat_cache, inherits = FALSE) | synonym_sources_rerun) {
      LUTs <- list()
      for (source_name in sources) {
        lut_path <- switch(source_name,
                           "NatureServe" = system.file("extdata", "NatureServe.csv", package = "BulkCAT"),
                           "SEINet"      = system.file("extdata", "SEINet.csv", package = "BulkCAT"),
                           "USDA"        = system.file("extdata", "USDA.csv", package = "BulkCAT"),
                           "WCVP"        = system.file("extdata", "WCVP.csv", package = "BulkCAT"),
                           NULL
        )
        if (!is.null(lut_path) && nzchar(lut_path)) {
          cat("Loading built in lookup table from: ", source_name, "\n")
          LUTs[[source_name]] <- utils::read.csv(lut_path, stringsAsFactors = FALSE)
          # Filter out trivial rows and keep only relevant names
          LUTs[[source_name]] <- LUTs[[source_name]] %>%
            dplyr::filter(
              inputName != outputName,
              inputName %in% checklist[[checklist_name_col]] |
                outputName %in% checklist[[checklist_name_col]]
            )
          #  Swap input/output if inputName is in checklist
          LUTs[[source_name]] <- LUTs[[source_name]] %>%
            dplyr::mutate(
              swap_flag = inputName %in% checklist[[checklist_name_col]],
              input_orig = inputName,          # store original inputName
              inputName  = dplyr::if_else(swap_flag, outputName, inputName),
              outputName = dplyr::if_else(swap_flag, input_orig, outputName)
            ) %>%
            dplyr::select(-swap_flag, -input_orig)
          write.csv(LUTs[[source_name]], paste0("filtered_built_ins_", source_name, ".csv"), row.names=FALSE)

        }
      }
      assign("builtin_LUTs", LUTs, envir = .bulkcat_cache)
      cat("Saved builtin_LUTs to cache!\n")
    } else {
      cat("Using cached builtin_LUTs\n")
    }

    return(.bulkcat_cache$builtin_LUTs)
  }
  builtin_LUTs <- get_builtin_LUTs(synonym_sources, synonym_sources_rerun)
  # ------------------------------------------
  # 1️⃣ Initialize acceptedName and translationSource
  # ------------------------------------------
  input_df <- dplyr::mutate(input_df,
                            acceptedName = NA_character_,
                            translationSource = NA_character_)

  if (missing(checklist) || is.null(checklist)) {
    message("No checklist supplied. Using NatureServe Network Tracheophyta checklist...")
    checklist <- utils::read.csv(system.file("extdata", "NatureServe.csv", package = "BulkCAT"), stringsAsFactors = FALSE)
    checklist_name_col <- "outputName"
  }

  # ------------------------------------------
  # 2️⃣ Direct matches to checklist
  # ------------------------------------------
    print("Finding direct matches to checklist...")
    direct_idx <- input_df[[name_col]] %in% checklist[[checklist_name_col]]
    input_df$acceptedName[direct_idx] <- input_df[[name_col]][direct_idx]
    input_df$translationSource[direct_idx] <- "Direct"


  # ------------------------------------------
  # 3️⃣ Function to process one LUT
  # ------------------------------------------
    process_LUT <- function(LUT, df, checklist, source_name = "LUT", name_col = "scientificName", ssp=FALSE) {

        # Filter out trivial rows and keep only relevant names
        LUT <- LUT %>%
          dplyr::filter(
            inputName != outputName,
            inputName %in% checklist[[checklist_name_col]] |
              outputName %in% checklist[[checklist_name_col]]
          )

        # 2️⃣ Swap input/output if inputName is in checklist
        LUT <- LUT %>%
          dplyr::mutate(
            swap_flag = inputName %in% checklist[[checklist_name_col]],
            input_orig = inputName,          # store original inputName
            inputName  = dplyr::if_else(swap_flag, outputName, inputName),
            outputName = dplyr::if_else(swap_flag, input_orig, outputName)
          ) %>%
          dplyr::select(-swap_flag, -input_orig)

        if (source_name == "NatureServe"){
          write.csv(LUT, "NatureServe_test.csv", row.names=FALSE)
        }

        # 3️⃣ Keep unique inputName, prioritizing outputName in checklist
        LUT <- LUT %>%
          dplyr::mutate(in_checklist = outputName %in% checklist[[checklist_name_col]]) %>%
          dplyr::group_by(inputName) %>%
          dplyr::slice_max(in_checklist, n = 1, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::select(-in_checklist)

        # 4️⃣ Join to main df
        if (ssp)
        {
          df <- df %>%
            dplyr::left_join(LUT, by = setNames("inputName", "binomialName")) %>%
            dplyr::mutate(
              acceptedName = ifelse(is.na(acceptedName) & !is.na(outputName),
                                    outputName, acceptedName),
              translationSource = ifelse(is.na(translationSource) & !is.na(outputName),
                                         source_name, translationSource)
            ) %>%
            dplyr::select(-outputName)
        }
        else{
        df <- df %>%
          dplyr::left_join(LUT, by = setNames("inputName", name_col)) %>%
          dplyr::mutate(
            acceptedName = ifelse(is.na(acceptedName) & !is.na(outputName),
                                  outputName, acceptedName),
            translationSource = ifelse(is.na(translationSource) & !is.na(outputName),
                                       source_name, translationSource)
          ) %>%
          dplyr::select(-outputName)
        }



      return(df)
    }

  # ------------------------------------------
  # 4️⃣ Run synonyms function
  # ------------------------------------------
  run_synonyms <- function(input_df, name_col) {

    # 4a apply user-provided LUTs in order
    print("Applying user supplied LUTs...")
    if (length(synonym_LUTs) > 0) {
      for (i in seq_along(synonym_LUTs)) {
        LUT <- synonym_LUTs[[i]]
        source_name <- ifelse(i <= length(synonym_sources),
                              synonym_sources[i],
                              paste0("LUT", i))
        input_df <- process_LUT(LUT, input_df, checklist, source_name, name_col)
      }
    }

    # 4b apply built-in LUTs (already loaded)
    print("Applying built-in LUTs...")

    if (length(builtin_LUTs) > 0) {
      for (source_name in names(builtin_LUTs)) {
        print(paste0("Applying built-in LUTs for ", source_name))
        input_df <- process_LUT(
          builtin_LUTs[[source_name]],
          input_df,
          checklist,
          source_name,
          name_col
        )
      }
    }

    return(input_df)
  }

  # Run synonyms on full names
  print("############# Running synonyms...")
  input_df <- run_synonyms(input_df, name_col = name_col)

  # ------------------------------------------
  # 5️⃣ Repeat the process for binomial names only (ssp_mods)
  # ------------------------------------------
  if (ssp_mods) {
    input_df <- input_df %>%
      dplyr::mutate(binomialName = stringr::word(.data[[name_col]], 1, 2))

    cat("Finding direct matches for binomial names...\n")

    # Direct updates for binomial matches
    direct_idx <- is.na(input_df$acceptedName) &
      input_df$binomialName %in% checklist[[checklist_name_col]]

    input_df$acceptedName[direct_idx] <- input_df$binomialName[direct_idx]
    input_df$translationSource[direct_idx] <- "Direct binomial"

    # create subset of those still needing matching
    subset_df <- dplyr::filter(input_df, is.na(acceptedName))

    print("############# Running binomial synonyms...")
    subset_df <- run_synonyms(subset_df, name_col = "binomialName")
    subset_LUT <- subset_df %>%
      dplyr::filter(!is.na(acceptedName)) %>%
      dplyr::mutate(
        inputName = .data[[name_col]],
        outputName = acceptedName
      ) %>%
      dplyr::select(inputName, outputName) %>%
      dplyr::distinct()

    input_df <- process_LUT(subset_LUT, input_df, checklist, source_name = "binomial", name_col = name_col)
  }

###################################
  # --------------------------------------------
    #fuzzy matching using rWCVP
    # ------------------------------------------
    wcvp_run <- function(input_df, name_col, rerun = FALSE, is_binomial = FALSE){
    dir.create("rWCVP", showWarnings = FALSE)
    if (!(fuzzy && requireNamespace("rWCVP", quietly = TRUE))) {
      if (!requireNamespace("rWCVP", quietly = TRUE)) {
        warning("rWCVP not installed; skipping WCVP translation.")
      }
      return(input_df)
    }

    # ---------------------------------------
    # Decide cache keys
    # ---------------------------------------

    cache_tag <- if (is_binomial) "binomial" else "full"

    exact_key <- paste0("rWCVP/wcvp_exact_", cache_tag, "_", name_col)
    fuzzy_key <- paste0("rWCVP/wcvp_fuzzy_", cache_tag, "_", name_col)

    if (is_binomial) {
      unique_df <- input_df %>%
        dplyr::filter(is.na(acceptedName)) %>%
        dplyr::distinct(binomialName) %>%
        dplyr::rename(query = binomialName)
    } else {
      unique_df <- input_df %>%
        dplyr::filter(is.na(acceptedName)) %>%
        dplyr::distinct(.data[[name_col]]) %>%
        dplyr::rename(query = .data[[name_col]])
    }

    if (nrow(unique_df) == 0) {
      return(input_df)
    }

    # ---------------------------------------
    # Try cache for fuzzy results
    # ---------------------------------------

    cache_file <- paste0(fuzzy_key, ".rds")

    if (file.exists(cache_file) && !rerun) {
      wcvp_fuzzy_LUT <- readRDS(cache_file)
      cat("Found cached rWCVP results at: ", cache_file, "loading stored results.", "\n")

    } else {
      cat("Cached rWCVP results not found or rerun selected, running rWCVP for: ", fuzzy_key)

      wcvp_matches <- rWCVP::wcvp_match_names(
        unique_df,
        name_col = "query",
        fuzzy = TRUE
      )
      # ---------------------------------------
      # Build FUZZY LUT
      # ---------------------------------------
      # Assume builtin_LUTs[["WCVP"]] exists and is preloaded
      wcvp_fuzzy_LUT <- wcvp_matches %>%
        dplyr::left_join(
          builtin_LUTs[["WCVP"]],
          by = c("wcvp_accepted_id" = "accepted_plant_name_id")) %>%
          dplyr::mutate(
          inputName = query,
          outputName = outputName
        ) %>%
        dplyr::select(inputName, outputName) %>%
        dplyr::filter(inputName != outputName)

      # 2️⃣ Ensure one output per inputName, preferring checklist matches
      wcvp_fuzzy_LUT <- wcvp_fuzzy_LUT %>%
        dplyr::filter(!is.na(outputName)) %>%
        dplyr::group_by(inputName) %>%
        dplyr::arrange(dplyr::desc(outputName %in% checklist[[checklist_name_col]])) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

      # ---------------------------------------
      # Save to cache
      # ---------------------------------------

      saveRDS(wcvp_fuzzy_LUT, cache_file)

    }

      # ---------------------------------------
      # Also write CSVs (for review)
      # ---------------------------------------
      if (is_binomial) {
        utils::write.csv(wcvp_fuzzy_LUT, "rWCVP/rWCVP_fuzzy_binomial.csv", row.names = FALSE)
      } else {
        utils::write.csv(wcvp_fuzzy_LUT, "rWCVP/rWCVP_fuzzy.csv", row.names = FALSE)
      }

    # ---------------------------------------
    # Apply LUTs
    # ---------------------------------------
    if (is_binomial){
      input_df <- process_LUT(wcvp_fuzzy_LUT, input_df, checklist, source_name = "rWCVP_binomial", name_col = "binomialName")
    }
    else{

      input_df <- process_LUT(wcvp_fuzzy_LUT, input_df, checklist, source_name = "rWCVP", name_col = name_col)
    }

    return(input_df)
  }
  #######################################

  if (fuzzy)
  {
    print("Running WCVP")
    input_df <- wcvp_run(input_df, name_col = name_col, rerun = wcvp_rerun, is_binomial = FALSE)
    if (ssp_mods) {
      subset_df <- dplyr::filter(input_df, is.na(acceptedName)) %>%
        dplyr::mutate(
          binomialName = stringr::word(.data[[name_col]], 1, 2)
        )
      print("Running WCVP binomial")
      input_df <- wcvp_run(input_df, name_col = name_col, rerun = wcvp_rerun, is_binomial = TRUE)
    }
  }

  return(input_df)
}



