library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

#########################
# Code for building lookup tables - not needed for package.
###########################################
ns_data = read_excel("../Synonyms/allPlants_NS_20260127.xlsx")

ns_lookup_df <- ns_data %>%
  filter(!is.na(Synonyms), Synonyms != "") %>%
  separate_rows(Synonyms, sep = ",\\s*") %>%
  mutate(
    inputName = str_trim(Synonyms),
    outputName = scientificName
  ) %>%
  select(inputName, outputName) %>%
  distinct()

write.csv(ns_lookup_df, "../Synonyms/allPlants_NS_LUT.csv")

wcvp_names = read.csv("../Synonyms/wcvp_names_clean.csv")

wcvp_lookup_df <- wcvp_names %>%
  # Keep only synonyms
  filter(taxon_status == "Synonym") %>%

  # Join to same table to get accepted name
  left_join(
    wcvp_names %>%
      select(
        plant_name_id,
        outputName = taxon_name
      ),
    by = c("accepted_plant_name_id" = "plant_name_id")
  ) %>%

  # Build final mapping
  transmute(
    inputName  = taxon_name,
    outputName = outputName
  ) %>%
  distinct()

SEINet_names <- read.csv("../Synonyms/symbiota_lookup.csv")

SEINet_lookup <- SEINet_names %>%
  filter(str_count(inputName, "\\S+") >= 2)

usda_names <- read.csv("../Synonyms/usda_plant_codes_20250107.csv")

accepted <- usda_names %>% filter(symbol == synonymSymbol) %>%
  select(
    acceptedSymbol = symbol,
    acceptedName   = scientificName
  )

usda_lookup_df <- usda_names %>%
  # Keep only true synonyms
  filter(symbol != synonymSymbol) %>%

  # Join: synonymSymbol (accepted code) -> symbol (accepted row)
  left_join(
    accepted,
    by = c("symbol" = "acceptedSymbol")
  ) %>%

  # Build final mapping
  transmute(
    inputName  = scientificName,   # the synonym name
    outputName = acceptedName      # the accepted name
  ) %>%
  distinct()

################################################
## synonymize function:
######################
# Sketch of function:
# synonym_LUTs is a list of user-supplied dataframes. Tables must be have two columns called "inputName" and "outputName". Multiple synonym dataframes can be passed to the function.

# add column to input_df that is the acceptedName
# for all rows in input_df where name_col == checklist_name_col, update acceptedName to name_col value.
# First use synonym LUTs supplied by the user in the order provided.
# filter user defined LUTs to where synonym_LUTs[i]$outputName != synonym_LUTs[i]$inputName
#if synonym_LUTs[i]$inputName is in checklist[[checklist_name_col]], switch synonym_LUTs[i]$inputName and synonym_LUTs[i]$outputName for that row
# filter user defined LUTs to synonym_LUTs[i]$outputName is in checklist[[checklist_name_col]].
# if acceptedName is NA and translation is available for input_df[[name_col]], update acceptedName

# next use synonym_sources to do the same thing in the order provided.

# for any names that could not be translated, use the rWCVP package and check for any other translations that meet the requirements defined above.
 Sketch
#############################



synonymize <- function(input_df,
                       name_col = "scientificName",
                       checklist = NA,
                       checklist_name_col = "SNAME",
                       synonym_LUTs = list(),
                       synonym_sources = c("NatureServe", "SEINet", "USDA", "WCVP"),
                       fuzzy = TRUE,
                       ssp_mods = FALSE) {

  # 1️⃣ Initialize acceptedName and translationSource
  input_df <- input_df %>%
    mutate(
      acceptedName = NA_character_,
      translationSource = NA_character_
    )
  if (missing(checklist) || is.null(checklist)){
    print("No checklist supplied. Using NatureServe Network Tracheophyta checklist...")
    checklist = read.csv("Synonyms/NatureServe.csv")
    checklist_name_col = "outputName"
  }
  # 2️⃣ Direct matches to checklist
  if (!is.na(checklist)) {
    # Rows where input name is already in checklist
    direct_idx <- input_df[[name_col]] %in% checklist[[checklist_name_col]]
    input_df$acceptedName[direct_idx] <- input_df[[name_col]][direct_idx]
    input_df$translationSource[direct_idx] <- "Direct"
  }

  # 3️⃣ Function to process one LUT
  process_LUT <- function(LUT, df, checklist, source_name = "LUT", name_col = "scientificName") {

    # Filter trivial rows
    LUT <- LUT %>% filter(inputName != outputName)

    # Swap input/output if inputName is in checklist
    if (!is.na(checklist)) {
      LUT <- LUT %>%
        mutate(
          swap_flag = inputName %in% checklist[[checklist_name_col]]
        ) %>%
        mutate(
          tmp = ifelse(swap_flag, inputName, outputName),
          outputName = ifelse(swap_flag, outputName, tmp),
          inputName  = ifelse(swap_flag, tmp, inputName)
        ) %>%
        select(-tmp, -swap_flag)
    }

    # Filter LUT to outputName in checklist if provided
    if (!is.na(checklist)) {
      LUT <- LUT %>% filter(outputName %in% checklist[[checklist_name_col]])
    }

    # Left join LUT to df
    df <- df %>%
      left_join(LUT, by = setNames("inputName", name_col)) %>%
      mutate(
        # Only update acceptedName if not yet assigned
        acceptedName = ifelse(is.na(acceptedName) & !is.na(outputName), outputName, acceptedName),
        translationSource = ifelse(is.na(translationSource) & !is.na(outputName), source_name, translationSource)
      ) %>%
      select(-outputName)

    return(df)
  }

run_synonyms <- function(input_df, name_col){


  # 4 apply user provided LUTs in order
  if (length(synonym_LUTs) > 0) {
    for (i in seq_along(synonym_LUTs)) {
      LUT <- synonym_LUTs[[i]]
      source_name <- ifelse(i <= length(synonym_sources), synonym_sources[i], paste0("LUT", i))
      input_df <- process_LUT(LUT, input_df, checklist, source_name, name_col)
    }
  }

  # 5 apply built-in LUTs
  if (length(synonym_sources) > 0) {
    for (source_name in synonym_sources) {

      LUT <- switch(source_name,
                    "NatureServe" = read.csv("Synonyms/NatureServe.csv"),
                    "SEINet"      = read.csv("Synonyms/SEINet.csv"),
                    "USDA"        = read.csv("Synonyms/USDA.csv"),
                    "WCVP"        = read.csv("Synonyms/WCVP.csv"),
                    NULL)
      if (!is.null(LUT)) {
        input_df <- process_LUT(LUT, input_df, checklist, source_name, name_col)
      }
    }
  }
  ##################################
  # fuzzy matching
  ###########################
  # --- WCVP translation LUTs ---
  if (fuzzy && requireNamespace("rWCVP", quietly = TRUE)) {

    # 1️⃣ Run WCVP name matching
    wcvp_matches <- rWCVP::wcvp_match_names(input_df[[name_col]], fuzzy = TRUE)

    # 2️⃣ Build Exact match LUT
    wcvp_exact_LUT <- wcvp_matches %>%
      dplyr::filter(grepl("^Exact", match_type)) %>%
      dplyr::mutate(
        inputName = .data[[name_col]],
        outputName = wcvp_name,
        translationSource = "rWCVP"
      ) %>%
      dplyr::select(inputName, outputName, translationSource) %>%
      dplyr::filter(inputName != outputName)

    # 3️⃣ Build Fuzzy match LUT
    wcvp_fuzzy_LUT <- wcvp_matches %>%
      dplyr::filter(grepl("^Fuzzy", match_type)) %>%
      dplyr::mutate(
        inputName = .data[[name_col]],
        outputName = wcvp_name,
        translationSource = "Fuzzy"
      ) %>%
      dplyr::select(inputName, outputName, translationSource) %>%
      dplyr::filter(inputName != outputName)

    # 4️⃣ Apply LUTs using existing process_LUT
    input_df <- process_LUT(wcvp_exact_LUT, input_df, checklist, source_name = "rWCVP")
    input_df <- process_LUT(wcvp_fuzzy_LUT, input_df, checklist, source_name = "Fuzzy")

  } else if (!requireNamespace("rWCVP", quietly = TRUE)) {
    warning("rWCVP not installed; skipping WCVP translation.")
  }
  return(input_df)
}
  input_df <- run_synonyms(input_df, name_col = name_col)

  ###############################
  # repeat the whole process with binomials only
  #################################
  if (ssp_mods)
  {
    subset_df <- input_df %>% filter(is.na(acceptedName)) %>%
      mutate(
        binomialName = word(.data[[name_col]], 1, 2)
      )
    subset_df <- subset_df %>% run_synonyms(name_col = "binomialName")

    subset_LUT <- subset_df %>%
      filter(!is.na(acceptedName)) %>%
      mutate(
        inputName = .data[[name_col]],   # original value
        outputName = acceptedName
      ) %>%
      select(inputName, outputName) %>%
      distinct()

    input_df <- process_LUT(subset_LUT, input_df, checklist, source_name = "binomial")
  }

  return(input_df)
}




