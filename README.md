# BulkCAT

Multispecies NatureServe element rarity rank calculator

Created by Clark Hollenberg at the Colorado Natural Heritage Program, October 2024

---

## **Introduction**

BulkCAT calculates conservation status metrics for multiple species using user-supplied point occurrence data with WGS84 latitude and longitude coordinates. 
The package implements a bulk-processing version of the "Conservation status Assessment Tool" (CAT) methodology from IUCN and NatureServe. 
For each species, it computes extent of occurrence (EOO in km²), area of occupancy (AOO, number of 2x2 km grid cells), number of hypothetical Element Occurrences (EOs), and an overall rarity rank (subnational rank - SRank). 
This allows efficient multi-species assessments at scales not feasible with interactive tools like GeoCAT or RARECAT. Calculated values may differ from RARECAT by ~1% for large datasets.

---

## **Installation**

You can install BulkCAT directly from GitHub using the `remotes` package:

```r
# Install remotes if you don't already have it
install.packages("remotes")

# Install BulkCAT from GitHub
remotes::install_github("chollenb-cnhp/BulkCAT")
```

After installation, load the package:

```r
library(BulkCAT)
```

---

## **Input data for ranking**

The input `.csv` should include at least three columns:

* `scientificName` – species name
* `decimalLatitude` – decimal degrees in WGS84
* `decimalLongitude` – decimal degrees in WGS84

You can change the names of these columns, but if you do, they must also be modified in the function call.

There is an option to include a polygon layer as a shapefile for calculation of AOO, which may be helpful for plant communities or large occurrences best defined by a polygon.

For vascular plants, occurrence data were downloaded from SEINet and iNaturalist research grade. Note that GBIF downloads may have obscured locations for some species. A SEINet login with special permissions allows use of the most accurate locations available.

Scientific names were translated to SNAMEs used in Biotics using an iterative approach with Biotics synonyms, rWCVP, Symbiota, GNAME/SNAME mapping, and infraspecific epithet dropping to reach \~99% coverage. Deduplication was performed conservatively to avoid deleting unique records, prioritizing duplicates based on herbarium record counts for the sample region (Colorado).

---

## **Ranking Rules**

This bulk calculator considers only rank factors in the rarity category.
It uses:

* **Range Extent (EOO)**
* **Area of Occupancy (AOO)** – double weighted
* **Number of hypothetical EOs**

Rules are based on the **RULES sheet** from the Element Rank Estimator Excel macro workbook (NatureServe): [link](https://www.natureserve.org/products/conservation-rank-calculator/download).
Bin divisions for AOO were adjusted to match the expected logic of the calculator.

---

## **Infraspecific ranking**

* Binomial species without infraspecific epithet are clustered with any infraspecific taxa for ranking.

  * Example: *Abies lasiocarpa* includes points from *Abies lasiocarpa var. bifolia*.
* Trinomial names are ranked using only exact matches.

---

## **Estimating Number of EOs**

* Plant occurrences are clustered using a **1 km separation distance**.
* Invertebrates or other taxa might use a **5 km separation distance**.
* This is adjustable via the function arguments or at the start of the script.

---

## **Example usage**

```r
# Load package
library(BulkCAT)

# Load example input CSV (if included in the package)
csv_path <- "path/to/input_file.csv"
input_df <- read.csv(csv_path)

# Remove duplicates
input_df <- deduplicate(
  input_df = input_df,
  cols = c("recordedBy", "recordNumber", "scientificName", "eventDate")
  institution_col = "institutionCode"
)

# Run BulkCAT
results <- run_bulkCAT(
  input_df = input_df,
  sname = "scientificName",
  lat = "decimalLatitude",
  lon = "decimalLongitude",
  eo_separation = 1000,  # meters
  grid_size = 2000,      # meters
  community = FALSE,
  poly_layer = NULL
)

# View results
head(results)

# export to csv
output_path <- "path/to/output.csv
results.to_csv(output_path)

```


---
