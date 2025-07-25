# BulkCAT
Multispecies NatureServe rarity rank calculator
Created by Clark Hollenberg at the Colorado Natural Heritage Program, October 2024

**Introduction**

This script imports a .csv file with rows containing multi-species point occurrence data and WGS84 latitude, longitude coordinates. The output is a xlsx file containing extent of occurrence (EOO), area of occupancy with count of 2 x 2km grid cells (AOO), number of hypothetical Element Occurrences (EOs), and calculated rarity rank. This enables multi-species calculations which are not currently possible on GeoCAT. The projection coordinate system matches the .prj file from the EOO toolbox provided by IUCN. Calculated values match those of the toolbox, but may differ from GeoCAT by 1% at larger scales.

**Input data for ranking**

The input csv should include at least three columns:

•	SNAME – species name

•	decimalLatitude – decimal degrees in WGS84

•	decimalLongitude – decimal degrees in WGS84

You can change the names of these columns, but if you do they must also be modified in the R script.
I downloaded statewide occurrence data from SEINet and iNaturalist research grade. Note that occurrence data downloaded directly from GBIF may have locations obscured for certain species. A SEINet login with special permissions allows us to use the most accurate locations available. I used many different tools to clean the source data and translate the scientific names to SNAMEs used in Biotics. I won’t cover the specifics of that process here.
I provided a sample_occurrence_points_to_rank.csv file which is a slice of the iNaturalist and SEINet data that I prepared for Colorado. This can give you a sense of the fields I used and the formatting.

**Ranking Rules**

Since this bulk calculator is designed to work without manual biologist review, it only considers rank factors in the rarity category. I used Range Extent (EOO), Area of Occupancy (AOO), and Number of EOs to assign a preliminary rarity rank to each species. For the ranking rules, I used the RULES sheet within the Element Rank Estimator excel macro workbook from NatureServe at https://www.natureserve.org/products/conservation-rank-calculator/download. Note that AOO is double weighted. There appears to be an error in the Step2/3 AOO values. Based on the Calculator Form, A should be assigned to 1 4km cell and B should be assigned to 2 4km cells, and C should be assigned to 3-5 km cells. I used these bin divisions.

**Infraspecific ranking**

When ranking a binomial species without infraspecific epithet, I clustered any infraspecific taxa under the species for ranking purposes. Ex: The rank for Abies lasiocarpa would include occurrence points from Abies lasiocarpa and Abies lasiocarpa var. bifolia, whereas trinomial names like Abies lasiocarpa var. bifolia would be ranked only considering exact matches (excluding occurrences identified as only as Abies lasiocarpa).


**Estimating Number of EOs**

I used a 1km separation distance for clustering plants into hypothetical EOs. For invertebrates, we used a 5km distance. This can be adjusted at the start of the script based on the taxonomic group in question.


