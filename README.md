# hb-mig-env
Repository containing data and code for hummingbird migration project using eBird and environmental data

## Code Authors
Sarah R. Supp (supps@denison.edu), Laura Graham (l.graham@bham.ac.uk), and Frank LaSorte ()

## Collaborators
Sarah R. Supp, Laura Graham, Frank LaSorte, Catherine Graham... 

**Resources**
* [eBird Best Practices](http://strimas.com/ebird-best-practices/)
* [package auk](https://ropensci.org/blog/2018/08/07/auk/)
* [package dggridR](https://cran.r-project.org/web/packages/dggridR/vignettes/dggridR.html)
* [La Sorte and Fink 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12534)
* [Supp et al. 2015](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/ES14-00290.1)


## Data Files
Data files were queried from the EBD database in 2019 (www.ebird.org) for records 01 January 2008 - 31 December 2018. They are not provided here on GitHub as they are too large to upload and update. 

## Code Files
* `cleaned_species_data_from_full.r`: Uses `auk` package to read in data from EBD database, and retains a single record from each group_identifier, after accounting for potential discrepancies in how individuals recorded location (lat-lon).
* `create_effort_data_from_full.r`: Reads in entire eBird database from 2008-2018 to calculate total eBirder effort within a hexagonal grid, which will later be used to weight bird observations in calculating migration patterns.
* `hb-migration.r`: Reads in cleaned eBird datasets for each species and the eBird effort dataset. Calculates weighted daily population centroids for each species in each year, and summary tables for different migrtaion estimates -- start and end dates, peak latitude, population migration speed, error across all years combined, cross track progression, and cross track separation (sensu La Sorte and Fink 2017). 

## Documentation
* analysis_doc.Rmd file explains our approach and methods

