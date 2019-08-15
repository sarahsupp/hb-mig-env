# hb-mig-env
Repository containing data and code for hummingbird migration project using eBird and environmental data

## Code Authors
Sarah R. Supp (supps@denison.edu)
Laura Graham ()
Frank LaSorte ()

## Collaborators
Sarah R. Supp, Laura Graham, Frank LaSorte, Catherine Graham... 

**Resources**
* [eBird Best Practices](http://strimas.com/ebird-best-practices/)
* [package auk](https://ropensci.org/blog/2018/08/07/auk/)


## Data Files
* `cleaned_species_data_from_full.r`: Uses `auk` package to read in data from EBD database, and retains a single record from each group_identifier, after accounting for potential discrepancies in how individuals recorded location (lat-lon).
* `create_effort_data_from_full.r`: Reads in entire eBird database from 2008-2018 to calculate total eBirder effort within a hexagonal grid, which will later be used to weight bird observations in calculating migration patterns.

## Documentation
* analysis_doc.Rmd file explains our approach and methods

