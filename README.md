
# INDIRECT IDENTIFICATION OF RESISTANCE BIOMARKERS

Note that this is an archival repository for the exact code and data used in:
Ayestaran et al., Identification of Intrinsic Drug Resistance and Its Biomarkers in High-Throughput Pharmacogenomic and CRISPR Screens, Patterns (2020), https://doi.org/10.1016/j.patter.2020.100065

Further work on the analysis pipeline, or inclusion of more data will take place in the development repository: https://github.com/ia327/indirect_resistance_stat_framework
 

## R package dependencies
- **DRANOVA**: Source code provided in the `DRANOVA/` folder. To install from within R, set the working directory to where this file is and run:
		`install.packages(‘./DRANOVA/', repos = NULL, type = 'source’)`
- tidyverse
- scales
- qvalue
- beeswarm
- readxl
- methods
- shiny
- parallel



## INPUT FILES
Included in the `data/` folder, there is a collection of .xlsx and .tsv files.

To speed data input, a .RData file is created at the beginning of the pipeline containing all data used. To change the input data, the script `create_combined_data.R` should be modified.

Column names are used throughout the analysis for subsetting and manipulating data. If these are not maintained in the new input files several scripts should be changed accordingly.

In the case of drug response data, used column names are: `COSMIC_ID`, `DRUG_ID`, `LN_IC50`, `MAX_CONC`.


## RUNNING THE ANALYSIS
The bash script `indirect_resistance.sh` calls the R scripts that perform the different steps of the analysis.

These include:
1. The generation of the `combined_data.RData` file
2. Tissue-specific ANOVA analysis
3. Filtering of ANOVA output
4. Detection of resistant outlier cell lines using the variance change method.
5. Identification of putative resistance markers in these outlier resistant cohorts.


A choice between GDSC and CTRP datasets is done by specifying the flag `-g` or `--gdsc`, and `-c` or `--ctrp`, respectively.

The first 2 steps are the most computationally expensive, therefore `indirect_resistance.sh` won’t run them unless `-C` and `-A` arguments are included, respectively.

## OUTPUT
A number of files and plots are created over the analysis and are stored in the `results/` and `plots/` folders.

By order of creation:
- Subfolder `ANOVA_out`: contains the results for the tissue-specific ANOVAs.
- `sign_anova_out.csv`: Table containing all drug-CFE associations that fulfil our filtering criteria.
- Volcano plots: Show results of the ANOVA before filtering (`volcano_all.pdf` and `volcano_sign.pdf`) with respective barcharts showing the breakdown in different tissues.
- `all_found_outliers.csv`: Table with all resistant outliers and outlier groups identified.
- `unique_outliers_detected.csv`: Summarises the previous table grouping it by drug-CFE association.
- `outlier_cases_details.txt`: Detailed list of each outlier identified and description of the type of mutation, difference in copy-number or other co-occurring alteration. Sorted by tissue.
- `putative_res_markers.csv`: Table with outliers that were found to contain some differential marker that could explain the resistance. Final output of the analysis.

Some .RData objects are also created, and are used to faster transmit the data between different R scripts.
