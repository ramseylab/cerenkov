# genome-wide feature_extraction folder without source_data folder and eigen values

## Things to do before running

- The absolute directory of the file "bigWigAverageOverBed" must be adjusted in the sys_tool.py module (line 102).
- Make sure to have the source_data folder. It is available in the EC2 instance.
- Make sure to have the eigen values: mart_export_hg19_chr22_SNP.score
- Yao has already made the eigen values and it can be downloaded from his workstation (it is already present in the EC2 instance).

## How to run the feature_extraction 

- python3 must be used to run all the modules.
- The main module to run is listed as "cerenkov_bef_genome_wide.py".
- For the input, the parser requires the input to be the .txt file with the rsids (single rsid per row).
- For the output, the parser requires an already created .csv file.
