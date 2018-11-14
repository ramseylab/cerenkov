## Instructions

### 1. Download all the input files

```bash
./wget_input_files.sh
```

4 files will be downloaded. They are:

- `osu18_features1.1_cerenkov2(anonymous).rds`
- `osu18_intra_locus_dist(anonymous).rds`
- `osu18_snp_coordinates(anonymous).rds`
- `osu18_feature_names_categories.txt`

### 2. Run `xgboost` models for performance comparison

```bash
Rscript run_performance_models.R
```

- Performance records will be saved in `records_of_performance_CERENKOV2_1337.Rdata`

### 3. Analyze the performance of `xgboost` models

```bash
Rscript analyze_performance_records.R records_of_performance_CERENKOV2_1337.Rdata analysis_of_performance
```

- `analysis_of_performance` can be any folder name you specify.
- If this folder did not exist, the script will create one automatically.
- Subplots of **Figure 4** and p-values of **Table 1** will be generated in the folder.

### 4. Run `ranger` models for feature importance

```bash
Rscript run_feature_importance_models.R
```

- Impurity and permutation records will be saved in `records_of_feature_importance_CERENKOV2_1337.Rdata`

### 5. Analyze the feature importance

```bash
Rscript analyze_feature_importance_records.R records_of_feature_importance_CERENKOV2_1337.Rdata analysis_of_feature_importance
```

- `analysis_of_feature_importance` can be any folder name you specify.
- If this folder did not exist, the script will create one automatically.
- Plots of **Figure 5** will be generated in the folder.