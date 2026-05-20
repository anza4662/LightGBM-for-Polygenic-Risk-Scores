# LightGBM for Polygenic Risk Scores

This is a way of using plink2 .bed files to train, validate and test LightGBM models.

# Data
Use plink2 to reformat your data from .bed file to tabular data. Use the following command

```shell
plink2 --bfile <bfile> \
--export A \
--out <outfile> 
```

# Running the analysis
Run the *lgbmPRS.R* file with the path to the data file that was created above and a directory for the outputs.

```shell
Rscript lgbmPRS.R "<path_to_tabular_data>" "<path_to_outputs>"
```

# Output
...