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

- A file *validationResults.tab* that has the hyperparamters and the resulting AUC for the validation data set on the final boosting iteration *M*.
- A file *featureImportance.tab* that has feature importance information on the final model.
- A file *featureImportance.png* which is a plot of feature importance as measured by Gain restricted to top 20.

# Runtime 
Runtime depends on how many SNPs you are using and how many observations you have. If you have the time please feel free to make a parralellized version. :)

An example: Using 200k SNPs and 10k observations with 36 threads and 150GB of RAM takes approx 4 hours to run.
