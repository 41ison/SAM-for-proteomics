# Significance Analysis of Microarrays (SAM) applied to Proteomics

## Motivation
The Significance Analysis of Microarrays (SAM) was developed by Tusher et al. (2001) and is based on the t-test with a permutation-based approach to estimate the false discovery rate (FDR). The permutation-based FDR has a reduced error type II rate compared to the Benjamini-Hochberg correction (BH) and is currently implemented in [Perseus software](https://link.springer.com/protocol/10.1007/978-1-4939-7493-1_7) for proteomics data analysis. Check the S0 parameter in the statistics session of Perseus.

Read the original paper [here](https://www.pnas.org/content/98/9/5116).
>Tusher, V., Tibshirani, R. and Chu, G. (2001): Significance analysis of microarrays applied to the ionizing radiation response" PNAS 2001 98: 5116-5121, (Apr 24). http://www-stat.stanford.edu/~tibs/sam

:bulb:Follow the instructions to install the package [samr](https://github.com/MikeJSeo/SAM) and its dependencies.
The `samr` package is easy to use and has a shiny app that work well. This repositroy is only for help you to customize the results.

### Load the important libraries

```
library(tidyverse)
library(samr)
```

### Load the data as a matrix of abundance values
Check the class of the data and convert it to a matrix if it is not already.
Check the `class()` of the object. If it is a dataframe, coerce it to matrix with `as.matrix()` function.
calculate the average abundance of the proteins per condition and correlate them.
In this example, we are using a matrix of log2 abundance values containing two groups of samples: Control and Treated. The samples are named as `Control_1`, `Control_2`, `Control_3`, `Treated_1`, `Treated_2`, and `Treated_3`. 

```
mean_plot <- abundance_matrix %>%
    as.data.frame() %>%
    rownames_to_column("protein") %>%
    gather(-protein ,key = "Sample", value = "Intensity") %>%
    dplyr::mutate(
        condition = str_extract(Sample, "Control|Treated"),
        condition = factor(condition, levels = c("Control", "Treated"))
    ) %>%
    dplyr::group_by(protein, condition) %>%
    dplyr::summarise(mean_intensity = mean(Intensity)) %>%
    spread(key = condition, value = mean_intensity) %>%
    ggplot(aes(x = Control, y = Treated)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(
        x = "Control",
        y = "Treated"
    )

ggsave("mean_plot.png", 
    mean_plot, path = "plots",
    width = 10, height = 10, 
    units = "cm", dpi = 300)
```

Run the SAM analysis using the `samr` package. The function `SAM` requires the following arguments:
- `x`: a matrix of the data to be analyzed.
- `y`: a vector of the condition for the design. Similiar to the condition in the `limma` design matrix.
- `resp.type`: the type of response variable. In this case, we are using a two-class unpaired design.
- `geneid`: the row names of the matrix. So the proteins are identified.
- `S0`: the threshold for the significance of the test statistic.
- `nperms`: the number of permutations to be used in the analysis.
- `fdr.output`: the false discovery rate threshold for the analysis.
- `logged2`: if the data is log2 transformed. If not, do it.

```
# run the SAM analysis
sam_result <- samr::SAM(x = abundance_matrix, 
                y = condition_for_design,
                resp.type = "Two class unpaired",
                geneid = rownames(abundance_matrix),
                genenames = rownames(abundance_matrix),
                s0 = 0.1, 
                s0.perc = NULL, 
                nperms = 1000,
                center.arrays=FALSE,
                testStatistic = "standard",
                time.summary.type = "slope",
                regression.method = "standard", 
                return.x = TRUE, 
                knn.neighbors = 10,
                random.seed = NULL,
                fdr.output = 0.05,
                logged2 = TRUE,
                eigengene.number = 1)

# Scatter plot of the observed relative difference d(i) versus the expected relative difference dE(i). The solid line indicates the line for d(i) = dE(i), where the observed relative difference is identical to the expected relative difference.
samr::samr.plot(sam_result$samr.obj,
            del = 0.45,     # set the threshold for the significance of the test statistic Δ = 0.45
             min.foldchange = 1.5)
```

## Extract the significant proteins and reconstruct the SAM Q-Q plot using ggplot2

```
# extract the significant proteins
proteins_up <- data.frame(sam_result$siggenes.table$genes.up)
proteins_low <- data.frame(sam_result$siggenes.table$genes.lo)

# combine the significant proteins in a single dataframe
significant_proteins <- rbind(proteins_up, proteins_low)

# reconstructing the SAM Q-Q plot from the results
#  extract the observed relative difference d(i)
di <- data.frame(sam_result$samr.obj$tt) %>%
    rownames_to_column("protein") %>%
    dplyr::rename("di" = "sam_result.samr.obj.tt") %>%
    dplyr::arrange(di)

# extract the expected relative difference dE(i)
dEi <- data.frame(sam_result$samr.obj$evo) %>%
    dplyr::rename("dEi" = "sam_result.samr.obj.evo") %>%
    dplyr::arrange(dEi)

# combine in a single dataframe and add the status of the proteins with the significant proteins dataframe
sam_df <- cbind(di, dEi) %>%
    dplyr::left_join(significant_proteins, by = join_by("protein" == "Gene.ID")) %>%
    dplyr::mutate(
        status = case_when(
         Score.d. > 0 ~ "High",
         Score.d. < 0 ~ "Low",
         is.na(Score.d.) ~ "Not significant"
        ),
        status = factor(status, levels = c("Low", "Not significant", "High"))
    )

# reconstruct the SAM plot using ggplot2
sam_scatter_plot <- sam_df %>%
    ggplot(aes(x = dEi, y = di, color = status)) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c("#0d0887", "black", "firebrick1")) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 0.7))) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid") +
    geom_abline(intercept = 0.45, slope = 1, linetype = "dashed") +
    geom_abline(intercept = -0.45, slope = 1, linetype = "dashed") +
    labs(
        y = "Observed relative difference - d(i)",
        x = "Expected relative difference - dE(i)",
    ) +
    theme(legend.position = "bottom") +
    labs(color = NULL)

ggsave("sam_scatter_plot.png", 
    sam_scatter_plot, path = "plots",
    width = 10, height = 10, 
    units = "cm", dpi = 300)
```
