Run Cox proportional hazards regression on one continuous variable with
quantile cutoff
================

- <a href="#using-single-continous-variable"
  id="toc-using-single-continous-variable">Using single continous
  variable</a>
  - <a href="#user-defined-quantile-cutoffs-to-group-samples-run_coxph_1d"
    id="toc-user-defined-quantile-cutoffs-to-group-samples-run_coxph_1d">User
    defined quantile cutoff(s) to group samples
    <code>run_coxph_1d</code></a>
    - <a href="#with-one-cutoff" id="toc-with-one-cutoff">With one cutoff</a>
    - <a href="#with-both-lower-and-upper-quantile-cutoff"
      id="toc-with-both-lower-and-upper-quantile-cutoff">With both lower and
      upper quantile cutoff</a>
  - <a
    href="#scanning-a-gradient-of-quantile-cutoffs-to-find-optimal-cutoff"
    id="toc-scanning-a-gradient-of-quantile-cutoffs-to-find-optimal-cutoff">Scanning
    a gradient of quantile cutoff(s) to find optimal cutoff</a>
    - <a href="#with-one-cutoff-1" id="toc-with-one-cutoff-1">With one
      cutoff</a>
    - <a href="#with-both-lower-and-upper-quantile-cutoff-1"
      id="toc-with-both-lower-and-upper-quantile-cutoff-1">With both lower and
      upper quantile cutoff</a>
- <a href="#using-two-continous-variable"
  id="toc-using-two-continous-variable">Using two continous variable</a>
  - <a href="#user-defined-quantile-cutoffs-to-group-samples-run_coxph_2d"
    id="toc-user-defined-quantile-cutoffs-to-group-samples-run_coxph_2d">User
    defined quantile cutoff(s) to group samples
    <code>run_coxph_2d</code></a>
    - <a href="#with-one-cutoff-2" id="toc-with-one-cutoff-2">With one
      cutoff</a>
    - <a href="#with-both-lower-and-upper-quantile-cutoff-2"
      id="toc-with-both-lower-and-upper-quantile-cutoff-2">With both lower and
      upper quantile cutoff</a>
  - <a
    href="#scanning-a-gradient-of-quantile-cutoffs-to-find-optimal-cutoff-1"
    id="toc-scanning-a-gradient-of-quantile-cutoffs-to-find-optimal-cutoff-1">Scanning
    a gradient of quantile cutoff(s) to find optimal cutoff</a>
    - <a href="#with-one-cutoff-3" id="toc-with-one-cutoff-3">With one
      cutoff</a>
    - <a href="#with-both-lower-and-upper-quantile-cutoff-3"
      id="toc-with-both-lower-and-upper-quantile-cutoff-3">With both lower and
      upper quantile cutoff</a>

``` r
library(cutsurv)
library(tidyverse)
library(survival)
# library(skimr)
```

We will use HNSC data from the package for the examples below. To know
more about HNSC data, refer to `?HNSC`

``` r
data(HNSC)

HNSC
#> # A tibble: 502 × 23
#>    sample_id        subject_id   tumor_type sample_type  age_at_diagnosi… gender race 
#>    <chr>            <chr>        <chr>      <chr>                   <dbl> <chr>  <chr>
#>  1 TCGA-4P-AA8J-01A TCGA-4P-AA8J HNSC       Primary Tum…               66 male   Blac…
#>  2 TCGA-BA-4074-01A TCGA-BA-4074 HNSC       Primary Tum…               69 male   White
#>  3 TCGA-BA-4075-01A TCGA-BA-4075 HNSC       Primary Tum…               49 male   Blac…
#>  4 TCGA-BA-4076-01A TCGA-BA-4076 HNSC       Primary Tum…               39 male   White
#>  5 TCGA-BA-4077-01B TCGA-BA-4077 HNSC       Primary Tum…               45 female White
#>  6 TCGA-BA-4078-01A TCGA-BA-4078 HNSC       Primary Tum…               83 male   White
#>  7 TCGA-BA-5151-01A TCGA-BA-5151 HNSC       Primary Tum…               72 male   White
#>  8 TCGA-BA-5152-01A TCGA-BA-5152 HNSC       Primary Tum…               56 male   White
#>  9 TCGA-BA-5153-01A TCGA-BA-5153 HNSC       Primary Tum…               51 male   White
#> 10 TCGA-BA-5555-01A TCGA-BA-5555 HNSC       Primary Tum…               54 male   Blac…
#> # … with 492 more rows, and 16 more variables: SURVIVAL_os_days <dbl>,
#> #   SURVIVAL_os_event <dbl>, SURVIVAL_pfi_days <dbl>, SURVIVAL_pfi_event <dbl>,
#> #   cancer_grading <chr>, tnm_stage_pathologic_overall_stage <chr>, TNFRSF1B <dbl>,
#> #   `t_nk_cell:CD4+ T` <dbl>, `t_nk_cell:CD8+ T` <dbl>, `t_nk_cell:Treg` <dbl>,
#> #   `t_nk_cell:Cytotoxic` <dbl>, `t_nk_cell:NK` <dbl>, `t_nk_cell:Exhausted` <dbl>,
#> #   `myeloid:Anti-inflammatory` <dbl>, `myeloid:Pro-inflammatory` <dbl>,
#> #   subtype <chr>

# skimr::skim(HNSC)
```

# Using single continous variable

## User defined quantile cutoff(s) to group samples `run_coxph_1d`

Here we use overall survival response to regress against the expression
of *TNFRSF1B*, stratified by patient cancer grading.

### With one cutoff

Discretising patients by the median expression of *TNFRSF1B*, with
higher than median group as baseline.

``` r
res <- run_coxph_1d(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading), input_data = HNSC, q_col="TNFRSF1B", q = 0.5)

res[["test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>            coef exp(coef) se(coef)     z      p
#> grouplow 0.3168    1.3728   0.1401 2.262 0.0237
#> 
#> Likelihood ratio test=5.17  on 1 df, p=0.02296
#> n= 498, number of events= 217
```

Result also returns discretized input data which contains all the
variables in formula and can be used for other model fit in the future.

``` r
res[["data"]]
#> # A tibble: 498 × 5
#>    SURVIVAL_os_days SURVIVAL_os_event TNFRSF1B cancer_grading group
#>               <dbl>             <dbl>    <dbl> <chr>          <fct>
#>  1              102                 0     3.61 G2             high 
#>  2              462                 1     1.68 G3             low  
#>  3              283                 1     1.87 G2             low  
#>  4              415                 1     1.81 G2             low  
#>  5             1134                 1     3.88 G2             high 
#>  6              276                 1     2.74 G2             low  
#>  7              722                 0     3.51 G1             high 
#>  8             1288                 0     4.41 G2             high 
#>  9             1762                 1     4.74 G2             high 
#> 10              520                 0     3.17 G2             low  
#> # … with 488 more rows
```

### With both lower and upper quantile cutoff

- Expression of *TNFRSF1B* can also be binned to lower than quantile `q`
  group and higher than quantile `q_2` group.

- We can also switch the baseline for regression. The default from above
  is “high”.

``` r
res <- run_coxph_1d(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading), input_data = HNSC, q_col="TNFRSF1B", q = 0.4, q_2=0.6, baseline = "low")

res[["test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>              coef exp(coef) se(coef)      z     p
#> grouphigh -0.3268    0.7212   0.1550 -2.109 0.035
#> 
#> Likelihood ratio test=4.48  on 1 df, p=0.03423
#> n= 398, number of events= 178
```

For more helps, `?run_coxph_1d`.

## Scanning a gradient of quantile cutoff(s) to find optimal cutoff

### With one cutoff

`run_coxph_1d_optimal` performs Cox proportional hazards regression with
multiple quantile cutoff values and returns the best test result, best
data, best quantile cutoff, and a dataframe of the quantile cutoffs
ranked by the selected p-value type.

- The cutoff gradient can be defined by `low_q`, `high_q` and `breaks`
  or `by`, very similarly to `{seq}` function.
- The best cutoff are selected by choosing threshold with the lowest
  model p-value. There are several types of p-value for ranking. `"log"`
  means likelihood ratio p-value, `"sc"` means log rank score test
  p-value, and `"wald"` means wald test p-value. Default is `"log"`
- After discretion, the groups with lower than `min_group_size` will be
  removed from regression.
- Like `run_coxph_1d`, we can also choose `baseline` for the group
  level.

``` r
res <- run_coxph_1d_optimal(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B, input_data = HNSC, q_col="TNFRSF1B", by = 0.1)
```

The result include 4 parts:

- “best_test”: The coxph test object with the lowest p-value.
- “best_data”: The data frame used to fit the best test.
- “best_q”: The quantile cutoff value that produced the best test.
- “rank_q”: A data frame with the ranked quantile cutoff values and
  their corresponding p-values for the selected p-value type.

``` r
res[["best_test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>            coef exp(coef) se(coef)     z      p
#> grouplow 0.3603    1.4338   0.1426 2.528 0.0115
#> 
#> Likelihood ratio test=6.16  on 1 df, p=0.01306
#> n= 501, number of events= 218
res[["best_q"]]
#>   q 
#> 0.3
res[["rank_q"]]
#> # A tibble: 7 × 3
#>       q cutoff p.value.log
#>   <dbl>  <dbl>       <dbl>
#> 1   0.3   2.69      0.0131
#> 2   0.4   2.96      0.0136
#> 3   0.5   3.25      0.0146
#> 4   0.2   2.31      0.0633
#> 5   0.7   3.93      0.0771
#> 6   0.8   4.21      0.0895
#> 7   0.6   3.54      0.108
```

To learn more, `?run_coxph_1d_optimal`.

### With both lower and upper quantile cutoff

`run_coxph_1d2b_optimal` is similar to `run_coxph_1d_optimal` but with
both lower `q` and upper `q_2` quantile cutoffs to group samples. All
other parameters are the same with 2 additional parameters.

- Method for generating the search space for optimal cutoff values.
  Inspired from cSurvival paper[^1], there are two methods included:
  “median-anchored” and “exhaustive”.

> 1)  Median-anchored greedy search: we construct a 2D grid using
>     percentiles of both predictors. Next, we determine the starting
>     point for a greedy search by locating the minimum P-value computed
>     from testing each percentile in predictor B against the median
>     percentile in predictor A. Then, we test the nearest three
>     unexplored points; if a lower P-value is found, we move the search
>     to that newly found minimum P-value point and test the nearest
>     unexplored points until no lower P-value can be found.

> 2)  Exhaustive search: We construct a 2D grid using percentiles of
>     both predictors. Next, we determine the optimal percentile
>     combination by locating the minimum P-value computed from testing
>     each percentile in predictor B against each percentile in
>     predictor A.

- set_n: The number of quantitle combinations closest to the previous
  set of cutoff generating lowest p-value for the median-anchored
  method. Default is 3.

``` r
res <- run_coxph_1d2b_optimal(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B, input_data = HNSC, q_col="TNFRSF1B", low_q = 0.2, high_q = 0.8,by = 0.1)
```

``` r
res[["best_test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>            coef exp(coef) se(coef)    z       p
#> grouplow 0.4325    1.5412   0.1539 2.81 0.00495
#> 
#> Likelihood ratio test=7.74  on 1 df, p=0.005398
#> n= 402, number of events= 175
res[["best_q"]]
#>   q q_2 
#> 0.3 0.5
res[["rank_q"]]
#> # A tibble: 16 × 5
#>        q   q_2 cutoff cutoff_2 p.value.log
#>    <dbl> <dbl>  <dbl>    <dbl>       <dbl>
#>  1   0.3   0.5   2.51     3.52     0.00540
#>  2   0.4   0.5   2.86     3.38     0.00785
#>  3   0.3   0.4   2.61     3.11     0.00962
#>  4   0.3   0.6   2.37     4.00     0.0130 
#>  5   0.3   0.3   2.69     2.69     0.0131 
#>  6   0.4   0.4   2.96     2.96     0.0136 
#>  7   0.5   0.5   3.25     3.25     0.0146 
#>  8   0.2   0.5   2.01     3.71     0.0182 
#>  9   0.4   0.6   2.72     3.85     0.0200 
#> 10   0.2   0.4   2.10     3.29     0.0299 
#> 11   0.2   0.6   1.91     4.11     0.0331 
#> 12   0.2   0.3   2.23     2.86     0.0377 
#> 13   0.2   0.2   2.31     2.31     0.0633 
#> 14   0.7   0.7   3.93     3.93     0.0771 
#> 15   0.8   0.8   4.21     4.21     0.0895 
#> 16   0.6   0.6   3.54     3.54     0.108
```

For more information, `?run_coxph_1d2b_optimal`.

# Using two continous variable

## User defined quantile cutoff(s) to group samples `run_coxph_2d`

Here we use overall survival response to regress against the expression
of *TNFRSF1B* and the GSVA enrichment score of Treg `"t_nk_cell:Treg"`,
stratified by patient cancer grading and gender.

### With one cutoff

Discretising patients by the median expression of *TNFRSF1B* and median
enrichment score, with higher than median group as baseline. It will
group samples into 4 groups: “high:high”, “high:low”, “low:high”,
“low:low”.

``` r
res <- run_coxph_2d(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading, gender), input_data = HNSC, q1_col="TNFRSF1B", q2_col = "t_nk_cell:Treg", q1=0.5, q2 = 0.5)

res[["test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>                  coef exp(coef) se(coef)      z      p
#> grouphigh:low -0.4016    0.6693   0.2675 -1.501 0.1333
#> grouplow:high  0.0312    1.0317   0.2129  0.147 0.8835
#> grouplow:low   0.3029    1.3538   0.1608  1.883 0.0597
#> 
#> Likelihood ratio test=9.21  on 3 df, p=0.02667
#> n= 498, number of events= 217
```

Same as `run_coxph_1d`, result also returns discretized input data which
contains all the variables in formula and can be used for other model
fit in the future.

``` r
res[["data"]]
#> # A tibble: 498 × 7
#>    SURVIVAL_os_days SURVIVAL_os_event TNFRSF1B `t_nk_cell:Treg` cancer_grading gender
#>               <dbl>             <dbl>    <dbl>            <dbl> <chr>          <chr> 
#>  1              102                 0     3.61            0.665 G2             male  
#>  2              462                 1     1.68            0.540 G3             male  
#>  3              283                 1     1.87            0.274 G2             male  
#>  4              415                 1     1.81            0.596 G2             male  
#>  5             1134                 1     3.88            0.803 G2             female
#>  6              276                 1     2.74            0.641 G2             male  
#>  7              722                 0     3.51            0.666 G1             male  
#>  8             1288                 0     4.41            0.760 G2             male  
#>  9             1762                 1     4.74            0.554 G2             male  
#> 10              520                 0     3.17            0.607 G2             male  
#> # … with 488 more rows, and 1 more variable: group <fct>
```

We can also combine some groups into one, so that instead of regressing
against 4 groups, we only interested in whether either TNFR2 low and
Treg high enrichment will result in poor prognosis. Thus, “high:low”,
“low:low”, “high:high” will be integrated into one group `"other"`.

``` r
res <- run_coxph_2d(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading, gender), input_data = HNSC, q1_col="TNFRSF1B", q2_col = "t_nk_cell:Treg", q1=0.5, q2 = 0.5, combine_group = c("high:low", "low:low", "high:high"))

res[["test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>                  coef exp(coef) se(coef)     z    p
#> grouplow:high -0.0542    0.9472   0.1937 -0.28 0.78
#> 
#> Likelihood ratio test=0.08  on 1 df, p=0.7783
#> n= 498, number of events= 217
```

### With both lower and upper quantile cutoff

Same as `run_coxph_1d`, both variables can also be binned to lower than
quantile `q1` and/or `q2` and higher than quantile `q1_2` and/or `q2_2`.

``` r
res <- run_coxph_2d(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading, gender), input_data = HNSC, q1_col="TNFRSF1B", q2_col = "t_nk_cell:Treg", q1=0.4, q1_2 = 0.6, q2 = 0.5)

res[["test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>                  coef exp(coef) se(coef)      z      p
#> grouphigh:low -0.5270    0.5904   0.3201 -1.646 0.0997
#> grouplow:high  0.1124    1.1189   0.2393  0.470 0.6387
#> grouplow:low   0.2602    1.2971   0.1773  1.467 0.1424
#> 
#> Likelihood ratio test=7.81  on 3 df, p=0.05019
#> n= 398, number of events= 178
```

For more helps, `?run_coxph_2d`.

## Scanning a gradient of quantile cutoff(s) to find optimal cutoff

`run_coxph_2d_optimal` and `run_coxph_2d2d_optimal` perform coxph
survival analysis by converting two continuous variables to categorical
variables with optimal quantile cutoffs. The optimal cutoffs are
determined by exhaustively searching or using median-anchored method.
The optimal cutoffs are determined by minimizing the p-value of the
score test of the Cox proportional hazards model.

### With one cutoff

`run_coxph_2d_optimal` is to search combination of `q1` and `q2` for
both variables. `q1` and `q2` separate var1 and var2 into two groups
respectively. It is flexible with all the parameters defined in
`run_coxph_1d_optimal` seen above, but search 2-dimensional space.

Below is the example to evaluate, comparing to both TNFR2 and Treg
enrichment high, whether single high or both low result in poor or
better prognosis using “median-anchored” to search the optimal pair of
`q1` and `q2` cutoff.

``` r
res <- run_coxph_2d_optimal(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading), input_data = HNSC, q1_col="TNFRSF1B", q2_col = "t_nk_cell:Treg", by=0.1, combine_group = c("low:high", "high:low"), baseline="high:high", set_n = 5)
```

``` r
res[["best_test"]]
#> Call:
#> survival::coxph(formula = new_formula, data = tmp)
#> 
#>                 coef exp(coef) se(coef)      z      p
#> grouplow:low  0.2874    1.3330   0.1591  1.807 0.0708
#> groupother   -0.1615    0.8509   0.1835 -0.880 0.3789
#> 
#> Likelihood ratio test=7.03  on 2 df, p=0.02968
#> n= 498, number of events= 217
res[["best_q"]]
#>  q1  q2 
#> 0.5 0.5
res[["rank_q"]]
#> # A tibble: 12 × 5
#>       q1 cutoff1    q2 cutoff2 p.value.log
#>    <dbl>   <dbl> <dbl>   <dbl>       <dbl>
#>  1   0.5    3.25   0.5   0.630      0.0297
#>  2   0.4    2.96   0.6   0.668      0.0312
#>  3   0.6    3.54   0.5   0.630      0.0319
#>  4   0.5    3.25   0.6   0.668      0.0399
#>  5   0.4    2.96   0.5   0.630      0.0506
#>  6   0.3    2.70   0.5   0.630      0.0655
#>  7   0.4    2.96   0.4   0.581      0.0950
#>  8   0.2    2.34   0.5   0.630      0.100 
#>  9   0.5    3.25   0.4   0.581      0.109 
#> 10   0.6    3.54   0.4   0.581      0.164 
#> 11   0.7    3.93   0.5   0.630      0.209 
#> 12   0.8    4.21   0.5   0.630      0.302
```

For more information, `?run_coxph_2d_optimal`

### With both lower and upper quantile cutoff

Same as `run_coxph_1d2b_optimal`, `run_coxph_2d2b_optimal` search the
combination of `q1`, `q2`, `q1_2` and `q2_2`, binning samples by lower
and upper quantile cutoff.

Below is the example same as above, but with lower and upper quantile
cutoffs.

``` r
res <- run_coxph_2d2b_optimal(formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading), input_data = HNSC, q1_col="TNFRSF1B", q2_col = "t_nk_cell:Treg", by=0.1, combine_group = c("low:high", "high:low"), baseline="high:high", set_n = 5)
```

[^1]: Xuanjin Cheng, Yongxing Liu, Jiahe Wang, Yujie Chen, Andrew Gordon
    Robertson, Xuekui Zhang, Steven J M Jones, Stefan Taubert,
    cSurvival: a web resource for biomarker interactions in cancer
    outcomes and in cell lines, Briefings in Bioinformatics, Volume 23,
    Issue 3, May 2022, bbac090, <https://doi.org/10.1093/bib/bbac090>
