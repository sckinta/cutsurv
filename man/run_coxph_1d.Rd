% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_coxph_1d.R
\name{run_coxph_1d}
\alias{run_coxph_1d}
\title{Run Cox proportional hazards regression on one continuous variable with quantile cutoff}
\usage{
run_coxph_1d(
  formula,
  input_data,
  q_col,
  q = 0.5,
  q_2 = NULL,
  min_group_size = 10,
  baseline = NULL
)
}
\arguments{
\item{formula}{a formula specifying the model to be fit.}

\item{input_data}{a data frame containing the variables specified in the formula.}

\item{q_col}{a character string specifying the column name of the continuous variable to be used for grouping.}

\item{q}{a numeric value between 0 and 1 specifying the quantile threshold for dividing the continuous variable into two groups.}

\item{q_2}{an optional numeric value between q and 1 specifying the upper quantile threshold for dividing the continuous variable into three groups.}

\item{min_group_size}{an integer specifying the minimum number of observations required for each group.}

\item{baseline}{an optional character string specifying the baseline group for Cox model. Must be "high" or "low".}
}
\value{
a list containing the fitted Cox model object and the filtered input data.
}
\description{
This function fits a Cox proportional hazards model using a continuous variable. The continuous variable is converted to a categorical variable using quantiles, and the resulting groups are filtered to remove those with less than a specified minimum number of observations.
}
\examples{
data(HNSC)
run_coxph_1d(
   formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading),
   input_data = HNSC,
   q_col="TNFRSF1B",
   q = 0.5)

run_coxph_1d(
   formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading),
   input_data = HNSC,
   q_col="TNFRSF1B",
   q = 0.25, q_2 = 0.75,
   min_group_size = 5, baseline = "low")
}
