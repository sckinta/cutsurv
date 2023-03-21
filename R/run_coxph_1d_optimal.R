#' Find optimal quantile cutoff for Cox proportional hazards regression
#'
#' This function performs Cox proportional hazards regression with multiple
#' quantile cutoff values and returns the best test result, best data, best
#' quantile cutoff, and a dataframe of the quantile cutoffs ranked by the
#' selected p-value type.
#'
#' @param formula A formula object specifying the survival object and
#'   covariates to include in the regression model.
#' @param input_data A data frame containing the variables specified in the formula.
#' @param q_col A character string specifying the name of the first continuous variable in
#'   `input_data` to regress on in  survival `formula`
#' @param low_q A numeric value specifying the lowest quantile cutoff. Defaults to 0.2.
#' @param high_q A numeric value specifying the highest quantile cutoff. Defaults to 0.8.
#' @param breaks An integer specifying the number of quantile cutoffs to use
#'   if `by` is not specified. Ignored if `by` is specified. Defaults to NULL.
#' @param by A numeric value specifying the increment between quantile cutoffs
#'   to use if `breaks` is not specified. Ignored if `breaks` is specified.
#'   Defaults to NULL.
#' @param pvalue_type A character string specifying the type of p-value to use
#'   when ranking the quantile cutoffs. Must be one of "log", "sc", or "wald".
#'   Defaults to "log".
#' @param min_group_size An integer specifying the minimum number of
#'   observations in each group for the regression to be performed. Defaults
#'   to 10.
#' @param baseline An optional character string specifying the method for
#'   estimating the baseline hazard function. Must be one of "low",
#'   "high". Defaults to NULL, which uses "high".
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{best_test}{The coxph test object with the lowest p-value.}
#'     \item{best_data}{The data frame used to fit the best test.}
#'     \item{best_q}{The quantile cutoff value that produced the best test.}
#'     \item{rank_q}{A data frame with the ranked quantile cutoff values and
#'       their corresponding p-values for the selected p-value type.}
#'   }
#'
#' @export
#' @importFrom purrr map_dfr
#' @importFrom broom glance
#' @importFrom dplyr across mutate select_at arrange
#' @examples
#' data(HNSC)
#' run_coxph_1d_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B,
#'     input_data = HNSC,
#'     q_col="TNFRSF1B",
#'     by = 0.1)

run_coxph_1d_optimal <- function(formula, input_data, q_col, low_q=0.2, high_q=0.8, breaks=NULL, by=NULL,  pvalue_type=c("log", "sc", "wald")[1], min_group_size=10, baseline = NULL){

        stopifnot(!all(c(is.null(breaks), is.null(by))))
        if(is.null(by)){
                seq_num <- seq(low_q, high_q, length.out=breaks)
        }else if(is.null(breaks)){
                seq_num <- seq(low_q, high_q, by=by)
        }

        df <- purrr::map_dfr(
                seq_num,
                function(q){
                        res <- run_coxph_1d(formula, input_data, q_col, q=q, min_group_size = min_group_size, baseline = baseline)

                        if(!is.null(res$test)){

                                res$test %>%
                                        broom::glance() %>%
                                        dplyr::mutate(q=q) %>%
                                        dplyr::mutate(cutoff=stats::quantile(res$data[[q_col]], q))
                        }else{
                                NULL
                        }

                }
        ) %>%
                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                dplyr::select_at(c("q","cutoff", paste0("p.value.",pvalue_type)))

        q <- unlist(df[1,c("q")])

        res <- run_coxph_1d(formula, input_data, q_col, q=q, min_group_size = min_group_size, baseline = baseline)

        list(
                best_test = res$test,
                best_data = res$data,
                best_q = q,
                rank_q = df

        )
}
