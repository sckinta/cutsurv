#' Run Cox proportional hazards model with one continuous variable using optimal two-boundary cutoffs
#'
#' Similar to `run_coxph_1d_optimal`, the funciton performs coxph survival analysis by converting single continuous variable to categorical variable with optimal quantile cutoff. However, the optimal cutoff have two boundaries now, including upper and lower quantile cutoff. "median-anchored method" controls upper quantile cutoff and starts with 0.5.
#'
#' @param formula A formula object specifying the survival object and
#'   covariates to include in the regression model.
#' @param input_data A data frame containing the variables specified in the formula.
#' @param q_col A character string specifying the name of the first continuous variable in
#'   `input_data` to regress on in  survival `formula`
#' @param low_q the lower quantile for generating the search space for optimal cutoff values (default: 0.2)
#' @param high_q the upper quantile for generating the search space for optimal cutoff values (default: 0.8)
#' @param breaks the number of intervals to divide the search space into (if \code{by} is not specified)
#' @param by the width of each interval in the search space (if \code{breaks} is not specified)
#' @param pvalue_type the type of p-value to be used for ranking the results (default: "log")
#' @param min_group_size the minimum number of observations required for each group (default: 10)
#' @param baseline An optional character string specifying the method for
#'   estimating the baseline hazard function. Must be one of "low",
#'   "high". Defaults to NULL, which uses "high".
#' @param method the method for generating the search space for optimal cutoff values (default: "median-anchored")
#' @param set_n The number of quantitle combinations closest to the previous set of cutoff generating lowest p-value for the median-anchored method. Default is 3.
#'
#' @return A list with the following elements:
#' \item{best_test}{The output of the coxph function for the best combination of cutoffs.}
#' \item{best_data}{The input data for the best combination of cutoffs.}
#' \item{best_q}{The best combination of cutoffs.}
#' \item{rank_q}{A data frame with the p-value of the score test for each combination of cutoffs.}
#'
#' @export
#'
#' @examples
#' data(HNSC)
#' run_coxph_1d2b_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B,
#'     input_data = HNSC,
#'     q_col="TNFRSF1B",
#'     by = 0.1)
#'
#'
#' @importFrom broom glance
#' @importFrom dplyr filter select_at arrange mutate across anti_join bind_rows
#' @importFrom purrr pmap_dfr
#' @importFrom stats quantile
#' @importFrom tidyr crossing
#' @importFrom utils head

run_coxph_1d2b_optimal <- function(formula, input_data, q_col, low_q=0.2, high_q=0.8, breaks=NULL, by=NULL, pvalue_type=c("log", "sc", "wald")[1], min_group_size=10, baseline = NULL, method = c("median-anchored","exhaustive")[1], set_n = 3){

        stopifnot(!all(c(is.null(breaks), is.null(by))))
        if(is.null(by)){
                seq_num <- seq(low_q, high_q, length.out=breaks)
        }else if(is.null(breaks)){
                seq_num <- seq(low_q, high_q, by=by)
        }

        full_space <- tidyr::crossing(
                q=seq_num,
                q_2=seq_num
        ) %>%
                dplyr::filter(q_2 >= q)

        if(method=="exhaustive"){
                search_space <- full_space
        }else if(method=="median-anchored"){
                search_space <- tidyr::crossing(
                        q = seq_num,
                        q_2 = NA
                )
        }

        df <- purrr::pmap_dfr(
                search_space,
                function(q, q_2){
                        res <- run_coxph_1d(formula, input_data, q_col, q=q, q_2=q_2, min_group_size = min_group_size, baseline = baseline)

                        if(!is.null(res$test)){

                                res$test %>%
                                        broom::glance() %>%
                                        dplyr::mutate(q=q, q_2=q) %>%
                                        dplyr::mutate(cutoff=stats::quantile(res$data[[q_col]], q), cutoff_2=stats::quantile(res$data[[q_col]], q_2))
                        }else{
                                NULL
                        }

                }
        ) %>%
                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                dplyr::select_at(c("q","q_2","cutoff", "cutoff_2", paste0("p.value.",pvalue_type)))

        q <- unlist(df[1,c("q","q_2")])

        if(method %in% c("median-anchored")){
                while(nrow(dplyr::anti_join(full_space, df))!=0){
                        min_p <- unlist(df[1, paste0("p.value.",pvalue_type)])

                        df2 <- purrr::pmap_dfr(
                                .locate_closest_set(q, seq_num, explored_df=df) %>%
                                        dplyr::filter(q_2 > q) %>%
                                        utils::head(set_n),
                                function(q, q_2){
                                        res <- run_coxph_1d(formula, input_data, q_col, q=q, q_2=q_2, min_group_size = min_group_size, baseline = baseline)

                                        if(!is.null(res$test)){

                                                res$test %>%
                                                        broom::glance() %>%
                                                        dplyr::mutate(q=q, q_2=q_2) %>%
                                                        dplyr::mutate(cutoff=stats::quantile(res$data[[q_col]], q), cutoff_2=stats::quantile(res$data[[q_col]], q_2))
                                        }else{
                                                NULL
                                        }

                                }
                        ) %>%
                                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                                dplyr::select_at(c("q","q_2","cutoff", "cutoff_2", paste0("p.value.",pvalue_type)))

                        new_min_p <- unlist(df2[1, paste0("p.value.",pvalue_type)])
                        df <- dplyr::bind_rows(df, df2) %>% dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type)))
                        q <- unlist(df[1,c("q", "q_2")])
                        if(new_min_p >= min_p){
                                break
                        }

                }
        }

        res <- run_coxph_1d(formula, input_data, q_col, q=q[1], q_2=q[2], min_group_size = min_group_size, baseline = baseline)

        list(
                best_test = res$test,
                best_data = res$data,
                best_q = q,
                rank_q = df

        )
}
