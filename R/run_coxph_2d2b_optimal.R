#' Run Cox proportional hazards model with two continuous variables using optimal two-boundary cutoffs
#'
#' Similar to `run_coxph_2d_optimal`, the funciton performs coxph survival analysis by converting two continuous variables to categorical variables with optimal quantile cutoffs. However, the optimal cutoffs have two boundaries now, including upper and lower quantile cutoff. "median-anchored method" controls second variable and starts with single median cutoff.
#'
#' @param formula a formula specifying the survival object and covariates of interest
#' @param input_data a data.frame containing the survival object and covariates specified in \code{formula}
#' @param q1_col the name of the column in \code{input_data} containing the first continuous covariate
#' @param q2_col the name of the column in \code{input_data} containing the second continuous covariate
#' @param low_q the lower quantile for generating the search space for optimal cutoff values (default: 0.2)
#' @param high_q the upper quantile for generating the search space for optimal cutoff values (default: 0.8)
#' @param breaks the number of intervals to divide the search space into (if \code{by} is not specified)
#' @param by the width of each interval in the search space (if \code{breaks} is not specified)
#' @param pvalue_type the type of p-value to be used for ranking the results (default: "log")
#' @param min_group_size the minimum number of observations required for each group (default: 10)
#' @param combine_group A character vector specifying the group levels to be combined into "other". It must be from "high:high", "high:low", "low:high" and "low:low". Default is NULL, meaning no combination.
#' @param baseline A character string specifying the group level to be used as the baseline for the analysis.  It must be "high:high", "high:low", "low:high", "low:low" or 'other'. Default is NULL, meaning using "high:high" as baseline.
#' @param method the method for generating the search space for optimal cutoff values (default: "median-anchored")
#' @param set_n The number of quantitle combinations closest to the previous set of cutoff generating lowest p-value for the median-anchored method. Default is 3.
#'
#' @return A list with the following elements:
#' \item{best_test}{The output of the coxph function for the best combination of cutoffs.}
#' \item{best_data}{The input data for the best combination of cutoffs.}
#' \item{best_q}{The best combination of cutoffs.}
#' \item{rank_q}{A data frame with the p-value of the score test for each combination of cutoffs.}
#'
#' @examples
#' data(HNSC)
#' run_coxph_2d2b_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~
#'         TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading),
#'     input_data = HNSC,
#'     q1_col="TNFRSF1B",
#'     q2_col = "t_nk_cell:Treg",
#'     by=0.1)
#'
#' run_coxph_2d2b_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~
#'         TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading),
#'     input_data = HNSC,
#'     q1_col="TNFRSF1B",
#'     q2_col = "t_nk_cell:Treg",
#'     by=0.2,
#'     combine_group = c("low:high", "high:low"),
#'     method = "exhaustive")
#'
#' @importFrom tidyr crossing
#' @importFrom dplyr filter mutate distinct arrange select_at across anti_join bind_rows
#' @importFrom purrr pmap_dfr
#' @importFrom stats quantile
#' @importFrom broom glance
#' @importFrom utils head
#'
#' @export

run_coxph_2d2b_optimal <- function(formula, input_data, q1_col, q2_col, low_q=0.2, high_q=0.8, breaks=NULL, by=NULL, pvalue_type=c("log", "sc", "wald")[1], min_group_size=10, combine_group = NULL, baseline = NULL, method = c("median-anchored", "exhaustive")[1], set_n = 3){

        stopifnot(!all(c(is.null(breaks), is.null(by))))
        if(is.null(by)){
                seq_num <- seq(low_q, high_q, length.out=breaks)
        }else if(is.null(breaks)){
                seq_num <- seq(low_q, high_q, by=by)
        }

        full_space <- tidyr::crossing(
                q1 = seq_num,
                q1_2 = seq_num,
                q2 = seq_num,
                q2_2 = seq_num,
        ) %>%
                dplyr::filter(q1_2 > q1, q2_2 > q2)

        if(method=="exhaustive"){
                search_space <- full_space
        }else if(method=="median-anchored"){
                search_space <- full_space %>%
                        dplyr::mutate(q2_2=NA, q2 = 0.5) %>%
                        dplyr::distinct()
        }

        df <- purrr::pmap_dfr(
                search_space,
                function(q1, q1_2, q2, q2_2){
                        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q1, q1_2=q1_2, q2=q2, q2_2=q2_2, min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

                        if(!is.null(res$test)){

                                res$test %>%
                                        broom::glance() %>%
                                        dplyr::mutate(q1=q1, q1_2 = q1_2, q2=q2, q2_2=q2) %>%
                                        dplyr::mutate(cutoff1=stats::quantile(res$data[[q1_col]], q1),
                                               cutoff1_2=stats::quantile(res$data[[q1_col]], q1_2),
                                               cutoff2=stats::quantile(res$data[[q2_col]], q2),
                                               cutoff2_2=stats::quantile(res$data[[q2_col]], q2)
                                        )
                        }else{
                                NULL
                        }

                }
        ) %>%
                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                dplyr::select_at(c("q1","q1_2", "cutoff1", "cutoff1_2", "q2", "q2_2","cutoff2", "cutoff2_2", paste0("p.value.",pvalue_type)))

        q <- unlist(df[1,c("q1", "q1_2","q2", "q2_2")])

        # loop through "median-anchored" method
        if(method=="median-anchored"){
                while(nrow(dplyr::anti_join(full_space, df))!=0){
                        min_p <- unlist(df[1, paste0("p.value.",pvalue_type)])
                        df2 <- purrr::pmap_dfr(
                                .locate_closest_set(q, seq_num, explored_df=df) %>%
                                        dplyr::filter(q1_2 > q1, q2_2 > q2) %>%
                                        utils::head(set_n),
                                function(q1, q1_2, q2, q2_2){
                                        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q1, q1_2=q1_2, q2=q2, q2_2=q2_2, min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

                                        if(!is.null(res$test)){

                                                res$test %>%
                                                        broom::glance() %>%
                                                        dplyr::mutate(q1=q1, q1_2 = q1_2, q2=q2, q2_2=q2_2) %>%
                                                        dplyr::mutate(cutoff1=stats::quantile(res$data[[q1_col]], q1),
                                                               cutoff1_2=stats::quantile(res$data[[q1_col]], q1_2),
                                                               cutoff2=stats::quantile(res$data[[q2_col]], q2),
                                                               cutoff2_2=stats::quantile(res$data[[q2_col]], q2_2)
                                                        )
                                        }else{
                                                NULL
                                        }

                                }
                        ) %>%
                                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                                dplyr::select_at(c("q1","q1_2", "cutoff1", "cutoff1_2", "q2", "q2_2","cutoff2", "cutoff2_2", paste0("p.value.",pvalue_type)))
                        new_min_p <- unlist(df2[1, paste0("p.value.",pvalue_type)])
                        df <- dplyr::bind_rows(df, df2) %>% dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type)))
                        q <- unlist(df[1,c("q1", "q1_2","q2", "q2_2")])
                        if(new_min_p >= min_p){
                                break
                        }

                }
        }


        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q["q1"], q1_2=q["q1_2"], q2=q["q2"], q2_2=q["q2_2"], min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

        list(
                best_test = res$test,
                best_data = res$data,
                best_q = q,
                rank_q = df

        )
}
