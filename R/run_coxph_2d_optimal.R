#' Run Cox proportional hazards model with two continuous variables using optimal cutoffs
#'
#' Coxph survival analysis by converting two continuous variables to categorical variables with optimal quantile cutoffs. The optimal cutoffs are determined by exhaustively searching or using median-anchored method. The optimal cutoffs are determined by minimizing the p-value of the score test of the Cox proportional hazards model.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param input_data A data frame containing the variables specified in the formula.
#' @param q1_col A character string specifying the name of the first continuous variable.
#' @param q2_col A character string specifying the name of the second continuous variable.
#' @param low_q The lower limit of the quantile range used to generate the search space for optimal cutoffs.
#' @param high_q The upper limit of the quantile range used to generate the search space for optimal cutoffs.
#' @param breaks The number of quantiles used to generate the search space for optimal cutoffs. Only used if by is NULL.
#' @param by The distance between the quantiles used to generate the search space for optimal cutoffs. Only used if breaks is NULL.
#' @param pvalue_type The type of p-value used for ranking optimal cutoffs. Can be "log", "sc", or "wald".
#' @param min_group_size An integer specifying the minimum number of observations in each group. Default is 10
#' @param combine_group A character vector specifying the group levels to be combined into "other". It must be from "high:high", "high:low", "low:high" and "low:low". Default is NULL, meaning no combination.
#' @param baseline A character string specifying the group level to be used as the baseline for the analysis.  It must be "high:high", "high:low", "low:high", "low:low" or 'other'. Default is NULL, meaning using "high:high" as baseline.
#' @param method The method used to determine the optimal cutoffs. Can be "median-anchored" or "exhaustive".
#' @param set_n The number of quantitle combinations closest to the previous set of cutoff generating lowest p-value for the median-anchored method. Default is 3.
#'
#' @return A list with the following elements:
#' \item{best_test}{The output of the coxph function for the best combination of cutoffs.}
#' \item{best_data}{The input data for the best combination of cutoffs.}
#' \item{best_q}{The best combination of cutoffs.}
#' \item{rank_q}{A data frame with the p-value of the score test for each combination of cutoffs.}
#'
#' @importFrom broom glance
#' @importFrom dplyr select_at arrange mutate across anti_join bind_rows
#' @importFrom purrr pmap_dfr
#' @importFrom stats quantile
#' @importFrom tidyr crossing
#'
#' @export
#' @examples
#' data(HNSC)
#' run_coxph_2d_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~
#'         TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading),
#'     input_data = HNSC,
#'     q1_col="TNFRSF1B",
#'     q2_col = "t_nk_cell:Treg",
#'     breaks=10)
#'
#' run_coxph_2d_optimal(
#'     formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~
#'         TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading),
#'     input_data = HNSC,
#'     q1_col="TNFRSF1B",
#'     q2_col = "t_nk_cell:Treg",
#'     by=0.1,
#'     combine_group = c("low:high", "high:low"),
#'     baseline="other",
#'     set_n = 5)

run_coxph_2d_optimal <- function(formula, input_data, q1_col, q2_col, low_q=0.2, high_q=0.8, breaks=NULL, by=NULL, pvalue_type=c("log", "sc", "wald")[1], min_group_size=10, combine_group = NULL, baseline = NULL, method = c("median-anchored", "exhaustive")[1], set_n = 3){

        stopifnot(!all(c(is.null(breaks), is.null(by))))
        if(is.null(by)){
                seq_num <- seq(low_q, high_q, length.out=breaks)
        }else if(is.null(breaks)){
                seq_num <- seq(low_q, high_q, by=by)
        }

        full_space <- tidyr::crossing(
                q1 = seq_num,
                q2 = seq_num
        )
        if(method=="exhaustive"){
                search_space <- full_space
        }else if(method=="median-anchored"){
                search_space <- tidyr::crossing(
                        q1 = seq_num,
                        q2 = 0.5
                )
        }

        df <- purrr::pmap_dfr(
                search_space,
                function(q1, q2){
                        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q1, q2=q2, min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

                        if(!is.null(res$test)){

                                res$test %>%
                                        broom::glance() %>%
                                        dplyr::mutate(q1=q1, q2=q2) %>%
                                        dplyr::mutate(cutoff1=stats::quantile(res$data[[q1_col]], q1), cutoff2=stats::quantile(res$data[[q2_col]], q2))
                        }else{
                                NULL
                        }

                }
        ) %>%
                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                dplyr::select_at(c("q1","cutoff1", "q2","cutoff2", paste0("p.value.",pvalue_type)))

        q <- unlist(df[1,c("q1", "q2")])

        # loop through "median-anchored" method
        if(method=="median-anchored"){
                while(nrow(dplyr::anti_join(full_space, df))!=0){
                        min_p <- unlist(df[1, paste0("p.value.",pvalue_type)])
                        df2 <- purrr::pmap_dfr(
                                .locate_closest_set(q, seq_num, explored_df=df, set_n=set_n),
                                function(q1, q2){
                                        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q1, q2=q2, min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

                                        if(!is.null(res$test)){

                                                res$test %>%
                                                        broom::glance() %>%
                                                        dplyr::mutate(q1=q1, q2=q2) %>%
                                                        dplyr::mutate(cutoff1=stats::quantile(res$data[[q1_col]], q1), cutoff2=stats::quantile(res$data[[q2_col]], q2))
                                        }else{
                                                NULL
                                        }
                                }
                        ) %>%
                                dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type))) %>%
                                dplyr::select_at(c("q1","cutoff1", "q2","cutoff2", paste0("p.value.",pvalue_type)))
                        new_min_p <- unlist(df2[1, paste0("p.value.",pvalue_type)])
                        df <- dplyr::bind_rows(df, df2) %>% dplyr::arrange(dplyr::across(paste0("p.value.",pvalue_type)))
                        q <- unlist(df[1,c("q1", "q2")])
                        if(new_min_p >= min_p){
                                break
                        }

                }
        }


        res <- run_coxph_2d(formula, input_data, q1_col, q2_col, q1=q["q1"], q2=q["q2"], min_group_size=min_group_size, combine_group=combine_group, baseline=baseline)

        list(
                best_test = res$test,
                best_data = res$data,
                best_q = q,
                rank_q = df

        )
}



