#' Run Cox proportional hazards regression on one continuous variable with quantile cutoff
#'
#' This function fits a Cox proportional hazards model using a continuous variable. The continuous variable is converted to a categorical variable using quantiles, and the resulting groups are filtered to remove those with less than a specified minimum number of observations.
#'
#'
#' @param formula a formula specifying the model to be fit.
#' @param input_data a data frame containing the variables specified in the formula.
#' @param q_col a character string specifying the column name of the continuous variable to be used for grouping.
#' @param q a numeric value between 0 and 1 specifying the quantile threshold for dividing the continuous variable into two groups.
#' @param q_2 an optional numeric value between q and 1 specifying the upper quantile threshold for dividing the continuous variable into three groups.
#' @param min_group_size an integer specifying the minimum number of observations required for each group.
#' @param baseline an optional character string specifying the baseline group for Cox model. Must be "high" or "low".
#'
#' @return a list containing the fitted Cox model object and the filtered input data.
#'
#' @importFrom dplyr select_at mutate case_when semi_join filter ungroup count group_by
#' @importFrom rlang is_formula
#' @importFrom stats na.omit quantile terms as.formula
#' @importFrom survival coxph Surv strata
#'
#' @export
#'
#' @examples
#' data(HNSC)
#' run_coxph_1d(
#'    formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading),
#'    input_data = HNSC,
#'    q_col="TNFRSF1B",
#'    q = 0.5)
#'
#' run_coxph_1d(
#'    formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~ TNFRSF1B + strata(cancer_grading),
#'    input_data = HNSC,
#'    q_col="TNFRSF1B",
#'    q = 0.25, q_2 = 0.75,
#'    min_group_size = 5, baseline = "low")


run_coxph_1d <- function(formula, input_data, q_col, q=0.5, q_2=NULL, min_group_size = 10, baseline = NULL){
        # control input
        stopifnot(rlang::is_formula(formula))

        stopifnot("data.frame" %in% class(input_data))

        stopifnot(q_col %in% names(input_data))

        stopifnot(all(is.numeric(c(q,q_2))))

        if(!any(c(is.null(q_2), is.na(q_2)))){
                stopifnot(q_2 >= q)
        }


        if(!is.null(baseline)){
                stopifnot(baseline %in% c("high", "low"))
        }

        # assign group based on q
        tmp <- input_data %>%
                dplyr::select_at(all.vars(formula)) %>%
                stats::na.omit()

        if(!any(c(is.null(q_2), is.na(q_2)))){
                tmp <- tmp %>%
                        dplyr::mutate(group = dplyr::case_when(get(q_col) <= stats::quantile(get(q_col), q) ~ "low",
                                                 get(q_col) >= stats::quantile(get(q_col), q_2) ~ "high",
                                                 TRUE ~ as.character(NA))) %>%
                        stats::na.omit()
        }else{
                tmp <- tmp %>%
                        dplyr::mutate(group = ifelse(get(q_col) <= stats::quantile(get(q_col), q), "low", "high"))
        }




        # factor group
        group_levels <- c("high", "low")

        if(!is.null(baseline)){
                group_levels <- c(baseline, group_levels[group_levels!=baseline])
        }

        tmp <- tmp  %>%
                dplyr::mutate(group = factor(group, levels=group_levels)) %>%
                I()


        # filter by min_group_size
        tmp <- tmp %>%
                dplyr::semi_join(
                        tmp %>%
                                dplyr::group_by(group) %>%
                                dplyr::count() %>%
                                dplyr::ungroup() %>%
                                dplyr::filter(n >= min_group_size)
                ) %>%
                dplyr::mutate(group=droplevels(group))

        if(length(levels(tmp$group)) > 1){
                # make formula
                x_sides <- attr(stats::terms(formula), 'term.labels') # formula x side
                x_sides <- x_sides[!x_sides %in% c(q_col)]
                x_sides <- x_sides[!x_sides %in% gsub("^", "`", gsub("$","`",c(q_col)))] # in case some labels are using ```
                y_sides <- rownames(attr(stats::terms(formula), 'factors'))[1] # forumula y side
                if(length(x_sides)==0){
                        new_formula <- stats::as.formula(paste(y_sides, "~", "group"))
                }else{
                        new_formula <- stats::as.formula(paste(y_sides, "~", "group", "+", x_sides))
                }


                # perform test
                test <- survival::coxph(new_formula, data = tmp)
        }else{
                test <- NULL
                warning("There is only 1 group left")
        }

        list(test = test, data = tmp)
}
