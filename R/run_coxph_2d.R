#' Run Cox proportional hazards regression with two continuous variables
#'
#' Perform coxph survival analysis by converting two continuous variables to categorical variables with user defined quantile cutoff.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param input_data A data.frame containing the variables specified in the formula.
#' @param q1_col A character string specifying the name of the first continuous variable.
#' @param q2_col A character string specifying the name of the second continuous variable.
#' @param q1 A numeric value specifying the quantile cutoff for the first continuous variable.
#' @param q1_2 A numeric value specifying the upper quantile cutoff for the first continuous variable.
#' Set to NULL or NA if no upper cutoff is desired.
#' @param q2 A numeric value specifying the quantile cutoff for the second continuous variable.
#' @param q2_2 A numeric value specifying the upper quantile cutoff for the second continuous variable.
#' Set to NULL or NA if no upper cutoff is desired.
#' @param min_group_size An integer specifying the minimum number of observations in each group.
#' @param combine_group A character vector specifying the group levels to be combined into "other". It must be from "high:high", "high:low", "low:high" and "low:low".
#' @param baseline A character string specifying the group level to be used as the baseline for the analysis.  It must be "high:high", "high:low", "low:high", "low:low" or 'other'.
#'
#' @return A list with two elements: "test", which contains the results of the Cox proportional hazards regression,
#' and "data", which contains the data used for the analysis.
#'
#' @importFrom dplyr select_at mutate case_when select semi_join filter ungroup count group_by
#' @importFrom glue glue
#' @importFrom rlang is_formula
#' @importFrom stats na.omit quantile terms as.formula
#' @importFrom survival coxph Surv strata
#'
#' @export
#' @examples
#' data(HNSC)
#' run_coxph_2d(
#'    formula=Surv(SURVIVAL_os_days, SURVIVAL_os_event) ~
#'        TNFRSF1B + `t_nk_cell:Treg` + strata(cancer_grading),
#'    input_data = HNSC,
#'    q1_col="TNFRSF1B",
#'   q2_col = "t_nk_cell:Treg",
#'   q1=0.4, q1_2 = 0.6, q2 = 0.5)
#'

run_coxph_2d <- function(formula, input_data, q1_col, q2_col, q1=0.5, q1_2=NULL, q2=0.5, q2_2 = NULL, min_group_size = 10, combine_group = NULL, baseline = NULL){

        # control input
        stopifnot(rlang::is_formula(formula))

        stopifnot("data.frame" %in% class(input_data))

        stopifnot(all(c(q1_col, q2_col) %in% names(input_data)))

        stopifnot(all(is.numeric(c(q1,q2,q1_2, q2_2))))

        if(!any(c(is.null(q1_2), is.na(q1_2)))){
                stopifnot(q1_2 >= q1)
        }

        if(!any(c(is.null(q2_2), is.na(q2_2)))){
                stopifnot(q2_2 >= q2)
        }


        if(!is.null(combine_group)){
                stopifnot(all(combine_group %in% c("high:high", "high:low", "low:high", "low:low")))
        }

        if(!is.null(baseline)){
                stopifnot(baseline %in% c("high:high", "high:low", "low:high", "low:low", "other"))
        }

        # assign group based on q
        tmp <- input_data %>%
                dplyr::select_at(all.vars(formula)) %>%
                stats::na.omit()

        if(!any(c(is.null(q1_2), is.na(q1_2)))){
                tmp <- tmp %>%
                        dplyr::mutate(q1_cat = dplyr::case_when(get(q1_col) <= stats::quantile(get(q1_col), q1) ~ "low",
                                                  get(q1_col) >= stats::quantile(get(q1_col), q1_2) ~ "high",
                                                  TRUE ~ as.character(NA))) %>%
                        stats::na.omit()
        }else{
                tmp <- tmp %>%
                        dplyr::mutate(q1_cat = ifelse(get(q1_col) <= stats::quantile(get(q1_col), q1), "low", "high"))
        }

        if(!any(c(is.null(q2_2), is.na(q2_2)))){
                tmp <- tmp %>%
                        dplyr::mutate(q2_cat = dplyr::case_when(get(q2_col) <= stats::quantile(get(q2_col), q2) ~ "low",
                                                  get(q2_col) >= stats::quantile(get(q2_col), q2_2) ~ "high",
                                                  TRUE ~ as.character(NA))) %>%
                        stats::na.omit()
        }else{
                tmp <- tmp %>%
                        dplyr::mutate(q2_cat = ifelse(get(q2_col) <= stats::quantile(get(q2_col), q2), "low", "high"))
        }

        tmp <- tmp %>%
                dplyr::mutate(group = glue::glue("{q1_cat}:{q2_cat}")) %>%
                dplyr::select(-q1_cat, -q2_cat)

        # deal with combine_group
        group_levels <- c("high:high", "high:low", "low:high", "low:low")

        if(!is.null(combine_group)){
                tmp <- tmp %>%
                        dplyr::mutate(group = ifelse(group %in% combine_group, "other", group))
                group_levels <- unique(tmp$group)
        }

        # factor group
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
                dplyr::mutate(group = droplevels(group))

        if(length(levels(tmp$group)) > 1){
                # make formula
                x_sides <- attr(stats::terms(formula), 'term.labels') # formula x side
                x_sides <- x_sides[!x_sides %in% c(q1_col, q2_col)]
                x_sides <- x_sides[!x_sides %in% gsub("^", "`", gsub("$","`",c(q1_col, q2_col)))] # in case some labels are using ```
                y_sides <- rownames(attr(stats::terms(formula), 'factors'))[1] # forumula y side
                if(length(x_sides)==0){
                        new_formula <- stats::as.formula(paste(y_sides, "~", "group"))
                }else{
                        new_formula <- stats::as.formula(paste(y_sides, "~", "group +", x_sides))
                }


                # perform test
                test <- survival::coxph(new_formula, data = tmp)
        }else{
                test <- NULL
                warning("There is only 1 group left")
        }

        list(test = test, data = tmp)
}


