#' Locate the closest set of items to a given query
#'
#' @param q a named vector of query values
#' @param seq_num a sequence of integers indicating the possible set sizes
#' @param explored_df a data frame containing previously explored sets
#' @param set_n an optional integer indicating the number of sets to return
#'
#' @return a data frame of the closest set(s) of items to the query


.locate_closest_set <- function(q, seq_num, explored_df, set_n = NULL){
        stopifnot(!is.null(names(q)))
        tmp <- do.call("crossing", purrr::map(1:length(q), ~seq_num))
        names(tmp) <- names(q)
        tmp_m <- as.matrix(tmp)
        tmp_m2 <- tmp_m - matrix(rep(q,each=nrow(tmp_m)),nrow=nrow(tmp_m))

        tmp <- tmp |>
                dplyr::mutate(dist2optimal = rowSums(tmp_m2^2)) |>
                dplyr::arrange(dist2optimal) |>
                dplyr::anti_join(explored_df)

        if(is.null(set_n)){
                tmp |>
                        dplyr::select_at(names(q))
        }else{
                tmp |>
                        utils::head(set_n) |>
                        dplyr::select_at(names(q))
        }

}
