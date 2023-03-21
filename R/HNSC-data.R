#' Head & Neck Squamous Cell Carcinomas (HNSC) survival data
#'
#' Data is derived from the primary tumor in HNSC from TCGA. The gene expression of *TNFRSF1B* is from RNA-seq FPKQ. The cell type value (from `t_nk_cell:CD4+ T` to `myeloid:Pro-inflammatory`) is calculated from canonical cell markers using `GSVA` to represent the enrichment of given cell type/state in individual tumor types.
#'
#'
#' @docType data
#'
#' @usage data(HNSC)
#'
#' @format a data.frame with 23 columns representing clinical data and
#'
#' @references Hoadley KA, Yau C, Hinoue T, Wolf DM, Lazar AJ, Drill E, Shen R, Taylor AM, Cherniack AD, Thorsson V, Akbani R, Bowlby R, Wong CK, Wiznerowicz M, Sanchez-Vega F, Robertson AG, Schneider BG, Lawrence MS, Noushmehr H, Malta TM; Cancer Genome Atlas Network; Stuart JM, Benz CC, Laird PW. Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer. Cell. 2018 Apr 5;173(2):291-304.e6. doi: 10.1016/j.cell.2018.03.022. PMID: 29625048; PMCID: PMC5957518.
#'
#' @source TCGA-HNSC https://www.cbioportal.org/study/summary?id=hnsc_tcga_pan_can_atlas_2018
#'
#' @examples
#' data(HNSC)

"HNSC"
