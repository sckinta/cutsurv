#' Extract TCGA clinical and omics data
#'
#' @param tumor_type A character of TCGA study abbreviation. Refer to https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations. One exception is colorectal cancer is combination of COAD and READ, labeled "COADREAD".
#' @param gene_names A character vector to define hugoGeneSymbol if type is not "clinical". Default is NULL
#' @param type A character to define type of data to extract. must be one of "clinical", "expression", "CNA call", "CNA value", "mutation", "methylation" and "protein". Default is "clinical"
#' @param short_clinical A Boolean vector to define whether shorten the clinical data to include just patient_id, sample_type, sample_id, os_months, os_status, age, race, sex.
#'
#' @return a data.frame for given data type
#'
#' @export
#'
#' @examples
#' get_tcga(tumor_type="BRCA") # get clinical
#' get_tcga(tumor_type="BRCA",
#'     gene_names=c("BRCA1", "BRCA2"),
#'     type="expression") # get gene expression for BRCA1 and BRCA2
#'
#' @importFrom cBioPortalData cBioPortal clinicalData getDataByGenes
#' @importFrom dplyr as_tibble select
#' @importFrom janitor clean_names
#'
get_tcga <- function(tumor_type, gene_names=NULL, type = c("clinical", "expression", "CNA call", "CNA value", "mutation", "methylation", "protein")[1],  short_clinical=T){

        stopifnot(tumor_type %in% c('ACC','BLCA','BRCA','CESC','CHOL','COADREAD','DLBC','ESCA','GBM','HNSC','LAML','KIRC','KICH','LGG','LIHC','KIRP','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','SKCM','SARC','STAD','TGCT','THYM','THCA','UCEC','UCS','UVM'))

        # cat(gene_names)

        cbio <- cBioPortalData::cBioPortal()

        study_id <- tolower(gsub("$","_tcga", tumor_type))

        if(type == "clinical"){
                out <- cBioPortalData::clinicalData(cbio, study_id) |>
                        dplyr::as_tibble() |>
                        janitor::clean_names()

                if(short_clinical) {
                        out <- out |>
                                dplyr::select(patient_id, sample_type, sample_id, os_months, os_status, age, race, sex)
                }

                return(out)
        }else{
                stopifnot(!is.null(gene_names))
        }

        if(type == "expression"){
                molecularprofile_id = paste0(tolower(tumor_type), "_tcga_rna_seq_v2_mrna")
        }else if(type == "CNA call"){
                molecularprofile_id = paste0(tolower(tumor_type), "_tcga_gistic")
        }else if(type == "CNA value"){
                if(tumor_type=="THCA"){
                        molecularprofile_id = NULL
                }else{
                        molecularprofile_id = paste0(tolower(tumor_type), "_tcga_linear_CNA")
                }

        }else if(type == "mutation"){
                molecularprofile_id = paste0(tolower(tumor_type), "_tcga_mutations")
        }else if(type == "methylation"){
                if(tumor_type=="OV"){
                        molecularprofile_id = paste0(tolower(tumor_type), "_tcga_methylation_hm27")
                }else{
                        molecularprofile_id = paste0(tolower(tumor_type), "_tcga_methylation_hm450")
                }
        }else if(type == "protein"){
                if(tumor_type %in% c("LAML", "LIHC", "UVM")){
                        molecularprofile_id = NULL
                }else{
                        molecularprofile_id = paste0(tolower(tumor_type), "_tcga_rppa")
                }
        }


        if(!is.null(molecularprofile_id)){
                data <- cBioPortalData::getDataByGenes(api = cbio,
                                       by = "hugoGeneSymbol",
                                       studyId = study_id,
                                       genes = gene_names,
                                       molecularProfileIds = molecularprofile_id
                ) # here I download both BRCA1 and BRCA2 expression.

                gene_expr <- data[[molecularprofile_id]]

        }else{
                gene_expr <- NULL
        }

        return(gene_expr)




}
