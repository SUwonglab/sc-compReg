"compute_lambda" <- function(peak_O_path,
                             p_symbol_path,
                             e_symbol_path,
                             e_matrix_path,
                             pe_path,
                             k_components,
                             ...){
    UseMethod("compute_lambda")
}
