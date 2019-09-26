project_eqtl_to_flash_factor <- function(eqtl_mat, factor_mat, full = F) {
  # input:
  #   eqtl_mat: n-by-k matrix, where each row is a vector of eqtl effect size across k tissues (note that NaN entry will be treated as zero!)
  #   factor_mat: k-by-m matrix, where each column is a flash factor
  #   full: indicator variable telling whether returns the full projection onto all factors or not
  # output: a list with one or three (if full is TRUE) objects
  #    1. best_membership: a data.frame with two columns
  #       membership: a length-n vector, indicating the column index of the closest flash factor of each eqtl vector
  #       fraction_of_explained: a length-n vector, where each entry f that 1 - f = \| eqtl_residual \|_2^2 / \| eqtl \|_2^2 indicating to what extend closest flash factor captures the cross-tissue pattern of eqtl vector
  #    2. full_projection: n-by-m matrix with each entry ij to be the 2-norm of eqtl i projected onto factor j
  #    3. full_variance_explained: n-by-m matrix with each entry ij to be f value for eqtl i projected onto factor j
  
  eqtl_mat[is.nan(eqtl_mat)] <- 0
  eqtl_mat[is.na(eqtl_mat)] <- 0
  norm2 <- apply(factor_mat, 2, norm, type = '2')
  factor_mat <- factor_mat %*% diag(1 / norm2)
  projection <- eqtl_mat %*% factor_mat
  membership <- unlist(apply(abs(projection), 1, which.max))
  inner_product <- apply(abs(projection), 1, max)
  explained_var <- inner_product ^ 2
  original_var <- rowSums(eqtl_mat ^ 2)
  fraction_of_explain <- explained_var / original_var
  df <- data.frame(membership = membership, fraction_of_explain = fraction_of_explain)
  if(isTRUE(full)) {
    obj <- list(best_membership = df, full_projection = abs(projection), full_variance_explained = sweep(projection ^ 2, 1, original_var, "/"))
    return(obj)
  } else {
    return(list(best_membership = df))
  }
}

# function to obtain membership from projection result
get_membership_from_pve_mat <- function(mat, threshold = 0.2) {
  midx <- apply(mat, 1, function(x){
    idxmax <- which.max(x)
    if(length(idxmax) != 1) {
      return(0)
    } else {
      return(idxmax)
    }})
  mval <- apply(mat, 1, function(x){
    idxmax <- which.max(x)
    if(length(idxmax) != 1) {
      return(-1)
    } else {
      return(max(x))
    }})
  mem <- midx
  mem[mval < threshold] <- NA
  return(mem)
}

