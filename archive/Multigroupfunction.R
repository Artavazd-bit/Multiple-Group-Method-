multigroup <- function(data , model, latent1, latent2){
  model_df <- lavaanify(model_est)
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  
  S <- cor(data[,c(unlist(listind1), unlist(listind2))])
  Sstar <- S
  
  corlist <- list(cor(data[,unlist(listind1)]), cor(data[,unlist(listind2)]) ) 
  communalities <- lapply(corlist, function(corlist) {
    sapply(1:nrow(corlist), function(x) {sum(combn(corlist[x,-x], m = 2, FUN = prod))/sum(corlist[-x,-x][upper.tri(corlist[-x,-x])])})
  })
  
  diag(Sstar) <- unlist(communalities)
  
  pattern_matrix <- S[,c(1,2)] * 0 
  colnames(pattern_matrix) <- c(latent1, latent2)
  pattern_matrix[unlist(listind1),latent1] <- 1 
  pattern_matrix[unlist(listind2),latent2] <- 1 
  pattern_matrix <- t(pattern_matrix)
  
  Fzero <- pattern_matrix %*% Sstar 
  Lzero <- Fzero %*% t(pattern_matrix)
  
  Dinvsqrt <- diag(1/sqrt(diag(diag(diag(Lzero)))))
  
  L <- Dinvsqrt %*% Lzero %*% Dinvsqrt
  
  Fmatrix <- t(Fzero) %*% Dinvsqrt 
  
  return(list(L, Fmatrix))
}

multigroup(data, model_est, "xi_1", "xi_2")
