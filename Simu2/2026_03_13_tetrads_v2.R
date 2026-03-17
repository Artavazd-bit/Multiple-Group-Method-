library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)
library(numDeriv)

source("2026_01_20_setup_v2.R")

tetrad <- function(data, model, latent1 = NULL, latent2 = NULL, scale = FALSE){
  model_df <- lavaanify(model)
  
  if(is.null(latent1) & is.null(latent2) ){
    latentvars <- unique(model_df$lhs[model_df$op == "=~"])
    listind1 <- list(model_df$rhs[model_df$lhs == latentvars[1] & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latentvars[2] & model_df$op == "=~"])
    if(length(latentvars) != 2){
      error("please select two latent variables")
    }
  }else{
    listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  }
  all_indicators <- unlist(list(listind1, listind2)) 
  
  subset_data <- data[, all_indicators]
  if(scale == FALSE){
    cor_subset_data <- cov(subset_data)
  } else { 
    cor_subset_data <- cor(subset_data)
  }
  
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row 
                  %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind2)] <- "het"
  
  mono <- cor_values[cor_values$type %in% c("mono1", "mono2"), ]
  cor_sym <- cor_values %>%
    bind_rows(
      cor_values %>% rename(col = row, row = col)
    )
  
  out <- mono %>% mutate(id = row_number()) %>% 
    cross_join(., ., suffix = c("_1", "_2")) %>% 
    filter(id_1 < id_2) %>% 
    filter(type_1 != type_2) %>% 
    left_join(., cor_sym, by = c("row_1" = "row", "col_2" = "col")) %>% 
    left_join(., cor_sym, by = c("col_1" = "col", "row_2" = "row"))
  
  out$tetrad <- out$val_1 * out$val_2 - out$val.x*out$val.y
  
  out <- out %>% select(- c("type.x", "type.y", "id_1", "id_2"))
  
  return(out)
}

tetrad_inner <- function(x, cor_subset_data, listind1, listind2){
  cor_subset_data[lower.tri(cor_subset_data)] <- x
  cor_subset_data[upper.tri(cor_subset_data)] <- t(cor_subset_data)[upper.tri(cor_subset_data)]
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row 
                  %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind2)] <- "het"
  
  mono <- cor_values[cor_values$type %in% c("mono1", "mono2"), ]
  cor_sym <- cor_values %>%
    bind_rows(
      cor_values %>% rename(col = row, row = col)
    )
  
  out <- mono %>% mutate(id = row_number()) %>% 
    cross_join(., ., suffix = c("_1", "_2")) %>% 
    filter(id_1 < id_2) %>% 
    filter(type_1 != type_2) %>% 
    left_join(., cor_sym, by = c("row_1" = "row", "col_2" = "col")) %>% 
    left_join(., cor_sym, by = c("col_1" = "col", "row_2" = "row"))
  
  out$tetrad <- out$val_1 * out$val_2 - out$val.x*out$val.y
 
 return(out$tetrad)
}


tetrad_outer <- function(data, model, scale = FALSE, latent1 = NULL, latent2 = NULL){
  model_df <- lavaanify(model)
  
  if(is.null(latent1) & is.null(latent2) ){
    latentvars <- unique(model_df$lhs[model_df$op == "=~"])
    listind1 <- list(model_df$rhs[model_df$lhs == latentvars[1] & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latentvars[2] & model_df$op == "=~"])
    if(length(latentvars) != 2){
      error("please select two latent variables")
    }
  }else{
    listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  }
  all_indicators <- unlist(list(listind1, listind2)) 
  
  subset_data <- data[, all_indicators]
  if(scale == FALSE){
    cor_subset_data <- cov(subset_data)
  } else { 
    cor_subset_data <- cor(subset_data)
  }
  values_lt <- cor_subset_data[lower.tri(cor_subset_data)]
  tetrads <- tetrad_inner(values_lt, cor_subset_data, listind1, listind2)
  grad_tetrads <- jacobian(tetrad_inner, x = values_lt, cor_subset_data = cor_subset_data, 
                           listind1 = listind1, listind2 = listind2)
  
  return(list(tetrads = tetrads, grad_tetrads = grad_tetrads))
}

out <- tetrad_outer(data_1, model_unconstrained)

covcov <- calcovcov(data_1)

covtetrads <-  out$grad_tetrads %*% (covcov / nrow(data_1) )  %*% t(out$grad_tetrads)
T <- nrow(data_1)* t(out$tetrads) %*% solve(covtetrads) %*% out$tetrads

pchisq(T, 9, lower.tail = FALSE)
