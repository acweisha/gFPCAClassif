#####
#
#Wrapping Functions for the gFPCAClassif package
#Author: Anthony Weishampel
#Date Updated: 10/12/2021
#
######



###
# wrapping function for single level fpca
#' @curves N x m matrix of binary data
#' @Ys_train N long vector of responses
#' @static_covariates N x Q dataframe of covariates
#' @pve Proportion of variation explained
#' @return list of information required to build the model and predict new groups
#' @export
###

gsFPCA <- function(X_dat_s, Ys_train, static_covariates = NA, pve = 0.95, k = NA, J = 14){

  D = dim(X_dat_s)[2]
  N = dim(X_dat_s)[1]
  tt=seq(0,1, len=D)


  if(!is.na(static_covariates)){
    if((N != (dim(static_covariates)[1])) ){
      stop("Dimensions of Covariates and Binary Curves do not match")
    }
  }
  if(N != length(Ys_train)){
    stop("Dimensions of Covariates and Binary Curves do not match")
  }

  ##
  #Step 1 of the proposed method
  ##
  vec = matrix(1:(N), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, X_dat_s, tt, k=J))))

  ##
  #Step 2 of the proposed method
  ##
  fpca.cur2 = refund::fpca.face(smoothed_x, pve = pve, p=3, m=2, knots = J) #lambda selected via grid search optim, #p=degree of splines
  get_multiplier = 1/D
  fpca.cur = fpca.cur2
  #correct eigenfunctions
  fpca.cur$efunctions = fpca.cur2$efunctions/sqrt(get_multiplier)
  #correct eigenvalues
  fpca.cur$evalues = fpca.cur2$evalues*get_multiplier
  #correct scores
  fpca.cur$scores = fpca.cur2$scores*sqrt(get_multiplier)

  ##
  #STEP 3:
  ##
  fit = list(mu = fpca.cur$mu,
             evalues = fpca.cur$evalues,
             efunctions = fpca.cur$efunctions)

  mu_t_hat = fit$mu
  eigen_vals1 = fit$evalues
  eigen_funcs1 = fit$efunctions

  #data frame used in bayesglm
  dta = data.frame(index = rep(tt, N),
                   value = c(t(X_dat_s)),
                   id = rep(1:N, each = D))

  npc = length(eigen_vals1)
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }

  #assign names to data frame
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  #repeat mean function in data frame once per user
  dta$mu = rep(mu_t_hat , N)

  #get formula for glm
  glm_structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm_structure , sep="")
  #set scale for the glm
  prior_scales_test = eigen_vals1

  #Estimate the Scores for the training set
  vec = matrix(1:N, ncol = 1)
  #vec = matrix(vec[users_to_keep_train,], ncol = 1)
  scores_train = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))

  #Step 3 for Testing Data
  #Get the scores for the testing set
  return_vals = list( )

  return_vals$scores_train = scores_train
  return_vals$eigen_funcs = eigen_funcs1
  return_vals$eigen_vals = eigen_vals1
  return_vals$static_covariates = static_covariates
  return_vals$classes = Ys_train

  return(gsFPCA.model = return_vals)

}




###
# wrapping function for single level fpca
#' @curves N x m matrix of binary data
#' @static_covariates N x Q dataframe of covariates
#' @pve Proportion of varation explained
#' @return A matrix of the infile
#' @export
###

gsFPCA_predict <- function(gsFPCA.model, X_dat_s_test, static_covariates_test = NA){

  D = dim(X_dat_s_test)[2]
  N_test = dim(X_dat_s_test)[1]
  tt=seq(0,1, len=D)

  if(!is.na(static_covariates_test)){
    if((N_test != (dim(static_covariates_test)[1]))){
      stop("Dimensions of Covariates and Binary Curves do not match")
    }
  }

  scores_train = gsFPCA.model$scores_train
  eigen_funcs1 = gsFPCA.model$eigen_funcs
  eigen_vals1 = gsFPCA.model$eigen_vals
  static_covariates = gsFPCA.model$static_covariates
  Ys_train  = gsFPCA.model$classes

  #if vector
  if(is.null(dim(eigen_funcs1))){
    matrix(eigen_funcs1, ncol = 1)
  }

  if(D != dim(eigen_funcs1)[1]){
    stop("Dimensions of new curves do not match eigenfunctions length")
  }

  #just like before define data frame
  dta = data.frame(index = rep(tt, N_test),
                   value = c(t(X_dat_s_test)),
                   id = rep(1:N_test, each = D))

  npc = length(eigen_vals1)

  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_test))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  dta$mu = rep(mu_t_hat , N_test)

  glm_structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm_structure , sep="")

  vec = matrix(1:N_test, ncol = 1)
  scores_test = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))

  #step 4
  #get propability of being in each group
  #Ys_train = Classes_train
  prior_g = c(table(Ys_train)/length(Ys_train))
  #run non parametric bayes classifier

  if(is.na(covariates_train)[1]){
    guess = nb_updated_grid_scores_only(scores_train,
                                        Ys_train,
                                        prior_g, scores_test,
                                        min.h = 0.3, max.h = 1.5)
  }else{

    numeric_cols = which(sapply(static_covariates, is.numeric))

    cur.mat = data.matrix(static_covariates[,numeric_cols])
    scores_train2 = cbind(scores_train, cur.mat)

    cur.mat = data.matrix(static_covariates_test[,numeric_cols])
    scores_test2 = cbind(scores_test, cur.mat)

    #need to update the categorical data

    cat_covariates_train  = static_covariates[,-numeric_cols]
    cat_covariates_test  = static_covariates_test[,-numeric_cols]


    #need to update for categorical variables

    guess = nb_updated_grid_scores_only(scores_train,
                                        Ys_train,
                                        prior_g, scores_test,
                                        min.h = 0.3, max.h = 1.5)

  }


  return(new_groups = guess)

}
