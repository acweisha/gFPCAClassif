#####
#
#Wrapping Functions for the gFPCAClassif package
#Author: Anthony Weishampel
#Date Updated: 11/29/2021
#
######



#' Function for modeling single-level fpca
#' @name gsFPCA
#' @param X_dat_s N x m matrix of binary data
#' @param Ys N long vector of groups
#' @param covariates N x Q data frame of covariates
#' @param pve Proportion of variation explained
#' @param Kb number of basis functions used when estimating latent trajectories
#' @param bs0 type of basis functions used when estimating latent trajectories
#' @param num_knots number of knots in basis functions
#' @return list of information required to build the model and predict new groups
gsFPCA <- function(X_dat_s, Ys, covariates = NA, pve = 0.95, Kb = 10, num_knots = 10, bs0 = "cr"){

  D = dim(X_dat_s)[2]
  N = dim(X_dat_s)[1]
  tt=seq(0,1, len=D)
  k= Kb
  Ys_train = Ys
  J = num_knots
  static_covariates = covariates

  if(!is.na(static_covariates)[1]){
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
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, X_dat_s, tt, k=k, bs0 = bs0))))

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
  return_vals$mu_t = mu_t_hat

  return(gsFPCA.model = return_vals)

}




#' Function for predicting new groups of based on single realizations of binary-valued functional data
#' @name gsFPCA_predict
#' @param gsFPCA.model trained gsFPCA model
#' @param X_dat_s_new N_new x m matrix of binary data
#' @param covariates_new N_new x Q data frame of covariates
#' @return Predicted new groups for the N_new users
gsFPCA_predict <- function(gsFPCA.model, X_dat_s_new, covariates_new = NA){

  X_dat_s_test = X_dat_s_new
  D = dim(X_dat_s_test)[2]
  N_test = dim(X_dat_s_test)[1]
  tt=seq(0,1, len=D)
  static_covariates_test = covariates_new


  if(!is.na(static_covariates_test)[1]){
    if((N_test != (dim(static_covariates_test)[1]))){
      stop("Dimensions of Covariates and Binary Curves do not match")
    }
  }

  scores_train = gsFPCA.model$scores_train
  eigen_funcs1 = gsFPCA.model$eigen_funcs
  eigen_vals1 = gsFPCA.model$eigen_vals
  static_covariates = gsFPCA.model$static_covariates
  Ys_train  = gsFPCA.model$classes
  mu_t_hat = gsFPCA.model$mu_t

  #if vector
  if(is.null(dim(eigen_funcs1))){
    matrix(eigen_funcs1, ncol = 1)
  }

  if(D != dim(eigen_funcs1)[1]){
    stop("Dimensions of new curves do not match eigenfunctions length")
  }

  prior_scales_test = eigen_vals1


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

  if(is.na(static_covariates)[1]){
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

    if(dim(cat_covariates_train)[2]==0){
      cat_covariates_train  = NA
      cat_covariates_test  = NA
      #need to update for categorical variables
      guess = nb_updated_grid_scores_only(scores_train2,
                                              Ys_train,
                                              prior_g, scores_test2,
                                              min.h = 0.3, max.h = 1.5)

    }else{
      #need to update for categorical variables
      guess = nb_updated_grid_scores_cat_only(scores_train2, cat_covariates_train,
                                              Ys_train,
                                              prior_g, scores_test2, cat_covariates_test,
                                              min.h = 0.3, max.h = 1.5)
    }




  }


  return(new_groups = guess)

}










#'  Function for fitting the gmFPCA
#' @name gMFPCA
#' @param X_dat_m N*J x m matrix of binary data
#' @param Ys N long vector of responses
#' @param J Number of realizations per individuals
#' @param N number of individuals
#' @param gAR variable to include generalized autoregressive structure
#' @param covariates N x Q Data frame of covariates
#' @param pve1 Proportion of variation explained in the first level eigenfunction
#' @param pve2 Proportion of variation explained in the first level eigenfunction
#' @param Kb number of basis functions in smoothing covariance functions
#' @param bs0 type of basis functions when smoothing covariance functions
#' @param approximation choice of which approximation to use
#' @param gar_covariates covariates to include in the gAR model
#' @param q Lag in the generalized Auto regressive model
#' @return list of information required to build the gmFPCA model and predict new groups
gMFPCA <- function(X_dat_m, Ys, J, N, covariates = NA, gAR = F, pve1 = 0.95,
                   pve2 = 0.95, Kb = 5, q = 3, approximation = "linear", gar_covariates = NA, bs0 = "cr"){


  #X_dat_m = t(matrix(t(X_dat_train), ncol = N*J))

  Ys_train = Ys
  static_covariates = covariates
  k = Kb
  Js = J

  D = dim(X_dat_m)[2]

  if(length(Js)>1){
    N = length(Js)
    #set check if Js and N do not match up
  }else{
    N = dim(X_dat_m)[1]/Js
    Js = rep(Js, N)
  }

  tt=seq(0,1, len=D)


  if(!is.na(static_covariates)[1]){
    if((N != (dim(static_covariates)[1])) ){
      stop("Dimensions of Covariates and Binary Curves do not match")
    }
  }
  if(N != length(Ys_train)){
    stop("Dimensions of Covariates and Binary Curves do not match")
  }

  J = Js[1]

  posting_days = 1-(rowSums(X_dat_m)==0)
  s_mat_train = t(matrix(as.numeric(matrix(posting_days, nrow = J)), nrow = J))
  Js_s_train = rowSums(s_mat_train)

  #Get the parsimonious distribution
  if(approximation == "linear"){
    cur.train = multilevel_linear_fpca(X_dat_m, J,
                                              pve1 = pve1, pve2 = pve2, k = k, bs0 = bs0)
  }else{
    cur.train = multilevel_exponential_fpca(X_dat_m, J,
                                             pve1 = pve1, pve2 = pve2, k = k, bs0 = bs0)
  }

  mu_t_hat = cur.train$mu_hat
  eigen_vals1 = cur.train$eigen_vals1
  eigen_funcs1 = cur.train$eigen_funcs1
  eigen_vals2 = cur.train$eigen_vals2
  eigen_funcs2 = cur.train$eigen_funcs2

  posting_days = (rowSums(X_dat_m)>1)
  s_mat_hat_train = t(matrix(as.numeric(matrix(posting_days, nrow = J)), nrow = J))

  scores_train = estimate_scores(X_dat_m, s_mat = s_mat_train, I=N,  J=J,
                                 eigen_vals1, eigen_vals2,
                                 eigen_funcs1, eigen_funcs2, mu_t_hat)

  return_vals = list( )

  return_vals$scores_train = scores_train
  return_vals$eigen_funcs1 = eigen_funcs1
  return_vals$eigen_vals1 = eigen_vals1
  return_vals$eigen_funcs2 = eigen_funcs2
  return_vals$eigen_vals2 = eigen_vals2
  return_vals$static_covariates = static_covariates
  return_vals$classes = Ys_train
  return_vals$mu_t = mu_t_hat
  return_vals$gAR = gAR
  return_vals$J = J

  if(gAR){

    gar_models_ls = list()

    ng = length(unique(Ys_train))

    for(l in 1:ng){

      gar_models_ls[[l]] = fit_ajs_model(l, q, s_mat_hat_train, classes = Ys_train, static_train = gar_covariates)

    }

    return_vals$gar_models_ls = gar_models_ls
    return_vals$s_mat_train = s_mat_train
    return_vals$q = q
    return_vals$gar_covariates = gar_covariates

  }

  return(gsFPCA.model = return_vals)

}








#' Function for predicting the groups for new gmFPCA Classifier
#' @name gmFPCA_predict
#' @param X_dat_m_new N_new x m matrix of binary data
#' @param gmFPCA.model Trained gmFPCA
#' @param covariates_new N_new x Q data frame of covariates
#' @param gar_covariates_new N_new X Q_g data frame of covariates for gAR model
#' @return predicted group values for the N_new users
gmFPCA_predict <- function(gmFPCA.model, X_dat_m_new, covariates_new = NA, gar_covariates_new = NA){

  scores_train = gmFPCA.model$scores_train
  eigen_funcs1 = gmFPCA.model$eigen_funcs1
  eigen_vals1 = gmFPCA.model$eigen_vals1
  eigen_funcs2 = gmFPCA.model$eigen_funcs2
  eigen_vals2 = gmFPCA.model$eigen_vals2
  static_covariates = gmFPCA.model$static_covariates
  Ys_train  = gmFPCA.model$classes
  mu_t_hat = gmFPCA.model$mu_t
  gAR = gmFPCA.model$gAR
  J = gmFPCA.model$J
  static_covariates_test = covariates_new
  X_dat_m_test = X_dat_m_new
  gar_covariates_test = gar_covariates_new

  D = dim(X_dat_m_new)[2]
  #N_test = dim(X_dat_new)[1]
  #X_dat_m_test = t(matrix(t(X_dat_new), nrow = N*J))

  X_dat_m_test = X_dat_m_new

  if(gAR){

    gar_models_ls = gmFPCA.model$gar_models_ls
    s_mat_train  = gmFPCA.model$s_mat_train
    gar_covariates = gmFPCA.model$gar_covariates
    q = gmFPCA.model$q

  }

  #if vector
  if(is.null(dim(eigen_funcs1))){
    matrix(eigen_funcs1, ncol = 1)
  }

  if(D != dim(eigen_funcs1)[1]){
    stop("Dimensions of new curves do not match eigenfunctions length")
  }

  D = dim(X_dat_m_test)[2]
  N_test = dim(X_dat_m_test)[1]/J
  tt=seq(0,1, len=D)


  if(!is.na(static_covariates_test)[1]){
    if((N_test != (dim(static_covariates_test)[1]))){
      stop("Dimensions of Covariates and Binary Curves do not match")
    }
  }

  J_test = J
  posting_days = 1-(rowSums(X_dat_m_test)==0)
  s_mat_test = t(matrix(as.numeric(matrix(posting_days, nrow = J_test)), nrow = J_test))
  Js_s_test = rowSums(s_mat_test)

  #X_dat_s_test = t(matrix(c(t(X_dat_m_test)), ncol = N_test))

  #estimate scores testing set
  scores_test=estimate_scores(X_dat_m_test, s_mat = s_mat_test, I=N_test, J=J_test,
                              eigen_vals1, eigen_vals2,
                              eigen_funcs1, eigen_funcs2, mu_t_hat)

  #step 4
  prior_g = c(table(Ys_train)/length(Ys_train))


  if(!gAR){
    if(is.na(static_covariates)[1]){
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

      guess = nb_updated_grid_scores_cat_only(scores_train2, cat_covariates_train,
                                              Ys_train,
                                              prior_g, scores_test2, cat_covariates_test,
                                              min.h = 0.3, max.h = 1.5)

    }
  }else{

    if(is.na(static_covariates)[1]){
      guess = nb_updated_grid(scores = scores_train, classes = Ys_train,
                              prior_g = c(table(Ys_train)/length(Ys_train)),
                              scores_test =  scores_test,
                              s_mat_hat_test =  s_mat_test,
                              s_mat_hat_train =  s_mat_train,
                              P_max = q,
                              static_train = gar_covariates,
                              static_test = gar_covariates_test)
    }else{


      numeric_cols = which(sapply(static_covariates, is.numeric))

      cur.mat = data.matrix(static_covariates[,numeric_cols])
      scores_train2 = cbind(scores_train, cur.mat)

      cur.mat = data.matrix(static_covariates_test[,numeric_cols])
      scores_test2 = cbind(scores_test, cur.mat)

      #need to update the categorical data

      cat_covariates_train  = static_covariates[,-numeric_cols]
      cat_covariates_test  = static_covariates_test[,-numeric_cols]



      if(dim(cat_covariates_train)[2]==0){
        cat_covariates_train  = NA
        cat_covariates_test  = NA
        #need to update for categorical variables
        guess = nb_updated_grid(scores = scores_train2, classes = Ys_train,
                                prior_g = c(table(Ys_train)/length(Ys_train)),
                                scores_test =  scores_test2,
                                s_mat_hat_test =  s_mat_test,
                                s_mat_hat_train =  s_mat_train,
                                P_max = q,
                                static_train = gar_covariates,
                                static_test = gar_covariates_test)

      }else{
        #need to update for categorical variables
        guess = nb_updated_grid_cat(scores = scores_train2, classes = Ys_train,
                                    cat_covariates_train, cat_covariates_test,
                                    prior_g = c(table(Ys_train)/length(Ys_train)),
                                    scores_test =  scores_test2,
                                    s_mat_hat_test =  s_mat_test,
                                    s_mat_hat_train =  s_mat_train,
                                    P_max = q,
                                    static_train = gar_covariates,
                                    static_test = gar_covariates_test)
      }



      #update for categorical variables



        }

    }

  return(new_groups = guess)

}







#' The binary-valued functional data summarizing 2-weeks of posting data for the 400 accounts in the training set
#'
#' @name X_dat_train
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords Binary-Valued Functional Data
NULL


#' The binary-valued functional data summarizing 2-weeks of posting data for the 100 accounts in the testing set
#'
#' @name X_dat_test
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords Binary-Valued Functional Data
NULL


#' The multilevel binary-valued functional data summarizing 2-weeks of posting data for the 400 accounts in the training set
#'
#' @name X_dat_m_train
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords Binary-Valued Functional Data
NULL

#' The multilevel binary-valued functional data summarizing 2-weeks of posting data for the 100 accounts in the testing set
#'
#' @name X_dat_m_test
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords Binary-Valued Functional Data
NULL

#'  Data about the 400 accounts in the training set
#' @description For each individual (id), we have the group of the account (group) bot (group = 2) and genuine (group = 1). Additionally for each individual, we have the number of accounts which the individual follows (friends\_count), the number of accounts which follow the individual (followers\_count), and the date when the account was created (created\_at).
#' @name acc_data_train
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords  Functional Data
NULL

#'  Data about the 100 accounts in the testing set
#' @description For each individual (id), we have the group of the account (group) bot (group = 2) and genuine (group = 1). Additionally for each individual, we have the number of accounts which the individual follows (friends\_count), the number of accounts which follow the individual (followers\_count), and the date when the account was created (created\_at).
#' @name acc_data_test
#' @docType data
#' @author Anthony Weishampel \email{acweisha@@ncsu.edu}
#' @keywords  Functional Data
NULL
