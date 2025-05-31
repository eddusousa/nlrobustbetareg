
library(betareg)

###############################################################################################################################

robustbetareg <- function (formula, 
                           data, 
                           functions = c("exp","log"), 
                           parameters = NULL,
                           alpha, 
                           type = c("LSMLE", "LMDPDE", "SMLE", "MDPDE", "MLE"), 
                           link = c("logit", "probit", "cloglog", "cauchit", "loglog"), 
                           link.phi = NULL, 
                           se.bootstrap = FALSE,
                           simulation = FALSE,
                           wald.test.bootstrap = FALSE,
                           control = robustbetareg.control(...), 
                           model = TRUE, 
                           ...)
{
  if( (type!="MLE")&(!is.null(alpha)) ){
    if( alpha==0 ){
      type <- "MLE"
      message("\nWARNING: Alpha = ",alpha," means MLE estimators. Using 'type=MLE' instead!")
    }
  }
  
  # Checking if type is MLE and alpha is not zero
  if( (type=="MLE")&(!is.null(alpha)) ){
    if(alpha!=0){
      alpha <- 0
      message("\nWARNING: For ",type," estimator, alpha = 0 will be used.")
    }
  }
  
  
  # Checking if se.bootstrap=TRUE for non linear models
  if( (!is.null(parameters)) & (se.bootstrap == FALSE) & (simulation==FALSE) ){
    message("\nWARNING: For non linear models, only standard errors via bootstrap are avaible! Using 'se.bootstrap = TRUE' instead!\n")
    se.bootstrap <- TRUE
  } 
  
  # Checking if wald.test.bootstrap=TRUE and se.bootstrap=FALSE
  if( (wald.test.bootstrap) & (se.bootstrap == FALSE) ){
    message("\nWARNING: For type-Wald bootstrap test you must use stabdard errors via bootstrap too! Using 'se.bootstrap = TRUE' instead!\n")
    se.bootstrap <- TRUE
  } 
  
  cl = match.call()
  type = match.arg(type)
  ocontrol = control
  if (missing(data)) {
    data = environment(formula)
  }
  if (!missing(alpha)) {
    control$alpha.optimal = FALSE
  }
  if (missing(alpha)) {
    alpha = NULL
  }

  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels = TRUE
  formula = Formula::as.Formula(formula)
  oformula = as.formula(formula)
  if (length(formula)[2L] < 2L) {
    formula = Formula::as.Formula(formula(formula), ~1)
    simple_formula = TRUE
  }
  else {
    if (length(formula)[2L] > 2L) {
      formula = Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf1 = model.frame(formula, data = data)
  y = model.response(mf1)
  x = model.matrix(formula, data = mf1, rhs = 1L)
  z = model.matrix(formula, data = mf1, rhs = 2L)
  if (simple_formula) {
    colnames(z)[1] = "(Phi)"
  }
  mt = terms(formula, data = data)
  mtX = terms(formula, data = data, rhs = 1L)
  mtZ = delete.response(terms(formula, data = data, rhs = 2L))
  if (length(y) < 1) {
    stop("empty model")
  }
  if (!(min(y) > 0 & max(y) < 1)) {
    stop("invalid dependent variable, all observations must be in (0, 1)")
  }
  if (!is.null(control$start) & (ncol(x) + ncol(z)) != length(control$start)) {
    stop("Invalid initial starting point")
  }
  if (!is.null(alpha)) {
    if (alpha < 0 || alpha > 1) {
      stop("invalid tuning constant, the value must be in [0, 1)")
    }
  }
  if (!is.null(link.phi)) {
    if (link.phi == "identity" & !simple_formula) {
      link.phi = "log"
      warning("Non suitable precision link function, log link used instead")
    }
  }
  else {
    link.phi = if (simple_formula) {
      "identity"
    }
    else "log"
  }
  
  link = match.arg(link)
  linkobj = set.link(link.mu = link, link.phi = link.phi)
  if ( (is.null(control$start)) & !is.null((parameters)) ) {
    # using a beta nonlinear regression for obtain de starting values for the optim funcion
    control$start = starting.points(formula = formula,
                                    data = data,
                                    y=y,
                                    x=x,
                                    z=z,
                                    type = type,
                                    parameters = parameters,
                                    functions = functions,
                                    link = link,
                                    link.phi = link.phi)

    # est.mle = suppressWarnings(betareg(oformula, data, link = link,
    #                                    link.phi = link.phi))
    # control$start = c(est.mle$coefficients$mean, est.mle$coefficients$precision)
    
  }else if( (is.null(control$start)) & is.null((parameters)) ){
    
    est.mle = suppressWarnings(betareg(oformula, data, link = link,
                                       link.phi = link.phi))
    control$start = c(est.mle$coefficients$mean, est.mle$coefficients$precision)
    
  }
  if (simple_formula) {
    control$start[length(control$start)] = linkobj$linkfun.phi$linkfun(control$start[length(control$start)])
  }
  ##############################################################################
  check.alpha.optimal <- FALSE
  if(type!="MLE"){
    if( (is.null(alpha))) {
      ########################################
      require(compiler)
      
      compiled_select_tuning_constant <- cmpfun(select_tuning_constant)
      
      enableJIT(3)
      ########################################
      control$alpha.optimal = FALSE
      output <- compiled_select_tuning_constant(y=y,
                                                x=x,
                                                z=z,
                                                type=type,
                                                link=link,
                                                link.phi=link.phi,
                                                start_theta=control$start,
                                                se.bootstrap=se.bootstrap,
                                                wald.test.bootstrap=FALSE,
                                                functions=functions,
                                                parameters=parameters,
                                                data=data,
                                                formula=formula)
      alpha <- output$alpha.optimal
      check.alpha.optimal <- TRUE
      if(alpha==0){
        type <- "MLE"
      }
    }
  }
  

  ##############################################################################
  if (type == "LMDPDE") {
    result = LMDPDE.fit(y, x, z, alpha = alpha, link = link,  type = type, se.bootstrap = se.bootstrap,
                        link.phi = link.phi, control = control)
  }
  ###########################################################################################################
  # nonlinearity avaible
  if (type == "LSMLE") {
    result = LSMLE.fit(y, x, z, alpha = alpha, link = link, type = type, se.bootstrap = se.bootstrap,
                       wald.test.bootstrap = wald.test.bootstrap, functions = functions, 
                       parameters = parameters, link.phi = link.phi, control = control)
  }
  ###########################################################################################################
  # nonlinearity avaible
  if (type == "SMLE") {
    result = SMLE.fit(y, x, z, alpha = alpha, link = link, type = type, se.bootstrap = se.bootstrap,
                      wald.test.bootstrap = wald.test.bootstrap, functions = functions, 
                      parameters = parameters, link.phi = link.phi, control = control)
  }
  if (type == "MDPDE") {
    result = MDPDE.fit(y, x, z, alpha = alpha, link = link,  type = type, se.bootstrap = se.bootstrap,
                       link.phi = link.phi, control = control)
  }
  ###########################################################################################################
  # nonlinearity avaible
  if (type == "MLE") {
    result = MLE.fit(y, x, z, link = link, type = type, se.bootstrap = se.bootstrap,
                     wald.test.bootstrap = wald.test.bootstrap, functions = functions, 
                     parameters = parameters, link.phi = link.phi, control = control,
                     data = data, formula = formula)
  }
  ###########################################################################################################
  if(exists("output")){
      result$message = output$message
      result$enpv = output$enpv
  }

  ###########################################################################################################
  result$y = y
  
  if (simple_formula) {
    precision.name = names(result$coefficients$precision)
    result$coefficients$precision = linkobj$linkfun.phi$inv.link(z %*% 
                                                                   result$coefficients$precision)[1]
    names(result$coefficients$precision) = precision.name
  }
  if (model) {
    result$model = list(mean = x, precision = z)
  }
  result$Optimal.Tuning <- check.alpha.optimal
  result$terms = list(mean = mtX, precision = mtZ, full = mt)
  result$call = cl
  result$data = mf1
  result$formula = as.formula(formula)
  nonlinear.args <- list()
  nonlinear.args$functions <- functions
  nonlinear.args$parameters <- parameters
  result$nonlinear.args <- nonlinear.args
  names(result$coefficients$mean) <- colnames(x)
  names(result$coefficients$precision) <- colnames(z)

  if(any(is.na(result$vcov))){
    message(paste("INFO: The standard-Error is unvailable"))
    result$message <- "Standard-Error is unavailable"
  }
  if( (!any(control$start == c(result$coefficients$mean, result$coefficients$precision)))==FALSE ){
    message(paste("INFO: The optimization process failed and the final coeficients are equals to the starting points!"))
  }
  
  gc()
  
  return(result)
}

######################################################################################################################

set.link <- function(link.mu = "logit", link.phi = "log") 
{
  switch(link.mu, logit = {
    linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(log(mu) - log(1 - mu))
    }
    d.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return((mu - mu^2)^(-1))
    }
    d2.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return((2 * mu - 1)/(mu * (1 - mu))^2)
    }
    inv.link <- function(eta) {
      return(as.numeric(pmax(pmin(exp(eta - Rmpfr::log1pexp(eta)), 
                                  1 - .Machine$double.eps), .Machine$double.eps)))
    }
  }, probit = {
    linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(qnorm(mu))
    }
    d.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return((dnorm(qnorm(mu)))^(-1))
    }
    d2.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(-ddnorm(qnorm(mu))/(dnorm(qnorm(mu)))^3)
    }
    inv.link <- function(eta) {
      return(as.numeric(pmax(pmin(pnorm(eta), 1 - .Machine$double.eps), 
                             .Machine$double.eps)))
    }
  }, cloglog = {
    linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(log(-log(1 - mu)))
    }
    d.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return((-(1 - mu) * log(1 - mu))^(-1))
    }
    d2.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(-(log(1 - mu) + 1)/((1 - mu) * log(1 - mu))^2)
    }
    inv.link <- function(eta) {
      return(as.numeric(pmax(pmin(1 - exp(-exp(-eta)), 
                                  1 - .Machine$double.eps), .Machine$double.eps)))
    }
  }, cauchit = {
    linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(tan(pi * (mu - 0.5)))
    }
    d.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(pi * pracma::sec(pi * (mu - 0.5))^2)
    }
    d2.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(2 * pi * tan(pi * (mu - 0.5)) * pracma::sec(pi * 
                                                           (mu - 0.5))^2)
    }
    inv.link <- function(eta) {
      return(as.numeric(pmax(pmin(0.5 + atan(eta)/pi, 1 - 
                                    .Machine$double.eps), .Machine$double.eps)))
    }
  }, loglog = {
    linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(-log(-log(mu)))
    }
    d.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return(-(mu * log(mu))^(-1))
    }
    d2.linkfun <- function(mu) {
      mu <- pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
      return((log(mu) + 1)/(mu * log(mu))^2)
    }
    inv.link <- function(eta) {
      return(as.numeric(pmax(pmin(exp(-exp(-eta)), 1 - 
                                    .Machine$double.eps), .Machine$double.eps)))
    }
  }, stop(gettextf("%s link.mu not recognised", sQuote(link.mu)), 
          domain = NA))
  Linkfun.Mu = list(linkfun = linkfun, d.linkfun = d.linkfun, 
                    d2.linkfun = d2.linkfun, inv.link = inv.link)
  switch(link.phi, log = {
    linkfun.phi <- function(phi) {
      phi <- pmax(phi, .Machine$double.eps)
      return(log(phi))
    }
    d.linkfun.phi <- function(phi) {
      phi <- pmax(phi, .Machine$double.eps)
      return((phi)^(-1))
    }
    d2.linkfun.phi <- function(phi) {
      phi <- pmax(phi, .Machine$double.eps)
      return(-(phi^(-2)))
    }
    inv.link.phi <- function(eta) {
      return(pmax(as.numeric(exp(eta)), .Machine$double.eps))
    }
  }, identity = {
    linkfun.phi <- function(phi) {
      phi <- pmax(phi, .Machine$double.eps)
      return(phi)
    }
    d.linkfun.phi <- function(phi) {
      return(rep(1, length(phi)))
    }
    d2.linkfun.phi <- function(phi) {
      return(rep(0, length(phi)))
    }
    inv.link.phi <- function(eta) {
      return(as.numeric(eta))
    }
  }, sqrt = {
    linkfun.phi <- function(phi) {
      phi <- pmax(phi, .Machine$double.eps)
      return(sqrt(phi))
    }
    d.linkfun.phi <- function(phi) {
      return((2 * sqrt(phi))^(-1))
    }
    d2.linkfun.phi <- function(phi) {
      return(-0.25 * (phi^(-3/2)))
    }
    inv.link.phi <- function(eta) {
      return(pmax(as.numeric(eta^2), .Machine$double.eps))
    }
  }, stop(gettextf("%s link.phi not recognised", sQuote(link.phi)), 
          domain = NA))
  Linkfun.Phi = list(linkfun = linkfun.phi, d.linkfun = d.linkfun.phi, 
                     d2.linkfun = d2.linkfun.phi, inv.link = inv.link.phi)
  linkobj = structure(list(linkfun.mu = Linkfun.Mu, linkfun.phi = Linkfun.Phi), 
                      name.link.mu = link.mu, name.link.phi = link.phi, class = "link-rbr")
  return(linkobj)
}

######################################################################################################################

robustbetareg.control <- function(start = NULL, alpha.optimal = TRUE, tolerance = 0.001, 
          maxit = 5000, L = 0.02, M = 3, ...) 
{
  return(list(start = start,
              alpha.optimal = alpha.optimal,
              tolerance = tolerance,
              maxit = maxit,
              L = L,
              M = M,
              alpha.max = 0.5))
}

######################################################################################################################

sweighted2_res <- function(mu_hat, phi_hat, y, X, linkobj) 
{
  n = length(mu_hat)
  m = length(phi_hat)
  if (m == 1) {
    phi_hat = rep(phi_hat, n)
  }
  d.link.mu = linkobj$linkfun.mu$d.linkfun(mu_hat)
  y_star = log(y) - log(1 - y)
  mu_star = digamma(mu_hat * phi_hat) - digamma((1 - mu_hat) * 
                                                  phi_hat)
  V_star = trigamma(mu_hat * phi_hat) + trigamma((1 - mu_hat) * 
                                                   phi_hat)
  W.PHI = diag(x = phi_hat * V_star * ((d.link.mu)^(-2)))
  H = sqrt(W.PHI) %*% X %*% solve(t(X) %*% W.PHI %*% X) %*% 
    t(X) %*% sqrt(W.PHI)
  nu = V_star * (1 - diag(H))
  diff = (y_star - mu_star)
  ri = diff/sqrt(nu)
  return(ri)
}


######################################################################################################################
# generate simulated data

FunSample_betareg <- function(n, X, Z,
                              beta, gamma,
                              functions, parameters,
                              link.mu='logit', link.phi='log'){
  
  
  #####################################################################################################
  # Applying splines in mu predictor
  if(!is.null(parameters)){
    eta <- apply_spline(vector = beta,
                        matrix = X,
                        functions = functions, 
                        parameters = parameters)
  }else{
    eta = as.vector(X%*%beta)
  }
  #####################################################################################################
  
  vartheta <- as.vector(Z%*%gamma) #Preditor linear de sigma
  
  # Funções de ligação
  if(link.mu == "logit"){
    mu <- exp(eta)/(1.0+exp(eta)) #logit
  } else if (link.mu == "probit"){
    mu <- pnorm(eta) #probit
  } else if (link.mu == "cloglog"){
    mu <- 1.0 - exp(-exp(eta)) # cloglog
  } else if (link.mu == "log"){
    mu <- exp(eta) #log
  } else if (link.mu == "loglog"){
    mu <- exp(-exp(-eta)) #loglog
  } else if (link.mu == "cauchit"){
    mu <- (pi^(-1))*atan(eta) + 0.5 #cauchit 
  } else {
    print("Funçao de ligação informada para a média mu é inválida! Verifique.")
    break
  }
  
  if(link.phi == "log"){
    phi <- exp(vartheta) #log
  } else if (link.phi == "identify") {
    phi <- vartheta #identify
  } else if (link.phi == "sqrt") {
    phi <- vartheta^2 #sqrt
  } else {
    print("Funçao de ligação informada para a precisão phi é inválida! Verifique.")
    break
  }
  
  a <- mu*phi  #Parâmetro de forma 1 da dist. beta 
  b <- (1-mu)*phi
  
  y <- rbeta(n=n, shape1=a, shape2=b) #Geração dos números aleatórios 
  
  results <- list(y=y, 
                  X=X, 
                  Z=Z,
                  mu=mu, 
                  phi=phi)
  
  return(results)
  
}
######################################################################################################################
# Calculate std error with bootstrap

bootstrap.std.error <- function(N = 200,
                                n,
                                alpha,
                                type,
                                beta,
                                gamma,
                                functions, 
                                parameters,
                                start_theta,
                                data = data,
                                y,
                                x,
                                z,
                                link,
                                link.phi,
                                sample_function = FunSample_betareg){
  
  message(paste("INFO: Standard error being obtained via bootstrap for", type, "estimators. This may take a while..."))
  
  #Leitura das bibliotecas que farão a paralelização
  library(foreach)
  library(doParallel)
  
  # Objeto para armazenar os valores estimados dos betas
  beta.hat <- lapply(rep(list(0), length(beta)), list)
  
  # Objeto para armazenar os valores estimados dos gammas
  gamma.hat <- lapply(rep(list(0), length(gamma)), list)
  
  # Objeto para armazenar os erros padrões
  standard.error.beta  <- rep(0, length(beta))
  standard.error.gamma  <- rep(0, length(gamma))
  
  alg1 <- list(list())

  # Deteccao da quantidade de nucleos da CPU e registro do backend p/ paralelizacao 
  numCores <- detectCores()
  registerDoParallel(numCores)

  # Algoritmo para geracao das amostras com paralelização
  alg1 <- foreach(j = 1:N,
                  .verbose = FALSE,
                  .export = c("FunSample_betareg",
                              "apply_spline",
                              "Psi_SMLE",
                              "L_q",
                              "L_alpha_R",
                              "Psi_Beta_LSMLE",
                              "Psi_Beta_SMLE",
                              "Psi_Gamma_LSMLE",
                              "Psi_Gamma_SMLE",
                              "Psi_LSMLE",
                              "Loglik_MLE",
                              "Psi_MLE",
                              "apply_derivative",
                              "set.link",
                              "dEGB"),
                  .packages = "base") %dopar% {
                    
     
    # Geracao da amostra aleatoria
    FunSample_betareg(n=n, 
                      X=x, 
                      Z=z,
                      beta = beta,
                      gamma = gamma,
                      functions = functions,
                      parameters = parameters,
                      link.mu=link, 
                      link.phi=link.phi) -> dados
    
    # Calculo dos parametros
    #####################################################################################################
    if(type=="SMLE"){
      stats::optim(par = start_theta,
                   fn = L_q,
                   gr = Psi_SMLE,
                   y = dados$y,
                   X = as.matrix(x),
                   Z = as.matrix(z), 
                   alpha = alpha,
                   link_mu = link,
                   link_phi = link.phi,
                   functions = functions,
                   parameters = parameters,
                   apply_spline = apply_spline,
                   apply_derivative = apply_derivative,
                   control = list(fnscale = -1))
    }else if(type=="LSMLE"){
      stats::optim(par = start_theta,
                   fn = L_alpha_R, 
                   gr = Psi_LSMLE,
                   y = dados$y,
                   X = as.matrix(x),
                   Z = as.matrix(z), 
                   alpha = alpha,
                   link_mu = link,
                   link_phi = link.phi,
                   functions = functions,
                   parameters = parameters,
                   apply_spline = apply_spline,
                   apply_derivative = apply_derivative,
                   control = list(fnscale = -1, 
                                  maxit = 10000))
    }else if(type=="MLE"){
      stats::optim(par = start_theta,
                   fn = Loglik_MLE,
                   gr = Psi_MLE,
                   y = dados$y,
                   X = as.matrix(x),
                   Z = as.matrix(z), 
                   link_mu = link, 
                   link_phi = link.phi,
                   functions = functions, 
                   parameters = parameters,
                   #data = data,
                   apply_spline = apply_spline, 
                   apply_derivative = apply_derivative,
                   control = list(fnscale = -1))

    }
    

      #####################################################################################################
                  }

  closeAllConnections()
  
  # Armazena os valores das estimativas dos coeficientes
  for(k in 1:N){
    for (m in 1:length(beta)) {
      beta.hat[[m]][k] <- alg1[[k]]$par[[m]]
    }
    for (m in 1:length(gamma)) {
      gamma.hat[[m]][k] <- alg1[[k]]$par[[m+length(beta)]]
    }
  }
  
  # Calcula o erro padrão amostral
  bootstrap.std.error.beta <- c()
  for (m in 1:length(beta)) {
    bootstrap.std.error.beta[m] = sd(unlist(beta.hat[[m]]))
  }
  bootstrap.std.error.gamma <- c()
  for (m in 1:length(gamma)) {
    bootstrap.std.error.gamma[m] = sd(unlist(gamma.hat[[m]]))
  }

  std.bootstrap.theta <- c(bootstrap.std.error.beta, bootstrap.std.error.gamma)

  vcov.bootstrap.theta <- diag(std.bootstrap.theta^2)

  message(paste("INFO: Standard error calculation completed!"))
  

 return(list(std.bootstrap.theta = std.bootstrap.theta,
             vcov.bootstrap.theta = vcov.bootstrap.theta))
}

######################################################################################################################
# Especific summarys
summary_nlrobustbetareg <- function (object){
  
  if(is.null(object$boot.wald.type.test)){
    stop("Summary not available for this object!! Hypotesys type-Wald test via bootstrap not found!!")
  }
  #type <- match.arg(type, c("sweighted2", "pearson", 
  #                          "weighted", "sweighted", "sweighted.gamma", 
  #                          "sweighted2.gamma", "combined", "combined.projection"))
  #object$residuals = residuals(object, type = type)
  object$residuals.type <- "quantile"
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$precision)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- data.frame(cf, se, (cf/se)^2, rstatix::p_format(object$boot.wald.type.test$pvalue,
                                                        accuracy = 1/object$boot.wald.type.test$valid.w.stat) )
  colnames(cf) <- c("Estimate", "Std. Error", "W value", 
                    "Pr(>W)")
  #rownames(cf) <- rownames(object$vcov)
  cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], 
             precision = cf[seq.int(length.out = m) + k, , drop = FALSE])

  rownames(cf$mean) <- row.names(object$vcov)[1:length(object$coefficients$mean)]
  rownames(cf$precision) <- row.names(object$vcov)[(length(object$coefficients$mean)+1):(length(object$coefficients$mean)+length(object$coefficients$precision))]
  object$coefficients <- cf
  #class(object) <- "summary.robustbetareg"
  return(object)
}


######################################################################################################################
# Starting points
starting.points = function(formula,
                           data,
                           y,
                           x,
                           z,
                           type,
                           parameters,
                           functions,
                           link = link,
                           link.phi = link.phi){
  # Using a robust start value for the parameters associated with the mean submodel
  # Generating the formulae to use in this function 
  # browser()
  if(length(formula)[2L] > 1L){
    mu_elements <- unlist(strsplit(as.character(formula[[3]][[2]]), " "))[!(unlist(strsplit(as.character(formula[[3]][[2]]), " ")) %in% "+")]
  }else{
    mu_elements <- unlist(strsplit(as.character(formula[[3]]), " "))[!(unlist(strsplit(as.character(formula[[3]]), " ")) %in% "+")]
  }
  
  if(link == "logit"){
    form <- paste("log(",formula[[2]],
                  "/(1-",formula[[2]],"))"
                  ,formula[[1]])
  }else if(link == "probit"){
    form <- paste("qnorm(",formula[[2]],")"
                  ,formula[[1]])
  }else if(link == "cloglog"){
    form <- paste("log(-log(1-",formula[[2]],"))"
                  ,formula[[1]])
  }else if(link == "loglog"){
    form <- paste("-log(-log(",formula[[2]],"))"
                  ,formula[[1]])
  }else if(link == "cauchit"){
    form <- paste("tan(pi*(",formula[[2]],"-0.5))"
                  ,formula[[1]])
  }else{
    stop("\nLink function for mean not availble for starting points!")
  }
  beta_names = c()
  for(i in 1:(length(mu_elements)+1)){
    if( (i==1)&(!(i %in% parameters)) ){
      form <- paste(form, "beta", i-1)
      beta_names[i] <- gsub(" ","",paste("beta", i-1))
    }else if((i==1)&((i %in% parameters))){
      form <- paste(form, functions[match(i,parameters)], "(","beta", i-1, ")")
      beta_names[i] <- gsub(" ","",paste("beta", i-1))
    }
    if((i!=1)&((i %in% parameters))){
      form <- paste(form, "+",functions[match(i,parameters)], "(", "beta", i-1, "*", mu_elements[i-1], ")")
      beta_names[i] <- gsub(" ","",paste("beta", i-1))
    }else if((i!=1)&(!(i %in% parameters[i-1]))){
      form <- paste(form, "+" ,"beta", i-1, "*", mu_elements[i-1])
      beta_names[i] <- gsub(" ","",paste("beta", i-1))
    }
  }
  form_final_mu <- stats::formula(gsub(" ","",as.character(form)))
  
  # beta regression to use as a starting point to robust nonlinear regression
  
  beta_pre_start <- as.list(rep(1, length(beta_names)))
  
  names(beta_pre_start) <- beta_names
  
  suppressWarnings(betareg::betareg(formula = formula,
                                    data = data,
                                    link = link,
                                    link.phi = link.phi)) -> start_betareg
  
  beta_exception <- as.list(start_betareg$coefficients$mean)
  
  names(beta_exception) <- beta_names
  
  tryCatch(as.list(coef(suppressWarnings(nls(formula = form_final_mu,
                                             data = data,
                                             algorithm = "port",
                                             start = beta_pre_start,
                                             trace = FALSE)))), 
           error = function(e) {
             return(beta_exception)
           }) -> start_nls
  
  
  
  midStart <- unname(c(unlist(start_nls), start_betareg$coefficients$precision))
  
  if(type=="MLE"){
    
    thetaStart <- midStart
    
  }else{
    
    est.mle <- list()
    est.mle = tryCatch(MLE.fit(y=y, 
                               x=as.matrix(x), 
                               z=as.matrix(z), 
                               link = link, 
                               link.phi = link.phi,
                               type="MLE",
                               se.bootstrap = FALSE,
                               wald.test.bootstrap = FALSE, 
                               functions = functions, 
                               parameters = parameters,
                               data = data,
                               formula = formula,
                               control = robustbetareg.control(start = midStart)), 
                       error = function(e) {
                         est.mle$converged <- FALSE
                         return(midStart)
                       })
    
    thetaStart <- c(est.mle$coefficients$mean, est.mle$coefficients$precision)
  }
  
  return(thetaStart)
}

######################################################################################################################
# Hypotesis wald test via bootstrap

exp_bootstrap.wald.test <- function(B2 = 200, #First level (number of w statistic replicas)
                                B1 = 50, # Second level (for the bootstrap sd)
                                n,
                                alpha,
                                beta,
                                gamma,
                                X,
                                Y,
                                Z,
                                type,
                                formula,
                                data.names,
                                std.bootstrap.theta,
                                functions, 
                                parameters,
                                start_theta,
                                link,
                                link.phi
){
  
  closeAllConnections()
  
  start_time <- Sys.time()
  
  #Leitura das bibliotecas que farão a paralelização
  library(foreach)
  library(doParallel)
  
  message(paste("INFO: Now starting to generate simulated distribution for the type-wald test statistics. It's a good time to grab some coffee..."))
  
  # Tratamento do ojeto que carrega os erros padrões
  std.bootstrap.theta <- c(as.numeric(std.bootstrap.theta$se.mean),as.numeric(std.bootstrap.theta$se.precision))
  
  # Deteccao da quantidade de nucleos da CPU e registro do backend p/ paralelizacao 
  numCores <- detectCores()
  registerDoParallel(numCores)
  
  # Objeto para armazenar os valores estimados dos betas e gammas
  beta_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  gamma_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  # Objeto para armazenar os valores estimados dos erros-padrão para cada replica
  bootstrap.std.error.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  bootstrap.std.error.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  bootstrap.std.error.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  bootstrap.std.error.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Objeto para armazenar os valores estimados das estatísticas de Wald para cada replica
  w.stat.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  w.stat.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  w.stat.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  w.stat.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Passo 1: Estimativas dos coeficientes e erros-padrão calculados anteriormente
  # Esses entram na função a partir dos argumentos iniciais, sendo:
  # - theta (beta e gamma), as estimativa pontuais iniciais obtidas e
  # - ep(theta) obtido via bootstrap
  
  # Passo 2: Obter as estatísticas W - Serão a estatistica de teste utilizada
  w.stat.beta <- (c(beta/std.bootstrap.theta[1:length(beta)]))^2
  w.stat.gamma <- (c(gamma/std.bootstrap.theta[(length(beta)+1):(length(beta)+length(gamma))]))^2
  
  # Passo 3: Novo bootstrap para Obter B2=500 réplicas de y e, para cada uma, a estimativa de theta (beta e gamma).
  
  alg2 <- list(list(),
               list(),
               list())
  
  # Algoritmo para geracao das amostras sob H0, com paralelização
  alg2 <- foreach(j = 1:B2,
                  .verbose = FALSE,
                  .export = c("FunSample_betareg",
                              "apply_spline",
                              "Psi_SMLE",
                              "L_q",
                              "L_alpha_R", 
                              "Psi_LSMLE",
                              "Loglik_MLE",
                              "Psi_MLE",
                              "dEGB",
                              "apply_derivative",
                              "set.link",
                              "starting.points"),
                  .packages = c("base",
                                "stats",
                                "foreach",
                                "doParallel",
                                "betareg",
                                "Formula")) %dopar% {
                                  
                                  # Nesse laço são geradas, a partir do theta inicial, B2 réplicas de y e estimados os thetas de cada uma
                                  # São retornados o beta_b2 e o gamma_b2
                                  
                                  # Geracao da amostra aleatoria
                                  #set.seed(2000+j)
                                  FunSample_betareg(n=n, 
                                                    X=X, 
                                                    Z=Z,
                                                    beta = replicate(n=length(beta), expr=c(0)),
                                                    gamma = replicate(n=length(gamma), expr=c(0)),
                                                    functions = functions,
                                                    parameters = parameters,
                                                    link.mu=link, 
                                                    link.phi=link.phi) -> dados_B2
                                  
                                  dados_B2$y[dados_B2$y>0.999999999999999] <- 0.999999999999999
                                  
                                  # Calculo dos parametros
                                  # Novo chute inicial, pois está alterando o y a cada passada
                                  sim_df_B2 <- base::data.frame(dados_B2$y,
                                                                 dados_B2$X,
                                                                 dados_B2$Z)

                                  names(sim_df_B2) <- c(data.names[1],colnames(X),colnames(Z))
                                  attach(sim_df_B2)
                                  
                                  tryCatch(starting.points(formula = formula,
                                                           data = sim_df_B2,
                                                           y=dados_B2$y,
                                                           x=dados_B2$X,
                                                           z=dados_B2$Z,
                                                           type = type,
                                                           parameters = parameters,
                                                           functions = functions,
                                                           link = link,
                                                           link.phi = link.phi), 
                                           error = function(e) {
                                             return(replicate(n=length(c(beta,gamma)), expr=c(0)))
                                           }) -> start_theta_B2
                                  
                                  tryCatch(
                                    # Calculo dos parametros
                                    if(type=="SMLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = L_q,
                                                   gr = Psi_SMLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   alpha = alpha,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }else if(type=="LSMLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = L_alpha_R,
                                                   gr = Psi_LSMLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   alpha = alpha,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }else if(type=="MLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = Loglik_MLE,
                                                   gr = Psi_MLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }, 
                                  error = function(e) {
                                    return(start_theta_B2)
                                  }) -> theta_B2
                                  
                                  
                                  # Obtenção dos coeficientes
                                  beta_B2 <- as.numeric(theta_B2$par[1:length(beta)])
                                  gamma_B2 <- as.numeric(theta_B2$par[(length(beta)+1):(length(beta)+length(gamma))])
                                  #####################################################################################################
                                  alg3 <- foreach(j = 1:B1,
                                                  .verbose = FALSE,
                                                  .export = c("FunSample_betareg",
                                                              "apply_spline",
                                                              "Psi_SMLE",
                                                              "L_q",
                                                              "L_alpha_R", 
                                                              "Psi_LSMLE",
                                                              "dEGB",
                                                              "apply_derivative",
                                                              "set.link",
                                                              "theta_B2",
                                                              "starting.points",
                                                              "beta_B2",
                                                              "gamma_B2"),
                                                  .packages = c("base",
                                                                "stats",
                                                                "foreach",
                                                                "doParallel",
                                                                "betareg",
                                                                "Formula")) %dopar% {
                                                                  
                                                                  # Nesse laço são geradas, a partir do theta_B2 anterior, B1 réplicas de y e estimados os thetas de cada uma
                                                                  # São retornados o theta_B1
                                                                  
                                                                  # Geracao da amostra aleatoria
                                                                  FunSample_betareg(n=n,
                                                                                    X=X,
                                                                                    Z=Z,
                                                                                    beta = beta_B2,
                                                                                    gamma = gamma_B2,
                                                                                    functions = functions,
                                                                                    parameters = parameters,
                                                                                    link.mu=link,
                                                                                    link.phi=link.phi) -> dados_B1
                                                                  
                                                                  dados_B1$y[dados_B1$y>0.999999999999999] <- 0.999999999999999
                                                                  
                                                                  
                                                                  # Calculo dos parametros
                                                                  # Novo chute inicial, pois está alterando o y a cada passada
                                                                  sim_df_B1 <- base::data.frame(dados_B1$y,
                                                                                                dados_B1$X,
                                                                                                dados_B1$Z)
                                                                  
                                                                  names(sim_df_B1) <- c(data.names[1],colnames(X),colnames(Z))
                                                                  attach(sim_df_B1)
                                                                  
                                                                  tryCatch(starting.points(formula = formula,
                                                                                           data = sim_df_B1,
                                                                                           y=dados_B1$y,
                                                                                           x=dados_B1$X,
                                                                                           z=dados_B1$Z,
                                                                                           type = type,
                                                                                           parameters = parameters,
                                                                                           functions = functions,
                                                                                           link = link,
                                                                                           link.phi = link.phi), 
                                                                           error = function(e) {
                                                                             return(c(beta_B2, gamma_B2))
                                                                           }) -> start_theta_B1
                                                                  
                                                                  # Obtenção dos coeficientes
                                                                  if(type=="SMLE"){
                                                                    tryCatch(stats::optim(par = start_theta_B1,
                                                                                          fn = L_q,
                                                                                          gr = Psi_SMLE,
                                                                                          y = dados_B1$y,
                                                                                          X=X,
                                                                                          Z=Z,
                                                                                          alpha = alpha,
                                                                                          link_mu = link,
                                                                                          link_phi = link.phi,
                                                                                          functions = functions,
                                                                                          parameters = parameters,
                                                                                          apply_spline = apply_spline,
                                                                                          apply_derivative = apply_derivative,
                                                                                          control = list(fnscale = -1)), 
                                                                             error = function(e) {
                                                                               return(NULL)
                                                                             }) -> theta_B1
                                                                  }else if(type=="LSMLE"){
                                                                    tryCatch(stats::optim(par = start_theta_B1,
                                                                                          fn = L_alpha_R,
                                                                                          gr = Psi_LSMLE,
                                                                                          y = dados_B1$y,
                                                                                          X=X,
                                                                                          Z=Z,
                                                                                          alpha = alpha,
                                                                                          link_mu = link,
                                                                                          link_phi = link.phi,
                                                                                          functions = functions,
                                                                                          parameters = parameters,
                                                                                          apply_spline = apply_spline,
                                                                                          apply_derivative = apply_derivative,
                                                                                          control = list(fnscale = -1)), 
                                                                             error = function(e) {
                                                                               return(NULL)
                                                                             }) -> theta_B1

                                                                  }else if(type=="MLE"){
                                                                    tryCatch(stats::optim(par = start_theta_B1,
                                                                                          fn = Loglik_MLE,
                                                                                          gr = Psi_MLE,
                                                                                          y = dados_B1$y,
                                                                                          X=X,
                                                                                          Z=Z,
                                                                                          link_mu = link,
                                                                                          link_phi = link.phi,
                                                                                          functions = functions,
                                                                                          parameters = parameters,
                                                                                          apply_spline = apply_spline,
                                                                                          apply_derivative = apply_derivative,
                                                                                          control = list(fnscale = -1)), 
                                                                             error = function(e) {
                                                                               return(replicate(n=length(start_theta), expr=c(NULL)))
                                                                             }) -> theta_B1
                                                                  }
                                                                  
                                                                  list("theta_B1" = theta_B1$par)
                                                                }
                                  
                                  #####################################################################################################
                                  list("alg3" = alg3,
                                       "beta_B2" = beta_B2,
                                       "gamma_B2" = gamma_B2)
                                }
  
  closeAllConnections()
  
  # Armazena os valores das estimativas dos coeficientes
  for (j in 1:B2) {
    for (k in 1:B1){
      for (m in 1:length(beta)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m]]) ){
          beta_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m]]
        }else{
          beta_B1.hat[[j]][[m]][k] <- NA
        }
      }
      for (m in 1:length(gamma)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]) ){
          gamma_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]
        }else{
          gamma_B1.hat[[j]][[m]][k] <- NA
        }
      }
    }
  }
  # Calcula os erros-padrão amostrais das estimativas dos coeficientes
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1[[j]][[m]] = sd(beta_B1.hat[[j]][[m]])
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1[[j]][[m]] = sd(gamma_B1.hat[[j]][[m]])
    }
  }
  # Calcula as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1[[j]][[m]] = (alg2[[j]]$beta_B2[m]/sd(beta_B1.hat[[j]][[m]]))^2
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1[[j]][[m]] = (alg2[[j]]$gamma_B2[m]/sd(gamma_B1.hat[[j]][[m]]))^2
    }
  }
  
  # Agrupa por coeficiente os erros-padrão
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1_vec[[m]][[j]] =  bootstrap.std.error.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1_vec[[m]][[j]] =  bootstrap.std.error.gamma_B1[[j]][[m]]
    }
  }
  
  # Agrupa por coeficiente as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1_vec[[m]][[j]] = w.stat.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1_vec[[m]][[j]] = w.stat.gamma_B1[[j]][[m]]
    }
  }
  # Obtem os p-valores referentes às estatísticas de Wald
  pvalue_w.stat.beta.bootstrap <- c()
  pvalue_w.stat.gamma.bootstrap <- c()
  for (m in 1:length(beta)) {
    pvalue_w.stat.beta.bootstrap[m] <- mean(unlist(w.stat.beta_B1_vec[[m]]) > w.stat.beta[m], na.rm = TRUE)
  }
  for (m in 1:length(gamma)) {
    pvalue_w.stat.gamma.bootstrap[m] <- mean(unlist(w.stat.gamma_B1_vec[[m]]) > w.stat.gamma[m], na.rm = TRUE)
  }
  
  valid.w.stat <- sum(!is.na(unlist(w.stat.beta_B1_vec[[m]])))
  
  message(paste("INFO: The simulated distribution for the type-Wald test statistics is ready and the p-value has been calculated!"))
  
  end_time <- Sys.time()
  total.processing.time <- (end_time - start_time)
  message("\nINFO: Process of type-Wald test finished! Total time: ", total.processing.time)
  
  return(list("pvalue_w.stat.bootstrap" = c(pvalue_w.stat.beta.bootstrap, pvalue_w.stat.gamma.bootstrap),
              "valid.w.stat" = valid.w.stat,
              "B2" = B2,
              "B1" = B1)
         )

}

################################################################################
# Calculate wald type test estatistics simulated distribution via bootstrap 

bootstrap.wald.test <- function(B2 = 500, #First level (number of w statistic replicas)
                                B1 = 200, # Second level (for the bootstrap sd)
                                n,
                                alpha,
                                beta,
                                gamma,
                                X,
                                Y,
                                Z,
                                type,
                                formula,
                                data.names,
                                std.bootstrap.theta,
                                functions, 
                                parameters,
                                start_theta,
                                link,
                                link.phi
){
  
  closeAllConnections()
  
  start_time <- Sys.time()
  
  #Leitura das bibliotecas que farão a paralelização
  library(foreach)
  library(doParallel)
  
  message(paste("INFO: Now starting to generate simulated distribution for the type-wald test statistics. It's a good time to grab some coffee..."))
  
  # Tratamento do ojeto que carrega os erros padrões
  std.bootstrap.theta <- c(as.numeric(std.bootstrap.theta$se.mean),as.numeric(std.bootstrap.theta$se.precision))
  
  # Deteccao da quantidade de nucleos da CPU e registro do backend p/ paralelizacao 
  # numCores <- detectCores()
  # registerDoParallel(numCores)
  
  # Objeto para armazenar os valores estimados dos betas e gammas
  beta_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  gamma_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  # Objeto para armazenar os valores estimados dos erros-padrão para cada replica
  bootstrap.std.error.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  bootstrap.std.error.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  bootstrap.std.error.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  bootstrap.std.error.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Objeto para armazenar os valores estimados das estatísticas de Wald para cada replica
  w.stat.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  w.stat.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  w.stat.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  w.stat.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Passo 1: Estimativas dos coeficientes e erros-padrão calculados anteriormente
  # Esses entram na função a partir dos argumentos iniciais, sendo:
  # - theta (beta e gamma), as estimativa pontuais iniciais obtidas e
  # - ep(theta) obtido via bootstrap
  
  # Passo 2: Obter as estatísticas W - Serão a estatistica de teste utilizada
  w.stat.beta <- (c(beta/std.bootstrap.theta[1:length(beta)]))^2
  w.stat.gamma <- (c(gamma/std.bootstrap.theta[(length(beta)+1):(length(beta)+length(gamma))]))^2
  
  # Passo 3: Novo bootstrap para Obter B2=500 réplicas de y e, para cada uma, a estimativa de theta (beta e gamma).
  
  alg2 <- list(list(),
               list(),
               list())
  
  for (j in 1:B2) {
    
    # Deteccao da quantidade de nucleos da CPU e registro do backend p/ paralelizacao 
    numCores <- detectCores()
    registerDoParallel(numCores)
    
    message(paste("INFO: (TYPE WALD TEST): Replica",j))
    
    # Nesse laço são geradas, a partir do theta inicial, B2 réplicas de y e estimados os thetas de cada uma
    # São retornados o beta_b2 e o gamma_b2
    
    # Geracao da amostra aleatoria
    #set.seed(2000+j)
    FunSample_betareg(n=n, 
                      X=X, 
                      Z=Z,
                      beta = replicate(n=length(beta), expr=c(0)),
                      gamma = replicate(n=length(gamma), expr=c(0)),
                      functions = functions,
                      parameters = parameters,
                      link.mu=link, 
                      link.phi=link.phi) -> dados_B2
    
    dados_B2$y[dados_B2$y>0.999999999999999] <- 0.999999999999999
    
    # Calculo dos parametros
    # Novo chute inicial, pois está alterando o y a cada passada
    sim_df_B2 <- base::data.frame(dados_B2$y,
                                  dados_B2$X,
                                  dados_B2$Z)
    
    names(sim_df_B2) <- c(data.names[1],colnames(X),colnames(Z))
    # attach(sim_df_B2)
    
    tryCatch(starting.points(formula = formula,
                             data = sim_df_B2,
                             y=dados_B2$y,
                             x=dados_B2$X,
                             z=dados_B2$Z,
                             type = type,
                             parameters = parameters,
                             functions = functions,
                             link = link,
                             link.phi = link.phi), 
             error = function(e) {
               return(replicate(n=length(c(beta,gamma)), expr=c(0)))
             }) -> start_theta_B2
    
    tryCatch(
      # Calculo dos parametros
      if(type=="SMLE"){
        stats::optim(par = start_theta_B2,
                     fn = L_q,
                     gr = Psi_SMLE,
                     y = dados_B2$y,
                     X=X,
                     Z=Z,
                     alpha = alpha,
                     link_mu = link,
                     link_phi = link.phi,
                     functions = functions,
                     parameters = parameters,
                     apply_spline = apply_spline,
                     apply_derivative = apply_derivative,
                     control = list(fnscale = -1))
      }else if(type=="LSMLE"){
        stats::optim(par = start_theta_B2,
                     fn = L_alpha_R,
                     gr = Psi_LSMLE,
                     y = dados_B2$y,
                     X=X,
                     Z=Z,
                     alpha = alpha,
                     link_mu = link,
                     link_phi = link.phi,
                     functions = functions,
                     parameters = parameters,
                     apply_spline = apply_spline,
                     apply_derivative = apply_derivative,
                     control = list(fnscale = -1))
      }else if(type=="MLE"){
        stats::optim(par = start_theta_B2,
                     fn = Loglik_MLE,
                     gr = Psi_MLE,
                     y = dados_B2$y,
                     X=X,
                     Z=Z,
                     link_mu = link,
                     link_phi = link.phi,
                     functions = functions,
                     parameters = parameters,
                     apply_spline = apply_spline,
                     apply_derivative = apply_derivative,
                     control = list(fnscale = -1))
      }, 
      error = function(e) {
        return(start_theta_B2)
      }) -> theta_B2
    
    
    # Obtenção dos coeficientes
    beta_B2 <- as.numeric(theta_B2$par[1:length(beta)])
    gamma_B2 <- as.numeric(theta_B2$par[(length(beta)+1):(length(beta)+length(gamma))])
    #####################################################################################################
    alg3 <- foreach(j = 1:B1,
                    .verbose = FALSE,
                    .export = c("Loglik_MLE",
                                "Psi_MLE",
                                "FunSample_betareg",
                                "apply_spline",
                                "Psi_SMLE",
                                "L_q",
                                "L_alpha_R", 
                                "Psi_LSMLE",
                                "dEGB",
                                "apply_derivative",
                                "set.link",
                                "theta_B2",
                                "starting.points"),
                    .packages = c("base",
                                  "stats",
                                  "foreach",
                                  "doParallel",
                                  "betareg",
                                  "Formula")) %dopar% {
                                    
                                    # Nesse laço são geradas, a partir do theta_B2 anterior, B1 réplicas de y e estimados os thetas de cada uma
                                    # São retornados o theta_B1
                                    
                                    # Geracao da amostra aleatoria
                                    FunSample_betareg(n=n,
                                                      X=X,
                                                      Z=Z,
                                                      beta = beta_B2,
                                                      gamma = gamma_B2,
                                                      functions = functions,
                                                      parameters = parameters,
                                                      link.mu=link,
                                                      link.phi=link.phi) -> dados_B1
                                    
                                    dados_B1$y[dados_B1$y>0.999999999999999] <- 0.999999999999999
                                    
                                    
                                    # Calculo dos parametros
                                    # Novo chute inicial, pois está alterando o y a cada passada
                                    sim_df_B1 <- base::data.frame(dados_B1$y,
                                                                  dados_B1$X,
                                                                  dados_B1$Z)
                                    
                                    names(sim_df_B1) <- c(data.names[1],colnames(X),colnames(Z))
                                    attach(sim_df_B1)
                                    
                                    tryCatch(starting.points(formula = formula,
                                                             data = sim_df_B1,
                                                             y=dados_B1$y,
                                                             x=dados_B1$X,
                                                             z=dados_B1$Z,
                                                             type = type,
                                                             parameters = parameters,
                                                             functions = functions,
                                                             link = link,
                                                             link.phi = link.phi), 
                                             error = function(e) {
                                               return(c(beta_B2, gamma_B2))
                                             }) -> start_theta_B1
                                    
                                    # Obtenção dos coeficientes
                                    if(type=="SMLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = L_q,
                                                            gr = Psi_SMLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            alpha = alpha,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(NULL)
                                               }) -> theta_B1
                                    }else if(type=="LSMLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = L_alpha_R,
                                                            gr = Psi_LSMLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            alpha = alpha,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(NULL)
                                               }) -> theta_B1
                                      
                                    }else if(type=="MLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = Loglik_MLE,
                                                            gr = Psi_MLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(replicate(n=length(start_theta), expr=c(NULL)))
                                               }) -> theta_B1
                                    }
                                    
                                    list("theta_B1" = theta_B1$par)
                                  }
    
    #####################################################################################################
    alg2[[j]] <- list("alg3" = alg3,
                     "beta_B2" = beta_B2,
                     "gamma_B2" = gamma_B2)
    
    closeAllConnections()

    
  }
  
  closeAllConnections()

  
  # Armazena os valores das estimativas dos coeficientes
  for (j in 1:B2) {
    for (k in 1:B1){
      for (m in 1:length(beta)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m]]) ){
          beta_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m]]
        }else{
          beta_B1.hat[[j]][[m]][k] <- NA
        }
      }
      for (m in 1:length(gamma)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]) ){
          gamma_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]
        }else{
          gamma_B1.hat[[j]][[m]][k] <- NA
        }
      }
    }
  }
  # Calcula os erros-padrão amostrais das estimativas dos coeficientes
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1[[j]][[m]] = sd(beta_B1.hat[[j]][[m]])
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1[[j]][[m]] = sd(gamma_B1.hat[[j]][[m]])
    }
  }
  # Calcula as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1[[j]][[m]] = (alg2[[j]]$beta_B2[m]/sd(beta_B1.hat[[j]][[m]]))^2
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1[[j]][[m]] = (alg2[[j]]$gamma_B2[m]/sd(gamma_B1.hat[[j]][[m]]))^2
    }
  }
  
  # Agrupa por coeficiente os erros-padrão
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1_vec[[m]][[j]] =  bootstrap.std.error.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1_vec[[m]][[j]] =  bootstrap.std.error.gamma_B1[[j]][[m]]
    }
  }
  
  # Agrupa por coeficiente as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1_vec[[m]][[j]] = w.stat.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1_vec[[m]][[j]] = w.stat.gamma_B1[[j]][[m]]
    }
  }
  # Obtem os p-valores referentes às estatísticas de Wald
  pvalue_w.stat.beta.bootstrap <- c()
  pvalue_w.stat.gamma.bootstrap <- c()
  for (m in 1:length(beta)) {
    pvalue_w.stat.beta.bootstrap[m] <- mean(unlist(w.stat.beta_B1_vec[[m]]) > w.stat.beta[m], na.rm = TRUE)
  }
  for (m in 1:length(gamma)) {
    pvalue_w.stat.gamma.bootstrap[m] <- mean(unlist(w.stat.gamma_B1_vec[[m]]) > w.stat.gamma[m], na.rm = TRUE)
  }
  
  valid.w.stat <- sum(!is.na(unlist(w.stat.beta_B1_vec[[m]])))
  
  message(paste("INFO: The simulated distribution for the type-Wald test statistics is ready and the p-value has been calculated!"))
  
  end_time <- Sys.time()
  total.processing.time <- (end_time - start_time)
  message("\nINFO: Process of type-Wald test finished! Total time: ", total.processing.time)
  
  return(list("pvalue_w.stat.bootstrap" = c(pvalue_w.stat.beta.bootstrap, pvalue_w.stat.gamma.bootstrap),
              "valid.w.stat" = valid.w.stat,
              "B2" = B2,
              "B1" = B1)
  )
  
}

######################################################################################################################

exp2_bootstrap.wald.test <- function(B2 = 200, #First level (number of w statistic replicas)
                                B1 = 200, # Second level (for the bootstrap sd)
                                n,
                                alpha,
                                beta,
                                gamma,
                                X,
                                Y,
                                Z,
                                type,
                                formula,
                                data.names,
                                std.bootstrap.theta,
                                functions, 
                                parameters,
                                start_theta,
                                link,
                                link.phi
){
  
  closeAllConnections()
  
  start_time <- Sys.time()
  
  #Leitura das bibliotecas que farão a paralelização
  library(foreach)
  library(doParallel)
  
  message(paste("INFO: Now starting to generate simulated distribution for the type-wald test statistics. It's a good time to grab some coffee..."))
  
  # Tratamento do ojeto que carrega os erros padrões
  std.bootstrap.theta <- c(as.numeric(std.bootstrap.theta$se.mean),as.numeric(std.bootstrap.theta$se.precision))
  
  # Deteccao da quantidade de nucleos da CPU e registro do backend p/ paralelizacao 
  numCores <- detectCores()
  registerDoParallel(numCores)
  
  # Objeto para armazenar os valores estimados dos betas e gammas
  beta_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  gamma_B1.hat <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  # Objeto para armazenar os valores estimados dos erros-padrão para cada replica
  bootstrap.std.error.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  bootstrap.std.error.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  bootstrap.std.error.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  bootstrap.std.error.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Objeto para armazenar os valores estimados das estatísticas de Wald para cada replica
  w.stat.beta_B1 <- replicate(n=B2, expr=list(replicate(n=length(beta), expr = list(0))))
  w.stat.gamma_B1 <- replicate(n=B2, expr=list(replicate(n=length(gamma), expr = list(0))))
  
  w.stat.beta_B1_vec <- lapply(rep(list(0), length(beta)), list)
  w.stat.gamma_B1_vec <- lapply(rep(list(0), length(gamma)), list)
  
  # Passo 1: Estimativas dos coeficientes e erros-padrão calculados anteriormente
  # Esses entram na função a partir dos argumentos iniciais, sendo:
  # - theta (beta e gamma), as estimativa pontuais iniciais obtidas e
  # - ep(theta) obtido via bootstrap
  
  # Passo 2: Obter as estatísticas W - Serão a estatistica de teste utilizada
  w.stat.beta <- (c(beta/std.bootstrap.theta[1:length(beta)]))^2
  w.stat.gamma <- (c(gamma/std.bootstrap.theta[(length(beta)+1):(length(beta)+length(gamma))]))^2
  
  # Passo 3: Novo bootstrap para Obter B2=500 réplicas de y e, para cada uma, a estimativa de theta (beta e gamma).
  
  alg2 <- list(list(),
               list(),
               list())
  
  alg3 <- list()
  
  # Algoritmo para geracao das amostras sob H0, com paralelização
  alg2 <- foreach(j = 1:B2,
                  .verbose = FALSE,
                  .export = c("FunSample_betareg",
                              "apply_spline",
                              "Psi_SMLE",
                              "L_q",
                              "L_alpha_R", 
                              "Psi_LSMLE",
                              "Loglik_MLE",
                              "Psi_MLE",
                              "dEGB",
                              "apply_derivative",
                              "set.link",
                              "starting.points"),
                  .packages = c("base",
                                "stats",
                                "foreach",
                                "doParallel",
                                "betareg",
                                "Formula")) %dopar% {
                                  
                                  # Nesse laço são geradas, a partir do theta inicial, B2 réplicas de y e estimados os thetas de cada uma
                                  # São retornados o beta_b2 e o gamma_b2
                                  
                                  # Geracao da amostra aleatoria
                                  #set.seed(2000+j)
                                  FunSample_betareg(n=n, 
                                                    X=X, 
                                                    Z=Z,
                                                    beta = replicate(n=length(beta), expr=c(0)),
                                                    gamma = replicate(n=length(gamma), expr=c(0)),
                                                    functions = functions,
                                                    parameters = parameters,
                                                    link.mu=link, 
                                                    link.phi=link.phi) -> dados_B2
                                  
                                  dados_B2$y[dados_B2$y>0.999999999999999] <- 0.999999999999999
                                  
                                  # Calculo dos parametros
                                  # Novo chute inicial, pois está alterando o y a cada passada
                                  sim_df_B2 <- base::data.frame(dados_B2$y,
                                                                dados_B2$X,
                                                                dados_B2$Z)
                                  
                                  names(sim_df_B2) <- c(data.names[1],colnames(X),colnames(Z))
                                  attach(sim_df_B2)
                                  
                                  tryCatch(starting.points(formula = formula,
                                                           data = sim_df_B2,
                                                           y=dados_B2$y,
                                                           x=dados_B2$X,
                                                           z=dados_B2$Z,
                                                           type = type,
                                                           parameters = parameters,
                                                           functions = functions,
                                                           link = link,
                                                           link.phi = link.phi), 
                                           error = function(e) {
                                             return(replicate(n=length(c(beta,gamma)), expr=c(0)))
                                           }) -> start_theta_B2
                                  
                                  tryCatch(
                                    # Calculo dos parametros
                                    if(type=="SMLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = L_q,
                                                   gr = Psi_SMLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   alpha = alpha,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }else if(type=="LSMLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = L_alpha_R,
                                                   gr = Psi_LSMLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   alpha = alpha,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }else if(type=="MLE"){
                                      stats::optim(par = start_theta_B2,
                                                   fn = Loglik_MLE,
                                                   gr = Psi_MLE,
                                                   y = dados_B2$y,
                                                   X=X,
                                                   Z=Z,
                                                   link_mu = link,
                                                   link_phi = link.phi,
                                                   functions = functions,
                                                   parameters = parameters,
                                                   apply_spline = apply_spline,
                                                   apply_derivative = apply_derivative,
                                                   control = list(fnscale = -1))
                                    }, 
                                    error = function(e) {
                                      return(start_theta_B2)
                                    }) -> theta_B2
                                  
                                  
                                  # Obtenção dos coeficientes
                                  beta_B2 <- as.numeric(theta_B2$par[1:length(beta)])
                                  gamma_B2 <- as.numeric(theta_B2$par[(length(beta)+1):(length(beta)+length(gamma))])
                                  #####################################################################################################
                                  for (k in 1:B1) {
                                    
                                    # Nesse laço são geradas, a partir do theta_B2 anterior, B1 réplicas de y e estimados os thetas de cada uma
                                    # São retornados o theta_B1
                                    
                                    # Geracao da amostra aleatoria
                                    FunSample_betareg(n=n,
                                                      X=X,
                                                      Z=Z,
                                                      beta = beta_B2,
                                                      gamma = gamma_B2,
                                                      functions = functions,
                                                      parameters = parameters,
                                                      link.mu=link,
                                                      link.phi=link.phi) -> dados_B1
                                    
                                    dados_B1$y[dados_B1$y>0.999999999999999] <- 0.999999999999999
                                    
                                    
                                    # Calculo dos parametros
                                    # Novo chute inicial, pois está alterando o y a cada passada
                                    sim_df_B1 <- base::data.frame(dados_B1$y,
                                                                  dados_B1$X,
                                                                  dados_B1$Z)
                                    
                                    names(sim_df_B1) <- c(data.names[1],colnames(X),colnames(Z))
                                    attach(sim_df_B1)
                                    
                                    tryCatch(starting.points(formula = formula,
                                                             data = sim_df_B1,
                                                             y=dados_B1$y,
                                                             x=dados_B1$X,
                                                             z=dados_B1$Z,
                                                             type = type,
                                                             parameters = parameters,
                                                             functions = functions,
                                                             link = link,
                                                             link.phi = link.phi), 
                                             error = function(e) {
                                               return(c(beta_B2, gamma_B2))
                                             }) -> start_theta_B1
                                    
                                    # Obtenção dos coeficientes
                                    if(type=="SMLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = L_q,
                                                            gr = Psi_SMLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            alpha = alpha,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(NULL)
                                               }) -> theta_B1
                                    }else if(type=="LSMLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = L_alpha_R,
                                                            gr = Psi_LSMLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            alpha = alpha,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(NULL)
                                               }) -> theta_B1
                                      
                                    }else if(type=="MLE"){
                                      tryCatch(stats::optim(par = start_theta_B1,
                                                            fn = Loglik_MLE,
                                                            gr = Psi_MLE,
                                                            y = dados_B1$y,
                                                            X=X,
                                                            Z=Z,
                                                            link_mu = link,
                                                            link_phi = link.phi,
                                                            functions = functions,
                                                            parameters = parameters,
                                                            apply_spline = apply_spline,
                                                            apply_derivative = apply_derivative,
                                                            control = list(fnscale = -1)), 
                                               error = function(e) {
                                                 return(replicate(n=length(start_theta), expr=c(NULL)))
                                               }) -> theta_B1
                                      
                                  }
                                  alg3[[k]] <- list(theta_B1 = theta_B1$par)
                                  }
                                   
                                  #####################################################################################################
                                  list("alg3" = alg3,
                                       "beta_B2" = beta_B2,
                                       "gamma_B2" = gamma_B2)
                                }
  
  closeAllConnections()
  
  # print("AQUI")
  # browser()
  # Armazena os valores das estimativas dos coeficientes
  for (j in 1:B2) {
    for (k in 1:B1){
      for (m in 1:length(beta)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m]]) ){
          beta_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m]]
        }else{
          beta_B1.hat[[j]][[m]][k] <- NA
        }
      }
      for (m in 1:length(gamma)) {
        if( !is.null(alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]) ){
          gamma_B1.hat[[j]][[m]][k] <- alg2[[j]]$alg3[[k]]$theta_B1[[m+length(beta)]]
        }else{
          gamma_B1.hat[[j]][[m]][k] <- NA
        }
      }
    }
  }
  # Calcula os erros-padrão amostrais das estimativas dos coeficientes
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1[[j]][[m]] = sd(beta_B1.hat[[j]][[m]])
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1[[j]][[m]] = sd(gamma_B1.hat[[j]][[m]])
    }
  }
  # Calcula as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1[[j]][[m]] = (alg2[[j]]$beta_B2[m]/sd(beta_B1.hat[[j]][[m]]))^2
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1[[j]][[m]] = (alg2[[j]]$gamma_B2[m]/sd(gamma_B1.hat[[j]][[m]]))^2
    }
  }
  
  # Agrupa por coeficiente os erros-padrão
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      bootstrap.std.error.beta_B1_vec[[m]][[j]] =  bootstrap.std.error.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      bootstrap.std.error.gamma_B1_vec[[m]][[j]] =  bootstrap.std.error.gamma_B1[[j]][[m]]
    }
  }
  
  # Agrupa por coeficiente as estatísticas de Wald
  for (j in 1:B2) {
    for (m in 1:length(beta)) {
      w.stat.beta_B1_vec[[m]][[j]] = w.stat.beta_B1[[j]][[m]]
    }
    for (m in 1:length(gamma)) {
      w.stat.gamma_B1_vec[[m]][[j]] = w.stat.gamma_B1[[j]][[m]]
    }
  }
  # Obtem os p-valores referentes às estatísticas de Wald
  pvalue_w.stat.beta.bootstrap <- c()
  pvalue_w.stat.gamma.bootstrap <- c()
  for (m in 1:length(beta)) {
    pvalue_w.stat.beta.bootstrap[m] <- mean(unlist(w.stat.beta_B1_vec[[m]]) > w.stat.beta[m], na.rm = TRUE)
  }
  for (m in 1:length(gamma)) {
    pvalue_w.stat.gamma.bootstrap[m] <- mean(unlist(w.stat.gamma_B1_vec[[m]]) > w.stat.gamma[m], na.rm = TRUE)
  }
  
  valid.w.stat <- sum(!is.na(unlist(w.stat.beta_B1_vec[[m]])))
  
  message(paste("INFO: The simulated distribution for the type-Wald test statistics is ready and the p-value has been calculated!"))
  
  end_time <- Sys.time()
  total.processing.time <- (end_time - start_time)
  message("\nINFO: Process of type-Wald test finished! Total time: ", total.processing.time)
  
  return(list("pvalue_w.stat.bootstrap" = c(pvalue_w.stat.beta.bootstrap, pvalue_w.stat.gamma.bootstrap),
              "valid.w.stat" = valid.w.stat,
              "B2" = B2,
              "B1" = B1)
  )
  
}

######################################################################################################################
# Tuning constant selection algorithm

SQV <- function (zq, n, p)
{
  if (!is.null(zq) && dim(zq)[1] != 1) {
    # return(sqrt(apply((diff(zq, 1, 1))^2, 1, sum))/(p))
    return(sqrt(apply((diff(zq, 1, 1))^2, 1, sum)))
  }
  else {
    return(c(0))
  }
}

# SQV <- function (zq, n, p)
# {
#   if (!is.null(zq) && dim(zq)[1] != 1) {
#     return(sqrt(apply((diff(zq, 1, 1))^2, 1, sum))/(sqrt(n)*p))
#   }
#   else {
#     return(c(0))
#   }
# }

select_tuning_constant <- function(y,
                                   x,
                                   z,
                                   type,
                                   link,
                                   link.phi,
                                   start_theta,
                                   se.bootstrap=FALSE,
                                   wald.test.bootstrap,
                                   functions,
                                   parameters,
                                   data,
                                   formula,
                                   control=robustbetareg.control()){
  
  message("\nINFO: Starting to run the algorithm for the optimal tuning contant selection!")
  
  # Captura dos parãmetros a serem usados no cálculo
  M = control$M
  
  n = nrow(x)
  alpha.max = control$alpha.max
  unstable = TRUE
  stability = TRUE
  loop.number = 1
  output = list()
  
  # Passo 1#####################################################################
  M1 = 11
  q.grade = seq(from=1, to=0.8, by= -0.02)[1:M1]
  alpha.grade = 1-q.grade
  
  # Passo 2#####################################################################
  while (unstable) {
    
    message("\nINFO: Ciclo ",loop.number)
    
    # Criação/reset de objetos
    est.list <- list()
    est.par <- list()
    conv.list <- list()
    z_q.list <- list()
    sqv.vec <- NULL
    k = length(alpha.grade)
    
    for (k in 1:length(alpha.grade)) {
      message("\nINFO: Loop ",k)
      if(alpha.grade[k]==0){
        est.par = tryCatch(MLE.fit(y=y, 
                                   x=x, 
                                   z=z, 
                                   link = link, 
                                   link.phi = link.phi,
                                   type="MLE",
                                   se.bootstrap = FALSE,
                                   wald.test.bootstrap = wald.test.bootstrap, 
                                   functions = functions, 
                                   parameters = parameters,
                                   data = data,
                                   formula = formula,
                                   control = robustbetareg.control(start = start_theta)), 
                           error = function(e) {
                             est.par$converged <- FALSE
                             return(est.par)
                           })
      }else{
        if(type=="SMLE"){
          est.par = tryCatch(SMLE.fit(y=y, 
                                      x=x, 
                                      z=z, 
                                      alpha = alpha.grade[k], 
                                      link = link, 
                                      link.phi = link.phi,
                                      type=type,
                                      se.bootstrap = FALSE,
                                      wald.test.bootstrap = wald.test.bootstrap, 
                                      functions = functions, 
                                      parameters = parameters,
                                      control = robustbetareg.control(start = start_theta)), 
                             error = function(e) {
                               est.par$converged <- FALSE
                               return(est.par)
                             })
        }else if(type=="LSMLE"){
          est.par = tryCatch(LSMLE.fit(y=y, 
                                       x=x, 
                                       z=z, 
                                       alpha = alpha.grade[k], 
                                       link = link, 
                                       link.phi = link.phi,
                                       type=type,
                                       se.bootstrap = FALSE,
                                       wald.test.bootstrap = wald.test.bootstrap, 
                                       functions = functions, 
                                       parameters = parameters,
                                       control = robustbetareg.control(start = start_theta)), 
                             error = function(e) {
                               est.par$converged <- FALSE
                               return(est.par)
                             })
        }
      }
      est.list[[k]] <- est.par
      
      # Passo 3#####################################################################
      
      z_q.list[[k]] <- c(est.par$coefficients$mean,est.par$coefficients$precision)

      
      # z_q.list[[k]] <- ((c(est.par$coefficients$mean,est.par$coefficients$precision))/
      #                     (sqrt(n)*c(est.par$std.error$se.mean,est.par$std.error$se.precision))
      # 
      # )
      
      conv.list[[k]] <- est.par$converged
    }
    # browser()
    # z_q.matrix <- matrix(data=unlist(z_q.list),
    #                      nrow = k,
    #                      ncol = length(start_theta),
    #                      byrow = TRUE)
    
    z_q.matrix <- tryCatch(matrix(data=unlist(z_q.list),
                                  nrow = k,
                                  ncol = length(start_theta),
                                  byrow = TRUE), 
                           error = function(e) {
                             return(NULL)
                           })
    
    sqv.vec <- SQV(zq = z_q.matrix,
                   n=nrow(x),
                   p=length(start_theta))
    
    # Tolerance for diference between params vectors
    if(is.null(functions)){
      L = control$L
    }else{
      L = 2.1*median(sqv.vec)
    }
    
    # Passo 4#####################################################################
    
    check_tolerance <- sqv.vec < L
    
    check_convergence <- unlist(conv.list)
    
    if( (is.na(match(FALSE,check_tolerance))) & (is.na(match(FALSE,check_convergence))) ){
      alpha.optimal <- round(min(alpha.grade),2)
      unstable = FALSE
    }else{
      
      # alpha.start <- alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1]+0.02
      
      if(is.na(alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1])){
        alpha.start <- alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02
      }else if(is.na(alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02)){
        alpha.start <- alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1]
      }else{
        alpha.start <- max(alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1],
                           alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02)
      }
      
      
      if( ((alpha.start) >= 0.48) | (is.na(alpha.start)) ){
        alpha.optimal <- 0
        unstable = FALSE
        stability = FALSE
      }else{
        alpha.grade <- seq(from=alpha.start, by= 0.02, length.out = M)
      }
    }
    
    if(loop.number==1){
      output$enpv <- sqv.vec
    }

    loop.number = loop.number + 1
    
    # print(alpha.grade)
    # print(sqv.vec)
    # print(check_tolerance)

    # Passo 5#####################################################################
    # Repetição dos passos 2 a 4 até que a condição de estabilidade seja satisfeita
    # ou até chegar ao q maximo, que é alpha=0.2 (q=0.8)
  }
  
  output$alpha.optimal = alpha.optimal
  
  if(stability==FALSE){
    message("\nINFO: The data driven algorithm has not reached stability! Selecting Alpha optimal = ",alpha.optimal)
    output$message = "Lack of stability"
  }else{
    message("\nINFO: Selection finished! Alpha optimal = ",alpha.optimal)
    output$message = "Stability reached"
  }
  

  ##############################################################################
  # return(list("alpha.optimal" = alpha.optimal,
  #             "q.grade"=q.grade,
  #             "k"= k,
  #             "M1" = M1,
  #             "est.list" = est.list,
  #             "conv.list" = conv.list,
  #             "z_q.list" = z_q.list,
  #             "z_q.matrix" = z_q.matrix,
  #             "sqv.vec" = sqv.vec)
  # )
  
  return(output)
  
}

###################################
# using bootstrap

select_tuning_constant_std <- function(y,
                                   x,
                                   z,
                                   type,
                                   link,
                                   link.phi,
                                   start_theta,
                                   se.bootstrap,
                                   wald.test.bootstrap,
                                   functions,
                                   parameters,
                                   data,
                                   formula,
                                   control=robustbetareg.control()){
  
  message("\nINFO: Starting to run the algorithm for the optimal tuning contant selection!")
  
  # Captura dos parãmetros a serem usados no cálculo
  M = control$M
  if(is.null(functions)){
    L = control$L
  }else{
    L = 0.1
  }
  
  n = nrow(x)
  alpha.max = control$alpha.max
  unstable = TRUE
  stability = TRUE
  loop.number = 1
  output = list()
  
  # Passo 1#####################################################################
  M1 = 11
  q.grade = seq(from=1, to=0.8, by= -0.02)[1:M1]
  alpha.grade = 1-q.grade
  
  # Passo 2#####################################################################
  while (unstable) {
    
    message("\nINFO: Ciclo ",loop.number)
    
    # Criação/reset de objetos
    est.list <- list()
    est.par <- list()
    conv.list <- list()
    z_q.list <- list()
    sqv.vec <- NULL
    k = length(alpha.grade)
    
    for (k in 1:length(alpha.grade)) {
      message("\nINFO: Loop ",k)
      if(alpha.grade[k]==0){
        est.par = tryCatch(MLE.fit(y=y, 
                                   x=x, 
                                   z=z, 
                                   link = link, 
                                   link.phi = link.phi,
                                   type="MLE",
                                   se.bootstrap = se.bootstrap,
                                   wald.test.bootstrap = wald.test.bootstrap, 
                                   functions = functions, 
                                   parameters = parameters,
                                   data = data,
                                   formula = formula,
                                   control = robustbetareg.control(start = start_theta)), 
                           error = function(e) {
                             est.par$converged <- FALSE
                             return(est.par)
                           })
      }else{
        if(type=="SMLE"){
          est.par = tryCatch(SMLE.fit(y=y, 
                                      x=x, 
                                      z=z, 
                                      alpha = alpha.grade[k], 
                                      link = link, 
                                      link.phi = link.phi,
                                      type=type,
                                      se.bootstrap = se.bootstrap,
                                      wald.test.bootstrap = wald.test.bootstrap, 
                                      functions = functions, 
                                      parameters = parameters,
                                      control = robustbetareg.control(start = start_theta)), 
                             error = function(e) {
                               est.par$converged <- FALSE
                               return(est.par)
                             })
        }else if(type=="LSMLE"){
          est.par = tryCatch(LSMLE.fit(y=y, 
                                       x=x, 
                                       z=z, 
                                       alpha = alpha.grade[k], 
                                       link = link, 
                                       link.phi = link.phi,
                                       type=type,
                                       se.bootstrap = se.bootstrap,
                                       wald.test.bootstrap = wald.test.bootstrap, 
                                       functions = functions, 
                                       parameters = parameters,
                                       control = robustbetareg.control(start = start_theta)), 
                             error = function(e) {
                               est.par$converged <- FALSE
                               return(est.par)
                             })
        }
      }
      est.list[[k]] <- est.par
      
      # Passo 3#####################################################################
      
      # z_q.list[[k]] <- c(est.par$coefficients$mean,est.par$coefficients$precision)
      
      
      z_q.list[[k]] <- ((c(est.par$coefficients$mean,est.par$coefficients$precision))/
                          (sqrt(n)*c(est.par$std.error$se.mean,est.par$std.error$se.precision))
                        
      )
      
      conv.list[[k]] <- est.par$converged
    }
    
    z_q.matrix <- matrix(data=unlist(z_q.list),
                         nrow = k,
                         ncol = length(start_theta),
                         byrow = TRUE)
    
    sqv.vec <- SQV(zq = z_q.matrix,
                   n=nrow(x),
                   p=length(start_theta))
    
    # Passo 4#####################################################################
    
    check_tolerance <- sqv.vec < L
    
    check_convergence <- unlist(conv.list)
    
    if( (is.na(match(FALSE,check_tolerance))) & (is.na(match(FALSE,check_convergence))) ){
      alpha.optimal <- round(min(alpha.grade),2)
      unstable = FALSE
    }else{
      
      # alpha.start <- alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1]+0.02
      
      if(is.na(alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1])){
        alpha.start <- alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02
      }else if(is.na(alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02)){
        alpha.start <- alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1]
      }else{
        alpha.start <- max(alpha.grade[Position(isFALSE, check_tolerance, right = TRUE)+1],
                           alpha.grade[Position(isFALSE, check_convergence, right = TRUE)]+0.02)
      }
      
      
      if( ((alpha.start) >= 0.48) | (is.na(alpha.start)) ){
        alpha.optimal <- 0
        unstable = FALSE
        stability = FALSE
      }else{
        alpha.grade <- seq(from=alpha.start, by= 0.02, length.out = M)
      }
    }
    
    loop.number = loop.number + 1
    
    # print(alpha.grade)
    # print(sqv.vec)
    # print(check_tolerance)
    
    # Passo 5#####################################################################
    # Repetição dos passos 2 a 4 até que a condição de estabilidade seja satisfeita
    # ou até chegar ao q maximo, que é alpha=0.2 (q=0.8)
  }
  
  output$alpha.optimal = alpha.optimal
  
  if(stability==FALSE){
    message("\nINFO: The data driven algorithm has not reached stability! Selecting Alpha optimal = ",alpha.optimal)
    output$message = "Lack of stability"
  }else{
    message("\nINFO: Selection finished! Alpha optimal = ",alpha.optimal)
    output$message = "Stability reached"
  }
  
  
  ##############################################################################
  # return(list("alpha.optimal" = alpha.optimal,
  #             "q.grade"=q.grade,
  #             "k"= k,
  #             "M1" = M1,
  #             "est.list" = est.list,
  #             "conv.list" = conv.list,
  #             "z_q.list" = z_q.list,
  #             "z_q.matrix" = z_q.matrix,
  #             "sqv.vec" = sqv.vec)
  # )
  
  return(output)
  
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
# NonLinear Specific functions
######################################################################################################################
######################################################################################################################
######################################################################################################################


######################################################################################################################
power2 <- function(x){return(x**2)}
power3 <- function(x){return(x**3)}
power_par <- function(power,
                      base){
  return(base**power)
}

######################################################################################################################


apply_spline <- function(vector,
                         matrix,
                         functions,
                         parameters){
  
  #Return dot product of nonlinear predictor after apply an spline in one or more elements
  
  if( any(!(functions %in% c('log','exp','sqrt','power2','power3','power_par'))) ){
    stop("\nThe functions to be used in splines must be exp(), log(), sqrt(), power2(), power3() or power_par()")
  }
  
  if(length(functions)!=length(parameters)){
    stop("\nThe number of functions and parameters informed must be the same!",
         "\nNumber of functions:", length(functions),
         "\nNumber of paramenters:", length(parameters))
  }
  
  if(!is.null(parameters)){
    if( (any(parameters<1)) | (any(parameters>length(vector))) ){
      stop(paste("The 'parameters' argument must have all elements between 1 and", length(vector)))
    }
  }
  
  # browser()
  if( isFALSE(any(functions%in%c('power_par','power2','power3'))) ){
    if( (length(functions)==1) ){
      ch_functions <- sapply(c(functions), match.fun)
    }else{
      ch_functions <- sapply(functions, match.fun)  
    }
  }
  
  
  nl_product <- matrix(ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:length(vector)) {
    if(i %in% parameters){
      if( functions[[match(i,parameters)]] == 'power_par' ){
        # nl_product[,i] = power_par(power = vector[parameters],
        #                            base = matrix[,parameters])
        nl_product[,i] = matrix[,parameters]**vector[parameters]
                                   
        # nl_product[,i] = matrix[,parameters]**vector[parameters]
                                  
      }else{
        nl_product[,i] = sapply(matrix[,i]*vector[i], ch_functions[[match(i,parameters)]])
      }
    }
    else{
      nl_product[,i] = matrix[,i]*vector[i]
    }
  }
  return(matrix(rowSums(nl_product)))
}

######################################################################################################################

apply_derivative <- function(vector,
                             x_matrix,
                             functions,
                             parameters){
  
  if(length(functions)==1){
    ch_functions <- sapply(c(functions), match.fun)
  }else{
    ch_functions <- sapply(functions, match.fun)  
  }
  
  nl_matrix = matrix(nrow = nrow(x_matrix), ncol = ncol(x_matrix))
  
  for (i in 1:dim(x_matrix)[2]) {
    
    if(i %in% parameters){
      
      if(functions[[match(i,parameters)]] == "exp"){
        
        nl_matrix[,i] = sapply( (vector[i]*x_matrix[,i]), ch_functions[[match(i,parameters)]])*x_matrix[,i]
        
      }else if(functions[[match(i,parameters)]] == "log"){
        
        nl_matrix[,i] = ( 1/(vector[i]*x_matrix[,i]) )*x_matrix[,i]
        
      }else if(functions[[match(i,parameters)]] == "sqrt"){
        
        nl_matrix[,i] = ( (1/2)*(vector[i]*x_matrix[,i])**(-1/2) )*x_matrix[,i]
        
      }else if(functions[[match(i,parameters)]] == "power2"){
        
        nl_matrix[,i] = ( 2*(vector[i]*x_matrix[,i]) )*x_matrix[,i]
        
      }else if(functions[[match(i,parameters)]] == "power3"){
        
        nl_matrix[,i] = ( 3*(vector[i]*x_matrix[,i])**2 )*x_matrix[,i]
        
      }else if(functions[[match(i,parameters)]] == "power_par"){
        
        nl_matrix[,i] = (x_matrix[,i]**vector[i])*log(x_matrix[,i])
        
      }else{
        stop(paste("Calculo de derivada só está disponível para as funções 'log','exp','sqrt','power2','power3','power_par'!"))
      }
      
    }
    else{
      nl_matrix[,i] = x_matrix[,i]
    }
    
  }
  colnames(nl_matrix) <- colnames(x_matrix)
  return(nl_matrix)
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
# SMLE Specific functions
######################################################################################################################
######################################################################################################################
######################################################################################################################

###############################################################################################################################
SMLE.fit <- function(y, x, z, alpha = NULL, link = "logit", link.phi = "log", type = type,
                     se.bootstrap = se.bootstrap, wald.test.bootstrap = wald.test.bootstrap,
                     functions = functions, parameters = parameters,
                     control = robustbetareg.control(...), ...) 
{
  result = theta = list()
  # if(is.null(alpha)){
  #   control$alpha.optimal <- TRUE
  # }
  alpha.optimal = control$alpha.optimal
  start_theta = control$start
  if (!is.null(alpha)) {
    alpha.optimal = FALSE
  }
  if (is.null(alpha) & alpha.optimal == FALSE) {
    alpha = 0
  }
  linkobj = set.link(link.mu = link, link.phi = link.phi)
  k = ncol(x)
  m = ncol(z)
  if (alpha.optimal) {
    return(Opt.Tuning.SMLE(y, x, z, link, link.phi, control))
  }
  if (is.null(start_theta)) {
    mle = tryCatch(suppressWarnings(betareg.fit(x, y, z, 
                                                link = link, link.phi = link.phi)), error = function(e) NULL)
    start_theta = as.numeric(do.call("c", mle$coefficients))
    theta$x = start_theta
    theta$converged = mle$converged
  }
  q = 1 - alpha
  check = TRUE
  #browser()
  theta = tryCatch(optim(par = start_theta,
                         fn = L_q, 
                         gr = Psi_SMLE, 
                         y = y,
                         X = x,
                         Z = z, 
                         alpha = alpha, 
                         link_mu = link, 
                         link_phi = link.phi,
                         functions = functions, 
                         parameters = parameters,
                         apply_spline = apply_spline, 
                         apply_derivative = apply_derivative,
                         control = list(fnscale = -1)), error = function(e) {
                           theta$msg <- e$message
                           check <<- F
                           return(theta)
                         }, warning = function(w) {
                           theta <- suppressWarnings(optim(par = start_theta,
                                                           fn = L_q, 
                                                           gr = Psi_SMLE, 
                                                           y = y, 
                                                           X = x, 
                                                           Z = z, 
                                                           alpha = alpha, 
                                                           link_mu = link,
                                                           link_phi = link.phi, 
                                                           functions = functions,
                                                           parameters = parameters,
                                                           apply_spline = apply_spline, 
                                                           apply_derivative = apply_derivative,
                                                           control = list(fnscale = -1)))
                           return(theta)
                         })
  if (check) {
    if (theta$convergence == 0) {
      theta$converged = T
      theta$x = theta$par
    }
    else {
      theta$converged = F
      theta$x = start_theta
    }
  }
  else {
    theta$converged = F
    theta$x = start_theta
  }
  beta = theta$x[1:k]
  gamma = theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  #####################################################################################################
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta <- apply_spline(vector = beta,
                        matrix = x,
                        functions = functions, 
                        parameters = parameters)
  }else{
    eta = x %*% beta
  }
  #####################################################################################################
  
  mu_hat = linkobj$linkfun.mu$inv.link(eta)
  phi_hat = linkobj$linkfun.phi$inv.link(z %*% gamma)
  #####################################################################################################
  # VCOC original
  # M.SMLE = tryCatch(SMLE_Cov_Matrix(mu_hat, phi_hat, x, z, 
  #                                   alpha = alpha, linkobj = linkobj), error = function(e) NULL)
  #####################################################################################################
  # Ajuste da variancia
  if(se.bootstrap == TRUE){
    # browser()
    M.SMLE = tryCatch(bootstrap.std.error(n = nrow(x),
                                          beta = beta,
                                          gamma = gamma,
                                          alpha = alpha,
                                          functions = functions,
                                          parameters = parameters,
                                          type = type,
                                          start_theta = start_theta,
                                          data = data,
                                          y = y,
                                          x = x,
                                          z = z,
                                          link = link,
                                          link.phi = link.phi), error = function(e) NULL)
  }else{
    M.SMLE = tryCatch(SMLE_Cov_Matrix(mu_hat, phi_hat, x, z, 
                                      alpha = alpha, linkobj = linkobj), error = function(e) NULL)
  }
  if (is.null(M.SMLE)) {
    vcov = std.error.SMLE = NaN
  }
  else {
    if(se.bootstrap == TRUE){
      vcov = M.SMLE$vcov.bootstrap.theta
      std.error.SMLE = M.SMLE$std.bootstrap.theta
    }else{
      vcov = M.SMLE$Cov
      std.error.SMLE = M.SMLE$Std.Error
    }
  }
  rownames(vcov) <- c(colnames(x), colnames(z))
  colnames(vcov) <- c(colnames(x), colnames(z))
  ####################################################################################################
  y_star = log(y) - log(1 - y)
  str1 = str2 = NULL
  result$coefficients = coefficients
  result$vcov = vcov
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0) {
    NA
  }
  else {
    cor(eta, linkobj$linkfun.mu$linkfun(y))^2
  }
  if (!is.null(theta$iter)) {
    result$iter = theta$iter
  }
  result$converged = (theta$converged & all(!is.na(std.error.SMLE)))
  result$converged = (theta$converged)
  if (m == 1) {
    phi_predict = phi_hat[1]
  }
  if (m != 1) {
    phi_predict = phi_hat
  }
  result$fitted.values = list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start = start_theta
  weights = (dbeta(y, mu_hat * phi_hat, (1 - mu_hat) * 
                     phi_hat))^(alpha)
  result$weights = weights
  result$Tuning = alpha  
  # result$residuals = sweighted2_res(mu_hat, phi_hat, y = y,
  #                                  X = x, linkobj = linkobj)

  result$residuals = nlrobustbetareg_residuals(mu = mu_hat,
                                               phi = phi_hat,
                                               y = y,
                                               wts = weights,
                                               type = c("quantile"))
  result$n = length(mu_hat)
  result$link = link
  result$link.phi = link.phi
  result$Optimal.Tuning = alpha.optimal
  result$pseudo.r.squared = pseudor2
  result$control = control
  result$call = match.call()
  result$method = "SMLE"
  if (any(is.na(std.error.SMLE))) {
    str1 = "Standard-Error is unvailable"
  }
  if (!is.null(theta$msg)) {
    str2 = theta$msg
  }
  if (!is.null(str1) || !is.null(str2)) {
    result$message = c(str1, str2)
  }

  ####################################################################################################
  if (!any(is.na(std.error.SMLE))) {
    names(std.error.SMLE) <- c(colnames(x), colnames(z))
    se.beta = std.error.SMLE[1:k]
    se.gamma = std.error.SMLE[1:m + k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error = coef.std.error
  }
  ####################################################################################################
  # including the especific summary

  if(wald.test.bootstrap){
    
    result.boot.wald.type.test <- bootstrap.wald.test(n = nrow(x),
                                                      alpha = alpha,
                                                      beta = beta,
                                                      gamma = gamma,
                                                      data.names = colnames(data),
                                                      X = x,
                                                      Y = y,
                                                      Z = z,
                                                      type=type,
                                                      formula = formula,
                                                      std.bootstrap.theta = coef.std.error,
                                                      functions = functions,
                                                      parameters = parameters,
                                                      start_theta = start_theta,
                                                      link = link,
                                                      link.phi = link.phi)
    
    result$boot.wald.type.test$pvalue = result.boot.wald.type.test$pvalue_w.stat.bootstrap
    result$boot.wald.type.test$valid.w.stat = result.boot.wald.type.test$valid.w.stat
    result$boot.wald.type.test$B2 = result.boot.wald.type.test$B2
    result$boot.wald.type.test$B1 = result.boot.wald.type.test$B1
    
  }

  class(result) = "robustbetareg"
  return(result)
}
###############################################################################################################################

Opt.Tuning.SMLE <- function(y, x, z, link, link.phi, control) 
{
  if (missing(control)) {
    control = robustbetareg.control()
  }
  control$alpha.optimal = FALSE
  SMLE.list = SMLE.par = list()
  zq.t = NULL
  alpha_tuning = seq(0, 0.3, 0.02)
  K = length(alpha_tuning)
  M1 = 11
  M = control$M
  L = control$L
  n = length(y)
  unstable = F
  sqv.unstable = T
  if (is.null(control$start)) {
    est.log.lik = tryCatch(suppressWarnings(betareg.fit(x, 
                                                        y, z, link = link, link.phi = link.phi)), error = function(e) NULL)
    if (is.null(est.log.lik)) {
      est.log.lik = tryCatch(suppressWarnings(betareg.fit(x, 
                                                          y, z, link = link, link.phi = link.phi, control = betareg.control(start = Initial.points(y, 
                                                                                                                                                   x, z)))), error = function(e) NULL)
    }
    if (!is.null(est.log.lik)) {
      control$start = do.call("c", est.log.lik$coefficients)
      names(control$start) = c(colnames(x), colnames(z))
    }
    else {
      control$start = Initial.points(y, x, z)
      names(control$start) = c(colnames(x), colnames(z))
    }
  }
  else {
    names(control$start) = c(colnames(x), colnames(z))
  }
  p = length(control$start)
  for (k in 1:M1) {
    SMLE.par = tryCatch(SMLE.fit(y, x, z, alpha = alpha_tuning[k], 
                                 link = link, link.phi = link.phi, control = control), 
                        error = function(e) {
                          SMLE.par$converged <- FALSE
                          return(SMLE.par)
                        })
    if (!SMLE.par$converged || is.null(SMLE.par) || any(is.na(do.call("c", 
                                                                      SMLE.par$coefficients)/do.call("c", SMLE.par$std.error))) || 
        is.null(do.call("c", SMLE.par$std.error))) {
      sqv.unstable = F
      unstable = T
      break
    }
    SMLE.list[[k]] <- SMLE.par
    zq.t = unname(rbind(zq.t, do.call("c", SMLE.par$coefficients)/do.call("c", 
                                                                          SMLE.par$std.error)))
  }
  sqv = as.numeric(SQV(zq.t, n, p))
  alpha.ind = max(0, which(sqv > L))
  if (alpha.ind == 0) {
    SMLE.par.star <- tryCatch(SMLE.list[[1]], error = function(e) {
      SMLE.par.star <- SMLE.par
      SMLE.par.star$message <- "The function cannot be evaluated on initial parameters"
      return(SMLE.par.star)
    })
    SMLE.par.star$sqv = sqv
    SMLE.par.star$Optimal.Tuning = TRUE
    return(SMLE.par.star)
  }
  if (unstable) {
    SMLE.par.star <- tryCatch(SMLE.list[[1]], error = function(e) {
      SMLE.par.star <- SMLE.par
      return(SMLE.par.star)
    })
    SMLE.par.star$sqv = sqv
    SMLE.par.star$Optimal.Tuning = TRUE
    SMLE.par.star$message = "Lack of stability"
    return(SMLE.par.star)
  }
  if (alpha.ind < 8) {
    SMLE.par.star <- SMLE.list[[alpha.ind + 1]]
    SMLE.par.star$sqv = sqv
    SMLE.par.star$Optimal.Tuning = TRUE
    rm(SMLE.list)
    return(SMLE.par.star)
  }
  reached = FALSE
  k = M1 + 1
  while (sqv.unstable & !reached) {
    SMLE.par = tryCatch(SMLE.fit(y, x, z, alpha = alpha_tuning[k], 
                                 link = link, link.phi = link.phi, control = control), 
                        error = function(e) {
                          SMLE.par$converged <- FALSE
                          return(SMLE.par)
                        })
    if (!SMLE.par$converged || any(is.na(do.call("c", 
                                                 SMLE.par$coefficients)/do.call("c", SMLE.par$std.error))) || 
        is.null(do.call("c", SMLE.par$std.error)) || 
        !SMLE.par$converged) {
      unstable = T
      break
    }
    SMLE.list[[k]] = SMLE.par
    zq.t = unname(rbind(zq.t, do.call("c", SMLE.par$coefficients)/do.call("c", 
                                                                          SMLE.par$std.error)))
    sqv = as.numeric(SQV(zq.t, n, p))
    sqv.test = sqv[(k - M):(k - 1)]
    if (all(sqv.test <= L)) {
      sqv.unstable = F
    }
    k = k + 1
    if (k >= K) {
      reached = TRUE
    }
  }
  if (reached) {
    k = suppressWarnings(max(1, min(which(zoo::rollapply(sqv < 
                                                           L, M, sum) == M))) + M + 1)
  }
  if (k >= K || unstable) {
    SMLE.par.star = SMLE.list[[1]]
    SMLE.par.star$sqv = sqv
    SMLE.par.star$Optimal.Tuning = TRUE
    SMLE.par.star$message = "Lack of stability"
  }
  else {
    SMLE.par.star = SMLE.list[[(k - 1 - M)]]
    SMLE.par.star$sqv = sqv
    SMLE.par.star$Optimal.Tuning = TRUE
  }
  return(SMLE.par.star)
}

#####################################################################################################################

L_q <- function(theta, y, X, Z, alpha, link_mu, link_phi,
                functions = functions, parameters = parameters,
                apply_spline = apply_spline, apply_derivative = apply_derivative) 
{
  q = 1 - alpha
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu_q = link.model$linkfun.mu$inv.link(eta_beta)
  phi_q = link.model$linkfun.phi$inv.link(Z %*% Gamma)
    phi <- (1/q) * (phi_q - 2) + 2
  mu <- ((1/q) * (mu_q * phi_q - 1) + 1)/phi
  if ((any(mu == 0) || any(mu == 1))) {
    mu = pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
  }
  a <- mu * phi
  b <- (1 - mu) * phi
  log_likS <- sum(dbeta(y, a, b, log = F)^(alpha))
  return(log_likS)
}

#####################################################################################################################

Psi_SMLE <- function (theta, y, X, Z, alpha, link_mu, link_phi,
                      functions = functions, parameters = parameters,
                      apply_spline = apply_spline, apply_derivative = apply_derivative) 
{
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predictr
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu_hat = link.model$linkfun.mu$inv.link(eta_beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z %*% Gamma)
  psi_beta = Psi_Beta_SMLE(mu_hat, phi_hat, y, X, Z, alpha, 
                           link_mu, link_phi)
  psi_gamma = Psi_Gamma_SMLE(mu_hat, phi_hat, y, X, Z, alpha, 
                             link_mu, link_phi)
  return(c(psi_beta, psi_gamma))
}

#####################################################################################################################

Psi_Beta_SMLE <- function (mu_hat, phi_hat, y, X, Z, alpha, link_mu, link_phi,
                           Beta = Beta, functions = functions, parameters = parameters,
                           apply_derivative = apply_derivative)
{
  q = 1 - alpha
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  phi_q = q * (phi_hat - 2) + 2
  mu_q = (q * (mu_hat * phi_hat - 1) + 1)/phi_q
  y_star = log(y) - log(1 - y)
  mu_star = digamma(mu_hat * phi_hat) - digamma((1 - mu_hat) * 
                                                  phi_hat)
  Tb = (link.model$linkfun.mu$d.linkfun(mu_q))^(-1)
  #####################################################################################################
  apply_derivative = apply_derivative
  # Applying new derivative if splines are used in beta vector
  if(is.null(parameters)){
    result = t(X) %*% diag(phi_q * Tb/q) %*% (y_star - mu_star)
  }else{
    J <- apply_derivative(vector = Beta,
                          x_matrix = X,
                          functions = functions,
                          parameters = parameters)
    colnames(J) <- colnames(X)
    result = t(J) %*% diag(phi_q * Tb/q) %*% (y_star - mu_star)
  }
  #####################################################################################################
  
  return(result)
}

#####################################################################################################################


Psi_Gamma_SMLE <- function(mu_hat, phi_hat, y, X, Z, alpha, link_mu, link_phi,
                           Gamma = Gamma, functions = functions, parameters = parameters,
                           apply_derivative = apply_derivative) 
{
  q = 1 - alpha
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  phi_q = q * (phi_hat - 2) + 2
  mu_q = (q * (mu_hat * phi_hat - 1) + 1)/phi_q
  y_dagger = log(1 - y)
  y_star = log(y) - y_dagger
  mu_star = digamma(mu_hat * phi_hat) - digamma((1 - mu_hat) * 
                                                  phi_hat)
  mu_dagger = digamma((1 - mu_hat) * phi_hat) - digamma(phi_hat)
  Tg = (link.model$linkfun.phi$d.linkfun(phi_q))^1
  result = t(Z) %*% diag(mu_q * Tg/q) %*% (mu_q * (y_star - 
                                                     mu_star) + (y_dagger - mu_dagger))
  return(result)
}

#####################################################################################################################

SMLE_Cov_Matrix <- function(muhat_q, phihat_q, X, Z, alpha, linkobj) 
{
  q = 1 - alpha
  ahat_q <- muhat_q * phihat_q
  bhat_q <- (1 - muhat_q) * phihat_q
  T_1_q = diag((linkobj$linkfun.mu$d.linkfun(muhat_q))^(-1))
  T_2_q = diag((linkobj$linkfun.phi$d.linkfun(phihat_q))^(-1))
  phihat_n <- (1/q) * (phihat_q - 2) + 2
  muhat_n <- ((1/q) * (muhat_q * phihat_q - 1) + 1)/phihat_n
  ahat_n <- muhat_n * phihat_n
  bhat_n <- (1 - muhat_n) * phihat_n
  phihat_2_q <- (2 - q) * (phihat_n - 2) + 2
  muhat_2_q <- ((2 - q) * (ahat_n - 1) + 1)/phihat_2_q
  a2_qhat <- muhat_2_q * phihat_2_q
  b2_qhat <- (1 - muhat_2_q) * phihat_2_q
  mustarhat_n <- digamma(ahat_n) - digamma(bhat_n)
  mudaggerhat_n <- digamma(bhat_n) - digamma(phihat_n)
  mustarhat_2_q <- digamma(a2_qhat) - digamma(b2_qhat)
  mudaggerhat_2_q <- digamma(b2_qhat) - digamma(phihat_2_q)
  muhat_d_2_q <- muhat_q * (mustarhat_2_q - mustarhat_n) + 
    mudaggerhat_2_q - mudaggerhat_n
  m_phiq <- diag(phihat_q)
  psi1_n <- trigamma(ahat_n)
  psi2_n <- trigamma(bhat_n)
  psi3_n <- trigamma(phihat_n)
  V_n <- diag(psi1_n + psi2_n)
  B1 <- diag(exp(q * (lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - 
                   (lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
  B2 <- diag(exp(lgamma(a2_qhat) + lgamma(b2_qhat) - lgamma(phihat_2_q) - 
                   (2 * (1 - q) * (lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) + 
                      lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
  C_q_0 <- diag(phihat_q * (muhat_q * psi1_n - (1 - muhat_q) * 
                              psi2_n))
  D_q_0 <- diag((muhat_q^2) * psi1_n + ((1 - muhat_q)^2) * 
                  psi2_n - psi3_n)
  psi1_2_q <- trigamma(a2_qhat)
  psi2_2_q <- trigamma(b2_qhat)
  psi3_2_q <- trigamma(phihat_2_q)
  V_2_q <- diag(psi1_2_q + psi2_2_q)
  C_q_2_q <- diag(phihat_q * (muhat_q * psi1_2_q - (1 - muhat_q) * 
                                psi2_2_q))
  D_q_2_q <- diag((muhat_q^2) * psi1_2_q + ((1 - muhat_q)^2) * 
                    psi2_2_q - psi3_2_q)
  M1 <- diag(mustarhat_2_q - mustarhat_n)
  M2 <- diag(muhat_d_2_q)
  M3 <- diag((mustarhat_2_q - mustarhat_n) * muhat_d_2_q)
  Jq_betabeta <- as.matrix(t(X) %*% B1 %*% (T_1_q^2) %*% (m_phiq^2) %*% 
                             V_n %*% X)
  Jq_betagamma <- as.matrix(t(X) %*% B1 %*% T_1_q %*% T_2_q %*% 
                              C_q_0 %*% Z)
  Jq_gammagamma <- as.matrix(t(Z) %*% B1 %*% (T_2_q^2) %*% 
                               D_q_0 %*% Z)
  Jq = -(q^(-1)) * rbind(cbind(Jq_betabeta, Jq_betagamma), 
                         cbind(t(Jq_betagamma), Jq_gammagamma))
  Kq_betabeta <- as.matrix(t(X) %*% B2 %*% (T_1_q^2) %*% (m_phiq^2) %*% 
                             (V_2_q + M1^2) %*% X)
  Kq_betagamma <- as.matrix(t(X) %*% B2 %*% T_1_q %*% T_2_q %*% 
                              m_phiq %*% ((C_q_2_q * (1/phihat_q)) + M3) %*% Z)
  Kq_gammagamma <- as.matrix(t(Z) %*% B2 %*% (T_2_q^2) %*% 
                               (D_q_2_q + M2^2) %*% Z)
  Kq = (q^(-2)) * rbind(cbind(Kq_betabeta, Kq_betagamma), cbind(t(Kq_betagamma), 
                                                                Kq_gammagamma))
  inv.Jq = tryCatch(solve(Jq), error = function(e) {
    e
  })
  if (!BBmisc::is.error(inv.Jq)) {
    Vq <- inv.Jq %*% Kq %*% t(inv.Jq)
  }
  else {
    inv.Jq = MASS::ginv(Jq)
    Vq <- inv.Jq %*% Kq %*% t(inv.Jq)
  }
  result = list()
  result$Lambda = Jq
  result$Sigma = Kq
  result$Cov = Vq
  result$Std.Error = suppressWarnings(t(sqrt(diag(Vq))))
  return(result)
}

######################################################################################################################


######################################################################################################################
######################################################################################################################
######################################################################################################################
# LSMLE Specific functions
######################################################################################################################
######################################################################################################################
######################################################################################################################

LSMLE.fit <- function (y, x, z, alpha = NULL, link = "logit", link.phi = "log", type = type,
                       se.bootstrap = se.bootstrap, wald.test.bootstrap = wald.test.bootstrap,
                       functions = functions, parameters = parameters,
                       control = robustbetareg.control(...), ...) 
{
  result = theta = list()
  alpha.optimal = control$alpha.optimal
  start_theta = control$start
  if (!is.null(alpha)) {
    alpha.optimal = FALSE
  }
  if (is.null(alpha) & alpha.optimal == FALSE) {
    alpha = 0
  }
  linkobj = set.link(link.mu = link, link.phi = link.phi)
  k = ncol(x)
  m = ncol(z)
  if (alpha.optimal) {
    return(Opt.Tuning.LSMLE(y, x, z, link, link.phi, control))
  }
  if (is.null(start_theta)) {
    mle = tryCatch(suppressWarnings(betareg.fit(x, y, z, 
                                                link = link, link.phi = link.phi)), error = function(e) NULL)
    start_theta = as.numeric(do.call("c", mle$coefficients))
    theta$x = start_theta
    theta$converged = mle$converged
  }
  q = 1 - alpha
  check = TRUE
  theta = tryCatch(optim(par = start_theta, 
                         fn = L_alpha_R, 
                         gr = Psi_LSMLE, 
                         y = y, 
                         X = x,
                         Z = z, 
                         alpha = alpha,
                         link_mu = link, 
                         link_phi = link.phi, 
                         functions = functions, 
                         parameters = parameters,
                         apply_spline = apply_spline, 
                         apply_derivative = apply_derivative,
                         control = list(fnscale = -1,
                                        maxit = 10000)), 
                   error = function(e) {
                     theta$msg <- e$message
                     check <<- F
                     return(theta)
                   })
  if (check) {
    if (theta$convergence == 0) {
      theta$converged = T
      theta$x = theta$par
    }
    else {
      theta$converged = F
      theta$x = start_theta
    }
  }
  else {
    theta$converged = F
    theta$x = start_theta
  }
  # browser()
  beta = theta$x[1:k]
  gamma = theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  #####################################################################################################
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta <- apply_spline(vector = beta,
                        matrix = x,
                        functions = functions, 
                        parameters = parameters)
  }else{
    eta = x %*% beta
  }
  #####################################################################################################  
  mu_hat = linkobj$linkfun.mu$inv.link(eta)
  phi_hat = linkobj$linkfun.phi$inv.link(z %*% gamma)
  # VCOC original
  #M.LSMLE = tryCatch(LSMLE_Cov_Matrix(mu_hat, phi_hat, x, z, 
  #                                    alpha = alpha, linkobj = linkobj), error = function(e) NULL)
  #####################################################################################################
  # Ajuste da variancia
  if(se.bootstrap == TRUE){
    
    M.LSMLE = tryCatch(bootstrap.std.error(n = nrow(x),
                                           beta = beta,
                                           gamma = gamma,
                                           alpha = alpha,
                                           functions = functions,
                                           parameters = parameters,
                                           type = type,
                                           start_theta = start_theta,
                                           data = data,
                                           y = y,
                                           x = x,
                                           z = z,
                                           link = link,
                                           link.phi = link.phi), error = function(e) NULL)
  }else{
    M.LSMLE = tryCatch(LSMLE_Cov_Matrix(mu_hat, phi_hat, x, z, 
                                        alpha = alpha, linkobj = linkobj), error = function(e) NULL)
  }
  
  if (is.null(M.LSMLE)) {
    vcov = std.error.LSMLE = NaN
  }
  else {
    if(se.bootstrap == TRUE){
      vcov = M.LSMLE$vcov.bootstrap.theta
      std.error.LSMLE = M.LSMLE$std.bootstrap.theta
    }else{
      vcov = M.LSMLE$Cov
      std.error.LSMLE = M.LSMLE$Std.Error
    }
  }
  rownames(vcov) <- c(colnames(x), colnames(z))
  colnames(vcov) <- c(colnames(x), colnames(z))
  ####################################################################################################
  y_star = log(y) - log(1 - y)
  str1 = str2 = NULL
  result$coefficients = coefficients
  result$vcov = vcov
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0) {
    NA
  }
  else {
    cor(eta, linkobj$linkfun.mu$linkfun(y))^2
  }
  if (!is.null(theta$iter)) {
    result$iter = theta$iter
  }
  result$converged = (theta$converged & all(!is.na(std.error.LSMLE)))
  if (m == 1) {
    phi_predict = phi_hat[1]
  }
  if (m != 1) {
    phi_predict = phi_hat
  }
  result$fitted.values = list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start = start_theta
  weights = (dEGB(y_star, mu_hat, phi_hat/q))^(alpha)
  result$weights = weights
  result$Tuning = alpha
  # result$residuals = sweighted2_res(mu_hat, phi_hat, y = y, 
  #                                   X = x, linkobj = linkobj)
  
  result$residuals = nlrobustbetareg_residuals(mu = mu_hat,
                                               phi = phi_hat,
                                               y = y,
                                               wts = weights,
                                               type = c("quantile"))
  result$n = length(mu_hat)
  result$link = link
  result$link.phi = link.phi
  result$Optimal.Tuning = alpha.optimal
  result$pseudo.r.squared = pseudor2
  result$control = control
  result$call = match.call()
  result$method = "LSMLE"
  if (any(is.na(std.error.LSMLE))) {
    str1 = "Standard-Error is unvailable"
  }
  if (!is.null(theta$msg)) {
    str2 = theta$msg
  }
  if (!is.null(str1) || !is.null(str2)) {
    result$message = c(str1, str2)
  }
  ####################################################################################################
  if (!any(is.na(std.error.LSMLE))) {
    names(std.error.LSMLE) <- c(colnames(x), colnames(z))
    se.beta = std.error.LSMLE[1:k]
    se.gamma = std.error.LSMLE[1:m + k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error = coef.std.error
  }
  ####################################################################################################
  # including the especific summary

  if(wald.test.bootstrap){
    
    result.boot.wald.type.test <- bootstrap.wald.test(n = nrow(x),
                                                      alpha = alpha,
                                                      beta = beta,
                                                      gamma = gamma,
                                                      data.names = colnames(data),
                                                      X = x,
                                                      Y = y,
                                                      Z = z,
                                                      type=type,
                                                      formula = formula,
                                                      std.bootstrap.theta = coef.std.error,
                                                      functions = functions,
                                                      parameters = parameters,
                                                      start_theta = start_theta,
                                                      link = link,
                                                      link.phi = link.phi)
    
    result$boot.wald.type.test$pvalue = result.boot.wald.type.test$pvalue_w.stat.bootstrap
    result$boot.wald.type.test$valid.w.stat = result.boot.wald.type.test$valid.w.stat
    result$boot.wald.type.test$B2 = result.boot.wald.type.test$B2
    result$boot.wald.type.test$B1 = result.boot.wald.type.test$B1
    
  }
  
  class(result) = "robustbetareg"
  return(result)
}

######################################################################################################################

LSMLE_Cov_Matrix <- function (mu, phi, X, Z, alpha, linkobj) 
{
  n = length(mu)
  q = 1 - alpha
  a.0 = mu * phi
  b.0 = (1 - mu) * phi
  phi.k0 = phi/(1 - alpha)
  a.k0 = mu * phi.k0
  b.k0 = (1 - mu) * phi.k0
  phi.k1 = phi * (1 + alpha)/(1 - alpha)
  a.k1 = mu * phi.k1
  b.k1 = (1 - mu) * phi.k1
  C.k0 = diag(exp(q * lbeta(a.k0, b.k0) - lbeta(a.0, b.0)))
  C.k1 = diag(exp(lbeta(a.k1, b.k1) - lbeta(a.0, b.0) - 2 * 
                    alpha * lbeta(a.k0, b.k0)))
  mu_star.k0 = digamma(a.k0) - digamma(b.k0)
  mu_dagger.k0 = digamma(b.k0) - digamma(phi.k0)
  mu_star.k1 = digamma(a.k1) - digamma(b.k1)
  mu_dagger.k1 = digamma(b.k1) - digamma(phi.k1)
  Tb = diag((linkobj$linkfun.mu$d.linkfun(mu))^(-1))
  Tg = diag((linkobj$linkfun.phi$d.linkfun(phi))^(-1))
  Phi = diag(phi)
  Q.inv = diag(n)/q
  Lambda_mu_mu = diag(trigamma(a.k0) + trigamma(b.k0))
  Lambda_mu_phi = diag(mu * trigamma(a.k0) - (1 - mu) * trigamma(b.k0))
  Lambda_phi_phi = diag(mu^2 * trigamma(a.k0) + (1 - mu)^2 * 
                          trigamma(b.k0) - trigamma(phi.k0))
  Lambda_beta_beta = t(X) %*% Tb %*% Q.inv %*% (Phi^2) %*% 
    C.k0 %*% Lambda_mu_mu %*% Tb %*% X
  Lambda_beta_gamma = t(X) %*% Tb %*% Q.inv %*% Phi %*% C.k0 %*% 
    Lambda_mu_phi %*% Tg %*% Z
  Lambda_gamma_gamma = t(Z) %*% Tg %*% Q.inv %*% C.k0 %*% Lambda_phi_phi %*% 
    Tg %*% Z
  Lambda = -1 * rbind(cbind(Lambda_beta_beta, Lambda_beta_gamma), 
                      cbind(t(Lambda_beta_gamma), Lambda_gamma_gamma))
  Sigma_mu_mu = diag(trigamma(a.k1) + trigamma(b.k1) + (mu_star.k1 - 
                                                          mu_star.k0)^2)
  Sigma_mu_phi = diag(mu * trigamma(a.k1) - (1 - mu) * trigamma(b.k1) + 
                        mu * (mu_star.k1 - mu_star.k0)^2 + (mu_star.k1 - mu_star.k0) * 
                        (mu_dagger.k1 - mu_dagger.k0))
  Sigma_phi_phi = diag(((mu * (mu_star.k1 - mu_star.k0) + (mu_dagger.k1 - 
                                                             mu_dagger.k0))^2) + (mu^2) * trigamma(a.k1) + ((1 - mu)^2) * 
                         trigamma(b.k1) - trigamma(phi.k1))
  Sigma_beta_beta = t(X) %*% Tb %*% (Q.inv^2) %*% (Phi^2) %*% 
    C.k1 %*% Sigma_mu_mu %*% Tb %*% X
  Sigma_beta_gamma = t(X) %*% Tb %*% (Q.inv^2) %*% Phi %*% 
    C.k1 %*% Sigma_mu_phi %*% Tg %*% Z
  Sigma_gamma_gamma = t(Z) %*% Tg %*% (Q.inv^2) %*% C.k1 %*% 
    Sigma_phi_phi %*% Tg %*% Z
  Sigma = rbind(cbind(Sigma_beta_beta, Sigma_beta_gamma), cbind(t(Sigma_beta_gamma), 
                                                                Sigma_gamma_gamma))
  inv.Lambda = tryCatch(solve(Lambda), error = function(e) {
    e
  })
  if (!BBmisc::is.error(inv.Lambda)) {
    V = n * inv.Lambda %*% Sigma %*% t(inv.Lambda)
  }
  else {
    inv.Lambda = MASS::ginv(Lambda)
    V = n * inv.Lambda %*% Sigma %*% t(inv.Lambda)
  }
  result = list()
  result$Lambda = Lambda
  result$Sigma = Sigma
  result$Cov = V/n
  result$K = t(X) %*% Tb %*% Q.inv %*% (Phi^2) %*% C.k0 %*% 
    Lambda_mu_mu %*% Tb %*% X
  result$Std.Error = suppressWarnings(c(sqrt(diag(V/n))))
  return(result)
}

######################################################################################################################
Psi_LSMLE <- function (theta, y, X, Z, alpha, link_mu, link_phi,
                       functions = functions, parameters = parameters,
                       apply_spline = apply_spline, apply_derivative = apply_derivative) 
{
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predict
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu_hat = link.model$linkfun.mu$inv.link(eta_beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z %*% Gamma)
  psi_beta = Psi_Beta_LSMLE(mu_hat, phi_hat, y, X, Z, alpha, 
                            link_mu, link_phi)
  psi_gamma = Psi_Gamma_LSMLE(mu_hat, phi_hat, y, X, Z, alpha, 
                              link_mu, link_phi)
  return(c(psi_beta, psi_gamma))
}

######################################################################################################################
Psi_Beta_SMLE
Psi_Beta_LSMLE <- function (mu_hat, phi_hat, y, X, Z, alpha, link_mu, link_phi,
                            Beta = Beta, functions = functions, parameters = parameters,
                            apply_derivative = apply_derivative) 
{
  q = 1 - alpha
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  phi_q = phi_hat/q
  aq = mu_hat * phi_q
  bq = (1 - mu_hat) * phi_q
  y_star = log(y) - log(1 - y)
  mu_star = digamma(aq) - digamma(bq)
  Tb = (link.model$link.mu$d.linkfun(mu_hat))^(-1)
  f_q_star = (dEGB(y_star, mu_hat, phi_q))^(alpha)
  #####################################################################################################
  apply_derivative = apply_derivative
  # Applying new derivative if splines are used in beta vector
  if(is.null(parameters)){
    result = t(X) %*% diag(phi_q * Tb * f_q_star) %*% (y_star - 
                                                         mu_star)
  }else{
    J <- apply_derivative(vector = Beta,
                          x_matrix = X,
                          functions = functions,
                          parameters = parameters)
    colnames(J) <- colnames(X)
    result = t(J) %*% diag(phi_q * Tb * f_q_star) %*% (y_star - 
                                                         mu_star)
  }
  #####################################################################################################
  return(result)
}

######################################################################################################################

Psi_Gamma_LSMLE <- function (mu_hat, phi_hat, y, X, Z, alpha, link_mu, link_phi) 
{
  q = 1 - alpha
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  phi_q = phi_hat/q
  aq = mu_hat * phi_q
  bq = (1 - mu_hat) * phi_q
  y_dagger = log(1 - y)
  y_star = log(y) - y_dagger
  mu_star = digamma(aq) - digamma(bq)
  mu_dagger = digamma(bq) - digamma(phi_q)
  Tg = (link.model$linkfun.phi$d.linkfun(phi_hat))^1
  f_q_star = (dEGB(y_star, mu_hat, phi_q))^(alpha)
  eta = mu_hat * (y_star - mu_star) + (y_dagger - mu_dagger)
  result = t(Z) %*% diag(Tg %*% f_q_star/q) %*% eta
  return(result)
}
######################################################################################################################

Opt.Tuning.LSMLE <- function (y, x, z, link, link.phi, control) 
{
  if (missing(control)) {
    control = robustbetareg.control()
  }
  control$alpha.optimal = FALSE
  LSMLE.list = LSMLE.par = list()
  zq.t = NULL
  alpha_tuning = seq(0, 0.6, 0.02)
  K = length(alpha_tuning)
  M1 = 11
  M = control$M
  L = control$L
  n = length(y)
  unstable = F
  sqv.unstable = T
  if (is.null(control$start)) {
    est.log.lik = tryCatch(suppressWarnings(betareg.fit(x, 
                                                        y, z, link = link, link.phi = link.phi)), error = function(e) NULL)
    if (is.null(est.log.lik)) {
      est.log.lik = tryCatch(suppressWarnings(betareg.fit(x, 
                                                          y, z, link = link, link.phi = link.phi, control = betareg.control(start = Initial.points(y, 
                                                                                                                                                   x, z)))), error = function(e) NULL)
    }
    if (!is.null(est.log.lik)) {
      control$start = do.call("c", est.log.lik$coefficients)
      names(control$start) = c(colnames(x), colnames(z))
    }
    else {
      control$start = Initial.points(y, x, z)
      names(control$start) = c(colnames(x), colnames(z))
    }
  }
  else {
    names(control$start) = c(colnames(x), colnames(z))
  }
  p = length(control$start)
  control.temp = control
  for (k in 1:M1) {
    LSMLE.par = tryCatch(LSMLE.fit(y, x, z, alpha = alpha_tuning[k], 
                                   link = link, link.phi = link.phi, control = control.temp), 
                         error = function(e) {
                           LSMLE.par$converged <- FALSE
                           return(LSMLE.par)
                         })
    if (LSMLE.par$converged) {
      control.temp$start = c(LSMLE.par$coefficients$mean, 
                             LSMLE.par$coefficients$precision)
    }
    if (!LSMLE.par$converged || is.null(LSMLE.par) || any(is.na(do.call("c", 
                                                                        LSMLE.par$coefficients)/do.call("c", LSMLE.par$std.error))) || 
        is.null(do.call("c", LSMLE.par$std.error))) {
      sqv.unstable = F
      unstable = T
      break
    }
    LSMLE.list[[k]] <- LSMLE.par
    zq.t = unname(rbind(zq.t, do.call("c", LSMLE.par$coefficients)/do.call("c", 
                                                                           LSMLE.par$std.error)))
  }
  sqv = as.numeric(SQV(zq.t, n, p))
  alpha.ind = max(0, which(sqv > L))
  if (alpha.ind == 0) {
    LSMLE.par.star <- tryCatch(LSMLE.list[[1]], error = function(e) {
      LSMLE.par.star <- LSMLE.par
      LSMLE.par.star$message <- "The function cannot be evaluated on initial parameters"
      return(LSMLE.par.star)
    })
    LSMLE.par.star$sqv = sqv
    LSMLE.par.star$Optimal.Tuning = TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  if (unstable) {
    LSMLE.par.star <- tryCatch(LSMLE.list[[1]], error = function(e) {
      LSMLE.par.star <- LSMLE.par
      return(LSMLE.par.star)
    })
    LSMLE.par.star$sqv = sqv
    LSMLE.par.star$Optimal.Tuning = TRUE
    LSMLE.par.star$message = "Lack of stability"
    return(LSMLE.par.star)
  }
  if (alpha.ind < 8) {
    LSMLE.par.star <- LSMLE.list[[alpha.ind + 1]]
    LSMLE.par.star$sqv = sqv
    LSMLE.par.star$Optimal.Tuning = TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  reached = FALSE
  k = M1 + 1
  while (sqv.unstable & !reached) {
    LSMLE.par = tryCatch(LSMLE.fit(y, x, z, alpha = alpha_tuning[k], 
                                   link = link, link.phi = link.phi, control = control.temp), 
                         error = function(e) {
                           LSMLE.par$converged <- FALSE
                           return(LSMLE.par)
                         })
    if (LSMLE.par$converged) {
      control.temp$start = c(LSMLE.par$coefficients$mean, 
                             LSMLE.par$coefficients$precision)
    }
    if (any(is.na(do.call("c", LSMLE.par$coefficients)/do.call("c", 
                                                               LSMLE.par$std.error))) || is.null(do.call("c", 
                                                                                                         LSMLE.par$std.error))) {
      unstable = T
      break
    }
    LSMLE.list[[k]] = LSMLE.par
    zq.t = unname(rbind(zq.t, do.call("c", LSMLE.par$coefficients)/do.call("c", 
                                                                           LSMLE.par$std.error)))
    sqv = as.numeric(SQV(zq.t, n, p))
    sqv.test = sqv[(k - M):(k - 1)]
    if (all(sqv.test <= L)) {
      sqv.unstable = F
    }
    k = k + 1
    if (k >= K) {
      reached = TRUE
    }
  }
  if (reached) {
    k = suppressWarnings(max(1, min(which(zoo::rollapply(sqv < 
                                                           L, M, sum) == M))) + M + 1)
  }
  if (k >= K || unstable) {
    LSMLE.par.star = LSMLE.list[[1]]
    LSMLE.par.star$sqv = sqv
    LSMLE.par.star$Optimal.Tuning = TRUE
    LSMLE.par.star$message = "Lack of stability"
  }
  else {
    LSMLE.par.star = LSMLE.list[[(k - 1 - M)]]
    LSMLE.par.star$sqv = sqv
    LSMLE.par.star$Optimal.Tuning = TRUE
  }
  return(LSMLE.par.star)
}

######################################################################################################################

dEGB <- function (y_star, mu, phi, log = FALSE) 
{
  if (any((-abs(2 * mu - 1) + 1) <= 0)) {
    return(warning("'mu' parameter must be within unit interval"))
  }
  if (any(phi <= 0)) {
    return(warning("'phi' parameter must be a positive value"))
  }
  a0 = mu * phi
  b0 = (1 - mu) * phi
  if (log == F) {
    k = tryCatch(exp(-(lgamma(a0) + lgamma(b0) - lgamma(a0 + 
                                                          b0) + b0 * y_star + phi * Rmpfr::log1pexp(-y_star))), 
                 error = function(e) {
                   stop("Error")
                 })
  }
  if (log == T) {
    k = tryCatch(-(lgamma(a0) + lgamma(b0) - lgamma(a0 + 
                                                      b0) + b0 * y_star + phi * Rmpfr::log1pexp(-y_star)), 
                 error = function(e) {
                   stop("Error")
                 })
  }
  return(k)
}

######################################################################################################################

L_alpha_R <- function (theta, y, X, Z, alpha, link_mu, link_phi,
                       functions = functions, parameters = parameters,
                       apply_spline = apply_spline, apply_derivative = apply_derivative) 
{
  q = 1 - alpha
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu_hat = link.model$linkfun.mu$inv.link(eta_beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z %*% Gamma)
  phi_q = phi_hat/q
  y_star = log(y) - log(1 - y)
  f_q_star = dEGB(y_star, mu_hat, phi_q)
  if (alpha == 0) {
    L_q = sum(log(f_q_star))
  }
  else {
    L_q = sum((f_q_star^(alpha) - 1)/alpha)
  }
  return(L_q)
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
# Maximum Likelyhood Estimators (MLE) especific functions
######################################################################################################################
######################################################################################################################
######################################################################################################################

MLE.fit <- function(y, x, z, link = "logit", link.phi = "log", type = type,
                    se.bootstrap = se.bootstrap, wald.test.bootstrap = wald.test.bootstrap,
                    functions = functions, parameters = parameters, data = data, formula = formula,
                    control = robustbetareg.control(...), ...) 
{

  result = theta = list()
  #alpha.optimal = control$alpha.optimal
  start_theta = control$start
  alpha = 0
  #if (!is.null(alpha)) {
  #  alpha.optimal = FALSE
  #}
  #if (is.null(alpha) & alpha.optimal == FALSE) {
  #  alpha = 0
  #}
  linkobj = set.link(link.mu = link, link.phi = link.phi)
  k = ncol(x)
  m = ncol(z)
  #if (alpha.optimal) {
  #  return(Opt.Tuning.SMLE(y, x, z, link, link.phi, control))
  #}
  if (is.null(start_theta)) {
    mle = tryCatch(suppressWarnings(betareg.fit(x, y, z, 
                                                link = link, link.phi = link.phi)), error = function(e) NULL)
    start_theta = as.numeric(do.call("c", mle$coefficients))
    theta$x = start_theta
    theta$converged = mle$converged
  }
  #q = 1 - alpha
  # browser()
  check = TRUE
  if(!is.null(parameters)){
    theta = tryCatch(optim(par = start_theta,
                           fn = Loglik_MLE,
                           gr = Psi_MLE,
                           y = y,
                           X = as.matrix(x),
                           Z = as.matrix(z), 
                           link_mu = link, 
                           link_phi = link.phi,
                           functions = functions, 
                           parameters = parameters,
                           #data = data,
                           apply_spline = apply_spline, 
                           apply_derivative = apply_derivative,
                           control = list(fnscale = -1)), error = function(e) {
                             theta$msg <- e$message
                             check <<- F
                             return(theta)
                           }, warning = function(w) {
                             theta <- suppressWarnings(optim(par = start_theta,
                                                             fn = Loglik_MLE,
                                                             gr = Psi_MLE,
                                                             y = y, 
                                                             X = as.matrix(x),
                                                             Z = as.matrix(z), 
                                                             link_mu = link,
                                                             link_phi = link.phi, 
                                                             functions = functions,
                                                             parameters = parameters,
                                                             #data = data,
                                                             apply_spline = apply_spline, 
                                                             apply_derivative = apply_derivative,
                                                             control = list(fnscale = -1)))
                             return(theta)
                           })
  }else{
    betareg::betareg(formula = formula, 
                     data = data,
                     link = link, 
                     link.phi = link.phi,
                     model = TRUE) -> linear_fit
    
    theta$par = linear_fit$optim$par
    theta$convergence = linear_fit$optim$convergence
    theta$message = linear_fit$optim$message

  }

  if (check) {
    if (theta$convergence == 0) {
      theta$converged = T
      theta$x = theta$par
    }
    else {
      theta$converged = F
      theta$x = start_theta
    }
  }
  else {
    theta$converged = F
    theta$x = start_theta
  }
  beta = theta$x[1:k]
  gamma = theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  #####################################################################################################
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta <- apply_spline(vector = beta,
                        matrix = x,
                        functions = functions, 
                        parameters = parameters)
  }else{
    eta = x %*% beta
  }
  #####################################################################################################
  
  mu_hat = linkobj$linkfun.mu$inv.link(eta)
  phi_hat = linkobj$linkfun.phi$inv.link(z %*% gamma)
  #####################################################################################################
  # VCOC original
  # M.MLE = tryCatch(SMLE_Cov_Matrix(mu_hat, phi_hat, x, z, 
  #                                   alpha = alpha, linkobj = linkobj), error = function(e) NULL)
  #####################################################################################################
  # Mle para extrair o erro padrão assintótico
  # mle_fit = tryCatch(suppressWarnings(betareg.fit(x, y, z, 
  #                                                 link = link, link.phi = link.phi)), error = function(e) NULL)
  # mle_fit.coefs = as.numeric(do.call("c", mle_fit$coefficients))

  # Ajuste da variancia
  if(se.bootstrap == TRUE){

    M.MLE = tryCatch(bootstrap.std.error(n = nrow(x),
                                         beta = beta,
                                         gamma = gamma,
                                         alpha=NULL,
                                         functions = functions,
                                         parameters = parameters,
                                         type = type,
                                         start_theta = start_theta,
                                         data=data,
                                         y = y,
                                         x = x,
                                         z = z,
                                         link = link,
                                         link.phi = link.phi), 
                     error = function(e) NULL)

  }else{
    betareg::betareg(formula = formula, 
                     data = data,
                     link = link, 
                     link.phi = link.phi,
                     model = TRUE) -> linear_fit
    M.MLE = list(Std.Error = sqrt(diag(linear_fit$vcov)),
                 Cov = linear_fit$vcov)
  }
  if (is.null(M.MLE)) {
    vcov = std.error.MLE = NaN
  }
  else {
    if(se.bootstrap == TRUE){
      vcov = M.MLE$vcov.bootstrap.theta
      std.error.MLE = M.MLE$std.bootstrap.theta
      }else{
      vcov = M.MLE$Cov
      std.error.MLE = M.MLE$Std.Error
    }
  }
  if( (!is.null(M.MLE)) & (se.bootstrap) ){
    rownames(vcov) <- c(colnames(x), colnames(z))
    colnames(vcov) <- c(colnames(x), colnames(z))
  }
  
  ####################################################################################################
  y_star = log(y) - log(1 - y)
  str1 = str2 = NULL
  result$coefficients = coefficients
  result$vcov = vcov
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0) {
    NA
  }
  else {
    cor(eta, linkobj$linkfun.mu$linkfun(y))^2
  }
  if (!is.null(theta$iter)) {
    result$iter = theta$iter
  }
  result$converged = (theta$converged & all(!is.na(std.error.MLE)))
  result$converged = (theta$converged)
  if (m == 1) {
    phi_predict = phi_hat[1]
  }
  if (m != 1) {
    phi_predict = phi_hat
  }
  result$fitted.values = list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start = start_theta
  weights = rep(1, length(y))
  result$weights = weights
  result$Tuning = 0
  # result$residuals = sweighted2_res(mu_hat, phi_hat, y = y, 
  #                                   X = x, linkobj = linkobj)
  
  result$residuals = nlrobustbetareg_residuals(mu = mu_hat,
                                               phi = phi_hat,
                                               y = y,
                                               wts = weights,
                                               type = c("quantile"))
  result$n = length(mu_hat)
  result$link = link
  result$link.phi = link.phi
  result$Optimal.Tuning = 0
  result$pseudo.r.squared = pseudor2
  result$control = control
  result$call = match.call()
  result$method = "MLE"
  if (any(is.na(std.error.MLE))) {
    str1 = "Standard-Error is unvailable"
  }
  if (!is.null(theta$msg)) {
    str2 = theta$msg
  }
  if (!is.null(str1) || !is.null(str2)) {
    result$message = c(str1, str2)
  }
  ####################################################################################################
  if (!any(is.na(std.error.MLE))) {
    names(std.error.MLE) <- c(colnames(x), colnames(z))
    se.beta = std.error.MLE[1:k]
    se.gamma = std.error.MLE[1:m + k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error = coef.std.error
  }
  ####################################################################################################
  # including the especific summary

  if(wald.test.bootstrap){
    
    result.boot.wald.type.test <- bootstrap.wald.test(n = nrow(x),
                                                      alpha = alpha,
                                                      beta = beta,
                                                      gamma = gamma,
                                                      data.names = colnames(data),
                                                      X = x,
                                                      Y = y,
                                                      Z = z,
                                                      type=type,
                                                      formula = formula,
                                                      std.bootstrap.theta = coef.std.error,
                                                      functions = functions,
                                                      parameters = parameters,
                                                      start_theta = start_theta,
                                                      link = link,
                                                      link.phi = link.phi)
    
    result$boot.wald.type.test$pvalue = result.boot.wald.type.test$pvalue_w.stat.bootstrap
    result$boot.wald.type.test$valid.w.stat = result.boot.wald.type.test$valid.w.stat
    result$boot.wald.type.test$B2 = result.boot.wald.type.test$B2
    result$boot.wald.type.test$B1 = result.boot.wald.type.test$B1

  }
  
  class(result) = "robustbetareg"
  return(result)
}

######################################################################################################################
Psi_MLE <- function (theta, y, X, Z, link_mu, link_phi,
                      functions = functions, parameters = parameters,
                      apply_spline = apply_spline, apply_derivative = apply_derivative) 
{
 
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predict
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu_hat = link.model$linkfun.mu$inv.link(eta_beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z %*% Gamma)
  psi_beta = Psi_Beta_MLE(mu_hat, phi_hat, y, X, Z,
                           link_mu, link_phi)
  psi_gamma = Psi_Gamma_MLE(mu_hat, phi_hat, y, X, Z,
                             link_mu, link_phi)
  return(c(psi_beta, psi_gamma))
}

######################################################################################################################
#ok
Psi_Beta_MLE <- function (mu_hat, phi_hat, y, X, Z, link_mu, link_phi,
                          Beta = Beta, functions = functions, parameters = parameters,
                          apply_derivative = apply_derivative)
{
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  y_star = log(y) - log(1 - y)
  mu_star = base::digamma(mu_hat*phi_hat) - base::digamma((1 - mu_hat)*phi_hat)
  Tb = (link.model$linkfun.mu$d.linkfun(mu_hat))^(-1)
  #####################################################################################################
  apply_derivative = apply_derivative
  # Applying new derivative if splines are used in beta vector
  if(is.null(parameters)){
    result = t(X) %*% diag(phi_hat * Tb) %*% (y_star - mu_star)
  }else{
    J <- apply_derivative(vector = Beta,
                          x_matrix = X,
                          functions = functions,
                          parameters = parameters)
    colnames(J) <- colnames(X)
    result = t(J) %*% diag(phi_hat * Tb) %*% (y_star - mu_star)
  }
  #####################################################################################################
  
  return(result)
}

#####################################################################################################################

#ok
Psi_Gamma_MLE <- function(mu_hat, phi_hat, y, X, Z, link_mu, link_phi,
                          Gamma = Gamma, functions = functions, parameters = parameters,
                          apply_derivative = apply_derivative) 
{
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  y_dagger = log(1 - y)
  y_star = log(y) - y_dagger
  mu_star = base::digamma(mu_hat*phi_hat)-base::digamma((1-mu_hat)*phi_hat)
  mu_dagger = base::digamma((1 - mu_hat) * phi_hat) - base::digamma(phi_hat)
  Tg = (link.model$linkfun.phi$d.linkfun(phi_q))^1
  
  result = t(Z) %*% diag(Tg) %*% (mu_hat*(y_star-mu_star) + (y_dagger-mu_dagger))
  
  return(result)
}

#####################################################################################################################

Loglik_MLE <- function(theta, y, X, Z, link_mu, link_phi,
                       functions = functions, parameters = parameters,
                       apply_spline = apply_spline, 
                       apply_derivative = apply_derivative) 
{
  #q = 1 - alpha
  k = ncol(X)
  m = ncol(Z)
  link.model = set.link(link.mu = link_mu, link.phi = link_phi)
  Beta = theta[1:k]
  #####################################################################################################
  apply_spline = apply_spline
  apply_derivative = apply_derivative
  # Applying splines in beta predictor
  if(!is.null(parameters)){
    eta_beta <- apply_spline(vector = Beta, 
                             matrix = X,
                             functions = functions, 
                             parameters = parameters)
  }else{
    eta_beta <- X %*% Beta
  }
  #####################################################################################################
  Gamma = theta[(k + 1):(k + m)]
  mu = link.model$linkfun.mu$inv.link(eta_beta)
  phi = link.model$linkfun.phi$inv.link(Z %*% Gamma)
  if ((any(mu == 0) || any(mu == 1))) {
    mu = pmax(pmin(mu, 1 - .Machine$double.eps), .Machine$double.eps)
  }
  a <- mu * phi
  b <- (1 - mu) * phi
  log_likS <- sum(dbeta(y, a, b, log = T))
  return(log_likS)
}

######################################################################################################################
  
######################################################################################################################

#####################################################################################################################

######################################################################################################################

######################################################################################################################


######################################################################################################################
######################################################################################################################
######################################################################################################################
#Simulation functions

################################################################################
# generate simulated parametric samples

FunSample_cont_betareg <- function(n,
                                   X,
                                   Z,
                                   beta,
                                   gamma,
                                   functions,
                                   parameters,
                                   contamination = 0.05,
                                   scenario,
                                   link.mu='logit',
                                   link.phi='log'){
  
  
  #####################################################################################################
  # Applying splines in mu predictor
  if(!is.null(parameters)){
    eta <- apply_spline(vector = beta,
                        matrix = X,
                        functions = functions, 
                        parameters = parameters)
  }else{
    eta = as.vector(X%*%beta)
  }
  #####################################################################################################
  
  vartheta <- as.vector(Z%*%gamma) #Preditor linear de sigma
  
  # Funções de ligação
  if(link.mu == "logit"){
    mu <- exp(eta)/(1.0+exp(eta)) #logit
  } else if (link.mu == "probit"){
    mu <- pnorm(eta) #probit
  } else if (link.mu == "cloglog"){
    mu <- 1.0 - exp(-exp(eta)) # cloglog
  } else if (link.mu == "log"){
    mu <- exp(eta) #log
  } else if (link.mu == "loglog"){
    mu <- exp(-exp(-eta)) #loglog
  } else if (link.mu == "cauchit"){
    mu <- (pi^(-1))*atan(eta) + 0.5 #cauchit 
  } else {
    print("Funçao de ligação informada para a média mu é inválida! Verifique.")
    break
  }
  
  if(link.phi == "log"){
    phi <- exp(vartheta) #log
  } else if (link.phi == "identify") {
    phi <- vartheta #identify
  } else if (link.phi == "sqrt") {
    phi <- vartheta^2 #sqrt
  } else {
    print("Funçao de ligação informada para a precisão phi é inválida! Verifique.")
    break
  }
  
  # Amostra original
  a <- mu*phi  #Parâmetro de forma 1 da dist. beta 
  b <- (1-mu)*phi
  
  y <- rbeta(n=n, shape1=a, shape2=b) #Geração dos números aleatórios 
  
  df_n_cont <- data.frame(resposta = y, 
                          X = X,
                          Z = Z)
  
  df_cont <- df_n_cont
  
  ##############################################################################
  # Amostra para extração da contaminação


  if(scenario %in% c(2, 5) ){
    
    # Mistura das amostras para obtenção da amostra contaminada
    pc <- (contamination*n)#percentage of contamination
    # Retorna os indices dos 5 maiores valores
    ind_c <- sort(mu, index.return=TRUE, decreasing =TRUE)$ix[1:pc] 
    
  }else if(scenario %in% c(1, 4) ){
    
    # Mistura das amostras para obtenção da amostra contaminada
    pc <- (contamination*n)#percentage of contamination
    # Retorna os indices dos 5 maiores valores
    ind_c <- sort(mu, index.return=TRUE, decreasing =FALSE)$ix[1:pc] 
    
  }else if(scenario==3){# 3

    pc1 <- ((contamination/2)*n) #percentage 1 of contamination
    pc2 <- ((contamination/2)*n) #percentage 2 of contamination
    ind1_C <- sort(mu, index.return=TRUE, decreasing =TRUE)$ix[1:pc2] 
    ind2_C <- sort(mu, index.return=TRUE, decreasing =FALSE)$ix[1:pc1]
    
  }else{
    stop("\nSomente estão previstos os cenários 1 a 5! Reinforme!")
  }
  
  
  if(scenario==1){# VAR OK
    mu_pc <- (2.0 - 1.5*mu[ind_c])/2.0 
  }else if(scenario==2){# VAR OK
    mu_pc <- (2.0 - 1.5*mu[ind_c])/2.0 
  }else if(scenario==3){# VAR 

    #a_c1 <- 0.2; a_c2 <- 2.4; # simulation
    
    a_c1 <- 0.1; a_c2 <- 15.0; # Aplication

    c_c1 <- mu[ind1_C]/(1.0 - mu[ind1_C]);
    c_c2 <- mu[ind2_C]/(1.0 - mu[ind2_C]);
    mu_pc1 <- (a_c1*c_c1)/(1.0 + a_c1*c_c1) 
    mu_pc2 <- (a_c2*c_c2)/(1.0 + a_c2*c_c2) 
    
  }else if(scenario==4){# CTE OK
    mu_pc <- (2.0 - 1.7*mu[ind_c])/2.0 
  }else if(scenario==5){# CTE OK
    mu_pc <- (2.0 - 2.0*mu[ind_c])/2.0 
  }
  
  if(scenario %in% c(1,2,4,5)){
    
    yaux_c <- rbeta(pc, mu_pc*phi[ind_c], (1.0- mu_pc)*phi[ind_c])  
    y_c <- y
    y_c[ind_c] <- yaux_c 
    
    
  }else if(scenario==3){
    
    yaux1_c <- rbeta(pc1, mu_pc1*phi[ind1_C], (1.0 - mu_pc1)*phi[ind1_C]) 
    yaux2_c <- rbeta(pc2, mu_pc2*phi[ind2_C], (1.0 - mu_pc2)*phi[ind2_C]) 
    y_c <- y
    y_c[ind1_C] <- yaux1_c
    y_c[ind2_C] <- yaux2_c
    
  }
  
  df_cont$resposta <- y_c
  
  ##############################################################################
  results <- list(y=y, 
                  X=X,
                  Z=Z,
                  mu=mu, 
                  phi=phi, 
                  df_n_cont=df_n_cont,
                  y_cont=df_cont$resposta,
                  df_cont=df_cont)
  
  return(results)
}

#####################################################################################
# Simulation algorithm

#******************Algoritmo para Simulação*************************************
simulate_nlrobust <- function(N=20,
                              n=c(40, 80,160,320),
                              functions = c("exp"),
                              parameters = c(2),
                              contamination=0.05,
                              scenario,
                              run_MLE_not_contaminated,
                              run_MLE_contaminated,
                              run_SMLE_not_contaminated,
                              run_SMLE_contaminated,
                              run_LSMLE_not_contaminated,
                              run_LSMLE_contaminated,
                              link.mu="logit",
                              link.phi="log"){
  
  start_time <- Sys.time()
  
  # Atribuindo os valores reais dos parametros conforme cenário simulado
  if(scenario==1){
    
    beta = c(-1.8, 0.8)
    gamma = c(2.5, 3.0)
    
  }else if(scenario==2){
    
    beta = c(-1.7, 1.2)
    gamma = c(2.5, 3.5)
    
  }else if(scenario==3){
    
    beta = c(-1.0, 1.0)
    gamma = c(6.5)
    
  }else if(scenario==4){
    
    beta = c(-1.0, -1.4)
    gamma = c(6.0)
    
  }else if(scenario==5){
    
    beta = c(-1.7, 1.2)
    gamma = c(4.7)
    
  }
  
  message("\nINFO: Starting the sumulation process...")
  
  # Criação dos objetos para recebimento dos valores diversos
  MLE <- SMLE <- LSMLE <- list(not.contaminated=list(alg1=list(list(),list(),list(),list()),
                                                     beta.hat=list(list(list(), list()), 
                                                                   list(list(), list()),
                                                                   list(list(), list()),
                                                                   list(list(), list())),
                                                     gamma.hat=list(list(list(), list()), 
                                                                    list(list(), list()),
                                                                    list(list(), list()),
                                                                    list(list(), list())),
                                                     vies.beta.hat=list(list(numeric(), numeric()),
                                                                        list(numeric(), numeric()),
                                                                        list(numeric(), numeric()),
                                                                        list(numeric(), numeric())),
                                                     eqm.beta.hat=list(list(list(), list()), 
                                                                       list(list(), list()),
                                                                       list(list(), list()),
                                                                       list(list(), list())),
                                                     vies.gamma.hat=list(list(numeric(), numeric()),
                                                                         list(numeric(), numeric()),
                                                                         list(numeric(), numeric()),
                                                                         list(numeric(), numeric())),
                                                     eqm.gamma.hat=list(list(list(), list()), 
                                                                        list(list(), list()),
                                                                        list(list(), list()),
                                                                        list(list(), list()))),
                               contaminated=list(alg1=list(list(),list(),list(),list()),
                                                 beta.hat=list(list(list(), list()), 
                                                               list(list(), list()),
                                                               list(list(), list()),
                                                               list(list(), list())),
                                                 gamma.hat=list(list(list(), list()), 
                                                                list(list(), list()),
                                                                list(list(), list()),
                                                                list(list(), list())),
                                                 vies.beta.hat=list(list(numeric(), numeric()),
                                                                    list(numeric(), numeric()),
                                                                    list(numeric(), numeric()),
                                                                    list(numeric(), numeric())),
                                                 eqm.beta.hat=list(list(list(), list()), 
                                                                   list(list(), list()),
                                                                   list(list(), list()),
                                                                   list(list(), list())),
                                                 vies.gamma.hat=list(list(numeric(), numeric()),
                                                                     list(numeric(), numeric()),
                                                                     list(numeric(), numeric()),
                                                                     list(numeric(), numeric())),
                                                 eqm.gamma.hat=list(list(list(), list()), 
                                                                    list(list(), list()),
                                                                    list(list(), list()),
                                                                    list(list(), list())))
  )
  
  ENPV <- list(SMLE=list(not.contaminated=list(matrix(ncol = 10,
                                                      nrow = N),
                                               matrix(ncol = 10,
                                                      nrow = N),
                                               matrix(ncol = 10,
                                                      nrow = N),
                                               matrix(ncol = 10,
                                                      nrow = N)
                                               ),
                         contaminated=list(matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N)
                                          )
                         ),
               LSMLE=list(not.contaminated=list(matrix(ncol = 10,
                                                       nrow = N),
                                                matrix(ncol = 10,
                                                       nrow = N),
                                                matrix(ncol = 10,
                                                       nrow = N),
                                                matrix(ncol = 10,
                                                       nrow = N)
                         ),
                         contaminated=list(matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N),
                                           matrix(ncol = 10,
                                                  nrow = N)
                                           )
                         )
               )
  

  #ENPV$SMLE$not.contaminated[[2]][3,] <- c(1,2,3,4,5,6,7,8,9,10)


  # Dados aleatórios para utilização como covariáveis (mantidas fixas ao longo das replicas)
  x1 <- runif(min(n))
  z1 <- x1
  
  for(l in 1:length(n)){
    
    X <- matrix(c(rep(1,min(n)), x1), ncol=2, byrow=F); #regressor matrix for the mean submodel
    
    if(length(gamma)>1){
      Z <- matrix(c(rep(1,min(n)), z1), ncol=2, byrow=F); #regressor matrix for the precision submodel
    }else{
      Z <- matrix(c(rep(1,min(n))), ncol=1, byrow=F); #regressor matrix for the precision submodel
    }
    
    
    if(n[l]==80){
      X <- rbind(X,X)
      Z <- rbind(Z,Z)
    }else if(n[l]==160){
      X <- rbind(X,X,X,X)
      Z <- rbind(Z,Z,Z,Z)
    }else if(n[l]==320){
      X <- rbind(X,X,X,X,X,X,X,X)
      Z <- rbind(Z,Z,Z,Z,Z,Z,Z,Z)
    }
    
    message("\nINFO: Sample size: ",n[l])
    
    # Estrutura de repetição principal para o experimento de monte carlo
    for (j in 1:N) {
      
      message("\nINFO: Monte Carlo replica: ",j)
      
      FunSample_cont_betareg(n=n[l], 
                             X=X,
                             Z=Z,
                             beta = beta,
                             gamma = gamma,
                             functions = functions,
                             parameters = parameters,
                             contamination = contamination,
                             scenario = scenario,
                             link.mu=link.mu, 
                             link.phi=link.phi) -> dados_cont_sim
      
      len_form <- length(substr(names(dados_cont_sim$df_cont), start = 1, stop = 1))
      
      subs_form <- substr(names(dados_cont_sim$df_cont), start = 1, stop = 1)
      
      for (i in 2:(len_form)) {
        if(i==2){form <- "resposta ~ "}
        if(subs_form[i]=="X"){
          if(sum(subs_form=="X")>1){
            if(i==2){form <- paste(form, "X",".",i)}
            if( (i>2) & (i<=sum(subs_form=="X"))){form <- paste(form, "+","X",".",i)}
          }
        }else if(subs_form[i]=="Z"){
          if(sum(subs_form=="Z")>1){
            if(i== (sum(subs_form=="X")+3) ){form <- paste(form, "|","Z",".",i-1-(sum(subs_form=="X")))}
            if(i> (sum(subs_form=="X")+3) ){form <- paste(form, "+","Z",".",i-1-(sum(subs_form=="X")))}
          }else if(sum(subs_form=="Z")==1){form <- paste(form, "|","1")}
        }
        form_final <- Formula::as.Formula(gsub(" ","",form))
      }
      ##########################################################################
      # SMLE - Contaminated
      if(run_SMLE_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'SMLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.smle.cont
        SMLE$contaminated$alg1[[l]][[j]] <- c(fit.smle.cont$coefficients$mean, 
                                              fit.smle.cont$coefficients$precision,
                                              fit.smle.cont$Tuning)
        if(is.null(fit.smle.cont$enpv)){
          print('SMLE-NULL')
          ENPV$SMLE$contaminated[[l]][j,] <- c(0,0,0,0,0,0,0,0,0,0)
        }else{
          ENPV$SMLE$contaminated[[l]][j,] <- fit.smle.cont$enpv
        }
        
      }
      # SMLE - Not contaminated
      if(run_SMLE_not_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_n_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'SMLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.smle.n.cont
        SMLE$not.contaminated$alg1[[l]][[j]] <- c(fit.smle.n.cont$coefficients$mean, 
                                                  fit.smle.n.cont$coefficients$precision,
                                                  fit.smle.n.cont$Tuning)
        
        ENPV$SMLE$not.contaminated[[l]][j,] <- fit.smle.n.cont$enpv

      }
      ###############################################################################################
      ##########################################################################
      # LSMLE - Contaminated
      if(run_LSMLE_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'LSMLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.lsmle.cont
        LSMLE$contaminated$alg1[[l]][[j]] <- c(fit.lsmle.cont$coefficients$mean, 
                                               fit.lsmle.cont$coefficients$precision,
                                               fit.lsmle.cont$Tuning)
        if(is.null(fit.lsmle.cont$enpv)){
          ENPV$LSMLE$contaminated[[l]][j,] <- c(0,0,0,0,0,0,0,0,0,0)
        }else{
          ENPV$LSMLE$contaminated[[l]][j,] <- fit.lsmle.cont$enpv
        }
        

      }
      # LSMLE - Not contaminated
      if(run_LSMLE_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_n_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'LSMLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.lsmle.n.cont
        LSMLE$not.contaminated$alg1[[l]][[j]] <- c(fit.lsmle.n.cont$coefficients$mean, 
                                                   fit.lsmle.n.cont$coefficients$precision,
                                                   fit.lsmle.n.cont$Tuning)
        ENPV$LSMLE$not.contaminated[[l]][j,] <- fit.lsmle.n.cont$enpv

      }
      ###############################################################################################
      ##########################################################################
      # MLE - Contaminated
      if(run_MLE_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'MLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.mle.cont
        MLE$contaminated$alg1[[l]][[j]] <- c(fit.mle.cont$coefficients$mean, 
                                             fit.mle.cont$coefficients$precision,
                                             fit.mle.cont$Tuning)
      }
      ##########################################################################
      # MLE - Not contaminated
      if(run_MLE_not_contaminated){
        tryCatch(robustbetareg(formula = form_final,
                               data = dados_cont_sim$df_n_cont,
                               alpha = NULL,
                               functions = functions, 
                               parameters = parameters,
                               link = link.mu,
                               link.phi = link.phi,
                               type = 'MLE',
                               se.bootstrap = FALSE,
                               simulation = TRUE,
                               wald.test.bootstrap = FALSE),
                 error = function(e) {
                   return(NULL)
                 }) -> fit.mle.n.cont
        MLE$not.contaminated$alg1[[l]][[j]] <- c(fit.mle.n.cont$coefficients$mean, 
                                                 fit.mle.n.cont$coefficients$precision,
                                                 fit.mle.n.cont$Tuning)
      }
      
      
      
      
    }
    
    # Armazenando os valores separados por tamanho de amostra a cada passada
    ##############################################################################################################
    ##############################################################################
    #MLE - Not contaminated
    if(run_MLE_not_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          MLE$not.contaminated$beta.hat[[l]][[m]][k] <- MLE$not.contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          MLE$not.contaminated$gamma.hat[[l]][[m]][k] <- MLE$not.contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        MLE$not.contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(MLE$not.contaminated$beta.hat[[l]][[m]])) - beta[m])
        MLE$not.contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(MLE$not.contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        MLE$not.contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(MLE$not.contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        MLE$not.contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(MLE$not.contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    
    ##############################################################################
    #MLE - Contaminated
    if(run_MLE_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          MLE$contaminated$beta.hat[[l]][[m]][k] <- MLE$contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          MLE$contaminated$gamma.hat[[l]][[m]][k] <- MLE$contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        MLE$contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(MLE$contaminated$beta.hat[[l]][[m]])) - beta[m])
        MLE$contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(MLE$contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        MLE$contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(MLE$contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        MLE$contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(MLE$contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    ##############################################################################################################
    ##############################################################################
    #SMLE - Not contaminated
    if(run_SMLE_not_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          SMLE$not.contaminated$beta.hat[[l]][[m]][k] <- SMLE$not.contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          SMLE$not.contaminated$gamma.hat[[l]][[m]][k] <- SMLE$not.contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        SMLE$not.contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(SMLE$not.contaminated$beta.hat[[l]][[m]])) - beta[m])
        SMLE$not.contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(SMLE$not.contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        SMLE$not.contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(SMLE$not.contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        SMLE$not.contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(SMLE$not.contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    
    ##############################################################################
    #SMLE - Contaminated
    if(run_SMLE_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          SMLE$contaminated$beta.hat[[l]][[m]][k] <- SMLE$contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          SMLE$contaminated$gamma.hat[[l]][[m]][k] <- SMLE$contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        SMLE$contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(SMLE$contaminated$beta.hat[[l]][[m]])) - beta[m])
        SMLE$contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(SMLE$contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        SMLE$contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(SMLE$contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        SMLE$contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(SMLE$contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    ##############################################################################################################
    ##############################################################################
    #LSMLE - Not contaminated
    if(run_LSMLE_not_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          LSMLE$not.contaminated$beta.hat[[l]][[m]][k] <- LSMLE$not.contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          LSMLE$not.contaminated$gamma.hat[[l]][[m]][k] <- LSMLE$not.contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        LSMLE$not.contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(LSMLE$not.contaminated$beta.hat[[l]][[m]])) - beta[m])
        LSMLE$not.contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(LSMLE$not.contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        LSMLE$not.contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(LSMLE$not.contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        LSMLE$not.contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(LSMLE$not.contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    ##############################################################################
    #LSMLE - Contaminated
    if(run_LSMLE_contaminated){
      # Armazena os valores das estimativas dos coeficientes
      for(k in 1:N){
        for (m in 1:length(beta)) {
          LSMLE$contaminated$beta.hat[[l]][[m]][k] <- LSMLE$contaminated$alg1[[l]][[k]][[m]]
        }
        for (m in 1:length(gamma)) {
          LSMLE$contaminated$gamma.hat[[l]][[m]][k] <- LSMLE$contaminated$alg1[[l]][[k]][[m+length(beta)]]
        }
      }
      
      # Computando os Vieses medios relativos e erros quadráticos medios
      for (m in 1:length(beta)) {
        LSMLE$contaminated$vies.beta.hat[[l]][m] <- (mean(unlist(LSMLE$contaminated$beta.hat[[l]][[m]])) - beta[m])
        LSMLE$contaminated$eqm.beta.hat[[l]][[m]] <- mean( ( (unlist(LSMLE$contaminated$beta.hat[[l]][[m]]) - beta[m]) )^2 )
      }
      for (m in 1:length(gamma)) {
        LSMLE$contaminated$vies.gamma.hat[[l]][m] <- (mean(unlist(LSMLE$contaminated$gamma.hat[[l]][[m]])) - gamma[m])
        LSMLE$contaminated$eqm.gamma.hat[[l]][[m]] <- mean( ( (unlist(LSMLE$contaminated$gamma.hat[[l]][[m]]) - gamma[m]) )^2 )
      }
    }
    
  }
  
  ##############################################################################
  # Geração das tabelas para utilização nos gráficos
  
  ###################################################
  #MLE

  tryCatch(data.frame("40" = unlist(MLE$not.contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(MLE$not.contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(MLE$not.contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(MLE$not.contaminated$beta.hat[[4]][[1]]),
                      "40" = unlist(MLE$contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(MLE$contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(MLE$contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(MLE$contaminated$beta.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.mle.beta1
  

  tryCatch(data.frame("40" = unlist(MLE$not.contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(MLE$not.contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(MLE$not.contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(MLE$not.contaminated$beta.hat[[4]][[2]]),
                      "40" = unlist(MLE$contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(MLE$contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(MLE$contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(MLE$contaminated$beta.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.mle.beta2
  
  tryCatch(data.frame("40" = unlist(MLE$not.contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(MLE$not.contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(MLE$not.contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(MLE$not.contaminated$gamma.hat[[4]][[1]]),
                      "40" = unlist(MLE$contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(MLE$contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(MLE$contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(MLE$contaminated$gamma.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.mle.gamma1
  
  tryCatch(data.frame("40" = unlist(MLE$not.contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(MLE$not.contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(MLE$not.contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(MLE$not.contaminated$gamma.hat[[4]][[2]]),
                      "40" = unlist(MLE$contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(MLE$contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(MLE$contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(MLE$contaminated$gamma.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.mle.gamma2
  
  ###################################################
  # SMLE
  tryCatch(data.frame("40" = unlist(SMLE$not.contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(SMLE$not.contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(SMLE$not.contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(SMLE$not.contaminated$beta.hat[[4]][[1]]),
                      "40" = unlist(SMLE$contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(SMLE$contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(SMLE$contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(SMLE$contaminated$beta.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.smle.beta1
  
  tryCatch(data.frame("40" = unlist(SMLE$not.contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(SMLE$not.contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(SMLE$not.contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(SMLE$not.contaminated$beta.hat[[4]][[2]]),
                      "40" = unlist(SMLE$contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(SMLE$contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(SMLE$contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(SMLE$contaminated$beta.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.smle.beta2
  
  tryCatch(data.frame("40" = unlist(SMLE$not.contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(SMLE$not.contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(SMLE$not.contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(SMLE$not.contaminated$gamma.hat[[4]][[1]]),
                      "40" = unlist(SMLE$contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(SMLE$contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(SMLE$contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(SMLE$contaminated$gamma.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.smle.gamma1
  
  tryCatch(data.frame("40" = unlist(SMLE$not.contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(SMLE$not.contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(SMLE$not.contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(SMLE$not.contaminated$gamma.hat[[4]][[2]]),
                      "40" = unlist(SMLE$contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(SMLE$contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(SMLE$contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(SMLE$contaminated$gamma.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.smle.gamma2
  
  ###################################################
  # LSMLE
  tryCatch(data.frame("40" = unlist(LSMLE$not.contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(LSMLE$not.contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(LSMLE$not.contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(LSMLE$not.contaminated$beta.hat[[4]][[1]]),
                      "40" = unlist(LSMLE$contaminated$beta.hat[[1]][[1]]),
                      "80" = unlist(LSMLE$contaminated$beta.hat[[2]][[1]]),
                      "160" = unlist(LSMLE$contaminated$beta.hat[[3]][[1]]),
                      "320" = unlist(LSMLE$contaminated$beta.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.beta1
  
  tryCatch(data.frame("40" = unlist(LSMLE$not.contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(LSMLE$not.contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(LSMLE$not.contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(LSMLE$not.contaminated$beta.hat[[4]][[2]]),
                      "40" = unlist(LSMLE$contaminated$beta.hat[[1]][[2]]),
                      "80" = unlist(LSMLE$contaminated$beta.hat[[2]][[2]]),
                      "160" = unlist(LSMLE$contaminated$beta.hat[[3]][[2]]),
                      "320" = unlist(LSMLE$contaminated$beta.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.beta2
  
  tryCatch(data.frame("40" = unlist(LSMLE$not.contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(LSMLE$not.contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(LSMLE$not.contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(LSMLE$not.contaminated$gamma.hat[[4]][[1]]),
                      "40" = unlist(LSMLE$contaminated$gamma.hat[[1]][[1]]),
                      "80" = unlist(LSMLE$contaminated$gamma.hat[[2]][[1]]),
                      "160" = unlist(LSMLE$contaminated$gamma.hat[[3]][[1]]),
                      "320" = unlist(LSMLE$contaminated$gamma.hat[[4]][[1]])),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.gamma1
  
  tryCatch(data.frame("40" = unlist(LSMLE$not.contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(LSMLE$not.contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(LSMLE$not.contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(LSMLE$not.contaminated$gamma.hat[[4]][[2]]),
                      "40" = unlist(LSMLE$contaminated$gamma.hat[[1]][[2]]),
                      "80" = unlist(LSMLE$contaminated$gamma.hat[[2]][[2]]),
                      "160" = unlist(LSMLE$contaminated$gamma.hat[[3]][[2]]),
                      "320" = unlist(LSMLE$contaminated$gamma.hat[[4]][[2]])),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.gamma2
  
  ##############################################################################
  
  end_time <- Sys.time()
  
  total.processing.time <- (end_time - start_time)
  
  message("\nINFO: Process finished! Total time: ", total.processing.time)
  
  return(list(coefficients=list(beta=beta,
                                gamma=gamma),
              data=list(MLE=list(beta1=data.mle.beta1,
                                 beta2=data.mle.beta2,
                                 gamma1=data.mle.gamma1,
                                 gamma2=data.mle.gamma2),
                        SMLE=list(beta1=data.smle.beta1,
                                 beta2=data.smle.beta2,
                                 gamma1=data.smle.gamma1,
                                 gamma2=data.smle.gamma2),
                        LSMLE=list(beta1=data.lsmle.beta1,
                                   beta2=data.lsmle.beta2,
                                   gamma1=data.lsmle.gamma1,
                                   gamma2=data.lsmle.gamma2)),
              MLE=MLE,
              SMLE=SMLE,
              LSMLE=LSMLE,
              ENPV=ENPV,
              total.processing.time=total.processing.time
  ))
  
}

#####################################################################################
# Funçaõ para obtenção das rzões entre os TMSEs

extract_measures <- function(object,
                             sample_size){
  
  if(sample_size==40){
    index <- 1
  }else if(sample_size==80){
    index <- 2
  }else if(sample_size==160){
    index <- 3
  }else if(sample_size==320){
    index <- 4
  }else{
    stop("\nTamanho amostral informado não está disponível!.\nInforme 40, 80, 160 ou 320!")
  }
  
  # Vieses
  sum(unlist(c(object$MLE$contaminated$vies.beta.hat[[index]],
               object$MLE$contaminated$vies.gamma.hat[[index]]))) -> vies.cont.tmse.mle
  
  sum(unlist(c(object$SMLE$contaminated$vies.beta.hat[[index]],
               object$SMLE$contaminated$vies.gamma.hat[[index]]))) -> vies.cont.tmse.smle
  
  sum(unlist(c(object$LSMLE$contaminated$vies.beta.hat[[index]],
               object$LSMLE$contaminated$vies.gamma.hat[[index]]))) -> vies.cont.tmse.lsmle
  
  sum(unlist(c(object$MLE$not.contaminated$vies.beta.hat[[index]],
               object$MLE$not.contaminated$vies.gamma.hat[[index]]))) -> vies.n.cont.tmse.mle
  
  sum(unlist(c(object$SMLE$not.contaminated$vies.beta.hat[[index]],
               object$SMLE$not.contaminated$vies.gamma.hat[[index]]))) -> vies.n.cont.tmse.smle
  
  sum(unlist(c(object$LSMLE$not.contaminated$vies.beta.hat[[index]],
               object$LSMLE$not.contaminated$vies.gamma.hat[[index]]))) -> vies.n.cont.tmse.lsmle

  
  # eqms
  sum(unlist(c(object$MLE$contaminated$eqm.beta.hat[[index]],
               object$MLE$contaminated$eqm.gamma.hat[[index]]))) -> eqm.cont.tmse.mle
  
  sum(unlist(c(object$SMLE$contaminated$eqm.beta.hat[[index]],
               object$SMLE$contaminated$eqm.gamma.hat[[index]]))) -> eqm.cont.tmse.smle
  
  sum(unlist(c(object$LSMLE$contaminated$eqm.beta.hat[[index]],
               object$LSMLE$contaminated$eqm.gamma.hat[[index]]))) -> eqm.cont.tmse.lsmle
  
  sum(unlist(c(object$MLE$not.contaminated$eqm.beta.hat[[index]],
               object$MLE$not.contaminated$eqm.gamma.hat[[index]]))) -> eqm.n.cont.tmse.mle
  
  sum(unlist(c(object$SMLE$not.contaminated$eqm.beta.hat[[index]],
               object$SMLE$not.contaminated$eqm.gamma.hat[[index]]))) -> eqm.n.cont.tmse.smle
  
  sum(unlist(c(object$LSMLE$not.contaminated$eqm.beta.hat[[index]],
               object$LSMLE$not.contaminated$eqm.gamma.hat[[index]]))) -> eqm.n.cont.tmse.lsmle

  
  
  return(list(vies=list(contaminated=list("mle/smle" = vies.cont.tmse.mle/vies.cont.tmse.smle,
                                          "mle/lsmle" = vies.cont.tmse.mle/vies.cont.tmse.lsmle,
                                          "smle/lsmle" = vies.cont.tmse.smle/vies.cont.tmse.lsmle),
                        not.contaminated=list("mle/smle" = vies.n.cont.tmse.mle/vies.n.cont.tmse.smle,
                                              "mle/lsmle" = vies.n.cont.tmse.mle/vies.n.cont.tmse.lsmle,
                                              "smle/lsmle" = vies.n.cont.tmse.smle/vies.n.cont.tmse.lsmle)
                        ),
              eqm=list(contaminated=list("mle/smle" = eqm.cont.tmse.mle/eqm.cont.tmse.smle,
                                          "mle/lsmle" = eqm.cont.tmse.mle/eqm.cont.tmse.lsmle,
                                          "smle/lsmle" = eqm.cont.tmse.smle/eqm.cont.tmse.lsmle),
                        not.contaminated=list("mle/smle" = eqm.n.cont.tmse.mle/eqm.n.cont.tmse.smle,
                                              "mle/lsmle" = eqm.n.cont.tmse.mle/eqm.n.cont.tmse.lsmle,
                                              "smle/lsmle" = eqm.n.cont.tmse.smle/eqm.n.cont.tmse.lsmle)
              )
            )
  )

}

# Extrai as contantes e monta os dataframes
extract_tuning <- function(object,
                           return_q = FALSE){
  
  tuning.vec <- list(SMLE=list(contaminated=list(numeric(),
                                                 numeric(),
                                                 numeric(),
                                                 numeric()
  ),
  not.contaminated=list(numeric(),
                        numeric(),
                        numeric(),
                        numeric()
  )
  ),
  LSMLE=list(contaminated=list(numeric(),
                               numeric(),
                               numeric(),
                               numeric()
  ),
  not.contaminated=list(numeric(),
                        numeric(),
                        numeric(),
                        numeric()
  )
  )
  )
  
  for (j in 1:length(object$SMLE$contaminated$alg1)) {
    for(i in 1:length(object$SMLE$contaminated$alg1[[j]])){
      # Obtendo a posição da contante no vetor
      position_alpha <- length(object$SMLE$contaminated$alg1[[j]][[i]])
      # SMLE
      tryCatch(tuning.vec$SMLE$contaminated[[j]][i] <- object$SMLE$contaminated$alg1[[j]][[i]][[position_alpha]],
               error = function(e) {
                 return(NULL)
               })
      tryCatch(tuning.vec$SMLE$not.contaminated[[j]][i] <- object$SMLE$not.contaminated$alg1[[j]][[i]][[position_alpha]],
               error = function(e) {
                 return(NULL)
               })
      
      # LSMLE
      tryCatch(tuning.vec$LSMLE$contaminated[[j]][i] <- object$LSMLE$contaminated$alg1[[j]][[i]][[position_alpha]],
               error = function(e) {
                 return(NULL)
               })
      tryCatch(tuning.vec$LSMLE$not.contaminated[[j]][i] <- object$LSMLE$not.contaminated$alg1[[j]][[i]][[position_alpha]],
               error = function(e) {
                 return(NULL)
               })
    }
    if(return_q){
      # SMLE
      tuning.vec$SMLE$contaminated[[j]] <- 1-tuning.vec$SMLE$contaminated[[j]]
      tuning.vec$SMLE$not.contaminated[[j]] <- 1-tuning.vec$SMLE$not.contaminated[[j]]
      # LSMLE
      tuning.vec$LSMLE$contaminated[[j]] <- 1-tuning.vec$LSMLE$contaminated[[j]]
      tuning.vec$LSMLE$not.contaminated[[j]] <- 1-tuning.vec$LSMLE$not.contaminated[[j]]
    }
  }
  # Organiza todos o dados em um dtaframe por estimador
  # SMLE
  data.frame("40" = tuning.vec$SMLE$not.contaminated[[1]],
             "80" = tuning.vec$SMLE$not.contaminated[[2]],
             "160" = tuning.vec$SMLE$not.contaminated[[3]],
             "320" = tuning.vec$SMLE$not.contaminated[[4]],
             "40" = tuning.vec$SMLE$contaminated[[1]],
             "80" = tuning.vec$SMLE$contaminated[[2]],
             "160" = tuning.vec$SMLE$contaminated[[3]],
             "320" = tuning.vec$SMLE$contaminated[[4]]) -> data.smle.tuning
  
  # LSMLE
  data.frame("40" = tuning.vec$LSMLE$not.contaminated[[1]],
             "80" = tuning.vec$LSMLE$not.contaminated[[2]],
             "160" = tuning.vec$LSMLE$not.contaminated[[3]],
             "320" = tuning.vec$LSMLE$not.contaminated[[4]],
             "40" = tuning.vec$LSMLE$contaminated[[1]],
             "80" = tuning.vec$LSMLE$contaminated[[2]],
             "160" = tuning.vec$LSMLE$contaminated[[3]],
             "320" = tuning.vec$LSMLE$contaminated[[4]]) -> data.lsmle.tuning
  
  return(list(tuning.vec=tuning.vec,
              data=list(SMLE=data.smle.tuning,
                        LSMLE=data.lsmle.tuning))
  )
}


generate_graphs <- function(object,
                            max_len=100){
  
  ##############################################################################
  # Geração das tabelas para utilização nos gráficos
  
  ###################################################
  #MLE
  
  tryCatch(data.frame("40" = unlist(object$MLE$not.contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$MLE$not.contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$MLE$not.contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$MLE$not.contaminated$beta.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$MLE$contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$MLE$contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$MLE$contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$MLE$contaminated$beta.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.mle.beta1
  
  
  tryCatch(data.frame("40" = unlist(object$MLE$not.contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$MLE$not.contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$MLE$not.contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$MLE$not.contaminated$beta.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$MLE$contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$MLE$contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$MLE$contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$MLE$contaminated$beta.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.mle.beta2
  
  tryCatch(data.frame("40" = unlist(object$MLE$not.contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$MLE$not.contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$MLE$not.contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$MLE$not.contaminated$gamma.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$MLE$contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$MLE$contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$MLE$contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$MLE$contaminated$gamma.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.mle.gamma1
  
  tryCatch(data.frame("40" = unlist(object$MLE$not.contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$MLE$not.contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$MLE$not.contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$MLE$not.contaminated$gamma.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$MLE$contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$MLE$contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$MLE$contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$MLE$contaminated$gamma.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.mle.gamma2
  
  ###################################################
  # SMLE
  tryCatch(data.frame("40" = unlist(object$SMLE$not.contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$SMLE$not.contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$SMLE$not.contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$SMLE$not.contaminated$beta.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$SMLE$contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$SMLE$contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$SMLE$contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$SMLE$contaminated$beta.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.smle.beta1
  
  tryCatch(data.frame("40" = unlist(object$SMLE$not.contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$SMLE$not.contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$SMLE$not.contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$SMLE$not.contaminated$beta.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$SMLE$contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$SMLE$contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$SMLE$contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$SMLE$contaminated$beta.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.smle.beta2
  
  tryCatch(data.frame("40" = unlist(object$SMLE$not.contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$SMLE$not.contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$SMLE$not.contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$SMLE$not.contaminated$gamma.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$SMLE$contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$SMLE$contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$SMLE$contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$SMLE$contaminated$gamma.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.smle.gamma1
  
  tryCatch(data.frame("40" = unlist(object$SMLE$not.contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$SMLE$not.contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$SMLE$not.contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$SMLE$not.contaminated$gamma.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$SMLE$contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$SMLE$contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$SMLE$contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$SMLE$contaminated$gamma.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.smle.gamma2
  
  ###################################################
  # LSMLE
  tryCatch(data.frame("40" = unlist(object$LSMLE$not.contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$LSMLE$not.contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$LSMLE$not.contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$LSMLE$not.contaminated$beta.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$LSMLE$contaminated$beta.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$LSMLE$contaminated$beta.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$LSMLE$contaminated$beta.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$LSMLE$contaminated$beta.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.beta1
  
  tryCatch(data.frame("40" = unlist(object$LSMLE$not.contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$LSMLE$not.contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$LSMLE$not.contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$LSMLE$not.contaminated$beta.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$LSMLE$contaminated$beta.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$LSMLE$contaminated$beta.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$LSMLE$contaminated$beta.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$LSMLE$contaminated$beta.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.beta2
  
  tryCatch(data.frame("40" = unlist(object$LSMLE$not.contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$LSMLE$not.contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$LSMLE$not.contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$LSMLE$not.contaminated$gamma.hat[[4]][[1]])[1:max_len],
                      "40" = unlist(object$LSMLE$contaminated$gamma.hat[[1]][[1]])[1:max_len],
                      "80" = unlist(object$LSMLE$contaminated$gamma.hat[[2]][[1]])[1:max_len],
                      "160" = unlist(object$LSMLE$contaminated$gamma.hat[[3]][[1]])[1:max_len],
                      "320" = unlist(object$LSMLE$contaminated$gamma.hat[[4]][[1]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.gamma1
  
  tryCatch(data.frame("40" = unlist(object$LSMLE$not.contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$LSMLE$not.contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$LSMLE$not.contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$LSMLE$not.contaminated$gamma.hat[[4]][[2]])[1:max_len],
                      "40" = unlist(object$LSMLE$contaminated$gamma.hat[[1]][[2]])[1:max_len],
                      "80" = unlist(object$LSMLE$contaminated$gamma.hat[[2]][[2]])[1:max_len],
                      "160" = unlist(object$LSMLE$contaminated$gamma.hat[[3]][[2]])[1:max_len],
                      "320" = unlist(object$LSMLE$contaminated$gamma.hat[[4]][[2]])[1:max_len]),
           error = function(e) {
             return(NULL)
           }) -> data.lsmle.gamma2
  
  return(data=list(MLE=list(beta1=data.mle.beta1,
                            beta2=data.mle.beta2,
                            gamma1=data.mle.gamma1,
                            gamma2=data.mle.gamma2),
                   SMLE=list(beta1=data.smle.beta1,
                             beta2=data.smle.beta2,
                             gamma1=data.smle.gamma1,
                             gamma2=data.smle.gamma2),
                   LSMLE=list(beta1=data.lsmle.beta1,
                              beta2=data.lsmle.beta2,
                              gamma1=data.lsmle.gamma1,
                              gamma2=data.lsmle.gamma2)))
  
}

######################################################################################################
# Quantile residuals

# mu_hat, phi_hat, y = y, 
# X = x, linkobj = linkobj

nlrobustbetareg_residuals <- function (mu, 
                                       phi,
                                       y,
                                       wts,
                                       type = c("quantile")) 
{
  
  if (!(type %in% c("quantile"))) 
    stop(sprintf("only type 'quantile' is currently implemented!"))
  
  if (is.null(wts))
    wts <- 1
  
  pit <- pbetar(y, 
                mu = mu,
                phi = phi,
                log.p = FALSE)
  
  pit[pit==1] <- 0.9999999999999999
  
  #res <- sqrt(wts)*qnorm(pit)
  res <- qnorm(pit)
  
  return(res)
}

###############################################################################################

nlrobust_envelope <- function(object,
                              faixa.fixed = NULL,
                              title_comp = "") { 
  
    # Creating objects

  y <- object$y
  X <- object$model$mean
  Z <- object$model$precision
  theta <- c(object$coefficients$mean,object$coefficients$precision)
  linkmu <- object$link
  linkphi <- object$link.phi
  type <- as.character(object$method)
  formula <- object$formula
  data <- object$data
  alpha <- object$Tuning
  weights <- object$weigths
  functions <- object$nonlinear.args$functions
  parameters <- object$nonlinear.args$parameters

  #source("SMLE.r") 
  B <- 100; #number of replicates
  kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X)
  #***parameters for parametric bootstrap***#
  beta_p <- theta[1:kk1]
  gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	
  
  if(!is.null(parameters)){
    etahat <- apply_spline(vector = beta_p,
                           matrix = X,
                           functions = functions, 
                           parameters = parameters)
  }else{
    etahat <- X%*%beta_p
  }

  deltahat <- Z%*%gama_p 	
  
  #***************link functions for mean submodel**********************#
  if(linkmu == "logit") muhat <- exp(etahat)/(1.0+exp(etahat))
  if(linkmu == "probit") muhat <- pnorm(etahat) 
  if(linkmu == "cloglog") muhat <- 1.0 - exp(-exp(etahat)) 
  if(linkmu == "log") muhat <- exp(etahat) 
  if(linkmu == "loglog") muhat <- exp(-exp(etahat)) 
  if(linkmu == "cauchit") muhat <- (pi^(-1))*atan(etahat) + 0.5 
  #************************************************************#
  #***************link functions for precision submodel**********************#
  if(linkphi == "log") phihat <- exp(deltahat) 
  if(linkphi == "identify") phihat <- deltahat 
  if(linkphi == "sqrt") phihat <- deltahat^2
  #************************************************************#
  
  Menvelope_rp2 <- matrix(numeric(0),nrow=n,ncol=B)
  
  #------------> residuals for the observed sample<--------------#
  RP2 <- object$residuals
  #set.seed(c(1994,1991), kind="Marsaglia-Multicarry")
  for(j in 1:B){		
    # ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
    
    ygen <- rbetar(n, 
                   muhat,
                   phihat)
    while(any(round(ygen,5)==0|round(ygen,5)==1)){
      ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
    }
    starting.points(formula = formula,
                    data = data,
                    y=ygen,
                    x=X,
                    z=Z,
                    type = "MLE",
                    parameters = parameters,
                    functions = functions,
                    link = linkmu,
                    link.phi = linkphi) -> start_theta_mle
    
    starting.points(formula = formula,
                    data = data,
                    y=ygen,
                    x=X,
                    z=Z,
                    type = "SMLE",
                    parameters = parameters,
                    functions = functions,
                    link = linkmu,
                    link.phi = linkphi) -> start_theta_smle_lsmle
    
    if(type=="SMLE"){
      fit = tryCatch(SMLE.fit(y=ygen, 
                              x=X,
                              z=Z,
                              alpha = alpha, 
                              link = linkmu, 
                              link.phi = linkphi,
                              type=type,
                              se.bootstrap = FALSE,
                              wald.test.bootstrap = FALSE, 
                              functions = functions, 
                              parameters = parameters,
                              control = robustbetareg.control(start = start_theta_smle_lsmle)), 
                     error = function(e) {
                       return(NULL)
                     })
      mu <- fit$fitted.values$mu.predict
      if(length(fit$fitted.values$phi.predict)==1){
        phi <- rep(fit$fitted.values$phi.predict, length(mu))
      }else{
        phi = fit$fitted.values$phi.predict
      }
      
      RP2_b <- nlrobustbetareg_residuals(mu = mu,
                                         phi = phi,
                                         y = ygen,
                                         wts = fit$weights,
                                         type = c("quantile"))
    }else if(type=="LSMLE"){
      fit = tryCatch(LSMLE.fit(y=ygen, 
                               x=X,
                               z=Z,
                               alpha = alpha, 
                               link = linkmu, 
                               link.phi = linkphi,
                               type=type,
                               se.bootstrap = FALSE,
                               wald.test.bootstrap = FALSE, 
                               functions = functions, 
                               parameters = parameters,
                               control = robustbetareg.control(start = start_theta_smle_lsmle)), 
                     error = function(e) {
                       return(NULL)
                     })
      mu <- fit$fitted.values$mu.predict
      if(length(fit$fitted.values$phi.predict)==1){
        phi <- rep(fit$fitted.values$phi.predict, length(mu))
      }else{
        phi = fit$fitted.values$phi.predict
      }
      
      RP2_b <- nlrobustbetareg_residuals(mu = mu,
                                         phi = phi,
                                         y = ygen,
                                         wts = fit$weights,
                                         type = c("quantile"))
    }else if(type=="MLE"){
      
      fit = tryCatch(MLE.fit(y=ygen, 
                             x=X,
                             z=Z,
                             link = linkmu, 
                             link.phi = linkphi,
                             type=type,
                             se.bootstrap = FALSE,
                             wald.test.bootstrap = FALSE, 
                             functions = functions, 
                             parameters = parameters,
                             data = data,
                             formula = formula,
                             control = robustbetareg.control(start = start_theta_mle)), 
                     error = function(e) {
                       return(NULL)
                     })
      mu <- fit$fitted.values$mu.predict
      if(length(fit$fitted.values$phi.predict)==1){
        phi <- rep(fit$fitted.values$phi.predict, length(mu))
      }else{
        phi = fit$fitted.values$phi.predict
      }
      
      RP2_b <- nlrobustbetareg_residuals(mu = mu,
                                         phi = phi,
                                         y = ygen,
                                         wts = fit$weights,
                                         type = c("quantile"))
    }
    
    Menvelope_rp2[,j] = RP2_b
  }
  # browser()
  main.title = paste("Envelope para",type,"-", title_comp)
  faixa.fixed = faixa.fixed
  labels.fixed = NULL
  
  # browser()

  Menvelope_rp2 <- apply(Menvelope_rp2,2,sort)
  res_rp2 <- RP2
  res_min_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.025)));         
  res_mean_rp2 <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.5)));                              
  res_max_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.975)));           
  faixa <- range(res_rp2,res_min_rp2,res_max_rp2)
  if(is.vector(faixa.fixed)) faixa <- faixa.fixed
  if(is.vector(labels.fixed)) labels <- labels.fixed
  # par(mar=c(5.0,
  #           5.0,
  #           4.0,
  #           2.0))
  v <- qqnorm(res_rp2, 
              main=main.title, 
              xlab="Quantis normais",
              ylab="Residuos",
              ylim=faixa,
              pch=20, 
              cex=2.5,
              cex.lab=2.0, 
              cex.axis=2.0, 
              cex.main=2.5)
  # identify(v$x,
  #          v$y,
  #          labels,
  #          cex =1.3) #identify points in the plot
  identify(v$x, v$y, cex=1.3)
  par(new=T)
  #
  qqnorm(res_min_rp2,axes=F,main = "",xlab="",ylab="",type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_max_rp2,axes=F,main = "",xlab="",ylab="", type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_mean_rp2,axes=F,xlab="",main = "", ylab="", type="l",ylim=faixa,lty=2,lwd=2.0)
}#ends function


#######################

nlrobust_envelope_2 <- function(y, 
                              X, 
                              Z,
                              theta, 
                              linkmu="logit",
                              linkphi="log", 
                              SMLE=T, 
                              main.title = "Envelope", 
                              faixa.fixed = NULL, 
                              labels.fixed = NULL) { 
  #source("SMLE.r") 
  B <- 100; #number of replicates
  kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X)
  #***parameters for parametric bootstrap***#
  beta_p <- theta[1:kk1]
  gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
  etahat <- X%*%beta_p
  deltahat <- Z%*%gama_p 	
  
  #***************link functions for mean submodel**********************#
  if(linkmu == "logit") muhat <- exp(etahat)/(1.0+exp(etahat))
  if(linkmu == "probit") muhat <- pnorm(etahat) 
  if(linkmu == "cloglog") muhat <- 1.0 - exp(-exp(etahat)) 
  if(linkmu == "log") muhat <- exp(etahat) 
  if(linkmu == "loglog") muhat <- exp(-exp(etahat)) 
  if(linkmu == "cauchit") muhat <- (pi^(-1))*atan(etahat) + 0.5 
  #************************************************************#
  #***************link functions for precision submodel**********************#
  if(linkphi == "log") phihat <- exp(deltahat) 
  if(linkphi == "identify") phihat <- deltahat 
  if(linkphi == "sqrt") phihat <- deltahat^2
  #************************************************************#
  
  Menvelope_rp2 <- matrix(numeric(0),nrow=n,ncol=B)
  
  #------------> residuals for the observed sample<--------------#
  RP2 <- residuals_beta(y, X, Z, c(beta_p,gama_p), linkmu= linkmu, linkphi= linkphi)
  set.seed(c(1994,1991), kind="Marsaglia-Multicarry")
  
  for(j in 1:B){		
    ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
    while(any(round(ygen,5)==0|round(ygen,5)==1)){
      ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
    }
    if(SMLE==T){
      fit <- SMLE_BETA(y=ygen, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS", startV="CP", linkmu=linkmu, linkphi=linkphi)
    }
    else{
      fit <- SMLE_BETA(y=ygen, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu=linkmu, linkphi=linkphi)
    }
    
    RP2_b <- residuals_beta(ygen,X,Z,c(fit$beta,fit$gama), linkmu=linkmu, linkphi=linkphi)
    Menvelope_rp2[,j] = RP2_b
  }
  Menvelope_rp2 <- apply(Menvelope_rp2,2,sort);          
  res_rp2 <-    RP2;    
  res_min_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.05)));         
  res_mean_rp2 <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.5)));                              
  res_max_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.95)));           
  faixa <- range(res_rp2,res_min_rp2,res_max_rp2)
  if(is.vector(faixa.fixed)) faixa <- faixa.fixed
  if(is.vector(labels.fixed)) labels <- labels.fixed
  par(mar=c(5.0,5.0,4.0,2.0))
  v <- qqnorm(res_rp2, main=main.title, xlab="Normal quantiles", ylab="Residuals", ylim=faixa, pch=16, cex=1.5, cex.lab=2.0, cex.axis=1.5, cex.main=2.0)
  identify(v$x,v$y,labels,cex =1.3) #identify points in the plot
  #identify(v$x[c(15,16,72)],v$y[c(15,16,72)],cex=1.3,labels=c("15","16","72"), cex=1.3) #Only for the the firm cost data
  par(new=T)
  #
  qqnorm(res_min_rp2,axes=F,main = "",xlab="",ylab="",type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_max_rp2,axes=F,main = "",xlab="",ylab="", type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_mean_rp2,axes=F,xlab="",main = "", ylab="", type="l",ylim=faixa,lty=2,lwd=2.0)
}#ends function

################################################################################
# Generate simulated data

generate_simulated_data <- function(n,
                                    beta,
                                    gamma,
                                    functions,
                                    parameters,
                                    contamination,
                                    scenario,
                                    link.mu,
                                    link.phi){
  
  # Dados aleatórios para utilização como covariáveis
  for(i in 1:max(length(beta), length(gamma))) {
    a = runif(n)
    if(i==1){cov_vec = rep(1,n)}
    else{ cov_vec = append(cov_vec,a) }
  }
  
  X <- matrix(cov_vec[1:(n*length(beta))], ncol=length(beta), byrow=F)
  Z <- matrix(cov_vec[1:(n*length(gamma))], ncol=length(gamma), byrow=F)
  
  # FunSample_betareg(n=n,
  #                   X=X,
  #                   Z=Z,
  #                   functions = functions,
  #                   parameters = parameters,
  #                   beta = beta,
  #                   gamma = gamma,
  #                   link.mu=link.mu,
  #                   link.phi=link.phi) -> sim_dados
  
  FunSample_cont_betareg(n=n, 
                         X=X, 
                         Z=Z,
                         beta=beta,
                         gamma=gamma,
                         functions=functions,
                         parameters = parameters,
                         contamination = contamination,
                         scenario = scenario,
                         link.mu = link.mu, 
                         link.phi = link.phi) -> sim_dados

  sim_df <- sim_dados$df_n_cont
  
  attach(sim_df)
  
  len_form <- length(substr(names(sim_df), start = 1, stop = 1))
  
  subs_form <- substr(names(sim_df), start = 1, stop = 1)
  
  for (i in 2:(len_form)) {
    if(i==2){form <- "resposta ~ "}
    if(subs_form[i]=="X"){
      if(sum(subs_form=="X")>1){
        if(i==2){form <- paste(form, "X",".",i)}
        if( (i>2) & (i<=sum(subs_form=="X"))){form <- paste(form, "+","X",".",i)}
      }
    }else if(subs_form[i]=="Z"){
      if(sum(subs_form=="Z")>1){
        if(i== (sum(subs_form=="X")+3) ){form <- paste(form, "|","Z",".",i-1-(sum(subs_form=="X")))}
        if(i> (sum(subs_form=="X")+3) ){form <- paste(form, "+","Z",".",i-1-(sum(subs_form=="X")))}
      }else if(sum(subs_form=="Z")==1){form <- paste(form, "|","1")}
    }
    form_final <- Formula::as.Formula(gsub(" ","",form))
  }
  
  return(list(formula = form_final,
              data = list(n.cont=list(data = sim_dados$df_n_cont,
                                      resposta = sim_dados$y, 
                                      X = sim_dados$X, 
                                      Z = sim_dados$Z),
                          cont=list(data = sim_dados$df_cont,
                                      resposta = sim_dados$df_cont$resposta, 
                                      X = sim_dados$df_cont$X, 
                                      Z = sim_dados$df_cont$Z)),
              functions = functions,
              parameters = parameters,
              link.mu = link.mu,
              link.phi = link.phi))
  
}

################################################################################
# Run the aplication with dimulated data exériment

run_simdata_experiment <- function(sample_size = c(40,80,160,320)){
  for(i in sample_size){
    
    closeAllConnections()
    
    start_time <- Sys.time()

    set.seed(965)
    generate_simulated_data(n = i,
                            beta = c(-0.6, 0.8),
                            gamma = c(3.9),
                            functions = "power_par",
                            parameters = 2, 
                            contamination = 0.05,
                            scenario = 3,# Decresing=FALSE
                            link.mu = "logit",
                            link.phi = "log") -> sim_data
    
    plot(sim_data$data$cont$data$X.2,
         sim_data$data$cont$data$resposta,
         xlab="Variável explicativa", 
         ylab="Resposta", 
         pch=21,
         ylim = c(-0.1,1.1))
    
    robustbetareg(formula = sim_data$formula,
                  data = sim_data$data$n.cont$data,
                  alpha = NULL,
                  functions = sim_data$functions, 
                  parameters = sim_data$parameters,
                  link = sim_data$link.mu,
                  link.phi = sim_data$link.phi,
                  type = "MLE",
                  se.bootstrap = FALSE,
                  wald.test.bootstrap = TRUE
    ) -> sim_exp_mle_n_cont
    
    
    robustbetareg(formula = sim_data$formula,
                  data = sim_data$data$cont$data,
                  alpha = NULL,
                  functions = sim_data$functions, 
                  parameters = sim_data$parameters,
                  link = sim_data$link.mu,
                  link.phi = sim_data$link.phi,
                  type = "MLE",
                  se.bootstrap = FALSE,
                  wald.test.bootstrap = TRUE
    ) -> sim_exp_mle
    
    robustbetareg(formula = sim_data$formula,
                  data = sim_data$data$cont$data,
                  alpha = NULL,
                  functions = sim_data$functions, 
                  parameters = sim_data$parameters,
                  link = sim_data$link.mu,
                  link.phi = sim_data$link.phi,
                  type = "SMLE",
                  se.bootstrap = FALSE,
                  wald.test.bootstrap = TRUE
    ) -> sim_exp_smle
    
    robustbetareg(formula = sim_data$formula,
                  data = sim_data$data$cont$data,
                  alpha = NULL,
                  functions = sim_data$functions, 
                  parameters = sim_data$parameters,
                  link = sim_data$link.mu,
                  link.phi = sim_data$link.phi,
                  type = "LSMLE",
                  se.bootstrap = FALSE,
                  wald.test.bootstrap = TRUE
    ) -> sim_exp_lsmle
    
    sim_exp_mle_n_cont |> summary()
    sim_exp_mle |> summary()
    sim_exp_smle |> summary()
    sim_exp_lsmle |> summary()
    
    dir_result_mle_n_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/results/results_mle_n_cont_",i,".Rdata")
    dir_result_mle <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/results/results_mle_",i,".Rdata")
    dir_result_smle <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/results/results_smle_",i,".Rdata")
    dir_result_lsmle <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/results/results_lsmle_",i,".Rdata")
    
    saveRDS(sim_exp_mle_n_cont, file=dir_result_mle_n_cont)
    saveRDS(sim_exp_mle, file=dir_result_mle)
    saveRDS(sim_exp_smle, file=dir_result_smle)
    saveRDS(sim_exp_lsmle, file=dir_result_lsmle)
    
    dir_envelope_mle_n_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/envelope_mle_n_cont_",i,".png")
    dir_envelope_mle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/envelope_mle_cont_",i,".png")
    dir_envelope_smle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/envelope_smle_cont_",i,".png")
    dir_envelope_lsmle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/envelope_lsmle_cont_",i,".png")
    
    dir_pond_mle_n_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/pond_mle_n_cont_",i,".png")
    dir_pond_mle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/pond_mle_cont_",i,".png")
    dir_pond_smle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/pond_smle_cont_",i,".png")
    dir_pond_lsmle_cont <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/pond_lsmle_cont_",i,".png")
    
    dir_dispersao_curvas <- paste0("C:/Users/eddus/Documents/UnB/Mestrado/Dissertação/Dados/Analise de Dados/Aplicações/Aplicação com dados simulados/plots/dispersao_curvas_",i,".png")
    
    png(dir_envelope_mle_n_cont, width=878, height=718)
    nlrobust_envelope(object=sim_exp_mle_n_cont, faixa.fixed = c(-3.3,3.3), title_comp = "Sem contaminação")
    dev.off()
    
    png(dir_envelope_mle_cont, width=878, height=718)
    nlrobust_envelope(object=sim_exp_mle, faixa.fixed = c(-3.3,3.3), title_comp = "Com contaminação")
    dev.off()
    
    png(dir_envelope_smle_cont, width=878, height=718)
    nlrobust_envelope(object=sim_exp_smle, faixa.fixed = c(-3,3), title_comp = "Com contaminação")
    dev.off()
    
    png(dir_envelope_lsmle_cont, width=878, height=718)
    nlrobust_envelope(object=sim_exp_lsmle, faixa.fixed = c(-3,3), title_comp = "Com contaminação")
    dev.off()
    
    ################################################################################
    # Ponderações versus resíduos 
    
    png(dir_pond_mle_n_cont, width=878, height=718)
    plot(sim_exp_mle_n_cont$residuals,
         sim_exp_mle_n_cont$weights,
         main = "MLE - Sem contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=19,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_mle_cont, width=878, height=718)
    plot(sim_exp_mle$residuals,
         sim_exp_mle$weights,
         main = "MLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=19,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_smle_cont, width=878, height=718)
    plot(sim_exp_smle$residuals,
         sim_exp_smle$weights,
         main = "SMLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=19,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_lsmle_cont, width=878, height=718)
    plot(sim_exp_lsmle$residuals,
         sim_exp_lsmle$weights,
         main = "LSMLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=19,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    ################################################################################
    # Grafico de dispersao com a curva ajustada
    
    par(mfrow=c(1,1),
        cex.main=2.5,
        cex.lab=2.0,
        cex.sub=1.5,
        cex.axis=1.5)
    
    png(dir_dispersao_curvas, width=878, height=718)
    plot(sim_data$data$cont$data$X.2,
         sim_data$data$cont$data$resposta,
         xlab="Variável explicativa", 
         ylab="Resposta", 
         pch=19,
         ylim = c(-0.1,1.1)
    )
    abline(h = 0)
    curve(plogis(sim_exp_mle_n_cont$coefficients$mean[1]+x**(sim_exp_mle_n_cont$coefficients$mean[2])),
          add = T, lty=6, lwd=2, col="blue")
    curve(plogis(sim_exp_mle$coefficients$mean[1]+x**(sim_exp_mle$coefficients$mean[2])),
          add = T, lty=1, lwd=2, col="red")
    curve(plogis(sim_exp_smle$coefficients$mean[1]+x**(sim_exp_smle$coefficients$mean[2])),
          add = T, lty=4, lwd=2, col="darkgreen")
    curve(plogis(sim_exp_lsmle$coefficients$mean[1]+x**(sim_exp_lsmle$coefficients$mean[2])),
          add = T, lty=5, lwd=2, col="black")
    grid()
    legend(0.15, 1.05,c("MLE - Sem contaminação", 
                        "MLE - Com contaminação", 
                        "SMLE - Com contaminação",
                        "LSMLE - Com contaminação"),
           lty=c(6,1,4,5), 
           lwd=c(2,2,2,2),
           col = c("blue", "red", "darkgreen", "black"), 
           cex = 1.2)
    dev.off()
    
    end_time <- Sys.time()
    total.processing.time <- (end_time - start_time)
    message("\nINFO: Process finished! Total time: ", total.processing.time)
  }
}
