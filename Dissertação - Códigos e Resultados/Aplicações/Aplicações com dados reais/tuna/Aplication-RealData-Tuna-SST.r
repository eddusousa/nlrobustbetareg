
data_tuna <- read.table("dados/data_2000_Indian.txt",h=T) #original data

View(data_tuna)

attach(data_tuna)
tropicalP <-  trop/100 
index0<- 1:length(tropicalP)
index0[tropicalP==0] #no observation equal to zero
index1<- 1:length(tropicalP)
index1[tropicalP==1] #one observation equal to one: #46
tropicalP[c(46)] <- 0.999 #replacement value
Southern_Indian <- tropicalP[lat<0] #response variable
SST_Southern_Indian <- SST[lat<0] #covariate
SST <- log(SST[lat<0]) #covariate

#**************************************************#
#******Beta regression with constant precision*****#
#**************************************************#
y <- Southern_Indian 
X <- matrix(c(rep(1,77), SST_Southern_Indian, SST), ncol=3, byrow=F);
Z <- matrix(c(rep(1,77)), ncol=1, byrow=F);
#*************************Subsets without outliers******************************#
y46 <- y[-c(46)]; X46<- X[-c(46),];Z46 <- as.matrix(Z[-c(46),])
#********************************************************************************#

data.frame(y,
           X,
           Z) -> df_tuna_full

View(df_tuna_full)

data.frame(y46,
           X46,
           Z46) -> df_tuna_w46

View(df_tuna_w46)


########################################################################################################
# Full data

functions_exp <- c("power_par")
parameters_exp <- c(2)

robustbetareg(formula = y~X3|1,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'MLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = FALSE) -> fit_tuna_mle_full

robustbetareg(formula = y~X3|1,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'SMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = FALSE) -> fit_tuna_smle_full

fit_tuna_smle_full |> summary()

robustbetareg(formula = y~X2|1,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'LSMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = FALSE) -> fit_tuna_lsmle_full

fit_tuna_mle_full |> summary_nlrobustbetareg()
fit_tuna_mle_w46 |> summary_nlrobustbetareg()
fit_tuna_smle_full |> summary_nlrobustbetareg()
fit_tuna_lsmle_full |> summary_nlrobustbetareg()

################################################################################
# Envelope simulado

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_mle_full, 
                  faixa.fixed = c(-4.0,8.2), 
                  title_comp = "Completo")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_smle_full, 
                  faixa.fixed = c(-4.0,8.2), 
                  title_comp = "Completo")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_lsmle_full,
                  faixa.fixed = c(-4.0,8.2),
                  title_comp = "Completo")

################################################################################
#*********************SCATTER PLOTS***********************#
values3 <- seq(15,29,length.out=1000)
###########################################################
# FULL DATA
#***MLE fit full data***#
beta_mle <- fit_tuna_mle_full$coefficients$mean
gamma_mle <- fit_tuna_mle_full$coefficients$precision
eta_mle <- beta_mle[1] +  values3**beta_mle[2]
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
phi_mle <- exp(gamma_mle)

#*****SMLE fit full data*****#
beta_smle <- fit_tuna_smle_full$coefficients$mean
gamma_smle <- fit_tuna_smle_full$coefficients$precision
eta_smle <- beta_smle[1] +  values3**beta_smle[2]
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
phi_smle <- exp(gamma_smle)

#*****LSMLE fit full data*****#
beta_lsmle <- fit_tuna_lsmle_full$coefficients$mean
gamma_lsmle <- fit_tuna_lsmle_full$coefficients$precision
eta_lsmle <-beta_lsmle[1] +  values3**beta_lsmle[2]
mu_lsmle <- exp(eta_lsmle)/(1+exp(eta_lsmle))
phi_lsmle <- exp(gamma_lsmle)

#***SCATTER PLOT****#
par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
yoriginal <- y; 
yoriginal[46] <- 1
plot(SST_Southern_Indian, yoriginal, 
     pch=20, 
     xlab="SST",
     ylab="TTP",
     ylim=c(-0.1,1.1), 
     xlim=c(15,29),
     cex=3.0, 
     cex.lab=2.5, 
     cex.axis=1.5, 
     cex.main=2.0)
lines(values3, mu_mle, lwd=2,col=2,lty=1)
lines(values3,mu_smle, lwd=2,col=6,lty=4)
lines(values3,mu_lsmle, lwd=2,col=4,lty=3)
identify(SST_Southern_Indian, yoriginal, cex=1.3)
grid()
legend(16,1.0,c("MLE","SMLE","LSMLE"),
       col=c(2,6,4),lty=c(1,4,3), cex=1.5, lwd=c(2,2,2))
################################################################################
#*****************RESIDUALS AGAINST WEIGHTS******************#
#***FULL DATA***
#***MLE***
weights_mle <- fit_tuna_mle_full$weights
qres_mle <- fit_tuna_mle_full$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_mle,
     weights_mle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="MLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_mle,weights_mle,cex=1.3)

#***SMLE***
weights_smle <- fit_tuna_smle_full$weights
qres_smle <- fit_tuna_smle_full$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_smle,
     weights_smle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="SMLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_smle,weights_smle,cex=1.3)

#***LSMLE***
weights_lsmle <- fit_tuna_lsmle_full$weights
qres_lsmle <- fit_tuna_lsmle_full$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_lsmle,
     weights_lsmle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="LSMLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_lsmle,weights_lsmle,cex=1.3)

########################################################################################################
# Without observation 46

robustbetareg(formula = y46~X3|1,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'MLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_mle_w46

robustbetareg(formula = y46~X3|1,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'SMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_smle_w46

robustbetareg(formula = y46~X3|1,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'LSMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_lsmle_w46

fit_tuna_mle_w46 |> summary()
fit_tuna_smle_w46 |> summary()
fit_tuna_lsmle_w46 |> summary()

fit_tuna_mle_w46 |> summary_nlrobustbetareg()
fit_tuna_smle_w46 |> summary_nlrobustbetareg()
fit_tuna_lsmle_w46 |> summary_nlrobustbetareg()

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_mle_w46, 
                  faixa.fixed = c(-4.0,4.0), 
                  title_comp = "Sem obs. 46")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_smle_w46, 
                  faixa.fixed = c(-4.0,4.0), 
                  title_comp = "Sem obs. 46")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_lsmle_w46, 
                  faixa.fixed = c(-4.0,4.0), 
                  title_comp = "Sem obs. 46")

##############################################################################################
#*********************SCATTER PLOTS***********************#
# Without obs 46
#***MLE fit w46***#
beta_mle_w46 <- fit_tuna_mle_w46$coefficients$mean
gamma_mle_w46 <- fit_tuna_mle_w46$coefficients$precision
eta_mle_w46 <- beta_mle_w46[1] +  values3**beta_mle_w46[2]
mu_mle_w46 <- exp(eta_mle_w46)/(1+exp(eta_mle_w46))
phi_mle_w46 <- exp(gamma_mle_w46)

#*****SMLE fit w46*****#
beta_smle_w46 <- fit_tuna_smle_w46$coefficients$mean
gamma_smle_w46 <- fit_tuna_smle_w46$coefficients$precision
eta_smle_w46 <- beta_smle_w46[1] +  values3**beta_smle_w46[2]
mu_smle_w46 <- exp(eta_smle_w46)/(1+exp(eta_smle_w46))
phi_smle_w46 <- exp(gamma_smle_w46)

#*****LSMLE fit w46*****#
beta_lsmle_w46 <- fit_tuna_lsmle_w46$coefficients$mean
gamma_lsmle_w46 <- fit_tuna_lsmle_w46$coefficients$precision
eta_lsmle_w46 <- beta_lsmle_w46[1] +  values3**beta_lsmle_w46[2]
mu_lsmle_w46 <- exp(eta_lsmle_w46)/(1+exp(eta_lsmle_w46))
phi_lsmle_w46 <- exp(gamma_lsmle_w46)


#***SCATTER PLOT****#
par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
yoriginal <- y; 
yoriginal[46] <- 1
plot(SST_Southern_Indian, yoriginal, 
     pch=20, 
     xlab="SST",
     ylab="TTP",
     ylim=c(-0.1,1.1), 
     xlim=c(15,29),
     cex=3.0, 
     cex.lab=2.5, 
     cex.axis=1.5, 
     cex.main=2.0)
lines(values3, mu_mle_w46, lwd=2,col=2,lty=1)
lines(values3,mu_smle_w46, lwd=2,col=6,lty=4)
lines(values3,mu_lsmle_w46, lwd=2,col=4,lty=3)
identify(SST_Southern_Indian, yoriginal, cex=1.3)
grid()
legend(16,1.0,c("MLE s/ 46","SMLE s/ 46","LSMLE s/ 46"),
       col=c(2,6,4),lty=c(1,4,3), cex=1.5, lwd=c(2,2,2))


################################################################################

#*****************RESIDUALS AGAINST WEIGHTS******************#
#***WITHOUT 46 DATA***
#***MLE***
weights_mle <- fit_tuna_mle_w46$weights
qres_mle <- fit_tuna_mle_w46$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_mle,
     weights_mle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="MLE - Sem obs. 46",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_mle,weights_mle,cex=1.3)

#***SMLE***
weights_smle <- fit_tuna_smle_w46$weights
qres_smle <- fit_tuna_smle_w46$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_smle,
     weights_smle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="SMLE - Sem obs. 46",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_smle,weights_smle,cex=1.3)

#***LSMLE***
weights_lsmle <- fit_tuna_lsmle_w46$weights
qres_lsmle <- fit_tuna_lsmle_w46$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_lsmle,
     weights_lsmle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="LSMLE - Sem obs. 46",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_lsmle,weights_lsmle,cex=1.3)

########################################################################################################
# Full data with var precision

robustbetareg(formula = y~X3|X3,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'MLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_mle_full_prec_var

robustbetareg(formula = y~X3|X3,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'SMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_smle_full_prec_var

# alpha fixado
robustbetareg(formula = y~X3|X3,
              data = df_tuna_full,
              alpha = 0.04,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'SMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_smle_full_prec_var_2

robustbetareg(formula = y~X3|X3,
              data = df_tuna_full,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'LSMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_lsmle_full_prec_var

fit_tuna_mle_full_prec_var |> summary()
fit_tuna_smle_full_prec_var |> summary()
fit_tuna_smle_full_prec_var_2 |> summary()
fit_tuna_lsmle_full_prec_var |> summary()

fit_tuna_mle_full_prec_var |> summary_nlrobustbetareg()
fit_tuna_smle_full_prec_var |> summary_nlrobustbetareg()
fit_tuna_smle_full_prec_var_2 |> summary_nlrobustbetareg()
fit_tuna_lsmle_full_prec_var |> summary_nlrobustbetareg()

################################################################################
# Envelope simulado

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_mle_full_prec_var, 
                  faixa.fixed = c(-4.0,8.2), 
                  title_comp = "Completo")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_smle_full_prec_var, 
                  faixa.fixed = c(-4.0,8.2), 
                  title_comp = "Completo")

par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_lsmle_full_prec_var, 
                  faixa.fixed = c(-4.0,8.2), 
                  title_comp = "Completo")

################################################################################
#*********************SCATTER PLOTS***********************#
values3 <- seq(15,29,length.out=1000)

###########################################################
# Full data with carying precision
#***MLE fit full data with varying precision***#
beta_mle_full_prec_var <- fit_tuna_mle_full_prec_var$coefficients$mean
gamma_mle_full_prec_var <- fit_tuna_mle_full_prec_var$coefficients$precision
eta_mle_full_prec_var <- beta_mle_full_prec_var[1] +  values3**beta_mle_full_prec_var[2]
mu_mle_full_prec_var <- exp(eta_mle_full_prec_var)/(1+exp(eta_mle_full_prec_var))
phi_mle_full_prec_var <- exp(gamma_mle_full_prec_var)

#***SMLE fit full data with varying precision***#
beta_smle_full_prec_var <- fit_tuna_smle_full_prec_var$coefficients$mean
gamma_smle_full_prec_var <- fit_tuna_smle_full_prec_var$coefficients$precision
eta_smle_full_prec_var <- beta_smle_full_prec_var[1] +  values3**beta_smle_full_prec_var[2]
mu_smle_full_prec_var <- exp(eta_smle_full_prec_var)/(1+exp(eta_smle_full_prec_var))
phi_smle_full_prec_var <- exp(gamma_smle_full_prec_var)

#***LSMLE fit full data with varying precision***#
beta_lsmle_full_prec_var <- fit_tuna_lsmle_full_prec_var$coefficients$mean
gamma_lsmle_full_prec_var <- fit_tuna_lsmle_full_prec_var$coefficients$precision
eta_lsmle_full_prec_var <- beta_lsmle_full_prec_var[1] +  values3**beta_lsmle_full_prec_var[2]
mu_lsmle_full_prec_var <- exp(eta_lsmle_full_prec_var)/(1+exp(eta_lsmle_full_prec_var))
phi_lsmle_full_prec_var <- exp(gamma_lsmle_full_prec_var)

#***SCATTER PLOT****#
par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
yoriginal <- y; 
yoriginal[46] <- 1
plot(SST_Southern_Indian, yoriginal, 
     pch=20, 
     xlab="SST",
     ylab="TTP",
     ylim=c(-0.1,1.1), 
     xlim=c(15,29),
     cex=3.0, 
     cex.lab=2.5, 
     cex.axis=1.5, 
     cex.main=2.0)
lines(values3, mu_mle_full_prec_var, lwd=2,col=2,lty=1)
lines(values3,mu_smle_full_prec_var, lwd=2,col=6,lty=4)
lines(values3,mu_lsmle_full_prec_var, lwd=2,col=4,lty=3)
identify(SST_Southern_Indian, yoriginal, cex=1.3)
grid()
legend(16,1.0,c("MLE","SMLE","LSMLE"),
       col=c(2,6,4),lty=c(1,4,3), cex=1.5, lwd=c(2,2,2))

################################################################################
#*****************RESIDUALS AGAINST WEIGHTS******************#
#***FULL DATA***
#***MLE***
weights_mle <- fit_tuna_mle_full_prec_var$weights
qres_mle <- fit_tuna_mle_full_prec_var$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_mle,
     weights_mle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="MLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,2.4)
)
grid()
identify(qres_mle,weights_mle,cex=1.3)

#***SMLE***
weights_smle <- fit_tuna_smle_full_prec_var$weights
qres_smle <- fit_tuna_smle_full_prec_var$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_smle,
     weights_smle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="SMLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,2.4)
)
grid()
identify(qres_smle,weights_smle,cex=1.3)

#***LSMLE***
weights_lsmle <- fit_tuna_lsmle_full_prec_var$weights
qres_lsmle <- fit_tuna_lsmle_full_prec_var$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_lsmle,
     weights_lsmle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="LSMLE - Dados completos",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,1.4)
)
grid()
identify(qres_lsmle,weights_lsmle,cex=1.3)

################################################################################
################################################################################
################################################################################
# Sem a observação 46

robustbetareg(formula = y46~X3|X3,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'MLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = TRUE) -> fit_tuna_mle_w46_prec_var

robustbetareg(formula = y46~X3|X3,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'SMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = FALSE) -> summary()

robustbetareg(formula = y46~X3|X3,
              data = df_tuna_w46,
              alpha = NULL,
              functions = functions_exp, 
              parameters = parameters_exp,
              link = "logit",
              link.phi = "log",
              type = 'LSMLE',
              se.bootstrap = FALSE,
              simulation = FALSE,
              wald.test.bootstrap = FALSE) -> summary()

################################################################################
fit_tuna_mle_w46_prec_var |> summary()

fit_tuna_mle_w46_prec_var |> summary_nlrobustbetareg()
################################################################################
par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
nlrobust_envelope(object=fit_tuna_mle_w46_prec_var, 
                  faixa.fixed = c(-4.0,4.0), 
                  title_comp = "Sem obs. 46")

################################################################################
#*********************SCATTER PLOTS***********************#
values <- seq(2.7,3.4,length.out=1000)
values2 <- seq(1.0,1.22,length.out=1000)
values3 <- seq(15,29,length.out=1000)

###########################################################
# Withoput observation 46 data with varying precision
#***MLE fit full data with varying precision***#
beta_mle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$mean
gamma_mle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$precision
eta_mle_w46_prec_var <- beta_mle_w46_prec_var[1] +  values3**beta_mle_w46_prec_var[2]
mu_mle_w46_prec_var <- exp(eta_mle_w46_prec_var)/(1+exp(eta_mle_w46_prec_var))
phi_mle_w46_prec_var <- exp(gamma_mle_w46_prec_var)

#***SMLE fit w46 data with varying precision***#
beta_smle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$mean
gamma_smle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$precision
eta_smle_w46_prec_var <- beta_smle_w46_prec_var[1] +  values3**beta_smle_w46_prec_var[2]
mu_smle_w46_prec_var <- exp(eta_smle_w46_prec_var)/(1+exp(eta_smle_w46_prec_var))
phi_smle_w46_prec_var <- exp(gamma_smle_w46_prec_var)

#***LSMLE fit w46 data with varying precision***#
beta_lsmle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$mean
gamma_lsmle_w46_prec_var <- fit_tuna_mle_w46_prec_var$coefficients$precision
eta_lsmle_w46_prec_var <- beta_lsmle_w46_prec_var[1] +  values3**beta_lsmle_w46_prec_var[2]
mu_lsmle_w46_prec_var <- exp(eta_lsmle_w46_prec_var)/(1+exp(eta_lsmle_w46_prec_var))
phi_lsmle_w46_prec_var <- exp(gamma_lsmle_w46_prec_var)

#***SCATTER PLOT****#
par(mfrow=c(1,1),
    cex.main=3.0,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
yoriginal <- y; 
yoriginal[46] <- 1
plot(SST_Southern_Indian, yoriginal, 
     pch=20, 
     xlab="SST",
     ylab="TTP",
     ylim=c(-0.1,1.1), 
     xlim=c(15,29),
     cex=3.0, 
     cex.lab=2.5, 
     cex.axis=1.5, 
     cex.main=2.0)
lines(values3, mu_mle_w46_prec_var, lwd=2,col=2,lty=1)
lines(values3,mu_smle_w46_prec_var, lwd=2,col=6,lty=4)
lines(values3,mu_lsmle_w46_prec_var, lwd=2,col=4,lty=3)
identify(SST_Southern_Indian, yoriginal, cex=1.3)
grid()
legend(16,1.0,c("MLE s/ 46","SMLE s/46","LSMLE s/46"),
       col=c(2,6,4),lty=c(1,4,3), cex=1.5, lwd=c(2,2,2))

################################################################################
#*****************RESIDUALS AGAINST WEIGHTS******************#
#***WITHOUT 46 DATA***
#***MLE***
weights_mle <- fit_tuna_mle_w46_prec_var$weights
qres_mle <- fit_tuna_mle_w46_prec_var$residuals
par(mfrow=c(1,1),
    cex.main=2.4,
    cex.lab=2.5,
    cex.sub=2.5,
    cex.axis=1.5,
    mgp=c(2.3, 0.5, 0))
plot(qres_mle,
     weights_mle,
     xlab="Resíduos", 
     ylab="Ponderações",
     main="MLE - Sem obs. 46",
     cex = 3.5,
     cex.lab = 1.9,
     pch=20,
     xlim=c(-4,9),
     ylim=c(-0.1,2.4)
)
grid()
identify(qres_mle,weights_mle,cex=1.3)
