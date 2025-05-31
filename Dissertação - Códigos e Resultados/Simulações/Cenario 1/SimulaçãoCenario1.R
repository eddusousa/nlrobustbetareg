


################################################################################################
# Simulação


require(compiler)

compiled_simulate_nlrobust <- cmpfun(simulate_nlrobust)

enableJIT(3)

is.compile <- function(func)
{
  # this function lets us know if a function has been byte-coded or not
  #If you have a better idea for how to do this - please let me know...
  if(class(func) != "function") stop("You need to enter a function")
  last_2_lines <- tail(capture.output(func),2)
  any(grepl("bytecode:", last_2_lines)) # returns TRUE if it finds the text "bytecode:" in any of the last two lines of the function's print
}

is.compile(simulate_nlrobust)

is.compile(compiled_simulate_nlrobust)

se1 = 158 ; se2 = 2 #random seed
set.seed(c(se1,se2), kind="Marsaglia-Multicarry") #To ensure repeatability of the experiment
simulation_results_scenario2 <- compiled_simulate_nlrobust(N=1000,
                                                           n=c(40,80,160,320),
                                                           functions = c("exp"),
                                                           parameters = c(2),
                                                           contamination=0.05,
                                                           scenario = 2,
                                                           run_MLE_not_contaminated=TRUE,
                                                           run_MLE_contaminated=TRUE,
                                                           run_SMLE_not_contaminated=TRUE,
                                                           run_SMLE_contaminated=TRUE,
                                                           run_LSMLE_not_contaminated=TRUE,
                                                           run_LSMLE_contaminated=TRUE,
                                                           link.mu="logit",
                                                           link.phi="log")



################################################################################
# Graficos Dissertação
par(mfrow=c(1,1),
    cex.main=2.5,
    cex.lab=2.0,
    cex.sub=1.5,
    cex.axis=1.5)

###############################################################################
# MLE

#######################
# beta 1

boxplot(simulation_results_scenario2$data$MLE$beta1,
        main=expression( paste("MLE for ", beta[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-2.05, -1.15), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -1.98,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -1.98,"Presença de contaminação",bty="n",cex=2.3)

#######################
# beta 2
boxplot(simulation_results_scenario2$data$MLE$beta2,
        main=expression( paste("MLE for ", beta[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.60, 1.40), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.66,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.66,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 1

boxplot(simulation_results_scenario2$data$MLE$gamma1,
        main=expression( paste("MLE for ", gamma[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.0, 6.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.50,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.50,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 2

boxplot(simulation_results_scenario2$data$MLE$gamma2,
        main=expression( paste("MLE for ", gamma[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-6.0, 8.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -5.00,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -5.00,"Presença de contaminação",bty="n",cex=2.3)

################################################################################
# SMLE

#######################
# beta 1

boxplot(simulation_results_scenario2$data$SMLE$beta1,
        main=expression( paste("SMLE for ", beta[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-2.05, -1.15), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -1.98,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -1.98,"Presença de contaminação",bty="n",cex=2.3)

#######################
# beta 2

boxplot(simulation_results_scenario2$data$SMLE$beta2,
        main=expression( paste("SMLE for ", beta[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.60, 1.40), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.66,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.66,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 1

boxplot(simulation_results_scenario2$data$SMLE$gamma1,
        main=expression( paste("SMLE for ", gamma[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.0, 6.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.5,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.5,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 2

boxplot(simulation_results_scenario2$data$SMLE$gamma2,
        main=expression( paste("SMLE for ", gamma[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-6.0, 8.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -5.00,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -5.00,"Presença de contaminação",bty="n",cex=2.3)


################################################################################
# LSMLE

#######################
# beta 1

boxplot(simulation_results_scenario2$data$LSMLE$beta1,
        main=expression( paste("LSMLE for ", beta[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-2.05, -1.15), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -1.98,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -1.98,"Presença de contaminação",bty="n",cex=2.3)
dev.off()

#######################
# beta 2

boxplot(simulation_results_scenario2$data$LSMLE$beta2,
        main=expression( paste("LSMLE for ", beta[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.60, 1.40), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$beta[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.66,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.66,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 1

boxplot(simulation_results_scenario2$data$LSMLE$gamma1,
        main=expression( paste("LSMLE for ", gamma[1])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.0, 6.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[1], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.50,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.50,"Presença de contaminação",bty="n",cex=2.3)

#######################
# gamma 2

boxplot(simulation_results_scenario2$data$LSMLE$gamma2,
        main=expression( paste("LSMLE for ", gamma[2])),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(-6.0, 8.5), 
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(h = simulation_results_scenario2$coefficients$gamma[2], col = "red", lwd=2, lty=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, -5.00,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, -5.00,"Presença de contaminação",bty="n",cex=2.3)


################################################################################
# Extração alphas

extract_tuning(object = simulation_results_scenario2,
               return_q = TRUE) -> constantes_cenario2

#######################
# SMLE

boxplot(constantes_cenario2$data$SMLE,
        main=expression(paste("Valor ótimo de q para o SMLE")),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.5, 1),
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.54,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.54,"Presença de contaminação",bty="n",cex=2.3)

#######################
# LSMLE

boxplot(constantes_cenario2$data$LSMLE,
        main=expression(paste("Valor ótimo de q para o LSMLE")),
        xlab="Tamanho amostral",
        ylab="",
        names=c(40,80,160,320,40,80,160,320),
        outline=T,
        ylim=c(0.5, 1),
        cex=2.5,
        cex.lab=3.0,
        cex.sub = 1.5,
        cex.axis=2.0,
        cex.main=3.5, 
        boxlwd = 2,
        pch=16, 
        col = NULL)
abline(v=4.5,lty=3,lwd=2)
legend(0.0, 0.54,"Ausência de contaminação",bty="n",cex=2.3)
legend(4.3, 0.54,"Presença de contaminação",bty="n",cex=2.3)


################################################################################
# Tabelas Dissertação


################################################################################


extract_measures(object = simulation_results_scenario2,
                 sample_size = 40) -> measures40

measures40$vies$not.contaminated
measures40$vies$contaminated
measures40$eqm$not.contaminated
measures40$eqm$contaminated


extract_measures(object = simulation_results_scenario2,
                 sample_size = 80) -> measures80

measures80$vies$not.contaminated
measures80$vies$contaminated
measures80$eqm$not.contaminated
measures80$eqm$contaminated

extract_measures(object = simulation_results_scenario2,
                 sample_size = 160) -> measures160

measures160$vies$not.contaminated
measures160$vies$contaminated
measures160$eqm$not.contaminated
measures160$eqm$contaminated

extract_measures(object = simulation_results_scenario2,
                 sample_size = 320) -> measures320

measures320$vies$not.contaminated
measures320$vies$contaminated
measures320$eqm$not.contaminated
measures320$eqm$contaminated

