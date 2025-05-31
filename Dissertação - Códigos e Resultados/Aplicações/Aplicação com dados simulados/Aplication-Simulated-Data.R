


#################################################################################
# Function to run the experiment

run_simdata_experiment <- function(sample_size = c(40,80,160,320)){
  for(i in sample_size){
    
    closeAllConnections()
    
    start_time <- Sys.time()
    
    # set.seed(965) # Original
    set.seed(566) # Melhor
    generate_simulated_data(n = i,
                            beta = c(-0.6, 0.8),
                            gamma = c(3.9),
                            functions = "power_par",
                            parameters = 2, 
                            contamination = 0.05,
                            scenario = 3,# Decresing=FALSE
                            link.mu = "logit",
                            link.phi = "log") -> sim_data
    
    dir_scatterplot <- paste0("paste_here_your_directory_path",i,".png")
    
    png(dir_scatterplot, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    comparison <- (sim_data$data$cont$data$resposta == sim_data$data$n.cont$data$resposta)
    col_index <- which(cumsum(comparison) & rev(cumsum(rev(comparison))) & !comparison)
    plot(sim_data$data$cont$data$X.2,
         sim_data$data$cont$data$resposta,
         main = paste("Amostra de", i, "observações"),
         xlab="Variável explicativa", 
         ylab="Resposta", 
         cex = 5.0,
         pch=20,
         cex.lab = 2.5,
         ylim = c(-0.05,1.05))
    points(sim_data$data$cont$data$X.2[col_index],
           sim_data$data$cont$data$resposta[col_index],
           col = "red",
           cex = 5.0,
           pch = 20)
    grid()
    dev.off()
    
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
    
    dir_result_mle_n_cont <- paste0("paste_here_your_directory_path",i,".Rdata")
    dir_result_mle <- paste0("paste_here_your_directory_path",i,".Rdata")
    dir_result_smle <- paste0("paste_here_your_directory_path",i,".Rdata")
    dir_result_lsmle <- paste0("paste_here_your_directory_path",i,".Rdata")
    
    saveRDS(sim_exp_mle_n_cont, file=dir_result_mle_n_cont)
    saveRDS(sim_exp_mle, file=dir_result_mle)
    saveRDS(sim_exp_smle, file=dir_result_smle)
    saveRDS(sim_exp_lsmle, file=dir_result_lsmle)
    
    dir_envelope_mle_n_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_envelope_mle_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_envelope_smle_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_envelope_lsmle_cont <- paste0("paste_here_your_directory_path",i,".png")
    
    dir_pond_mle_n_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_pond_mle_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_pond_smle_cont <- paste0("paste_here_your_directory_path",i,".png")
    dir_pond_lsmle_cont <- paste0("paste_here_your_directory_path",i,".png")
    
    dir_dispersao_curvas <- paste0("paste_here_your_directory_path",i,".png")
    
    png(dir_envelope_mle_n_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    nlrobust_envelope(object=sim_exp_mle_n_cont,
                      faixa.fixed = c(-10.0,8.5),
                      title_comp = "Sem contaminação")
    dev.off()

    png(dir_envelope_mle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    nlrobust_envelope(object=sim_exp_mle,
                      faixa.fixed = c(-10.0,8.5),
                      title_comp = "Com contaminação")
    dev.off()

    png(dir_envelope_smle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    nlrobust_envelope(object=sim_exp_smle,
                      faixa.fixed = c(-10.0,8.5),
                      title_comp = "Com contaminação")
    dev.off()

    png(dir_envelope_lsmle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    nlrobust_envelope(object=sim_exp_lsmle,
                      faixa.fixed = c(-10.0,8.5),
                      title_comp = "Com contaminação")
    dev.off()
    
    ################################################################################
    # Ponderações versus Resíduos 
    
    png(dir_pond_mle_n_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    plot(sim_exp_mle_n_cont$residuals,
         sim_exp_mle_n_cont$weights,
         main = "MLE - Sem contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         cex = 3.5,
         cex.lab = 2.5,
         pch=20,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_mle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    plot(sim_exp_mle$residuals,
         sim_exp_mle$weights,
         main = "MLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=20,
         cex = 3.5,
         cex.lab = 2.5,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_smle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    plot(sim_exp_smle$residuals,
         sim_exp_smle$weights,
         main = "SMLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=20,
         cex = 3.5,
         cex.lab = 2.5,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    png(dir_pond_lsmle_cont, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    plot(sim_exp_lsmle$residuals,
         sim_exp_lsmle$weights,
         main = "LSMLE - Com contaminação",
         xlab="Resíduos", 
         ylab="Ponderações", 
         pch=20,
         cex = 3.5,
         cex.lab = 2.5,
         xlim = c(-9,9),
         ylim = c(-0.1,1.3)
    )
    abline(h = 0)
    grid()
    dev.off()
    
    ################################################################################
    # Grafico de dispersao com a curva ajustada
    
    png(dir_dispersao_curvas, width=878, height=718)
    par(mfrow=c(1,1),
        cex.main=3.0,
        cex.lab=2.5,
        cex.sub=2.5,
        cex.axis=1.5,
        mgp=c(2.3, 0.5, 0))
    plot(sim_data$data$cont$data$X.2,
         sim_data$data$cont$data$resposta,
         main = paste("Amostra de", i, "observações"),
         xlab="Variável explicativa", 
         ylab="Resposta", 
         pch=20,
         cex = 3.0,
         cex.lab = 2.5,
         ylim = c(-0.1,1.1)
    )
    abline(h = 0)
    curve(plogis(sim_exp_mle_n_cont$coefficients$mean[1]+x**(sim_exp_mle_n_cont$coefficients$mean[2])),
          add = T, lty=6, lwd=3, col="blue")
    curve(plogis(sim_exp_mle$coefficients$mean[1]+x**(sim_exp_mle$coefficients$mean[2])),
          add = T, lty=1, lwd=3, col="red")
    curve(plogis(sim_exp_smle$coefficients$mean[1]+x**(sim_exp_smle$coefficients$mean[2])),
          add = T, lty=4, lwd=3, col="darkgreen")
    curve(plogis(sim_exp_lsmle$coefficients$mean[1]+x**(sim_exp_lsmle$coefficients$mean[2])),
          add = T, lty=5, lwd=3, col="black")
    grid()
    legend(0.15, 1.05,c("MLE - Sem contaminação", 
                        "MLE - Com contaminação", 
                        "SMLE - Com contaminação",
                        "LSMLE - Com contaminação"),
           lty=c(6,1,4,5), 
           lwd=c(2,2,2,2),
           col = c("blue", "red", "darkgreen", "black"), 
           cex = 2.0)
    dev.off()
    
    end_time <- Sys.time()
    total.processing.time <- (end_time - start_time)
    message("\nINFO: Process finished! Total time: ", total.processing.time)
  }
}

################################################################################
# Runing the experiment

run_simdata_experiment(sample_size = c(40,80,160,320))


##############
# Summary

results_mle_n_cont_40 |> summary_nlrobustbetareg()
results_mle_n_cont_80 |> summary_nlrobustbetareg()
results_mle_n_cont_160 |> summary_nlrobustbetareg()
results_mle_n_cont_320 |> summary_nlrobustbetareg()

results_mle_40 |> summary_nlrobustbetareg()
results_mle_80 |> summary_nlrobustbetareg()
results_mle_160 |> summary_nlrobustbetareg()
results_mle_320 |> summary_nlrobustbetareg()


results_smle_40 |> summary_nlrobustbetareg()
results_smle_80 |> summary_nlrobustbetareg()
results_smle_160 |> summary_nlrobustbetareg()
results_smle_320 |> summary_nlrobustbetareg()

results_lsmle_40 |> summary_nlrobustbetareg()
results_lsmle_80 |> summary_nlrobustbetareg()
results_lsmle_160 |> summary_nlrobustbetareg()
results_lsmle_320 |> summary_nlrobustbetareg()


