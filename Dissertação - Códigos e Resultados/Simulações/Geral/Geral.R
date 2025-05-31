

#################################################################################
# Function to run the experiment


run_scatter_simulation <- function(n,
                                   beta,
                                   gamma,
                                   functions,
                                   parameters,
                                   scenario,
                                   contamination,
                                   link.mu,
                                   link.phi){

  generate_simulated_data(n = n,
                          beta = beta,
                          gamma = gamma,
                          functions = functions,
                          parameters = parameters, 
                          contamination = contamination,
                          scenario = scenario,
                          link.mu = link.mu,
                          link.phi = link.phi) -> sim_data
  
  dir_scatterplot <- paste0("insert_here_you_directory_path",scenario,".png")
  png(dir_scatterplot, width=878, height=718)
  comparison <- (sim_data$data$cont$data$resposta == sim_data$data$n.cont$data$resposta)
  col_index <- which(cumsum(comparison) & rev(cumsum(rev(comparison))) & !comparison)
  par(mfrow=c(1,1),
      cex.main=3.0,
      cex.lab=2.5,
      cex.sub=2.5,
      cex.axis=1.5,
      mgp=c(2.3, 0.5, 0))
  plot(sim_data$data$cont$data$X.2,
       sim_data$data$cont$data$resposta,
       main = paste("Cenário", scenario-1),
       xlab="Variável explicativa", 
       ylab="Resposta", 
       cex = 5.0,
       pch=20,
       cex.lab = 2.5,
       ylim = c(0.00,1.05))
  points(sim_data$data$cont$data$X.2[col_index],
         sim_data$data$cont$data$resposta[col_index],
         col = "red",
         cex = 5.0,
         pch = 20)
  grid()
  dev.off()

}

# Scenario 2
run_scatter_simulation(n=80,
                       beta = c(-1.7, 1.2),
                       gamma = c(2.5, 3.5),
                       functions = c("exp"),
                       parameters = 2, 
                       contamination = 0.05,
                       scenario = 2,
                       link.mu="logit",
                       link.phi="log")

# Scenario 3
run_scatter_simulation(n=80,
                       beta = c(-1.0, 1.0),
                       gamma = c(6.5),
                       functions = c("exp"),
                       parameters = 2, 
                       contamination = 0.05,
                       scenario = 3,
                       link.mu="logit",
                       link.phi="log")

# Scenario 4
run_scatter_simulation(n=80,
                       beta = c(-1.0, -1.4),
                       gamma = c(6.0),
                       functions = c("exp"),
                       parameters = 2, 
                       contamination = 0.05,
                       scenario = 4,
                       link.mu="logit",
                       link.phi="log")

# Scenario 5
run_scatter_simulation(n=80,
                       beta = c(-1.7, 1.2),
                       gamma = c(4.7),
                       functions = c("exp"),
                       parameters = 2, 
                       contamination = 0.05,
                       scenario = 5,
                       link.mu="logit",
                       link.phi="log")

