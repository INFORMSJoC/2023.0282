#generate 100 data from DGP.r
data <- readRDS("simu3_data_p140.rds") 

data_B <- data$beta
data_Z <- data$Z
data_Y <- data$Y
p = 140

results_file <- "every_iteration_p140.csv"

write.csv(data.frame(Iteration = integer(0), Method = character(0), Q_Value = numeric(0),
                     FDP = numeric(0), TDP = numeric(0), FDP_Plus = numeric(0), TDP_Plus = numeric(0)),
          file = results_file, row.names = FALSE)


fdp_results <- list()
tdp_results <- list()
fdp_plus_results <- list()
tdp_plus_results <- list()

q_values <- c(0.1, 0.2, 0.3)

iter = length(data_B)

for (ii in 1:iter) {
  
  Z <- data_Z[[ii]]
  y <- data_Y[[ii]]
  beta0 <- data_B[[ii]]
  
  St <- which(beta0 != 0)
  n <- dim(Z)[1]; p <- dim(Z)[2]
  
  #beta_dp <- DP(Z, y, tau)
  beta_dp_de <- DP_de(Z, y, tau)
  beta_rank <- RANK(Z, y)
  pval_bhq <- BHq(Z, y, tau)
  
  cat("Iteration:", ii, "\n")
  
  for (q in q_values) {
    reco_dp_de <- fdRnew(beta_dp_de, St, q)
    reco_rank <- fdRnew(beta_rank, St, q)
    reco_bhq <- bhq(pval_bhq, St, q)
    reco_bhq$f1 <- 0
    reco_bhq$t1 <- 0

    methods_list <- list(
      "dp" = reco_dp_de,
      "rank" = reco_rank,
      "bhq" = reco_bhq
    )
    
    for (method_name in names(methods_list)) {
      reco <- methods_list[[method_name]]
      
      fdp_results[[paste(method_name, q, sep = "_")]][ii] <- reco$f1
      tdp_results[[paste(method_name, q, sep = "_")]][ii] <- reco$t1
      fdp_plus_results[[paste(method_name, q, sep = "_")]][ii] <- reco$f2
      tdp_plus_results[[paste(method_name, q, sep = "_")]][ii] <- reco$t2
      

      iteration_results <- data.frame(
        Iteration = ii,
        Method = method_name,
        Q_Value = q,
        FDP = reco$f1,
        TDP = reco$t1,
        FDP_Plus = reco$f2,
        TDP_Plus = reco$t2
      )
     
      write.table(iteration_results, results_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    }
    
}
  
}

wr_data <- read.csv(results_file)

mean_values <- aggregate(cbind(FDP, FDP_Plus, TDP, TDP_Plus) ~ Method + Q_Value, data = wr_data, FUN = mean)
filtered_data <- mean_values[order(mean_values$Method, mean_values$Q_Value), ]


filtered_data$p <- p
filtered_data$fdr <- filtered_data$FDP
filtered_data$power <- filtered_data$TDP
filtered_data$Method <- gsub("bhq", "BHq", filtered_data$Method)
filtered_data$Method <- gsub("rank", "Rank", filtered_data$Method)
filtered_data$Method <- gsub("dp", "DP", filtered_data$Method)

bhq_plus <- filtered_data[filtered_data$Method == "BHq", ]
bhq_plus$Method <- "BHq"
bhq_plus$fdr <- bhq_plus$FDP_Plus
bhq_plus$power <- bhq_plus$TDP_Plus

rank_plus <- filtered_data[filtered_data$Method == "Rank", ]
rank_plus$Method <- "Rank+"
rank_plus$fdr <- rank_plus$FDP_Plus
rank_plus$power <- rank_plus$TDP_Plus

dp_plus <- filtered_data[filtered_data$Method == "DP", ]
dp_plus$Method <- "DP+"
dp_plus$fdr <- dp_plus$FDP_Plus
dp_plus$power <- dp_plus$TDP_Plus

filtered_data <- filtered_data[filtered_data$Method != "BHq", ]
final_data <- rbind(filtered_data[, c("Method", "Q_Value", "p", "fdr", "power")],
                    bhq_plus[, c("Method", "Q_Value", "p", "fdr", "power")],
                    rank_plus[, c("Method", "Q_Value", "p", "fdr", "power")],
                    dp_plus[, c("Method", "Q_Value", "p", "fdr", "power")])

final_data <- final_data[order(final_data$Q_Value, final_data$Method), ]
rownames(final_data) <- NULL

write.csv(final_data, file ="simu3_ave_p140.csv", row.names = FALSE)
