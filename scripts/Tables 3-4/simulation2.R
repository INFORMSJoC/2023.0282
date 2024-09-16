#generate 100 data from DGP.r
data <- readRDS("simu2_data_p400.rds") 
#"simu2_data_p500.rds" "simu2_data_p600.rds"
data_B <- data$beta
data_Z <- data$Z
data_Y <- data$Y
p = 400 #500 600

results_file <- "every_iteration_p400.csv"

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
  
  beta_dp <- DP(Z, y, tau)
  beta_rank <- RANK(Z, y)
  pval_bhq <- BHq(Z, y, tau)
  
  
  for (q in q_values) {
    reco_dp <- fdRnew(beta_dp, St, q)
    reco_rank <- fdRnew(beta_rank, St, q)
    reco_bhq <- bhq(pval_bhq, St, q)
    #reco_dp <- reco_rank
    
    
    # dp
    fdp_results[[paste("dp", q, sep = "_")]][ii] <- reco_dp$f1
    tdp_results[[paste("dp", q, sep = "_")]][ii] <- reco_dp$t1
    fdp_plus_results[[paste("dp", q, sep = "_")]][ii] <- reco_dp$f2
    tdp_plus_results[[paste("dp", q, sep = "_")]][ii] <- reco_dp$t2
    
    #  rank
    fdp_results[[paste("rank", q, sep = "_")]][ii] <- reco_rank$f1
    tdp_results[[paste("rank", q, sep = "_")]][ii] <- reco_rank$t1
    fdp_plus_results[[paste("rank", q, sep = "_")]][ii] <- reco_rank$f2
    tdp_plus_results[[paste("rank", q, sep = "_")]][ii] <- reco_rank$t2
    
    #  bhq 
    fdp_results[[paste("bhq", q, sep = "_")]][ii] <- 0
    tdp_results[[paste("bhq", q, sep = "_")]][ii] <- 0
    fdp_plus_results[[paste("bhq", q, sep = "_")]][ii] <- reco_bhq$f2
    tdp_plus_results[[paste("bhq", q, sep = "_")]][ii] <- reco_bhq$t2
    
    iteration_results <- data.frame(
      Iteration = ii,
      Method = "dp",
      Q_Value = q,
      FDP = reco_dp$f1,
      TDP = reco_dp$t1,
      FDP_Plus = reco_dp$f2,
      TDP_Plus = reco_dp$t2
    )
    write.table(iteration_results, results_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    
    iteration_results <- data.frame(
      Iteration = ii,
      Method = "rank",
      Q_Value = q,
      FDP = reco_rank$f1,
      TDP = reco_rank$t1,
      FDP_Plus = reco_rank$f2,
      TDP_Plus = reco_rank$t2
    )
    write.table(iteration_results, results_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    
    iteration_results <- data.frame(
      Iteration = ii,
      Method = "bhq",
      Q_Value = q,
      FDP = 0,
      TDP = 0,
      FDP_Plus = reco_bhq$f2,
      TDP_Plus = reco_bhq$t2
    )
    write.table(iteration_results, results_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    
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

write.csv(final_data, file ="simu2_ave_p400.csv", row.names = FALSE)

