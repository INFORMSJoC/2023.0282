library(ggplot2)
library(gridExtra)

file_names <- c("simu3_ave_p140.csv", "simu3_ave_p180.csv", "simu3_ave_p220.csv", "simu3_ave_p260.csv", "simu3_ave_p300.csv")
p_values <- c(140, 180, 220, 260, 300)

combined_data <- data.frame()

for (i in 1:length(file_names)) {

  temp_data <- read.csv(file_names[i])
  temp_data$p <- p_values[i]
  combined_data <- rbind(combined_data, temp_data)
}


q_values <- c(0.1, 0.2, 0.3)
for (q in q_values) {
  q_value_filtered <- combined_data[combined_data$Q_Value == q, ]
  q_value_filtered$Method <- factor(q_value_filtered$Method, levels = c("BHq", "Rank", "Rank+", "DP", "DP+"))
  

  methods <- unique(q_value_filtered$Method)
  result <- data.frame(Method = methods)
  

  for (p in p_values) {
    current_p_data <- q_value_filtered[q_value_filtered$p == p, ]
    result <- cbind(result,
                    p = rep(p, nrow(result)),
                    fdr = current_p_data$fdr[match(result$Method, current_p_data$Method)],
                    power = current_p_data$power[match(result$Method, current_p_data$Method)])
  }
  result <- result[c(1,4,5,2,3),]
  output_file <- paste0("reshaped_q_", q, ".csv")
  write.csv(result, output_file, row.names = FALSE)
}

#plot
df1 <- read.csv("reshaped_q_0.1.csv")
df1 <- read.csv("reshaped_q_0.2.csv")
df1 <- read.csv("reshaped_q_0.3.csv")

qq <- 0.1
qq <- 0.2
qq <- 0.3

FDR <- numeric(0)
Power <- numeric(0)
for(i in 1:nrow(df1)){
  FDR <- append(FDR,  as.numeric(df1[i, c(3,6,9,12,15) ]))
  Power <- append(Power,  as.numeric(df1[i, c(4,7,10,13,16) ]))
}

data <- data.frame(
  Method = rep(c("BHq", "Rank", "Rank+", "DP", "DP+"), each = 5),
  Dimension = rep(c(140, 180, 220, 260, 300), times = 5),
  FDR = FDR,
  Power =  Power
)

colors <- c("darkviolet", "deepskyblue", "limegreen", "gold", "tomato")

# Create FDR plot
fdr_plot <- ggplot(data, aes(x = Dimension, y = FDR, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = qq, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Dimension", y = "FDR") +
  theme_minimal() +
  theme(legend.position = "bottom right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.grid.major = element_line(color = "gray", size = 0.1), # Add major grid lines
        panel.grid.minor = element_line(color = "gray", size = 0.1)
)

# Create Power plot
power_plot <- 
  ggplot(data, aes(x = Dimension, y = Power, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors) +
  labs(x = "Dimension", y = "Power") +
  theme_bw() +
  theme(legend.box.background = element_rect(colour = "gray", fill = NA),
        legend.justification = c(0, 0), 
        legend.position = c(0.75, 0.102),
        #legend.position = c(0.75, 0.082),
        legend.key.size = unit(0.4, "lines"),         
        legend.key.height = unit(0.3, "lines"),          
        legend.key.width = unit(0.6, "lines"), 
        legend.background = element_rect(fill = "white", colour = "gray"), 
        legend.key = element_rect(fill = "white", colour = "white"),        
        panel.border = element_rect(colour = "gray", fill=NA, size=1),
        panel.grid.major = element_line(color = "gray", size = 0.1),
        panel.grid.minor = element_line(color = "gray", size = 0.1), 
  ) +
  scale_y_continuous(limits = c(0.85, 1), breaks = seq(0.85, 1, by = 0.02))

grid.arrange(fdr_plot, power_plot, ncol = 2)

