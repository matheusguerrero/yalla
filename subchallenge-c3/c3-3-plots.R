
##############################################################################
# Create a scatter plots for the real and the simulated data (with size 10^4):
##############################################################################
library(readr)
library(ggplot2)
library(gridExtra)

# read the simualted data of size 10^4:
sims.new <- read.csv("sim.scatterplot.csv") 
X <- read.csv("Data/Coputopia.csv")
X=X[,-c(1,2)]
# Create data frames:
my_data1 <- data.frame(x = X[,1], y = X[,2])
my_data11 <- data.frame(x1 = sims.new[,1], y1 = sims.new[,2])
my_data2 <- data.frame(x = X[,2], y = X[,3])
my_data22 <- data.frame(x1 = sims.new[,2], y1 = sims.new[,3])
my_data3 <- data.frame(x = X[,1], y = X[,3])
my_data33 <- data.frame(x1 = sims.new[,1], y1 = sims.new[,3])


# Scatter plots:
fig1 <- ggplot(my_data1, aes(x = x, y = y,group=1)) +coord_fixed()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+ylim(range(X[,2],X[,1]))+xlim(range(X[,2],X[,1]))+
  geom_point(shape = 20, size = 0.9) +  # Plot original data in black
  geom_vline(xintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_hline(yintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_point(aes(x = x, y = y), shape = 20, size = 0.9, color = "black") +  
  geom_point(data = my_data11, aes(x = x1, y = y1,group=1), shape = 20, size = 0.9,
             color = "#56b1f7") +  # Plot subset data in red
  labs(title = "",
       x = expression(Y[1]),
       y = expression(Y[2])) + guides(alpha = "none", color = "none")
ggsave(fig1, filename =paste0("C3_1.pdf"), height = 7, width = 7, bg = "transparent" )

fig2 <- ggplot(my_data2, aes(x = x, y = y,group=1)) +coord_fixed()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+ylim(range(X[,2],X[,3]))+xlim(range(X[,2],X[,3]))+
  geom_point(shape = 20, size = 0.9) +  # Plot original data in black
  geom_vline(xintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_hline(yintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_point(aes(x = x, y = y), shape = 20, size = 0.9, color = "black") +  
  geom_point(data = my_data22, aes(x = x1, y = y1,group=1), shape = 20, size = 0.9,
             color = "#56b1f7") +  # Plot subset data in red
  labs(title = "",
       x = expression(Y[2]),
       y = expression(Y[3])) + guides(alpha = "none", color = "none")
ggsave(fig2, filename =paste0("C3_2.pdf"), height = 7, width = 7, bg = "transparent" )

fig3 <- ggplot(my_data3, aes(x = x, y = y,group=1)) +coord_fixed()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+ylim(range(X[,1],X[,3]))+xlim(range(X[,1],X[,3]))+
  geom_point(shape = 20, size = 0.9) +  # Plot original data in black
  geom_vline(xintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_hline(yintercept = -log(-log(0.95)), linetype = 1, color = "red", size = 0.3) +
  geom_point(aes(x = x, y = y), shape = 20, size = 0.9, color = "black") +  
  geom_point(data = my_data33, aes(x = x1, y = y1,group=1), shape = 20, size = 0.9,
             color = "#56b1f7") +  # Plot subset data in red
  labs(title = "",
       x = expression(Y[1]),
       y = expression(Y[3])) + guides(alpha = "none", color = "none")
ggsave(fig3, filename =paste0("C3_3.pdf"), height = 7, width = 7, bg = "transparent" )


#################################################
# Diagnostics plots for model fit (QQ-plots)
#################################################

# Calculate the min, max, sum for simulated and real data:
input.sim <- read.csv("sims.results95.csv") 
X <- read.csv("Data/Coputopia.csv")

#min.sim <- apply(input.sim, 1, min)
#max.sim <- apply(input.sim, 1, max)
sum.sim <- apply(input.sim, 1, sum)

# min.data <- apply(X[3:5], 1, min)
# max.data <- apply(X[3:5], 1, max)
sum.data <- apply(X[3:5], 1, sum)

# Calculate quantiles for simulated and real data:
# min.sim1 <- quantile(min.sim, probs = seq(0, 1,  0.01))
# min.data1 <- quantile(min.data, probs = seq(0, 1,  0.01))
# max.sim1 <- quantile(max.sim, probs = seq(0, 1, 0.01))
# max.data1 <- quantile(max.data, probs = seq(0, 1, 0.01))
sum.sim1 <- quantile(sum.sim, probs = seq(0, 1-1/length(sum.data), length=length(sum.data)))
sum.data1 <- quantile(sum.data, probs = seq(0, 1-1/length(sum.data), length=length(sum.data)))

# Create a data frame:
# min_data <- data.frame(x=min.sim1,y=min.data1)
# max_data <- data.frame(x=max.sim1,y=max.data1)
sum_data <- data.frame(x=sum.sim1,y=sum.data1)

# QQ plots:
# qq1 <- ggplot(min_data, aes(x = x, y = y)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.2) +
#   coord_fixed(ratio = 1) +
#   labs(title = "",
#        x = "Theoretical min",
#        y = "Observed min")
# 
# qq2 <- ggplot(max_data, aes(x = x, y = y)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed" , color = "red", size = 0.2) +
#   coord_fixed(ratio = 1) +
#   labs(title = "",
#        x = "Theoretical max",
#        y = "Observed max")

qq3 <- ggplot(sum_data, aes(x = x, y = y)) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+ylim(range(sum.sim1,sum.data1))+xlim(range(sum.sim1,sum.data1))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.2) +
  coord_fixed(ratio = 1) +
  labs(title = "",
       x = expression("Simulated" ~ R),
       y = expression("Empirical" ~ R))+ guides(alpha = "none", color = "none")
ggsave(qq3, filename =paste0("C3_agg.pdf"), height = 7, width = 7, bg = "transparent" )

# Combine the three QQ plots:
#combined_plot1 <- grid.arrange(qq1, qq2, qq3, ncol = 3) 
