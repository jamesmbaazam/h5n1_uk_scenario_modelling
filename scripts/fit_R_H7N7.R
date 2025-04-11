#Code adapated from Epichains "how to estimate R and k from cluster size data"
#https://github.com/epiverse-trace/howto/blob/3831f8c9f37b395f43a9385fb654798931046dcf/analyses/quantify_transmission/reproduction_number_cluster_size.qmd

# check whether {pak} is installed
if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/epiparameter")

# Load required packages
library(epichains)
library(MCMCpack)
library(epiparameter)
library(coda)
library(ggplot2)
library(dplyr)



# Define H7N7 clusters https://pmc.ncbi.nlm.nih.gov/articles/PMC337057/#abstract1
H7_clusters = c(rep(1,84),c(5,1))  

# Show summary table of frequencies
freq_df <- as.data.frame(table(H7_clusters)); names(freq_df) <- c("Cluster size", "Frequency")

# Create a table for the HTML document
knitr::kable(freq_df, caption = "Frequencies of MERS Clusters")

# Define likelihood function
lik_function <- function(param) {
  if (any(param <= 0)) return(-Inf) # Ensure positive parameters
  
  # Extract values of R and k
  r_val <- as.numeric(param[1])
  k_val <- as.numeric(param[2])
  
  # Define likelihood
  log_likelihood <- likelihood(
    chains = H7_clusters,
    statistic = "size",
    offspring_dist = rnbinom,
    size = k_val,
    mu = r_val
  )
  
  # Assume non-informative priors for R and k
  log_prior <- 0 # But could add informative priors here if required
  
  # Return log-posterior (log-likelihood + log-prior)
  return(log_likelihood + log_prior)
}

# Define number of MCMC iterations
n_iter <- 1e4

# Define 'burn in' period for fitting, to be discarded
n_burn <- 1e3

# Initial guess for c(R,k):
init_param <- c(R=0.1, k=1)

# Run MCMC to estimate parameters
result_mcmcpack <- MCMCmetrop1R(lik_function, 
                                theta.init = init_param, 
                                burnin = n_burn, 
                                mcmc = n_iter, 
                                thin = 1)

# Calculate effective sample size (i.e. measure of MCMC mixing)
ess_mcmcpack <- effectiveSize(result_mcmcpack)

# Plot posterior estimates
plot(result_mcmcpack)

# Define helper function to calculate median and 95% credible interval from data.frame of MCMC samples
get_param <- function(x){
  apply(x,2,function(y){val = signif(quantile(y,c(0.5,0.025,0.975)),3);
  val_text <- paste0(val[1]," (95%: CrI: ",val[2],"-",val[3],")")})
}

# Get posterior median and 95% CrI
posterior_estimates <- get_param(result_mcmcpack)

# Compile table
results_table <- data.frame(
  Package = "MCMCpack",
  Posterior_R = posterior_estimates[1],
  Posterior_k = posterior_estimates[2],
  ESS_R = ess_mcmcpack[1],
  ESS_k = ess_mcmcpack[2]
)

# Output the table with kable
knitr::kable(results_table, caption = "MCMC Comparison Table", align = 'c')


# Convert MCMC samples to a dataframe
posterior_df <- as.data.frame(result_mcmcpack)
colnames(posterior_df) <- c("R0", "k")  # Rename for clarity

# Calculate median and 95% CrI
R0_median <- median(posterior_df$R0)
R0_CrI <- quantile(posterior_df$R0, c(0.025, 0.975))

# Create density plot with ggplot2
R_1 <- ggplot(posterior_df, aes(x = R0)) +
  # Posterior density
  geom_density(fill = "white", alpha = 0.6, color = "black") +
  
  # Prior density: using exponential(1) as an example
  stat_function(fun = dexp, args = list(rate = 10), 
                aes(color = "Prior"), linetype = "dashed", size = 1.2) +
  
  # Posterior median
  geom_vline(xintercept = R0_median, 
             color = "black", linetype = "dotted", linewidth = 1.2) +
  
  # Labels and theme
  labs(
    title = "",
    x = expression(R),
    y = "Density",
    color = ""
  ) +
  theme_classic() +
  scale_x_continuous(
    limits = c(0, 0.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) + 
  #scale_color_manual(values = c("" = "red")) +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )

ggsave(
  filename = file.path("plots", "H7N7_R.png"), 
  plot = R_1, 
  width = 8, 
  height = 6, 
  dpi = 300
)
