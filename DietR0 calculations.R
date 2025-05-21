############################################################################
############### Sauer et al. Diet-SIR script           #####################
############### R code adapted by Erin L. Sauer from   #####################
############### Wolfram code developed by Jessica Hite #####################
############################################################################

# Load required library
library(ggplot2)

# Define parameters for Protein and Lipid diets
params_protein <- list(
  gamma1 = 0.048, gamma2 = 0.081,
  mu1 = 0.003, mu2 = 0.0015,
  nu1 = 0.25, nu2 = 0.125,
  delta = 0.01
)

params_lipid <- list(
  gamma1 = 0.033, gamma2 = 0.567,
  mu1 = 0.009, mu2 = 0.0045,
  nu1 = 0.25, nu2 = 0.125,
  delta = 0.01
)

# Function to compute R0 using the Next-Generation Matrix method
compute_R0 <- function(beta, S, params) {
  # Unpack parameters
  with(params, {
    V <- matrix(c(
      gamma1 + delta + mu1 + nu1, 0,
      0, gamma2 + delta + mu2 + nu2
    ), nrow = 2, byrow = TRUE)
    
    F <- matrix(c(
      beta * S, beta * S,
      0, 0
    ), nrow = 2, byrow = TRUE)
    
    K <- F %*% solve(V)
    eigenvalues <- eigen(K)$values
    Re(max(Mod(eigenvalues)))  # Return dominant eigenvalue (R0)
  })
}

# Compute R0 for fixed beta values
R0_protein <- compute_R0(0.00275, 100, params_protein)
R0_lipid   <- compute_R0(0.0055, 100, params_lipid)

# Print Results
R0_protein
R0_lipid

# Define range of beta values, assuming beta.protein > beta.lipid
beta_range <- seq(0.0028, 0.006875, by = 0.001)

# Compute R0 over beta range
R0_lipid_values   <- sapply(beta_range, function(b) compute_R0(b, 100, params_lipid))

# Create data frame for plotting
plot_data <- data.frame(
  beta = rep(beta_range, times = 2),
  R0 = c(R0_lipid_values),
  Diet = factor(rep(c("Variable β(Lipid)"), each = length(beta_range)))
)

protein.ref <- data.frame(
  beta = 0.00275,
  R0 = 0.8842444,
  Diet = factor("β(Protein)")
)

lipid.ref <- data.frame(
  beta = 0.0055,
  R0 = 1.821192,
  Diet = factor("β(Lipid)")
)

# Append rows
plot_data <- rbind(plot_data, protein.ref, lipid.ref)
str(plot_data)

# Plot 
ggplot(plot_data, aes(x = beta, y = R0, shape = Diet)) +
  geom_line(size = 1.2) +
  geom_point(size=5)+
  labs(
    x = "Transmission Rate (β)",
    y = "R0"
  ) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.35,0.87),
        legend.text= element_text(size=20),
        legend.title= element_blank()
  )+
  scale_colour_grey()+
  theme(plot.margin = margin(.5,1.5,.5,.5, "cm"))+
  scale_shape_manual(values = c("Variable β(Lipid)" = NA, "β(Protein)" = 15, "β(Lipid)" = 17))
