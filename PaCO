# Load required libraries
library(ape)
library(paco)

# Input: Host-Parasite association data
data <- data.frame(
  Host = c("fish2", "fish3", "fish4", "fish6", "fish8", "fish9", "fish10", 
           "fish11", "fish13", "fish14", "fish15", "fish16", "fish18", 
           "fish19", "fish20", "fish21", "fish23", "fish24", "fish25", 
           "fish26", "fish27", "fish28", "fish30", "fish31", "fish32", 
           "fish33", "fish34", "fish35", "fish36", "fish38", "fish39", 
           "fish43"),
  Parasite = c("Co2", "Co3", "Co4", "Co6", "Co8", "Co9", "Co10", "Co11", 
               "Co13", "Co14", "Co15", "Co16", "Co18", "Co19", "Co20", 
               "Co21", "Co23", "Co24", "Co25", "Co26", "Co27", "Co28", 
               "Co30", "Co31", "Co32", "Co33", "Co34", "Co35", "Co36", 
               "Co38", "Co39", "Co43")
)

# Create a list of unique hosts and parasites
hosts <- unique(data$Host)
parasites <- unique(data$Parasite)

# Initialize a binary association matrix
association_matrix <- matrix(0, nrow = length(hosts), ncol = length(parasites),
                             dimnames = list(hosts, parasites))

# Fill the matrix based on the associations
for (i in 1:nrow(data)) {
  association_matrix[data$Host[i], data$Parasite[i]] <- 1
}

# View the association matrix
print(association_matrix)
# Compute cophenetic distance matrices
D_host <- cophenetic(host_tree)
D_parasite <- cophenetic(parasite_tree)

# Prepare PACo data
paco_data <- prepare_paco_data(D_host, D_parasite, association_matrix)

# Continue with the PACo workflow...
paco_data <- add_pcoord(paco_data)
paco_result <- PACo(paco_data, nperm = 999, method = "r0")
summary(paco_result)
str(paco_result)
# Extract first two axes for host and parasite
host_pcoord <- paco_result$H_PCo[, 1:2]  # First two principal coordinates of hosts
parasite_pcoord <- paco_result$P_PCo[, 1:2]  # First two principal coordinates of parasites
# Plotting
plot(host_pcoord[, 1], host_pcoord[, 2], col = "blue", pch = 19, xlab = "PCo1", ylab = "PCo2", 
     main = "PACo: Host vs Parasite Principal Coordinates", xlim = range(c(host_pcoord[, 1], parasite_pcoord[, 1])),
     ylim = range(c(host_pcoord[, 2], parasite_pcoord[, 2])))
points(parasite_pcoord[, 1], parasite_pcoord[, 2], col = "red", pch = 19)

# Adding legend
legend("topright", legend = c("Host", "Parasite"), col = c("blue", "red"), pch = 19)
# Extract Procrustes coordinates for host and parasite
host_procrustes <- paco_result$proc$X[, 1:2]
parasite_procrustes <- paco_result$proc$Yrot[, 1:2]

# Plot Procrustes transformation
plot(host_procrustes[, 1], host_procrustes[, 2], col = "blue", pch = 19, xlab = "PCo1", ylab = "PCo2", 
     main = "PACo Procrustes: Host vs Parasite", xlim = range(c(host_procrustes[, 1], parasite_procrustes[, 1])),
     ylim = range(c(host_procrustes[, 2], parasite_procrustes[, 2])))
points(parasite_procrustes[, 1], parasite_procrustes[, 2], col = "red", pch = 19)

# Adding legend
legend("topright", legend = c("Host", "Parasite"), col = c("blue", "red"), pch = 19)
# Plot the PCoA or Procrustes results
plot(paco_result)

# Add labels to the points (for both hosts and parasites)
# Assuming the data is stored in the variables 'H' for hosts and 'P' for parasites

# Plot host points (blue) with labels
text(paco_result$H[, 1], paco_result$H[, 2], labels = rownames(paco_result$H), col = "blue", cex = 0.7)

# Plot parasite points (red) with labels
text(paco_result$P[, 1], paco_result$P[, 2], labels = rownames(paco_result$P), col = "red", cex = 0.7)


# Extract the host and parasite coordinates from the PACo result
host_coords <- paco_result$H_PCo  # Coordinates for hosts
parasite_coords <- paco_result$P_PCo  # Coordinates for parasites

# Create an empty plot (this will create the axes and space for the plot)
plot(host_coords[, 1], host_coords[, 2], 
     xlab = "PCoA Axis 1", ylab = "PCoA Axis 2", 
     pch = 19, col = "blue", xlim = range(c(host_coords[, 1], parasite_coords[, 1])), 
     ylim = range(c(host_coords[, 2], parasite_coords[, 2])))

# Add parasite points to the plot (in red)
points(parasite_coords[, 1], parasite_coords[, 2], pch = 19, col = "red")

# Add labels to the host points (blue)
text(host_coords[, 1], host_coords[, 2], labels = rownames(host_coords), col = "blue", cex = 0.7)

# Add labels to the parasite points (red)
text(parasite_coords[, 1], parasite_coords[, 2], labels = rownames(parasite_coords), col = "red", cex = 0.7)
