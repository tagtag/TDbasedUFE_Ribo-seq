#mixOmics
dim(Z)<- c(dim(Z)[1],18,3) #not scaled
library(mixOmics)

## --------------------------------------------------
# 1. Create a list of data (please insert your own data)
# --------------------------------------------------
# * Replace your_data_1, your_data_2... with your actual data frame names
X <- list(
  mRNA = t(Z[,,1]),  
  trans = t(Z[,,2]),
  protein = t(Z[,,3])   # Remove this line if you only have 2 datasets
)

# --------------------------------------------------
# 2. Create a list for variable selection (names will automatically match X here)
# --------------------------------------------------
list.keepX <- list(
  c(10, 10), # Number of variables to select for DataA
  c(10, 10), # Number of variables to select for DataB
  c(10, 10)  # Number of variables to select for DataC (Remove this line if you only have 2 datasets)
)
# * Forcibly copy the exact same names from X to keepX *
names(list.keepX) <- names(X) 

# --------------------------------------------------
# 3. Create the design matrix (names from X are also used automatically here)
# --------------------------------------------------
design <- matrix(1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# --------------------------------------------------
# 4. Run the model
# --------------------------------------------------
unsupervised.sgcca <- wrapper.sgcca(X = X, 
                                    design = design, 
                                    ncomp = 2, 
                                    keepX = list.keepX)

# If completed without errors, visualize!
plotIndiv(unsupervised.sgcca, ind.names = FALSE, legend = TRUE)


# -------------------------------------------------------------
# 1. Extract the representative values of the "3 omics blocks" for each component
# -------------------------------------------------------------
# 1st component
T1 <- cbind(unsupervised.sgcca3$variates[[1]][, 1], 
            unsupervised.sgcca3$variates[[2]][, 1], 
            unsupervised.sgcca3$variates[[3]][, 1])

# 2nd component
T2 <- cbind(unsupervised.sgcca3$variates[[1]][, 2], 
            unsupervised.sgcca3$variates[[2]][, 2], 
            unsupervised.sgcca3$variates[[3]][, 2])

# 3rd component
T3 <- cbind(unsupervised.sgcca3$variates[[1]][, 3], 
            unsupervised.sgcca3$variates[[2]][, 3], 
            unsupervised.sgcca3$variates[[3]][, 3])

# -------------------------------------------------------------
# 2. Calculate the connection weight vectors between blocks (length 3) using Singular Value Decomposition (SVD)
# -------------------------------------------------------------
# Calculate the weights (coefficients) of the 3 blocks when forming each component
u1k <- svd(T1)$v[, 1]
u2k <- svd(T2)$v[, 1]
u3k <- svd(T3)$v[, 1]

# -------------------------------------------------------------
# 3. Plot in the exact same format as the image
# -------------------------------------------------------------
# Split the screen vertically into 3 sections
pdf(file="mixOmics.pdf",width=5,height=10)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

# k on the x-axis is 1, 2, 3 (representing omics blocks 1, 2, and 3, respectively)
k_axis <- 1:3


# --- Plot for u1k ---
plot(k_axis, u1k, type = "h", lwd = 2, 
     ylab = "uk1", xlab = "k", 
     xlim = c(1, 3),cex.axis=2,cex.lab=1.5,ylim=c(0,1)) # Set X-axis from 1.0 to 3.0
abline(h = 0, col = "red", lty = 2)

# --- Plot for u2k ---
plot(k_axis, u2k, type = "h", lwd = 2, 
     ylab = "uk1", xlab = "k", 
     xlim = c(1, 3),cex.axis=2,cex.lab=1.5,ylim=c(0,1))
abline(h = 0, col = "red", lty = 2)

# --- Plot for u3k ---
plot(k_axis, u3k, type = "h", lwd = 2, 
     ylab = "uk3", xlab = "k", 
     xlim = c(1, 3),cex.axis=2,cex.lab=1.5,ylim=c(0,1))
abline(h = 0, col = "red", lty = 2)
# 描画設定を元に戻す
par(mfrow = c(1, 1))
dev.off()
