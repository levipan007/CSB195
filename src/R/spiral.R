# src/R/spiral.R
# Logarithmic spiral demo; writes plot to ./out/spiral.png

# Ensure output folder exists
if (!dir.exists("out")) dir.create("out", recursive = TRUE)

# Parameters for a logarithmic spiral r = a * e^(b*theta)
a <- 0.5
b <- 0.20
theta <- seq(0, 6*pi, length.out = 2000)
r <- a * exp(b * theta)

# Convert polar to Cartesian
x <- r * cos(theta)
y <- r * sin(theta)

# Plot and save
png("out/spiral.png", width = 1000, height = 1000, res = 150)
par(mai = c(0,0,0,0))
plot(x, y, type = "l", lwd = 3, axes = FALSE, xlab = "", ylab = "", asp = 1)
dev.off()

