library(ggplot2)

fn <- "../fields/edge-d40-idensp1-edgeInfo.bz2"
expected <- scan(fn, n=1, skip=0)
data <- read.table(fn, skip=1)
names(data) <- c("x", "rho", "gradrho", "gradcolor")
edgerho <- ggplot(data, aes(x = x, y = (gradrho / rho) / expected)) + geom_point(size = .1) + xlim(5, 7) + xlab("$x$ coordinate") + ylab("$\\textsc{EdgeDetect}()$")
ggsave("edgerho-ip1.svg", plot = edgerho, width = 3, height = 2)
edgecol <- ggplot(data, aes(x = x, y = gradcolor / expected)) + xlim(5, 7) + geom_point(size = .1) + xlab("$x$ coordinate") + ylab("$\\textsc{EdgeDetectColor}(c)$")
ggsave("edgecol-ip1.svg", plot = edgecol, width = 3, height = 2)
