######################################
##        Network Causal Tree       ##
######################################

## Graphs Causal Rules

rm(list=ls())
setwd("G:/Il mio Drive/Research/Networks/Results_simulations/v4")

par(xpd=TRUE)

################################################################
##     Plot of Two Rules with Different Number of Clusters    ##
################################################################

# 2 RULES

# COMPOSITE TREE
# 10 clusters
composite <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_composite_1000.csv")
singular_main <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_main_1000.csv")
singular_spil <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_spillover_1000.csv")

# Plot
layout(matrix(c(1,2,3,4,4,4), ncol=3, nrow=2, byrow=TRUE), heights=c(5, 2))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules/2,
     main = "10 Clusters",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 3,
     ylim=c(0,2))
lines(singular_main$correct_rules_main/2, type = "o", col = "blue", lty=2, lwd = 3)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=2, lwd = 3)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# 20 clusters
composite <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_composite_2000.csv")
singular_main <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_main_2000.csv")
singular_spil <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_spillover_2000.csv")

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules/2,
     main = "20 Clusters",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 3,
     ylim=c(0,2))
lines(singular_main$correct_rules_main/2, type = "o", col = "blue", lty=2, lwd = 3)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=2, lwd = 3)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# 30 clusters
composite <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_composite_3000.csv")
singular_main <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_main_3000.csv")
singular_spil <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_spillover_3000.csv")

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules/2,
     main = "30 Clusters",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 3,
     ylim=c(0,2))
lines(singular_main$correct_rules_main/2, type = "o", col = "blue", lty=2, lwd = 3)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=2, lwd = 3)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

op <- par(cex = 1.25)
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", 
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"),
       col=c("red", "blue", "forestgreen"), lty=c(1,2,2), cex=0.7, lwd = 3)
rm(op)

###########################################
##   Additional Trees (not in the paper) ##
###########################################

# TREATMENT EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# SPILLOVER EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(2*singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(2*singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(2*singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))


################################################################
##    Plot of Four Rules with Different Number of Clusters    ##
################################################################

rm(list=ls())
setwd("G:/Il mio Drive/Research/Networks/Results_simulations/v4")
composite <- read.csv("./main_text/4_causal_rules/two_main_spillover_four_effects_composite_3000.csv")
singular_main <- read.csv("./main_text/4_causal_rules/two_main_spillover_four_effects_singular_main_3000.csv")
singular_spil <- read.csv("./main_text/4_causal_rules/two_main_spillover_four_effects_singular_spil_3000.csv")

# Plot
layout(matrix(c(1,2,3,4,4,4), ncol=3, nrow=2, byrow=TRUE), heights=c(5, 2))

# COMPOSITE TREE
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Composite Tree",
     sub = "Four Leaves",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_main + singular_spil$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil + singular_main$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# TREATMENT EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_composite,
     main = "Treatment Effects Tree",
     sub = "Two Leaves",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# SPILLOVER EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_composite,
     main = "Spillover Effects Tree",
     sub = "Two Leaves",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))


op <- par(cex = 1.25)
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", 
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"),
       col=c("red", "blue", "forestgreen"), lty=c(1,2,2), cex=0.7, lwd = 3)
rm(op)


####################################
##      Plots for the Appendix    ##
####################################

# CORRELATED REGRESSORS

rm(list=ls())
setwd("G:/Il mio Drive/Research/Networks/Results_simulations/v4")
composite <- read.csv("./appendix/correlation/two_main_spillover_effects_composite_corr_0.25_3000.csv")
singular_main <- read.csv("./appendix/correlation/two_main_spillover_effects_singular_main_corr_0.25_3000.csv")
singular_spil <- read.csv("./appendix/correlation/two_main_spillover_effects_singular_spillover_corr_0.25_3000.csv")

# Plot
layout(matrix(c(1,2,3,3), ncol=2, nrow=2, byrow=TRUE), heights=c(5, 2))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correlation 0.25",
     sub = "Four Leaves",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)


composite <- read.csv("./appendix/correlation/two_main_spillover_effects_composite_corr_0.50_3000.csv")
singular_main <- read.csv("./appendix/correlation/two_main_spillover_effects_singular_main_corr_0.50_3000.csv")
singular_spil <- read.csv("./appendix/correlation/two_main_spillover_effects_singular_spillover_corr_0.50_3000.csv")

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correlation 0.50",
     sub = "Four Leaves",
     xlab = "Effect Size",
     ylab = "Number of Correctly Identified Leaves", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

op <- par(cex = 1.25)
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", 
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"),
       col=c("red", "blue", "forestgreen"), lty=c(1,2,2), cex=0.7, lwd = 3)
rm(op)
