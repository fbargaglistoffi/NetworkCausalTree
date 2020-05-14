######################################
##        Network Causal Tree       ##
######################################

## Graphs Causal Rules

rm(list=ls())
setwd("G:/Il mio Drive/Research/Networks/Results_simulations/v4")

par(xpd=TRUE)

##################
##     1000     ##
##################

# 2 RULES

# COMPOSITE TREE

composite <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_composite_1000.csv")
singular_main <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_main_1000.csv")
singular_spil <- read.csv("./main_text/2_causal_rules/two_main_spillover_effects_singular_spillover_1000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correctly Detected Rules (All Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil*2, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

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


# 4 RULES

composite <- read.csv("two_main_spillover_four_effects_composite_1500.csv")
singular_main <- read.csv("two_main_spillover_four_effects_singular_main_1500.csv")
singular_spil <- read.csv("two_main_spillover_four_effects_singular_spil_1500.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correctly Detected Rules (All Effects)",
     sub = "Four True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# TREATMENT EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# SPILLOVER EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

##################
##     3000     ##
##################

rm(list=ls())
setwd("G:\\Il mio Drive\\Research\\Networks\\Results_simulations\\v3\\3000")

par(xpd=TRUE)

######################################
##     Network-based Causal Tree    ##
######################################

# 2 RULES

singular_main <- read.csv("two_main_spillover_effects_singular_main_3000.csv")
singular_spil <- read.csv("two_main_spillover_effects_singular_spillover_3000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules (3000)",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

# 3 RULES

composite <- read.csv("two_main_spillover_three_effects_composite_3000.csv")
singular_main <- read.csv("two_main_spillover_three_effects_singular_main_3000.csv")
singular_spil <- read.csv("two_main_spillover_three_effects_singular_spillover_3000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correctly Detected Rules (All Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,3))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

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
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# SPILLOVER EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))


# 4 RULES

composite <- read.csv("two_main_spillover_four_effects_composite_3000.csv")
singular_main <- read.csv("two_main_spillover_four_effects_singular_main_3000.csv")
singular_spil <- read.csv("two_main_spillover_four_effects_singular_spil_3000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite$correct_rules,
     main = "Correctly Detected Rules (All Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,4))
lines(singular_main$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# TREATMENT EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))

# SPILLOVER EFFECT
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (main)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(composite$correct_rules), labels=c(seq(0.1,10.1,1)))


########################
##     Correlated     ##
########################


rm(list=ls())
setwd("G:\\Il mio Drive\\Research\\Networks\\Results_simulations\\v3\\Correlated")

par(xpd=TRUE)

######################################
##     Network-based Causal Tree    ##
######################################

##################
##     1500     ##
##################

# 2 RULES

singular_main <- read.csv("two_main_spillover_effects_singular_main_1500.csv")
singular_spil <- read.csv("two_main_spillover_effects_singular_spillover_1500.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules (Correlated Regressors)",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules (Correlated Regressors)", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

##################
##     3000     ##
##################

# 2 RULES

singular_main <- read.csv("two_main_spillover_effects_singular_main_3000.csv")
singular_spil <- read.csv("two_main_spillover_effects_singular_spillover_3000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules (Correlated Regressors)",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules (Correlated Regressors)", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))


########################
##      Homophily     ##
########################


rm(list=ls())
setwd("G:\\Il mio Drive\\Research\\Networks\\Results_simulations\\v3\\homophily")

par(xpd=TRUE)

######################################
##     Network-based Causal Tree    ##
######################################

##################
##     1500     ##
##################

# 2 RULES

singular_main <- read.csv("two_main_spillover_effects_singular_main_homophily_1500.csv")
singular_spil <- read.csv("two_main_spillover_effects_singular_spillover_homophily_1500.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules (Homophily)",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules (Homophily)", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

##################
##     3000     ##
##################

# 2 RULES

singular_main <- read.csv("two_main_spillover_effects_singular_main_homophily_3000.csv")
singular_spil <- read.csv("two_main_spillover_effects_singular_spillover_homophily_3000.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$correct_rules_composite,
     main = "Correctly Detected Rules (Treatment Effects)",
     sub = "Two True Rules (Homophily)",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_main$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_main$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$correct_rules_composite,
     main = "Correctly Detected Rules (Spillover Effects)",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules (Homophily)", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_spil$correct_rules_main, type = "o", col = "blue", lty=2, lwd = 2)
lines(singular_spil$correct_rules_spil, type = "o", col = "forestgreen", lty=3, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (treatment)", "NCT (spillover)"), 
       col=c("red", "blue", "forestgreen"), lty=1:3, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

###################
## PLOT VARIANCE ##
###################

main <- read.csv("G:\\Il mio Drive\\Research\\Networks\\Results_simulations\\v3\\1500\\two_main_spillover_effects_singular_main_1500.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_main$Rule.1.se.tau.10.00..,
     main = "Standard Errors (Main Effects)",
     #sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Standard Error", 
     xaxt='n',
     type = "o",
     col = "green",
     lwd = 2,
     ylim=c(0,4))
lines(main$Rule.1.se.tau.10.00.., type = "o", col = "orange", lty=2, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("ERGM (homophily)", "Erdos-Renyi"), 
       col=c("green", "orange"), lty=1:2, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

spil <- read.csv("G:\\Il mio Drive\\Research\\Networks\\Results_simulations\\v3\\1500\\two_main_spillover_effects_singular_spillover_1500.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(singular_spil$Rule.2.se.eta.01.00..,
     main = "Standard Errors (Spillover Effects)",
     #sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Standard Error", 
     xaxt='n',
     type = "o",
     col = "green",
     lwd = 2,
     ylim=c(0,2))
lines(spil$Rule.1.se.eta.01.00.., type = "o", col = "orange", lty=2, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("ERGM (homophily)", "Erdos-Renyi"), 
       col=c("green", "orange"), lty=1:2, cex=0.7, lwd = 2)
axis(1, at=1:length(singular_main$correct_rules_composite), labels=c(seq(0.1,10.1,1)))

