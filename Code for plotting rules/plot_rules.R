## Graphs Causal Rules

rm(list=ls())
setwd("G:\\Il mio Drive\\Research\\Networks\\Draft Costanza\\networks\\Tables")

par(xpd=TRUE)

######################################
##     Network-based Causal Tree    ##
######################################

# 2 RULES

composite_effect_size <- read.csv("netcausaltree_main_effects.csv")
singular_effect_size <- read.csv("one_netcausaltree_main_effect.csv")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(composite_effect_size$correct_rules,
     main = "Correctly Detected Rules",
     sub = "Two True Rules",
     xlab = "Effect Size",
     ylab = "Number of Correct Rules", 
     xaxt='n',
     type = "o",
     col = "red",
     lwd = 2,
     ylim=c(0,2))
lines(singular_effect_size$correct_rules, type = "o", col = "blue", lty=2, lwd = 2)
legend("topright", inset=c(-0.15,0.8),
       legend=c("NCT (composite)", "NCT (singular)"), 
       col=c("red", "blue"), lty=1:2, cex=0.7, lwd = 2)
axis(1, at=1:length(results_effect_size$correct_rules), labels=c(seq(0.1,10.1,1)))
