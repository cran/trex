
## Test script for optimalDesign function in trex package;

require(trex)


######### examples presented in Table 4, Schaid and Sinnwell (Human Genetics)

example1 <- optimalDesign(probCarrier=.01, or=5, alpha=.05, power=.9)

example2 <- optimalDesign(probCarrier=.05, or=5, alpha=.05, power=.9)

example1

example2

