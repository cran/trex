##
## Demo file for trex package, a truncated exact test for
## 2-stage case-control designs for genetic studies of rare variants
##


## load trex package like this:
require(trex)

args(trex) 

## The below examples are in the help file, runnable by:
##    R> example(trex)



## test like this:

tab11 <- cbind(c(2,1,8), c(46,46,98))
out11 <- trex(tab11, threshold=2)

print.trex(out11)



tab10 <- cbind(c(8, 0, 2), c(92, 0, 98))
out10 <- trex(tab10, threshold=2)

print(out10)


# try a different threshold:
out10.t4 <- trex(tab10, threshold=4)

print(out10.t4)



# this should get an error before the C function
# out11.fail <- trex(tab11, threshold=3)
# Error in trex(tab11, threshold = 3) : Stage 1 does not meet threshold.

#tab0 <- cbind(c(8,2,0), c(42, 48, 0))
#out0.fail <- trex(tab0, threshold=2)

