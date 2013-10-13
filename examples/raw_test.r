
## example with customer satisfaction analysis
## group comparison based on the segmentation variable "gender"

library(plspm2)

# load dataset satisfaction
data(satisfaction, package="plspm2")

# path matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0) 
LOY = c(1,0,0,0,1,0)
path_matrix = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# blocks of outer model
blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

# vector of modes (reflective indicators)
modes = rep("A", 6)

Data = satisfaction