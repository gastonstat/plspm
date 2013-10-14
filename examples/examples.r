
## example with customer satisfaction analysis
## group comparison based on the segmentation variable "gender"

library(plspm)

# load dataset satisfaction
data(satisfaction)

# path matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0) 
LOY = c(1,0,0,0,1,0)
sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# plot diagram of path matrix
innerplot(sat_path)

# blocks of outer model
sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)


# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply plspm
satpls = plspm(satisfaction, sat_path, sat_blocks, modes = sat_mod, 
                scaled = FALSE)

# plot diagram of the inner model
innerplot(satpls)

# plot loadings
outerplot(satpls, what="loadings")

# plot outer weights
outerplot(satpls, what="weights")
         
# rescale LVs
rescores = rescale(satpls)
head(rescores)

# permutation test with 100 permutations
group_perm = plspm.groups(satpls, satisfaction$gender, 
                            method="permutation", reps=100)
group_perm

# permutation test with 100 permutations
group_boot = plspm.groups(satpls, satisfaction$gender, 
                          method="bootstrap", reps=100)
group_boot
