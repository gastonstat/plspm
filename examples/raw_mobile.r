
## sourcing code
library(diagram)
library(turner)
setwd("/Users/Gaston/Documents/R_project_plspm/plspm2/R")
files = system('ls', intern=TRUE)
for (i in 1:length(files))
  source(files[i])
# =======================================================

# another dataset when sourcing files
setwd("/Users/Gaston/Documents/R_project_plspm/datasets")
mobile = read.csv('mobile.csv', row.names=1)


# path matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
COM = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COMP, LOY)

# plot diagram of path matrix
innerplot(sat_path)

# blocks of outer model
sat_blocks = list(1:5, 6:8, 9:15, 16:17, 18:20, 21, 22:24)

# vector of modes (reflective indicators)
sat_mod = rep("A", 7)

# apply plspm
satpls = plspm(mobile, sat_path, sat_blocks, modes = sat_mod, 
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
