plspm
============================

R package **plspm** exclusively dedicated to Partial Least Squares Path Modeling (PLS-PM) analysis. Versions later than 0.3.0 only contain functions related to Partial Least Squares Path Modeling. Other methods such as nipasl, pls regression, and pls canonical analysis are now in the brother package **plsdepot**.

## Requirements and Installation

To install the stable version of **plspm** from CRAN, run in your R console:
```r
install.packages("plspm")
```

To install the latest development version of **plspm** from github (using the package "devtools"), run in your R console:
```
# install.packages("devtools") 
library(devtools)
install_github('plspm', username='gastonstat')
```

## Example Usage with a Customer Satisfaction Model 
```
library(plspm)

# load dataset satisfaction
data(satisfaction)

# define inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0) 
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# define outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

# vector of modes
sat_mod = rep("A", 6)   ## reflective indicators

# apply plspm with bootstrap validation
satpls = plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
  
# summary of results
summary(satpls)

# plot inner model results
plot(satpls, what="inner")

# plot outer model loadings
plot(satpls, what="loadings")

# plot outer model weights
plot(satpls, what="weights")
```

More info at [www.gastonsanchez.com/plspm](http://www.gastonsanchez.com/plspm)

Links
-----
[plspm package github](http://github.com/gastonstat/plspm)

[plspm slides](http://www.gastonsanchez.com/plspm)

[PLS Modeling stuff](http://www.plsmodeling.com)


Author Contact
--------------
Gaston Sanchez (gaston.stat at gmail.com)