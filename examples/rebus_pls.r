## Example of REBUS PLS with simulated data
library(plspm)

# load simdata
data(simdata)
  
# Calculate global plspm
sim_path = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
dimnames(sim_path) = list(c("Price", "Quality", "Satisfaction"),
                           c("Price", "Quality", "Satisfaction"))
sim_blocks = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
sim_mod = c("A", "A", "A")  # reflective indicators
sim_global = plspm(simdata, sim_path, 
                   sim_blocks, modes=sim_mod)
sim_global

## Then compute cluster analysis on residuals of global model
sim_clus = res.clus(sim_global)

## To complete REBUS, run iterative algorithm
rebus_sim = it.reb(sim_global, sim_clus, nk=2, 
                  stop.crit=0.005, iter.max=100)

## You can also compute complete outputs 
## for local models by running:
local_rebus = local.models(sim_global, rebus_sim)

# Display plspm summary for first local model
summary(local_rebus$loc.model.1)

# apply rebus.test
sim_permu = rebus.test(sim_global, rebus_sim)
  
# inspect sim.rebus
sim_permu
sim_permu$test_1_2
  
# or equivalently
sim_permu[[1]]
