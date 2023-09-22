# PROTOCOL FOR: 
# Modeling effects of host behavioral defenses and parasite Allee effects on ectoparasite aggregation 

_Last updated: September 21, 2023_

## Authors: 

* Anna Willoughby

### Summary: 

Individual hosts serve as island habitats for ectoparasites like fleas that exploit their hosts for food and habitat, and between which parasites must disperse. Hosts defend themselves by limiting exposure to ectoparasites or reducing parasite loads via grooming. Despite the patchiness of host habitats and ubiquitousness of host behavioral defenses, ectoparasites are pervasive, reaching high infection prevalence across thousands of host species. The distributions of ectoparasite loads within populations, however, are rarely uniform and tend to follow the “80:20” rule: 80% of individuals are never infested or have low burdens, whereas 20% of individuals host high burdens and are responsible for the bulk of transmission. Processes that underlie the 80:20 rule remain elusive. While the variation in parasite burdens is often attributed to host trait heterogenity or age structure, there may be more fundamental processes at play. A majority of hosts being negligbly colonized by parasites may be due to parasite population Allee effect, where parasites fail to establish due to required density-dependent benefits. Alternatively, hosts maintaining low parasite populations could indicate effective regulation through behavioral defenses. Since ectoparasite populations on their hosts are small, we utilize a discrete, stochastic modeling approach to simulate trajectories of ectoparasite populations on their hosts. To investigate how parasite demographic stochasticity through the Allee effect and host density and defenses change ectoparasite burden distributions, we build a birth-death-immigration model. We will explore these factors under various scenarious of density-dependent changes to parasite birth rate (the Allee effect) or parasite carrying capacities and immigration rates. We will measure outcomes from each simulation from both the host and parasite perspective. 

### Research questions:
 1) Can inclusion of the Allee effect in ectoparasite population dynamics explain 80:20 distribution of parasite loads across hosts? 
 2) Can changes in host density explain 80:20 distribution of parasite loads across hosts? 
 3) How do these two drivers of population dynamics interact?

### Study design:

We will examine how parasite Allee effect and host density impact parasite population invasion success and stationary distribution sizes through simulated birth-death models. 

Table 1. Model components
| Parameter  | BD.Model.0 Value | BDI.Model.1 Value | BDI.Model.2 Value | BDI.Model.3 Value
|:-:|:-:|:-:|:-:|:-:|
| birth rate, $\lambda$  | 0.2 |  0.2 | 0.2 |  0.2 | 
| death rate, $\mu$ | 0.1 + logistic growth | 0.1 + logistic growth | 0.1 + logistic growth  | 0.1 + logistic growth | 
| Allee effect, $\theta$  | 0 | $\theta \in 1:K$ | 0 |  $\theta \in 1:K$ |
| immigration rate, $\iota$ | 0 |  0  | $\iota \in \iota_{min}:\iota_{max}$ ) |   $\iota \in \iota_{min}:\iota_{max}$ |
| carrying capacity, $K$ | K_{min} | K_{min} | $K \in K_{min}:K_{max}$ |  $K \in K_{min}:K_{max}$  |
| grooming rate, $\zeta$ | 0 | 0 | $\zeta \in \zeta_{min}:\zeta_{max}$ | $\zeta \in \zeta_{min}:\zeta_{max}$ |

Each model represents the population of parasites on one host. Each simulation will be run for 100 time points. 

Let $N(t)$ denote the random variable for the total parasite population size $N$ at time $t$. 

$$
\begin{aligned}
N_{t+1}=N_{t} + (\lambda {N_t}/({\theta + N_t}))N_t + \iota {N_t} - \mu {N_t}
\end{aligned}
$$

Apply this function with parameters from Table for each model. 
- Model 0. Birth-Death 
- Model 1. Birth-Death-Allee
- Model 2. Birth-Death-Immigration
- Model 3. Birth-Death-Alleee-Immigration


### Analysis 

_Host Population Infection Metrics_

* Calculate traditional infection metrics from host perspective (prevalence, mean abundance, mean intensity)
* How these values vary with Allee effect and host density changes

_Parasite Population Metrics_

* Calculate infection metrics from parasite perspective (total population size, crowding)
* How these values vary with Allee effect and host density changes
* Calculate stability? 

_Parasite Infrapopulation Distribution_

* From the stationary distribution, calculate aggregation metrics (mean intensity:variance intensity, $k$, crowding, skew, Pareto fraction)
* How these values relate to each other and vary with Allee effect and host density changes 

### Checklist: 

* Confirm model math 
* Decide on parameter values
* Write code for model runs
* Write code for outcomes
* Run model + calculate outcomes
* Make figures and perform statistical analysis
* Explore metapop? 

### Permits

### Important background papers: 
- Dennis B. Allee effects in stochastic populations. Oikos. 2002 Mar;96(3):389-401.
- Drake JM, Lodge DM. Allee effects, propagule pressure and the probability of establishment: risk analysis for biological invasions. Biological Invasions. 2006 Mar;8:365-75.

### CHANGE-LOG:
- (September 21,2023) Restructured protocol for population growth model instead of compartmental model. 
- (January 11,2023) Replaced "Food provisioning, behavioral defense, and host-ectoparasite dynamics" title added abstract
