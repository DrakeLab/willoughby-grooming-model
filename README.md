# PROTOCOL FOR: 
# Modeling effects of host behavioral defenses and parasite Allee effects on ectoparasite aggregation 

_Last updated: March 3, 2025_

## Authors: 

* Anna Willoughby

### Summary: 

Individual hosts serve as island habitats for ectoparasites like fleas that exploit their hosts for food and habitat, and between which parasites must disperse. Hosts defend themselves by limiting exposure to ectoparasites or reducing parasite loads by grooming. Despite the patchiness of host habitats and ubiquitousness of host behavioral defenses, ectoparasites are pervasive, reaching high infection prevalences across thousands of host species. The distributions of ectoparasite loads within populations, however, are rarely uniform and tend to follow the “80:20” rule: 80% of individuals are never infested or have low burdens, whereas 20% of individuals host high loads and are responsible for the bulk of transmission. Processes that underlie the 80:20 rule remain elusive. While the variation in parasite burdens is often attributed to host trait heterogeneity or age structure, more fundamental parasite population processes may be at play. A majority of hosts being negligibly colonized by parasites may be due to the Allee effect operating within the parasite population. Parasites fail to establish on hosts at low population levels due to density-dependent benefits of crowding (e.g., immune evasion). Alternatively, hosts maintaining low parasite populations could indicate effective regulation through behavioral defenses. Grooming frequently is a social activity among animals, with frequency increasing with host density. Therefore, while increases in host density could lead to increasing opportunity for parasites to immigrate onto a new host, there also may be an increased risk of mortality from allogrooming. As ectoparasite populations on their hosts are small (commonly <10 individuals), we utilize a discrete, stochastic modeling approach to simulate trajectories of ectoparasite populations on their hosts. To investigate how parasite demographic stochasticity through the Allee effect and host density-dependent defenses change ectoparasite burden distributions, we build a birth-death-immigration model. We will explore these factors under various scenarios of parasite birth rate (the Allee effect), immigration, and host-induced mortality rates. We will measure outcomes from each simulation from both the host and parasite perspective. 

### Research questions:
 1) Can the inclusion of the Allee effect in ectoparasite population dynamics explain 80:20 distribution of parasite loads across hosts? 
 2) Can host density changes to parasite immigration and mortality explain 80:20 distribution of parasite loads across hosts? 
 3) How do these two drivers of population dynamics interact?

### Study design:

We will examine how the parasite Allee effect and host density impact parasite population invasion success and stationary distribution sizes through simulated birth-death continuous time Markov chain (CTMC) models. 

Table 1. Model components
| Parameter  | BD.Model.0 Value | BDI.Model.1 Value | BDI.Model.2 Value | BDI.Model.3 Value
|:-:|:-:|:-:|:-:|:-:|
| birth rate, $\lambda$  | 0.2 |  0.2 | 0.2 |  0.2 | 
| death rate, $\mu$ | $d_{x}$ | $d_{x}$ | $d_{x}$   | $d_{x}$  | 
| Allee effect, $\theta$  | 0 | $\theta \in 1:K$ | 0 |  $\theta \in 1:K$ |
| immigration rate, $\iota$ | 0 |  0  | $\iota \in \iota_{min}:\iota_{max}$ ) |   $\iota \in \iota_{min}:\iota_{max}$ |
| carrying capacity, $K$ | $K$ | $K$ | $K$ |  $K$  |

Table 2. Different mortality based on grooming scenarios 
| Grooming Scenario |        Function        | Notes
|:-----------------:|:----------------------:|:------------------------|
| none | $d$(x) = 0 | only flea density dependence (need to find citations)
| programmed  | $d$(x) = $d_{0}$ | static grooming rate
| stimulus-driven |  $d$(x) = $d_{0}$ + $d_{1}x$ | positive linear relationship between grooming time and ectoparasite burden; may explore saturating function |
| costly | $d$(x) =  | left skewed gaussian distribution |


Each model represents the population of parasites on one host. Each simulation will be run until a quasi-stationary distribution is reached.  

We describe a stochastic population undergoing a continuous-time, discrete-state Markov Chain process. 

The birth rate will incorporate the Allee effect, $\theta$, and the mortality rate will incorporate density-dependence, . 

$$
\begin{aligned}
{\lambda_{x}} = b_{0} x \frac{x}{\theta + x}  \text{      and      } \mu_{x} = (d_{0} + d_{1}x)x 
\end{aligned}
$$

Together with immigration, $\iota$, which is not impacted by flea density, the deterministic rate equation is: 

$$
\begin{aligned}
\frac{dx}{dt} = b_{0} x \frac{x}{\theta + x}(t) - (d_{0} + d_{1}x)x(t) + \iota(t)
\end{aligned}
$$

Apply this function with parameters from Table 1 for each model. Relationships of $\iota$, $K$, $\zeta$ will have linear or non-linear relationships with host density (as a proxy for provisioning) that need to be determined. 
- Basic Model (0). Birth-Death with density dependence 

$$
\begin{aligned}
\frac{dx}{dt} = b_{0} x (t) - (d_{0} + d_{1}x)x(t) 
\end{aligned}
$$


- Allee Effect Model (1). Birth-Death-Allee

$$
\begin{aligned}
\frac{dx}{dt} = b_{0} x \frac{x}{\theta + x}(t) - (d_{0} + d_{1}x)x(t) 
\end{aligned}
$$

- Resource Model (2). Birth-Death-Immigration
  
$$
\begin{aligned}
\frac{dx}{dt} = b_{0}x(t) - (d_{0} + d_{1}x)x(t) + \iota(t)
\end{aligned}
$$
  
- Complex Model (3). Birth-Death-Allee-Immigration

$$
\begin{aligned}
\frac{dx}{dt} = b_{0}x(t) \frac{x}{\theta + x}(t) - (d_{0} + d_{1}x)x(t) + \iota(t)
\end{aligned}
$$

### Analysis 

_Host Population Infection Metrics_

* Calculate traditional infection metrics from host perspective (prevalence, mean abundance, mean intensity)
* From the stationary distribution, calculate aggregation metrics ($k$, mean intensity:variance intensity, skew, Pareto fraction). Is the 80:20 rule met? 
* How these values vary with Allee effect, immigration, and host-induced mortality changes

_Parasite Population Metrics_

* Calculate infection metrics from parasite perspective (total population size, crowding)
* How these values vary with Allee effect, immigration, and host-induced mortality changes
* Calculate stability? 

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
- Allen LJ, Allen EJ. A comparison of three different stochastic population models with regard to persistence time. Theoretical Population Biology. 2003 Dec 1;64(4):439-49.
[link](https://paperpile.com/app/p/74de65e5-f1a5-09d9-ad25-c40f4255334c)
- Dennis B. Allee effects in stochastic populations. Oikos. 2002 Mar;96(3):389-401. [link](https://paperpile.com/app/p/9b250b51-ad73-0bb8-a125-bf388cb4a48d)
- Drake JM, Lodge DM. Allee effects, propagule pressure and the probability of establishment: risk analysis for biological invasions. Biological Invasions. 2006 Mar;8:365-75. [link](https://paperpile.com/app/p/3bc8340a-08e5-0d17-ba8e-c389fd6b2fd4)
- May RM. Togetherness among schistosomes: its effects on the dynamics of the infection. Mathematical biosciences. 1977 Jan 1;35(3-4):301-43. [link](https://paperpile.com/app/p/674123bf-e4d1-0e4a-beb6-083ab0bba326)
- Sackett LC. Does the host matter? Variable influence of host traits on parasitism rates. International journal for parasitology. 2018 Jan 1;48(1):27-39. [link](https://paperpile.com/app/p/0d7f3b9f-da58-0d37-87ab-1b5ed9b5dab6)
  
### CHANGE-LOG:
- (September 25, 2023) Removed host body size as a factor; this eliminates changes in K 
- (September 21, 2023) Restructured protocol for population growth model instead of compartmental model. 
- (January 11, 2023) Replaced "Food provisioning, behavioral defense, and host-ectoparasite dynamics" title added abstract
