#a value has been taken as the default from this blog post: https://timchurches.github.io/blog/posts/2020-03-18-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-2/ 
#unless stated why it has been changed in the comments. 

simulate <- function(# control.icm params
  #type of modeL: Susceptible, Exposed, Infected, Quarintined, Hospitalised, Recovered, Fatalaties
  type = "SEIQHRF",
  #Number of days for simulation. 
  #Note that day 1 is for initialisation, day 2 is the first day of the simulation, 
  nsteps = 50, 
  #Number of simulations to run and then take an average.
  nsims = 8,
  #Number of CPU cores to use for parallel execution.
  ncores = 1,
  #Method for progression from E compartment to I. if FALSE, random draws from a Weibull distribution
  prog.rand = FALSE,
  #	Method for recovery transition from I, Q or H to R. If FALSE, random draws from a random draws from a Weibull distribution
  rec.rand = FALSE,
  #Method for case fatality transition from H to F. If FALSE, random sample with a sample fraction also given by fat.rate.base.  
  fat.rand = FALSE,
  #Method for self-isolation transition from I to Q. If TRUE, random binomial draws at quar.rate, if FALSE, random sample with a sample fraction also given by `quar.rate.
  quar.rand = FALSE,
  #Method for transition from I or Q to H. If TRUE, random binomial draws at hosp.rate, if FALSE, random sample with a sample fraction also given by `hosp.rate. 
  hosp.rand = FALSE,
  #Method for transition from H to R. If TRUE, random binomial draws at disch.rate, if FALSE, random sample with a sample fraction also given by disch.rate.
  disch.rand = FALSE,
  #name of the function to implement infection processes. Using the default.
  infection.FUN = infection.seiqhrf.icm,
  #name of the function to implement recovery processes. Using the default.
  recovery.FUN = progress.seiqhrf.icm,
  #Function that handles background demographics, specifically departures.
  departures.FUN = departures.seiqhrf.icm,
  #Function that handles background demographics, specifically arrivals (births and immigration). Using the default. 
  arrivals.FUN = arrivals.icm,
  #Utility function that collects prevalence and transition time data from each run and stores it away in the simulation result object. Use the default.
  get_prev.FUN = get_prev.seiqhrf.icm,
  # init.icm params
  #Initial number of *S compartment individuals in the simulated population. We take 400,000 as the essential people in Punjab. 
  s.num = 400000,
  #Initial number of E compartment individuals in the simulated population. Number of people exposed at the start.
  e.num=664, #come back to
  #Initial number of I compartment individuals in the simulated population.
  i.num = 31,
  #Initial number of Quarantined people. 
  q.num=30000,
  #Initial number of Hospitalised people.
  h.num=132,
  #Initial number of recovered people.
  r.num = 0,
  #Initial number of people who have died due to the virus.
  f.num = 1,
  # param.icm params
  #Probability of passing on infection at each exposure event for interactions between infectious people in the E compartment and susceptibles in S.
  #For now, kept same as inf.prob.i
  inf.prob.e = 0.05, #reduce introducing hygience 
  #The number of exposure events (acts) between infectious individuals in the E compartment and susceptible individuals in the S compartment, per day.
  act.rate.e = 20,
  #	Probability of passing on infection at each exposure event for interactions between infectious people in the I compartment and susceptibles in S.
  #Reducing inf.prob.i is equivalent to increasing hygiene measures, such as not putting hands in eyes, nose or moth, use of hand sanitisers, wearing masks by the infected, and so on.
  inf.prob.i = 0.05, 
  #	The number of exposure events (acts) between infectious individuals in the I compartment and susceptible individuals in the S compartment, per day.
  # Reducing act.rate.i is equivalent to increasing social distancing by people in the I compartment.
  #This assumes 20 people for now, with social distancing and optimising resources, this would reduce. 
  act.rate.i = 20, #social distancing
  #Probability of passing on infection at each exposure event for interactions between infectious people in the Q compartment and susceptibles in S. 
  #Note the default is lower than for inf.prob.i reflecting the greater care that self-isolated individuals will, on average, take regarding hygiene measures, such as wearing masks, to limit spread to others. 
  inf.prob.q = 0.02, 
  #The number of exposure events (acts) between infectious individuals in the Q compartment (isolated, self or otherwise) and susceptible individuals in the S compartment, per day.
  #This has been taken as people living in families in Punjab. In normal times, household size is 4.9. For now, lower due to precautions being taken. 
  act.rate.q = 2.5, 
  #Rate per day at which symptomatic (or tested positive), infected I compartment people enter self-isolation (Q compartment). Default for now. 
  quar.rate = 1/30, 
  #	Rate per day at which symptomatic (or tested positive), infected I compartment people or self-isolated Q compartment people enter the state of requiring hospital care -- that is, become serious cases.
  #Take date of first case in Punjab and where you are now, to get this rate. 
  hosp.rate = ?? #date of first case in Punjab and get the rate for this. 
  #Rate per day at which people needing hospitalisation recover.
  disch.rate = 1/15,
  # Rate per day at which people who are infected but asymptomatic (E compartment) progress to becoming symptomatic (or test-positive), the I compartment. 
  prog.rate = 1/10,
  #	Scale parameter for Weibull distribution for progression
  prog.dist.scale = 5,
  #Shape parameter for Weibull distribution for progression, see prog.rand for details. 
  prog.dist.shape = 1.5,
  #	Rate per day at which people who are infected and symptomatic (I compartment) recover, thus entering the R compartment.
  #This rates is between 10-24 days globally. Taken 20 for now. 
  rec.rate = 1/20,
  #	Scale parameter for Weibull distribution for recovery
  rec.dist.scale = 35,
  #Shape parameter for Weibull distribution for recovery
  rec.dist.shape = 1.5,
  #Baseline mortality rate per day for people needing hospitalisation (deaths due to the virus). 
  fat.rate.base = 1/50,
  #Number of available hospital beds for the modelled population. 
  hosp.cap = , #fill out exact number of beds
  #Mortality rate per day for people needing hospitalisation but who can't get into hospital due to the hospitals being full.
  #The default rate is twice that for those who do get into hospital.
  # Increasing hospital capacity would reduce this. 
  fat.rate.overcap = 1/25,
  #Time co-efficient for increasing mortality rate as time in the H compartment increases for each individual in it
  fat.tcoeff = 0.5,
  #Enables demographics, that is, arrivals and departures, to and from the simulated population.
  vital = TRUE,
  #Currently all arrivals go into the S compartment, the default is approximately the daily birth rate for Punjab. 
  a.rate = (15.2/365)/1000, 
  #This is the default
  a.prop.e = 0.01,
  #This is the default
  a.prop.i = 0.001,
  #This is the default
  a.prop.q = 0.01,
  #Background demographic departure (death not due to virus) rates. Defaults based on Punjab crude death rates.
  ds.rate = (6.2/365)/1000, 
  de.rate = (6.2/365)/1000, 
  #global death rate for Covid-19
  di.rate = (157/365)/1000,
  #global death rate for Covid-19
  dq.rate = (157/365)/1000,
  #The default is 7 (death rate), 20(for hospitalised). Based that, 157 (infected), in hospital, 449. 
  dh.rate = (449/365)/1000,
  dr.rate = (6.2/365)/1000,
  out="mean"
) {
  
  control <- control.icm(type = type, 
                         nsteps = nsteps, 
                         nsims = nsims,
                         ncores = ncores,
                         prog.rand = prog.rand,
                         rec.rand = rec.rand,
                         infection.FUN = infection.FUN,
                         recovery.FUN = recovery.FUN,
                         arrivals.FUN = arrivals.FUN,
                         departures.FUN = departures.FUN,
                         get_prev.FUN = get_prev.FUN)
  
  init <- init.icm(s.num = s.num,
                   e.num = e.num,
                   i.num = i.num,
                   q.num = q.num,
                   h.num = h.num,
                   r.num = r.num,
                   f.num = f.num)
  
  param <-  param.icm(inf.prob.e = inf.prob.e, 
                      act.rate.e = act.rate.e,
                      inf.prob.i = inf.prob.i, 
                      act.rate.i = act.rate.i,
                      inf.prob.q = inf.prob.q, 
                      act.rate.q = act.rate.q,                    
                      quar.rate = quar.rate,
                      hosp.rate = hosp.rate,
                      disch.rate = disch.rate,
                      prog.rate = prog.rate,
                      prog.dist.scale = prog.dist.scale,
                      prog.dist.shape = prog.dist.shape,
                      rec.rate = rec.rate,
                      rec.dist.scale = rec.dist.scale,
                      rec.dist.shape = rec.dist.shape,
                      fat.rate.base = fat.rate.base,
                      hosp.cap = hosp.cap,
                      fat.rate.overcap = fat.rate.overcap,
                      fat.tcoeff = fat.tcoeff,
                      vital = vital,
                      a.rate = a.rate, 
                      a.prop.e = a.prop.e,
                      a.prop.i = a.prop.i,
                      a.prop.q = a.prop.q,
                      ds.rate = ds.rate, 
                      de.rate = de.rate, 
                      di.rate = di.rate,
                      dq.rate = dq.rate,
                      dh.rate = dh.rate,
                      dr.rate = dr.rate)
  
  sim <- icm.seiqhrf(param, init, control)
  sim_df <- as.data.frame(sim, out=out)
  
  return(list(sim=sim, df=sim_df))
}
