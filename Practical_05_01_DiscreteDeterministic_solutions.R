##########################################################
# Discrete time deterministic models Section 1           #
##########################################################

#############################
# A. Recurrence equations
#############################

# A.1 Calculate the daily population density of a colony of liver flukes where
# the density is given by the recurrence equation with a roughly monthly 
# forcing term
#
# N(t+1) = a*N(t) + b*(1 - cos(2*pi*t/28))
#
# for a = 0.953 and b = 0.05, and N(0) = 0.03. Let t = [0, 365].
# Plot the solution

times_lf <- seq(0, 365)

parms_lf <- c(a = 0.953, b = 0.05)

f_lf <- function(t, y, parms){
    parms["a"]*y + parms["b"]*(1 - cos(2*pi*t/28))
}

N_lf <- matrix(data = NA, ncol = 1, nrow = length(times_lf))
N_lf[1] <- 0.03

for (i in 1:(length(times) - 1)){
    N_lf[i+1] <- f_lf(times_lf[i], N_lf[i], parms_lf)
}

plot(times_lf, N_lf, type = "l")

# A.2 Describe the behaviour of the system. Why does it behave this way?
#
# Answer: The density grows over time, with the oscillation term from the
# forcing eventually balancing out the growth. There is no steady state
# solution which can be reached for the model due to the forcing.

# A.3 Change the initial condition to be 2.1. Describe the behaviour of the 
# system, particularly comparing it to the previous answer. You may choose
# to do this by writing a function which has the following inputs:
#  - an initial condition
#  - a vector of times
#  - a function f(t,y,parms)
#  - a vector of parameters
# and makes a plot of the solution to the above system. Check that it works by
# passing in the values from the model above. 
# 
# Answer: The density decreases from the initial condition, indicating that
# there are too many liver flukes for the body to support. 

lf_plot <- function(y0, times, f, parms){
    
    ny <- length(y0)
    nt <- length(times)
    y <- data.frame(matrix(data = NA,
                           nrow = nt,
                           ncol = ny))
    # set names
    names(y) <- names(y0)
    
    # initial condition
    y[1, ] <- y0
    
    
    for (i in 2:nt){
        y[i, ] <- f(times[i-1], y[i-1,], parms)
    }
    
    nr <- floor(sqrt(ny))
    nc <- ceiling(ny/nr)
    par(mfrow = c(nr, nc))
    
    for (i in 1:ny){
        plot(y[,i] ~ times, xlab = "Time")
    }
    
}

lf_plot(y0 = c(N = 2.1), times = times_lf, f = f_lf, parms = parms_lf)

#############################
# B. SIR model
#############################

# B.1 Adapt the SIR model in the slides to incorporate logistic growth for new
# susceptibles. Assume that "recovery" from the disease is actually death and 
# so the new susceptibles are born by both susceptibles and infectious.
#
# Write the system in the form y(t+1) = y(t) + delta t * f(t,y,theta)
#
# Use these parameter values: r = 0.33, K = 1500, beta = 4e-4, gamma = 0.12
# Let the initial population be S(0) = 1499, I(0) = 1, indicating that the
# population is at carrying capacity and that 1 person contracts the infection

f_sir <- function(t, y, parms){
    S <- y[1]
    I <- y[2]
    R <- y[3]
    
    dS <- parms["r"]*(S+I)*(1 - (S+I)/parms["K"]) - 
        parms["beta"]*S*I
    
    dI <- parms["beta"]*S*I - parms["gamma"]*I
    
    dR <- parms["gamma"]*I
    
    return(c(dS, dI, dR))
    
}

parms_sir <- c(r = 0.33, K = 1500, beta = 0.0004, gamma = 0.12)

times_sir <- seq(0, 100)
N_sir <- length(times_sir)

dat_sir <- data.frame(S = numeric(N_sir),
                      I = numeric(N_sir),
                      R = numeric(N_sir))

dat_sir[1, ] <- c(1499, 1, 0)

for (i in 2:N_sir){
    dat_sir[i, ] <- dat_sir[i-1, ] + f_sir(i, dat_sir[i-1, ], parms_sir)
}

# B.2 Calculate N(t) = S(t) + I(t) the total number of alive individuals. Make
# a plot of S(t), I(t), R(t) and N(t)

dat_sir$N <- dat_sir$S + dat_sir$I

par(mfrow = c(2,2))
for (i in 1:ncol(dat_sir)){
    plot(dat_sir[,i] ~ times_sir, type = "p", 
         xlab = "Time (years)",
         ylab = names(dat_sir)[i])
}


# B.3 Discuss what happens to the population of S and I over time. Consider
# the parameters of the model, what they represent, and whether the assumptions
# they represent appear to be met.
#
# Answer: The presence of the infection leads to a long-term population which
# is below the carrying capacity, K, as births of new susceptibles occure at a
# rate slower than  gamma*I. The initial behaviour is that in the first ten
# or so years there are very few fatalities and the sum of the susceptible and
# infectious populations are approximately equal to S(0) + I(0) = K.

#############################
# C. Markov models
#############################

# C.1 Consider the languorem docens outbreak from the slides. Modify the system
# of equations used so that 20% of the stage 1 population recover back to a 
# susceptible state instead of progressing from stage 1 to stage 2. 

P <- matrix(data = c(0.99, 0.01, 0,
                     0.2,  0.55, 0.25,
                     0,    0,    1), 
            nrow = 3, byrow = T)

x <- matrix(c(X = 1, Y1 = 0, Y2 = 0), nrow = 1)

# solve the system
time <- 0:100
N_ld <- length(time)

dat <- data.frame(X  = numeric(length = N_ld),
                  Y1 = numeric(length = N_ld),
                  Y2 = numeric(length = N_ld))

dat[1, ] <- c(1, 0, 0)

for (i in 2:length(time)){
    dat[i, ] <- as.matrix(dat[i-1, ] ) %*% P
}

dat$time <- time

par(mfrow = c(3,1), mar = c(4,4,2,1)+0.1)
for (i in 1:3){
    plot(dat[,i] ~ dat[,4], type = "l", xlab="Time (years)",
         ylab = names(dat)[i])
}

# C.2 Describe the behaviour of the system over time. Does it appear to be
# trending towards a stable solution?
# 
# Answer: the Y1 population increases over time and appears to reach a steady
# value where the number of transitions from susceptible to stage 1 is balanced
# by those progressing from stage 1 to stage 2.

# C.3 Redo the simulation over a time scale of 100 years. What does it appear
# the long-term steady state is and why?
#
# Answer: Eventually all susceptibles will end up in stage 2 of the disease
# due to the third row of the matrix states that there is no way to recover
# from stage 2.
