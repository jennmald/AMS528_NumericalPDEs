## Jennefer Maldonado
## Homework #1, Problem #4
## 03/01/2021

## import statements
import math
import numpy
import matplotlib.pyplot as plt

## FORWARD TIME CENTERED SPACE ##
def FTCS(N, xmin, xmax, nu, a, tmax):
    # compute dx, dt, and initial t
    dx = (xmax - xmin)/(N-1.0)
    dt = min([(0.75*dx**2)/(2*nu),(0.75*dx)/a])
    t = 0.0
    # define x and u
    x = []
    u = []
    u_max = 100
    #sharp gaussian initial condition
    for i in range(1,N+1):
        x.append( xmin + (i-1)*dx )
        u.append( u_max * math.exp(-x[i-1]*x[i-1]*10.0) )
    #set initial cond    
    u0 = u
    v = u
    step = 0
    # loop through
    while t < tmax:
        u_n = []
        if t > tmax-dt:
            dt = tmax-t
        t = t+dt
        step = step+1
        for i in range(0,N-1):
            #FTcs
            u_n.append( u[i]-((a*dt)/(2*dx))*(u[i+1]-u[i-1])+ (nu*dt)/(dx**2)*(u[i+1]-2*u[i]+u[i-1]) )  
        #periodic boundary conditions
        u_n.insert(0, u_n[N-2])
        u_n.append( u_n[1] )
        u = u_n
    #return the final state
    return u_n


## FTCS for diffusion, upwind for advection ##
def FTCS_upwind(N, xmin, xmax, nu, a, tmax):
    # compute dx, dt, and initial t
    dx = (xmax - xmin)/(N-1.0)
    dt = min([(0.75*dx**2)/(2*nu),(0.75*dx)/a])
    t = 0.0
    # define x and u
    x = []
    u = []
    u_max = 100
    #sharp gaussian initial condition
    for i in range(1,N+1):
        x.append( xmin + (i-1)*dx )
        u.append( u_max * math.exp(-x[i-1]*x[i-1]*10.0) )
    #set initial cond    
    u0 = u
    v = u
    step = 0
    # loop through
    while t < tmax:
        u_n = []
        if t > tmax-dt:
            dt = tmax-t
        t = t+dt
        step = step+1
        for i in range(0,N-1):
            #FTcs
            u_n.append( u[i]-((a*dt)/dx)*(u[i]-u[i-1]) + (nu*dt)/(dx**2)*(u[i-1]-2*u[i]+u[i+1]) )  
        #periodic boundary conditions
        u_n.insert(0, u_n[N-2])
        u_n.append( u_n[1] )
        u = u_n
    #return the final state
    return u_n


## PROBLEM CONDITIONS ##
# number of time steps
N = 500
#values for x
xmin = -5.0
xmax = 15.0
#advection coefficient
a = 1
# values for t
tmax = 10

# diffusion coefficients
# values change to become smaller for analysis
nu_vals= [1, 0.5, 0.1, 0.05, 0.01, 0.001]

ftcs_values = []
upwind_values = []
for nu in nu_vals:
    ftcs_values.append(FTCS(N, xmin, xmax, nu, a, tmax))
    upwind_values.append(FTCS_upwind(N, xmin, xmax, nu, a, tmax))

plt.plot(ftcs_values[0])
plt.title('Figure 1: FTCS for nu = ' + str(nu_vals[0]))
plt.show()
plt.plot(ftcs_values[1])
plt.title('FTCS for nu = ' + str(nu_vals[1]))
plt.show()
plt.plot(ftcs_values[2])
plt.title('FTCS for nu = ' + str(nu_vals[2]))
plt.show()
plt.plot(ftcs_values[3])
plt.title('FTCS for nu = ' + str(nu_vals[3]))
plt.show()
plt.plot(ftcs_values[4])
plt.title('FTCS for nu = ' + str(nu_vals[4]))
plt.show()
plt.plot(ftcs_values[5])
plt.title('FTCS for nu = ' + str(nu_vals[5]))
plt.show()

plt.plot(upwind_values[0])
plt.title('FTCS upwind for nu = ' + str(nu_vals[0]))
plt.show()
plt.plot(upwind_values[1])
plt.title('FTCS upwind for nu = ' + str(nu_vals[1]))
plt.show()
plt.plot(upwind_values[2])
plt.title('FTCS upwind for nu = ' + str(nu_vals[2]))
plt.show()
plt.plot(upwind_values[3])
plt.title('FTCS upwind for nu = ' + str(nu_vals[3]))
plt.show()
plt.plot(upwind_values[4])
plt.title('FTCS upwind for nu = ' + str(nu_vals[4]))
plt.show()
plt.plot(upwind_values[5])
plt.title('FTCS upwind for nu = ' + str(nu_vals[5]))
plt.show()



