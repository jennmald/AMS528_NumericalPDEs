## Jennefer Maldonado
## Homework #1, Problem #2
## 03/01/2021

# Implement the Forward Time Centered Space scheme
# v_t = v_xx, x in (0,1)
#v(0,t) = v(1,t) = 0
# v(x,0) = sin(3pix)
# Exact Solution:
# v(x,t) = exp(-9pi^2t)sin(3pix)

#################################################
# To change norms or values of delta t
# Comment and uncomment lines:
# 29 & 30 for delta t
# 67 & 69 for norms
#################################################

## import statements
import math
import numpy


def FTCS(N, xmin, xmax, t, tmax):
    # delta x
    dx = (xmax-xmin)/(N-1.0)

    # CHANGE DELTA T HERE
    dt = dx**2/2
    #dt = dx**2/6

    # initial conditions and boundary conditions
    x=[] #x step values
    v=[] #v values
    for i in range(1,N+1):
        x.append(xmin+(i-1)*dx)  
        v.append(math.sin(3*math.pi*x[i-1]))

    # time step loop
    while t < tmax:
        v_n = []
        # update solution
        for i in range(1,N-1):
            # forward time centered space
            v_n.append(v[i] +(dt/dx**2)*(v[i-1]-2*v[i]+v[i+1]))
     
        # boundary conditions
        # initial
        v_n.insert(0, 0.0);
        # final
        v_n.append(0.0);
    
        # update state for next time step
        v = v_n
        # new time
        t = t + dt
    
    # compute analytic solution at t   
    actual = []
    for i in range(0,N):
        actual.append(math.exp(-9.0*(math.pi**2.0)*t)* math.sin(3.0*math.pi*x[i]))

    sub= []
    for i in range(0, len(v)):
        sub.append(v[i]-actual[i])
    # CHANGE NORM HERE (default is 2 norm)
    norm = numpy.linalg.norm(sub)
    # infinity norm
    #norm = numpy.linalg.norm(sub, numpy.inf)
    
    error=norm/math.sqrt(N)
    return [error, dx]



## GIVEN CONDITIONS FOR X AND T
# Interval for x creates the max and min values
xmin = 0
xmax = 1.0
#values of t
t = 0.0
tmax = 0.02

# VALUES OF N
N = [32,64,128,256]

# store the error and dx for each FTCS with different N
error = []
deltax = []
for n in N:
    [err, dx] = FTCS(n, xmin, xmax, t, tmax)
    error.append(err)
    deltax.append(dx)

# Display errors 
for i in range(0,len(N)):
    print('The error for N='+ str(N[i])+ ' is ' + str(error[i]))    
print()

# Compute order of convergence
for i in range(1, len(error)):
    err = error[i]/error[i-1]
    delx = deltax[i]/deltax[i-1]
    R = math.log(err,2)/math.log(delx)
    print('The convergence order for N=' + str(N[i]) +' is '+ str(R))
        
