## Jennefer Maldonado
## Homework #1, Problem #3
## 03/01/2021

# Implement the Forward Time Centered Space scheme
# with left boundary condition of Neumann Type
# v_t = v_xx, x in (0,1)
#v(0,t) = 3piexp(-9pi^2t)
#v(1,t) = 0
# v(x,0) = sin(3pix)
# Exact Solution:
# v(x,t) = exp(-9pi^2t)sin(3pix)

## import statements
import math
import numpy


def FTCS(N, xmin, xmax, t, tmax):
    # delta x & t
    dx = (xmax-xmin)/(N-1.0)
    dt = dx**2/2

    x=[] #x step values
    v=[] #v values
    u=[]
    for i in range(1,N+1):
        x.append(xmin+(i-1)*dx)  
        v.append(math.sin(3*math.pi*x[i-1]))
        u.append(3*math.pi*math.exp(-9*(math.pi**2)*t))
    # time step loop
    while t < tmax:
        # for approx sol
        v_n = []
        # update solution
        for i in range(1,N-1):
            # forward time centered space
            v_n.append(v[i] +(dt/dx**2)*(v[i-1]-2*v[i]+v[i+1]))
            ## CHANGE NEUMANN BOUNDARY SCHEME HERE (default is vertex 2nd order)
            # vertex grid: centering 2nd order difference
            u_0 = (1-2*(dt/dx**2)*u[0]+2*(dt/dx**2)*u[1])
            # vertex grid one sided 1st order differencing
            # cell centered grid
            #u_0 = u[1]
            
        # boundary conditions
        # initial
        v_n.insert(0, u_0)
        # final
        v_n.append(0.0)
    
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
N = 40

#error and dx for each FTCS with different N
[err, dx] = FTCS(N, xmin, xmax, t, tmax)

print(err)