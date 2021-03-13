## Jennefer Maldonado
## Homework #2, Problem #5
## 03/15/2021

# Implement the Backward Time Centered Space and Crank Nicolson Scheme scheme
# v_t = v_xx, x in (0,1)
#v(0,t) = v(1,t) = 0
# v(x,0) = sin(3pix)
# Exact Solution:
# v(x,t) = exp(-9pi^2t)sin(3pix)


## import statements
import math
import numpy as np

# Implement the Backward Time Centered Space Scheme
def BTCS(N, xmin, xmax, t, tmax):
    # delta x
    dx = (xmax-xmin)/(N-1.0)
    # delta t
    dt = 0.0004
    
    # define Q of size NxN
    Q = np.zeros((N+1, N+1))
    
    # define what r is for current dt and dt
    r = dt/dx**2    
    col = 1
    for row in range(1,N):
        Q[row][col-1] = -r
        Q[row][col] = 1+(2*r)
        Q[row][col+1] = -r
        col+=1
        
    # first row initialization
    Q[0][0] = 1
    # last row initialization
    Q[N][N] = 1

    # initial and boundary conditions
    u_neg1 = np.zeros(N+1)
    x = np.linspace(0,1, N+1)
    for i in range(0, N+1):
        u_neg1[i] = math.sin(3*math.pi*x[i]) 
    
    # define b vector size N
    b = np.zeros(N+1)
    
    # time step loop
    while t < tmax:  
        for i in range(0,N):
            b[i] = -u_neg1[i]
        b[0] = 0
        b[N] = 0
        # use numpy solver
        u = np.linalg.solve(Q, u_neg1)
        # update u
        u_neg1 = u
        # update time
        t = t + dt    
    # compute analytic solution at t   
    actual = []
    for i in range(0,N+1):
        actual.append(math.exp(-9.0*(math.pi**2.0)*t)* math.sin(3.0*math.pi*x[i]))
    sub= []
    for i in range(0, len(u)):
        sub.append(u_neg1[i]-actual[i])
    # CHANGE NORM HERE (default is 2 norm)
    norm_2 = np.linalg.norm(sub)
    # infinity norm
    norm_inf = np.linalg.norm(sub, np.inf)
    # compute error
    error_2 = norm_2/math.sqrt(N)
    error_inf = norm_inf/math.sqrt(N)
    return (error_2,error_inf)



## GIVEN CONDITIONS FOR X AND T
# Interval for x creates the max and min values
xmin = 0
xmax = 1.0
#values of t
t = 0.0
tmax = 0.04

# VALUES OF N
N = 64
[err2, errinf] = BTCS(N, xmin, xmax, t, tmax)
print("Backward Time Center Space for Diffusion Equation")
print(f'Error with 2-norm: {err2}')
print(f'Error with infinity norm: {errinf}')

def CN_implicit(N, xmin, xmax, t, tmax):
    # delta x
    dx = (xmax-xmin)/(N-1.0)
    # delta t
    dt = 0.0004
    
    # define Q of size NxN
    Q = np.zeros((N+1, N+1))
    
    # define what r is for current dt and dt
    r = dt/dx**2    
    col = 1
    for row in range(1,N):
        Q[row][col-1] = -0.5*r
        Q[row][col] = 1 + r
        Q[row][col+1] = -0.5*r
        col+=1
        
    # first row initialization
    Q[0][0] = 1
    #Q[0][1] = 0
    # last row initialization
    #Q[N][N-1] = 0
    Q[N][N] = 1

    # initial and boundary conditions
    u_neg1 = np.zeros(N+1)
    x = np.linspace(0,1, N+1)
    for i in range(0, N+1):
        u_neg1[i] = math.sin(3*math.pi*x[i]) 
    
    # define b vector size N
    b = np.zeros(N+1)
    # time step loop
    while t < tmax:  
        for i in range(0,N):
            b[i] = -u_neg1[i]
        b[0] = 0
        b[N] = 0
        u = np.linalg.solve(Q, u_neg1)
        # update u
        u_neg1 = u
        # update time
        t = t + dt    
    # compute analytic solution at t   
    actual = []
    for i in range(0,N+1):
        actual.append(math.exp(-9.0*(math.pi**2.0)*t)* math.sin(3.0*math.pi*x[i]))
    sub= []
    for i in range(0, len(u)):
        sub.append(u_neg1[i]-actual[i])
    # 2 norm 
    norm_2 = np.linalg.norm(sub)
    # infinity norm
    norm_inf = np.linalg.norm(sub, np.inf)
    
    error_2 = norm_2/math.sqrt(N)
    error_inf = norm_inf/math.sqrt(N)
    return (error_2,error_inf)



## GIVEN CONDITIONS FOR X AND T
# Interval for x creates the max and min values
xmin = 0
xmax = 1.0
#values of t
t = 0.0
tmax = 0.04
# VALUE OF N
N = 64

[err2, errinf] = CN_implicit(N, xmin, xmax, t, tmax)
print("Implicit Crank Nicolson for Diffusion Equation")
print(f'Error with 2-norm: {err2}')
print(f'Error with infinity norm: {errinf}')

def FTCS(N, xmin, xmax, t, tmax):
    # delta x
    dx = (xmax-xmin)/(N-1.0)
    dt = dx**2/2

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
        v_n.insert(0, 0.0)
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
    
    # 2 norm 
    norm_2 = np.linalg.norm(sub)
    # infinity norm
    norm_inf = np.linalg.norm(sub, np.inf)
    
    error_2 = norm_2/math.sqrt(N)
    error_inf = norm_inf/math.sqrt(N)
    return (error_2,error_inf)

## GIVEN CONDITIONS FOR X AND T
# Interval for x creates the max and min values
xmin = 0
xmax = 1.0
#values of t
t = 0.0
tmax = 0.04
# VALUE OF N
N = 64

[err2, errinf] = FTCS(N, xmin, xmax, t, tmax)
print("Explicit FTCS for Diffusion Equation")
print(f'Error with 2-norm: {err2}')
print(f'Error with infinity norm: {errinf}')