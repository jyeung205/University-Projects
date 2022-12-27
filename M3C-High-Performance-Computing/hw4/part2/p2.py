"""Final project, part 2"""
import numpy as np
import matplotlib.pyplot as plt
from m1 import bmodel as bm #assumes p2.f90 has been compiled with: f2py3 -c p2.f90 -m m1
import time
from scipy import optimize

def simulate_jacobi(n,input_num=(10000,1e-8),input_mod=(1,1,1,2,1.5),display=False):
    """ Solve contamination model equations with
        jacobi iteration.
        Input:
            input_num: 2-element tuple containing kmax (max number of iterations
                        and tol (convergence test parameter)
            input_mod: 5-element tuple containing g,k_bc,s0,r0,t0 --
                        g: bacteria death rate
                        k_bc: r=1 boundary condition parameter
                        s0,r0,t0: source function parameters
            display: if True, a contour plot showing the final concetration field is generated
        Output:
            C,deltac: Final concentration field and |max change in C| each iteration
    """
    #Set model parameters------

    kmax,tol = input_num
    g,k_bc,s0,r0,t0 = input_mod
    #-------------------------------------------
    #Set Numerical parameters
    Del = np.pi/(n+1)
    r = np.linspace(1,1+np.pi,n+2)
    t = np.linspace(0,np.pi,n+2) #theta
    tg,rg = np.meshgrid(t,r) # r-theta grid

    #Factors used in update equation
    rinv2 = 1.0/(rg*rg)
    fac = 1.0/(2 + 2*rinv2+Del*Del*g)
    facp = (1+0.5*Del/rg)*fac
    facm = (1-0.5*Del/rg)*fac
    fac2 = fac*rinv2

    #set initial condition/boundary conditions
    C = (np.sin(k_bc*tg)**2)*(np.pi+1.-rg)/np.pi

    #set source function, Sdel2 = S*del^2*fac
    Sdel2 = s0*np.exp(-20.*((rg-r0)**2+(tg-t0)**2))*(Del**2)*fac

    deltac = []
    Cnew = C.copy()

    #Jacobi iteration
    for k in range(kmax):
        #Compute Cnew
        Cnew[1:-1,1:-1] = Sdel2[1:-1,1:-1] + C[2:,1:-1]*facp[1:-1,1:-1] + C[:-2,1:-1]*facm[1:-1,1:-1] + (C[1:-1,:-2] + C[1:-1,2:])*fac2[1:-1,1:-1] #Jacobi update
        #Compute delta_p
        deltac += [np.max(np.abs(C-Cnew))]
        C[1:-1,1:-1] = Cnew[1:-1,1:-1]
        if k%1000==0: print("k,dcmax:",k,deltac[k])
        #check for convergence
        if deltac[k]<tol:
            print("Converged,k=%d,dc_max=%28.16f " %(k,deltac[k]))
            break

    deltac = deltac[:k+1]

    if display:
        plt.figure()
        plt.contour(t,r,C,50)
        plt.xlabel('theta')
        plt.ylabel('r')
        plt.title('Final concentration field')

    return C,deltac


def simulate(n,input_num=(10000,1e-8),input_mod=(1,1,1,2,1.5),display=True):
    """ Solve contamination model equations with
        OSI method, input/output same as in simulate_jacobi above
    """
    #Set model parameters------

    kmax,tol = input_num
    g,k_bc,s0,r0,t0 = input_mod

    #-------------------------------------------
    #Set Numerical parameters
    Del = np.pi/(n+1)
    r = np.linspace(1,1+np.pi,n+2)
    t = np.linspace(0,np.pi,n+2) #theta
    tg,rg = np.meshgrid(t,r) # r-theta grid

    #Factors used in update equation
    rinv2 = 1.0/(rg*rg)
    fac = 1.0/(2 + 2*rinv2+Del*Del*g)
    facp = (1+0.5*Del/rg)*fac
    facm = (1-0.5*Del/rg)*fac
    fac2 = fac*rinv2

    #set initial condition/boundary conditions
    C = (np.sin(k_bc*tg)**2)*(np.pi+1.-rg)/np.pi

    #set source function, Sdel2 = S*del^2*fac
    Sdel2 = s0*np.exp(-20.*((rg-r0)**2+(tg-t0)**2))*(Del**2)*fac

    deltac = []
    Cnew = C.copy()

    #Over-step iteration
    for k in range(kmax):
        for i in range(1,n+1):
            for j in range(1,n+1):
                Cnew[i,j] = -0.5*C[i,j] + 1.5*Sdel2[i,j] + 1.5*C[i+1,j]*facp[i,j] + 1.5*Cnew[i-1,j]*facm[i,j] + 1.5*(Cnew[i,j-1] + C[i,j+1])*fac2[i,j]
        #Compute delta_p
        deltac += [np.max(np.abs(C-Cnew))]
        C[1:-1,1:-1] = Cnew[1:-1,1:-1]
        if k%1000==0: print("k,dcmax:",k,deltac[k])
        #check for convergence
        if deltac[k]<tol:
            print("Converged,k=%d,dc_max=%28.16f " %(k,deltac[k]))
            break

    deltac = deltac[:k+1]

    if display:
        plt.figure()
        plt.contour(t,r,C,50)
        plt.xlabel('theta')
        plt.ylabel('r')
        plt.title('Final concentration field')
        plt.show()

    #C,deltac = None,None #Must be replaced
    return C,deltac


def simulate_jacobi4(n,tstar,input_num=(10000,1e-8),input_mod=(1,2,2,1+np.pi/2,1.5),display=True):
    """Boundary condition modified for question 4
    """
    #Set model parameters------

    kmax,tol = input_num
    g,k_bc,s0,r0,t0 = input_mod
    #-------------------------------------------
    #Set Numerical parameters
    Del = np.pi/(n+1)
    r = np.linspace(1,1+np.pi,n+2)
    t = np.linspace(0,np.pi,n+2) #theta
    tg,rg = np.meshgrid(t,r) # r-theta grid

    #Factors used in update equation
    rinv2 = 1.0/(rg*rg)
    fac = 1.0/(2 + 2*rinv2+Del*Del*g)
    facp = (1+0.5*Del/rg)*fac
    facm = (1-0.5*Del/rg)*fac
    fac2 = fac*rinv2

    #set initial condition/boundary conditions

    C = np.exp(-10*(tg-tstar)**2)*(np.sin(k_bc*tg)**2)*(np.pi+1.-rg)/np.pi

    #set source function, Sdel2 = S*del^2*fac
    Sdel2 = s0*np.exp(-20.*((rg-r0)**2+(tg-t0)**2))*(Del**2)*fac

    deltac = []
    Cnew = C.copy()

    #Jacobi iteration
    for k in range(kmax):
        #Compute Cnew
        Cnew[1:-1,1:-1] = Sdel2[1:-1,1:-1] + C[2:,1:-1]*facp[1:-1,1:-1] + C[:-2,1:-1]*facm[1:-1,1:-1] + (C[1:-1,:-2] + C[1:-1,2:])*fac2[1:-1,1:-1] #Jacobi update
        #Compute delta_p
        deltac += [np.max(np.abs(C-Cnew))]
        C[1:-1,1:-1] = Cnew[1:-1,1:-1]
        if k%1000==0: print("k,dcmax:",k,deltac[k])
        #check for convergence
        if deltac[k]<tol:
            print("Converged,k=%d,dc_max=%28.16f " %(k,deltac[k]))
            break

    deltac = deltac[:k+1]

    if display:
        plt.figure()
        plt.contour(t,r,C,50)
        plt.xlabel('theta')
        plt.ylabel('r')
        plt.title('Final concentration field')
        plt.show()

    return C,deltac


def performance():
    """
    Looking at p31, we can see that the time ratio of Fortran JAC against Fortran OSI stays the same even as n increases.
    This is demonstrated in p33, where the times seem to increase proportionally as n increases. Fortran OSI is
    faster than Fortran JAC which is expected as the Jacobi iteration is the most inefficient iterative solver.

    The time ratio of Python JAC against Python OSI increases as n increases. Looking at p32, we can see that the time for
    Python OSI increases much more than the time for Python JAC. This is expected as I could not vectorise the
    code for the OSI method and had to use loops.

    Looking at the ratio of Python JAC against Fortran JAC, we can see that the time ratio is originally very high,
    but it decreases as n increases.
    Looking at p34, we can see that the time ratio increases as n increases up to n=25. But it begins to decrease after n = 25.
    Python JAC is much slower than Fortran JAC.
    """
    times_osi, times_jac, times_forjac, times_forosi = list(), list(), list(), list()
    x_axis = list(range(5,50,10))
    for i in x_axis:
        t0 = time.time()
        simulate_jacobi(i)
        t1 = time.time()
        time_jac = t1-t0
        times_jac.append(time_jac)

        t0 = time.time()
        simulate(i)
        t1 = time.time()
        time_osi = t1-t0
        times_osi.append(time_osi)

        t0 = time.time()
        bm.simulate_jacobi(i)
        t1 = time.time()
        time_forjac = t1-t0
        times_forjac.append(time_forjac)

        t0 = time.time()
        bm.simulate(i)
        t1 = time.time()
        time_forosi = t1-t0
        times_forosi.append(time_forosi)

    plt.figure()
    plt.plot(x_axis,np.asarray(times_osi)/np.asarray(times_jac),label="Python OSI / Python JAC")
    plt.plot(x_axis,np.asarray(times_forjac)/np.asarray(times_forosi),label="Fortran JAC/ Fortran OSI")
    plt.plot(x_axis,np.asarray(times_jac)/np.asarray(times_forjac),label="Python JAC/ Fortran JAC")
    plt.title("Performance - Ratio walltimes of OSI and JAC")
    plt.xlabel('n')
    plt.ylabel('Time ratios')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(x_axis,times_osi,label="Python OSI")
    plt.plot(x_axis,times_jac,label="Python JAC")
    plt.title("Performance - Python OSI and JAC times")
    plt.xlabel('n')
    plt.ylabel('Time,s')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(x_axis,times_forosi,label="Fortran OSI")
    plt.plot(x_axis,times_forjac,label="Fortran JAC")
    plt.title("Performance - Fortran OSI and JAC times")
    plt.xlabel('n')
    plt.ylabel('Time,s')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(x_axis,np.asarray(times_osi)/np.asarray(times_forosi),label="Python JAC/ Fortran JAC")
    plt.title("Performance - Ratio walltimes Python JAC and Python OSI")
    plt.xlabel('n')
    plt.ylabel('Time ratios')
    plt.legend()
    plt.grid()
    plt.show()


    return None


def analyze():
    """The best choice of theta* is 0 or pi.
    I decided that the best theta* would be the value that gives me the lowest total concentration,
    so I plotted the sum of the concentration against theta*, using the parameters given in the question.

    Looking at figure p41, we can see that there is a local minimum at pi/2. However the minimum values are at 0 and pi.
    This makes sense the boundary condition shows that the concentration gets less as you go away from theta*,
    so you want to have theta* as far away as possible in order to minimise the concentration.

    Testing different parameters, I changed k to 1. Looking at figure p42, we can see that there is a maximum
    at theta* = pi/2 and that the minimums are at 0 and pi. The minimums will always be at 0 and pi.
    This reaffirms my conclusion that the best theta is 0 or pi.

    Looking at figure p43 and p44 we get similar results when we replace sum of concentration with max concentration.
    This reaffirms my conclusion.
    """

    Csums = list()
    x_axis = list(np.linspace(0,np.pi,50))
    for tstar in x_axis:
        C,_ = simulate_jacobi4(50,tstar,(1000,1e-8),(1,2,2,1+np.pi/2,1.5),display=False)
        Csum = np.sum(C)
        Csums.append(Csum)
    plt.figure()
    plt.plot(x_axis,Csums,label="k=2")
    plt.title("Sum of concentration against theta*")
    plt.xlabel('Theta*')
    plt.ylabel('Sum of concentration')
    plt.legend()
    plt.grid()
    plt.show()

    Csums = list()
    x_axis = list(np.linspace(0,np.pi,50))
    for tstar in x_axis:
        C,_ = simulate_jacobi4(50,tstar,(1000,1e-8),(1,1,2,1+np.pi/2,1.5),display=False)
        Csum = np.sum(C)
        Csums.append(Csum)
    plt.figure()
    plt.plot(x_axis,Csums,label="k=1")
    plt.title("Sum of concentration against theta*")
    plt.xlabel('Theta*')
    plt.ylabel('Sum of concentration')
    plt.legend()
    plt.grid()
    plt.show()

    Cmaxs = list()
    x_axis = list(np.linspace(0,np.pi,50))
    for tstar in x_axis:
        C,_ = simulate_jacobi4(50,tstar,(1000,1e-8),(1,2,2,1+np.pi/2,1.5),display=False)
        Cmax = np.max(C)
        Cmaxs.append(Cmax)
    plt.figure()
    plt.plot(x_axis,Cmaxs,label="k=2")
    plt.title("Max of concentration against theta*")
    plt.xlabel('Theta*')
    plt.ylabel('Max concentration at a point')
    plt.legend()
    plt.grid()
    plt.show()

    Cmaxs = list()
    x_axis = list(np.linspace(0,np.pi,50))
    for tstar in x_axis:
        C,_ = simulate_jacobi4(50,tstar,(1000,1e-8),(1,1,2,1+np.pi/2,1.5),display=False)
        Cmax = np.max(C)
        Cmaxs.append(Cmax)
    plt.figure()
    plt.plot(x_axis,Cmaxs,label="k=1")
    plt.title("Max of concentration against theta*")
    plt.xlabel('Theta*')
    plt.ylabel('Max concentration at a point')
    plt.legend()
    plt.grid()
    plt.show()



    return None


if __name__=='__main__':
    #Add code below to call performance and analyze
    #and generate figures you are submitting in
    #your repo.

    input=()
    output = performance()
    output = analyze()
