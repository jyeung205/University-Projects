"""M3C 2018 Homework 3 01203261 Jimmy Yeung
Contains five functions:
    plot_S: plots S matrix -- use if you like
    simulate2: Simulate tribal competition over m trials. Return: all s matrices at final time
        and fc at nt+1 times averaged across the m trials.
    performance: To be completed -- analyze and assess performance of python, fortran, and fortran+openmp simulation codes
    analyze: To be completed -- analyze influence of model parameter, g
    visualize: To be completed -- generate animation illustrating "interesting" tribal dynamics
"""
import numpy as np
import matplotlib.pyplot as plt
from m1 import tribes as tr #assumes that hw3_dev.f90 has been compiled with: f2py3 --f90flags='-fopenmp' -c hw3_dev.f90 -m m1 -lgomp
#May also use scipy and time modules as needed
import time
import matplotlib.animation as animation


def plot_S(S):
    """Simple function to create plot from input S matrix
    """
    ind_s0 = np.where(S==0) #C locations
    ind_s1 = np.where(S==1) #M locations
    plt.plot(ind_s0[1],ind_s0[0],'rs')
    plt.plot(ind_s1[1],ind_s1[0],'bs')
    plt.show()
    return None
#------------------

def plot(S):
    """Similar to plot_S. Created so that I can use it in my updatefig function in visualise.
    """
    ind_s0 = np.where(S==0) #C locations
    ind_s1 = np.where(S==1) #M locations
    plot = plt.plot(ind_s0[1],ind_s0[0],'rs', ind_s1[1],ind_s1[0],'bs')
    return plot

def simulate2(N,Nt,b,e,g,m):
    """Simulate m trials of C vs. M competition on N x N grid over
    Nt generations. b, e, and g are model parameters
    to be used in fitness calculations.
    Output: S: Status of each gridpoint at end of simulation, 0=M, 1=C
            fc_ave: fraction of villages which are C at all Nt+1 times
                    averaged over the m trials
    """
    #Set initial condition
    S  = np.ones((N,N,m),dtype=int) #Status of each gridpoint: 0=M, 1=C
    j = int((N-1)/2)
    S[j,j,:] = 0
    N2inv = 1./(N*N)

    fc_ave = np.zeros(Nt+1) #Fraction of points which are C
    fc_ave[0] = S.sum()

    #Initialize matrices
    NB = np.zeros((N,N,m),dtype=int) #Number of neighbors for each point
    NC = np.zeros((N,N,m),dtype=int) #Number of neighbors who are Cs
    S2 = np.zeros((N+2,N+2,m),dtype=int) #S + border of zeros
    F = np.zeros((N,N,m)) #Fitness matrix
    F2 = np.zeros((N+2,N+2,m)) #Fitness matrix + border of zeros
    A = np.ones((N,N,m)) #Fitness parameters, each of N^2 elements is 1 or b
    P = np.zeros((N,N,m)) #Probability matrix
    Pden = np.zeros((N,N,m))
    #---------------------

    #Calculate number of neighbors for each point
    NB[:,:,:] = 8
    NB[0,1:-1,:],NB[-1,1:-1,:],NB[1:-1,0,:],NB[1:-1,-1,:] = 5,5,5,5
    NB[0,0,:],NB[-1,-1,:],NB[0,-1,:],NB[-1,0,:] = 3,3,3,3
    NBinv = 1.0/NB
    #-------------

    #----Time marching-----
    for t in range(Nt):
        R = np.random.rand(N,N,m) #Random numbers used to update S every time step

        #Set up coefficients for fitness calculation
        A = np.ones((N,N,m))
        ind0 = np.where(S==0)
        A[ind0] = b

        #Add boundary of zeros to S
        S2[1:-1,1:-1,:] = S

        #Count number of C neighbors for each point
        NC = S2[:-2,:-2,:]+S2[:-2,1:-1,:]+S2[:-2,2:,:]+S2[1:-1,:-2,:] + S2[1:-1,2:,:] + S2[2:,:-2,:] + S2[2:,1:-1,:] + S2[2:,2:,:]

        #Calculate fitness matrix, F----
        F = NC*A
        F[ind0] = F[ind0] + (NB[ind0]-NC[ind0])*e
        F = F*NBinv
        #-----------

        #Calculate probability matrix, P-----
        F2[1:-1,1:-1,:]=F
        F2S2 = F2*S2
        #Total fitness of cooperators in community
        P = F2S2[:-2,:-2,:]+F2S2[:-2,1:-1,:]+F2S2[:-2,2:,:]+F2S2[1:-1,:-2,:] + F2S2[1:-1,1:-1,:] + F2S2[1:-1,2:,:] + F2S2[2:,:-2,:] + F2S2[2:,1:-1,:] + F2S2[2:,2:,:]

        #Total fitness of all members of community
        Pden = F2[:-2,:-2,:]+F2[:-2,1:-1,:]+F2[:-2,2:,:]+F2[1:-1,:-2,:] + F2[1:-1,1:-1,:] + F2[1:-1,2:,:] + F2[2:,:-2,:] + F2[2:,1:-1,:] + F2[2:,2:,:]

        P = (P/Pden)*g + 0.5*(1.0-g) #probability matrix
        #---------

        #Set new affiliations based on probability matrix and random numbers stored in R
        S[:,:,:] = 0
        S[R<=P] = 1

        fc_ave[t+1] = S.sum()
        #----Finish time marching-----

    fc_ave = fc_ave*N2inv/m

    return S,fc_ave
#------------------

def simulate3(N,Nt,b,e,g,m):
    "Simulate3 has been created to return an array of the fc at every year. This is used in visualize.
    "

    #Set initial condition
    S  = np.ones((N,N,m),dtype=int) #Status of each gridpoint: 0=M, 1=C
    j = int((N-1)/2)
    S[j,j,:] = 0
    N2inv = 1./(N*N)

    #Initialize matrices
    NB = np.zeros((N,N,m),dtype=int) #Number of neighbors for each point
    NC = np.zeros((N,N,m),dtype=int) #Number of neighbors who are Cs
    S2 = np.zeros((N+2,N+2,m),dtype=int) #S + border of zeros
    F = np.zeros((N,N,m)) #Fitness matrix
    F2 = np.zeros((N+2,N+2,m)) #Fitness matrix + border of zeros
    A = np.ones((N,N,m)) #Fitness parameters, each of N^2 elements is 1 or b
    P = np.zeros((N,N,m)) #Probability matrix
    Pden = np.zeros((N,N,m))
    #---------------------

    #Calculate number of neighbors for each point
    NB[:,:,:] = 8
    NB[0,1:-1,:],NB[-1,1:-1,:],NB[1:-1,0,:],NB[1:-1,-1,:] = 5,5,5,5
    NB[0,0,:],NB[-1,-1,:],NB[0,-1,:],NB[-1,0,:] = 3,3,3,3
    NBinv = 1.0/NB
    #-------------
    S_evo = np.zeros((N,N,Nt+1))
    S_evo[:,:,0] = S[:,:,0]

    #----Time marching-----
    for t in range(Nt):
        R = np.random.rand(N,N,m) #Random numbers used to update S every time step

        #Set up coefficients for fitness calculation
        A = np.ones((N,N,m))
        ind0 = np.where(S==0)
        A[ind0] = b

        #Add boundary of zeros to S
        S2[1:-1,1:-1,:] = S

        #Count number of C neighbors for each point
        NC = S2[:-2,:-2,:]+S2[:-2,1:-1,:]+S2[:-2,2:,:]+S2[1:-1,:-2,:] + S2[1:-1,2:,:] + S2[2:,:-2,:] + S2[2:,1:-1,:] + S2[2:,2:,:]

        #Calculate fitness matrix, F----
        F = NC*A
        F[ind0] = F[ind0] + (NB[ind0]-NC[ind0])*e
        F = F*NBinv
        #-----------

        #Calculate probability matrix, P-----
        F2[1:-1,1:-1,:]=F
        F2S2 = F2*S2
        #Total fitness of cooperators in community
        P = F2S2[:-2,:-2,:]+F2S2[:-2,1:-1,:]+F2S2[:-2,2:,:]+F2S2[1:-1,:-2,:] + F2S2[1:-1,1:-1,:] + F2S2[1:-1,2:,:] + F2S2[2:,:-2,:] + F2S2[2:,1:-1,:] + F2S2[2:,2:,:]

        #Total fitness of all members of community
        Pden = F2[:-2,:-2,:]+F2[:-2,1:-1,:]+F2[:-2,2:,:]+F2[1:-1,:-2,:] + F2[1:-1,1:-1,:] + F2[1:-1,2:,:] + F2[2:,:-2,:] + F2[2:,1:-1,:] + F2[2:,2:,:]

        P = (P/Pden)*g + 0.5*(1.0-g) #probability matrix
        #---------

        #Set new affiliations based on probability matrix and random numbers stored in R
        S[:,:,:] = 0
        S[R<=P] = 1

        S_evo[:,:,t+1] = S[:,:,0]

    return S_evo

def performance(input=(None),display=False):
    """I wanted to analyse the speedup of Fortran+OMP against Fortran as the no. of years, Nt, increases for a fixed m = 100.
    Looking at hw321, we can see that the speedup is between 1.6 and 1.85 for all years 10 to 200.
    This agrees with the fact that the speedup will be less than the number of threads used, which in this case is 2.

    In hw322, I have fixed m = 100 and plotted the time it takes for the functions to run against
    the number of years, Nt, it is iterated over.
    I have done this to see the relative performances of Python, Fortran and Fortran+OMP.
    Python is by far the slowest - takes multiple times longer than Fortran and Fortran+OMP.
    The time it takes for Fortran is just under twice as long as Fortran+OMP - just as expected from
    our calculated speedup.
    The time each function takes is proportional to the number of years iterated.
    """

    tr.tr_b = 1.1
    tr.tr_e = 0.01
    tr.tr_g = 0.95
    tr.numthreads = 2

    #Speedup of parallelised code
    speedups = list()
    x_axis = list(range(10,210,10))
    for i in x_axis:
        t0 = time.time()
        for i in range(20):
            tr.simulate2_f90(21,i,100)
        t1 = time.time()
        time_vec = t1-t0

        t0 = time.time()
        for i in range(20):
            tr.simulate2_omp(21,i,100)
        t1 = time.time()
        time_parallel = (t1-t0)

        speedup = time_vec/time_parallel
        speedups.append(speedup)
    plt.figure()
    plt.plot(x_axis,speedups)
    plt.title("Speedup of Fortran+OMP against Fortran - m = 100")
    plt.xlabel('Nt, years')
    plt.ylabel('Time, s')
    plt.show()


    #Relative performance of Python, Fortran and Fortran+OMP
    times_for, times_omp, times_py = list(), list(), list()
    x_axis = list(range(10,110,10))
    for nt in x_axis:
         t0 = time.time()
         for i in range(20):
            tr.simulate2_f90(21,nt,100)
         t1 = time.time()
         time_for = (t1-t0)/20
         times_for.append(time_for)

         t0 = time.time()
         for i in range(20):
            tr.simulate2_omp(21,nt,100)
         t1 = time.time()
         time_omp = (t1-t0)/20
         times_omp.append(time_omp)

         t0 = time.time()
         for i in range(20):
            simulate2(21,nt,1.1,0.01,0.95,100)
         t1 = time.time()
         time_py = (t1-t0)/20
         times_py.append(time_py)
    plt.figure()
    plt.plot(x_axis,times_for, label="Fortran")
    plt.plot(x_axis,times_omp, label="Fortran + OMP")
    plt.plot(x_axis,times_py, label="Python")
    plt.title("Performance - m = 100")
    plt.xlabel('Nt, years')
    plt.ylabel('Time, s')
    plt.legend()
    plt.show()

    return None

#---------------------------
def analyze():
    """To investigate the influence off g, I plotted three graphs of different values of g for a fixed value of b.
    I plotted the % of villages which are C against the no. of years.
    We know that a larger value of b means that C decreases at a quicker rate.

    For b = 1.01, 1.25 and 1.5, the effect of g was very similar.
    For g between 0.8 and 0.95, a larger g meant that C decreases at a slower rate,
    but the value that it converges at also becomes smaller.

    For g = 1, the pattern was completely different to the other values of g.
    g = 1 meant that the rate at which C decreased was very slow, and the value of which
    C converges depends on b. As b increases, then the value at which C
    converges decreases.
    """

    x=range(200+1)

    plt.figure(figsize=(12,4))
    for g in np.linspace(0.8,1,5):
        _, fc_ave = simulate2(51,200,1.01,0.01,g,50)
        plt.plot(x, fc_ave*100, label=g)
    plt.title('Evolution of C for different values of g, b = 1.01')
    plt.xlabel('Number of years past')
    plt.ylabel('% of villages which are C')
    plt.legend()
    plt.show()

    plt.figure(figsize=(12,4))
    for g in np.linspace(0.8,1,5):
        _, fc_ave = simulate2(51,200,1.25,0.01,g,50)
        plt.plot(x, fc_ave*100, label=g)
    plt.title('Evolution of C for different values of g, b = 1.25')
    plt.xlabel('Number of years past')
    plt.ylabel('% of villages which are C')
    plt.legend()
    plt.show()

    plt.figure(figsize=(12,4))
    for g in np.linspace(0.8,1,5):
        _, fc_ave = simulate2(51,200,1.5,0.01,g,50)
        plt.plot(x, fc_ave*100, label=g)
    plt.title('Evolution of C for different values of g, b = 1.5')
    plt.xlabel('Number of years past')
    plt.ylabel('% of villages which are C')
    plt.legend()
    plt.show()

    return None #Modify as needed

#------------------------------------


def visualize():
    """Generate an animation illustrating the evolution of
        villages during C vs M competition
    """
    fig = plt.figure()
    S = simulate3(21,100,1.1,0.01,0.95,1)

    def updatefig(i):
        plt = plot(S[:,:,i])
        return plt

    ani = animation.FuncAnimation(fig, updatefig,frames=100,repeat=False,blit=True)
    ani.save("hw3movie.mp4")

    return

#------------------------------------
if __name__ == '__main__':
    #Modify the code here so that it calls performance analyze and
    #generates the figures that you are submitting with your code

    output_a = performance()

    output_b = analyze()

    output_c = visualize()
