"""M345SC Homework 3, part 1
Jimmy Yeung 01203261
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import hann
import scipy
import time


def nwave(alpha,beta,Nx=256,Nt=801,T=200,display=False):
    """
    Question 1.1
    Simulate nonlinear wave model

    Input:
    alpha, beta: complex model parameters
    Nx: Number of grid points in x
    Nt: Number of time steps
    T: Timespan for simulation is [0,T]
    Display: Function creates contour plot of |g| when true

    Output:
    g: Complex Nt x Nx array containing solution
    """

    #generate grid
    L = 100
    x = np.linspace(0,L,Nx+1)
    x = x[:-1]

    def RHS(f,t,alpha,beta):
        """Computes dg/dt for model eqn.,
        f[:N] = Real(g), f[N:] = Imag(g)
        Called by odeint below
        """
        g = f[:Nx]+1j*f[Nx:]
        
        #add code here
        c = np.fft.fft(g)/Nx
        n = np.fft.fftshift(np.arange(-Nx/2,Nx/2))
        k = 2*np.pi*n/L
        d2g = Nx*np.fft.ifft(-(k**2)*c)
        #-----------
        dgdt = alpha*d2g + g - beta*g*g*g.conj()
        df = np.zeros(2*Nx)
        df[:Nx] = dgdt.real
        df[Nx:] = dgdt.imag
        return df

    #set initial condition
    g0 = np.random.rand(Nx)*0.1*hann(Nx)
    f0=np.zeros(2*Nx)
    f0[:Nx]=g0
    t = np.linspace(0,T,Nt)

    #compute solution
    f = odeint(RHS,f0,t,args=(alpha,beta))
    g = f[:,:Nx] + 1j*f[:,Nx:]

    if display:
        plt.figure()
        plt.contour(x,t,g.real)
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('Contours of Real(g)')

    return g


def analyze():
    """
    Question 1.2
    Add input/output as needed

    Discussion:
     
    In fig1 and fig2, I have plotted the amplitude of the Fourier coefficients 
    against Nt for case A and case B. For case A, the amplitudes of Nt = 801 is 
    and Nt = 200 are very similar. However, for case B then the amplitudes are
    similar near n=0 but for n near the boundaries then the amplitudes for Nt=200
    are slightly larger than Nt=801. 
    
    In fig3 and fig4, I have plotted the amplitude of the Fourier coefficients
    against Nx for case A and case B. For both case A and case B, we can clearly
    see that the amplitudes for Nx=100 is less than the amplitudes of Nx=256. 
    However, for case B the amplitudes are more similar for Nx=100 and Nx=256 - 
    the difference between the amplitudes between Nx=256 and Nx=100 are greater
    for case A. 
    
    In fig5 and fig6, I have plotted the amplitudes of the Fourier coefficients 
    against Nx for case A and case B. For both cases, the amplitudes of T=200 and
    T=100 are both very similar. 
    
    For all of the plots, the distribution of the amplitudes are
    very similar in the fact that the amplitudes are greatest at n=0 and 
    become smaller as n goes towards the boundaries. However, Case A is more chaotic 
    than case B as there is more fluctuation in the coefficients when Nt, Nx, T is varied.
    
    """
    Nx=256
    
    #Nt
    #caseA
    n = np.arange(-Nx/2,Nx/2)
    gA = nwave(1-2j,1+2j,Nt=801)
    gA1 = nwave(1-2j,1+2j,Nt=200)
    cA = np.fft.fft(gA[50:,])/Nx
    cA1 = np.fft.fft(gA1[50:,])/Nx
    
    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cA[-1,])),'.')
    plt.plot(n,np.fft.fftshift(np.abs(cA1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case A - Amplitude of Fourier Coefficients for difference values of Nt')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('Nt=801','Nt=200'))
    
    #caseB
    n = np.arange(-Nx/2,Nx/2)
    gB = nwave(1-1j,1+2j,Nt=801)
    gB1 = nwave(1-1j,1+2j,Nt=200)
    cB = np.fft.fft(gB[50:,])/Nx
    cB1 = np.fft.fft(gB1[50:,])/Nx
    
    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cB[-1,])),'.')
    plt.plot(n,np.fft.fftshift(np.abs(cB1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case B - Amplitude of Fourier Coefficients for difference values of Nt')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('Nt=801','Nt=200'))

    #Nx
    #caseA
    gA = nwave(1-2j,1+2j,Nx=256)
    gA1 = nwave(1-2j,1+2j,Nx=100)
    cA = np.fft.fft(gA[50:,])/100
    cA1 = np.fft.fft(gA1[50:,])/100
    n2 = np.arange(-100/2,100/2)

    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cA[-1,])),'.')
    plt.plot(n2,np.fft.fftshift(np.abs(cA1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case A - Amplitude of Fourier Coefficients for difference values of Nx')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('Nx=256','Nx=100'))
    
    #caseB
    gA = nwave(1-1j,1+2j,Nx=256)
    gA1 = nwave(1-1j,1+2j,Nx=100)
    cA = np.fft.fft(gA[50:,])/100
    cA1 = np.fft.fft(gA1[50:,])/100
    n2 = np.arange(-100/2,100/2)

    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cA[-1,])),'.')
    plt.plot(n2,np.fft.fftshift(np.abs(cA1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case B - Amplitude of Fourier Coefficients for difference values of Nx')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('Nx=256','Nx=100'))

    #T
    #caseA
    gA = nwave(1-2j,1+2j,T=200)
    gA1 = nwave(1-2j,1+2j,T=100)
    cA = np.fft.fft(gA[50:,])/Nx
    cA1 = np.fft.fft(gA1[50:,])/Nx

    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cA[-1,])),'.')
    plt.plot(n,np.fft.fftshift(np.abs(cA1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case A - Amplitude of Fourier Coefficients for difference values of T')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('T=200','T=100'))
    
    #caseB
    gB = nwave(1-1j,1+2j,T=200)
    gB1 = nwave(1-1j,1+2j,T=100)
    cB = np.fft.fft(gA[50:,])/Nx
    cB1 = np.fft.fft(gA1[50:,])/Nx

    plt.figure()
    plt.plot(n,np.fft.fftshift(np.abs(cB[-1,])),'.')
    plt.plot(n,np.fft.fftshift(np.abs(cB1[-1,])),'.')
    plt.yscale("log")
    plt.title('Jimmy Yeung \n analyze() \n Case B - Amplitude of Fourier Coefficients for difference values of T')
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend(('T=200','T=100'))


    return None


def wavediff():
    """
    Question 1.3
    Add input/output as needed

    Discussion:    
    
    In fig7, I have plotted the amplitude of the coefficients for the fourier 
    transform method and finite difference method. As we can see from the figure,
    the amplitudes of the coefficients are very similar for both methods.

    Looking at fig8, we can see that the absolute difference between the 
    amplitudes of fourier and FD are quite small throughout but are greatest at the edges.
    The FD method uses the one-sided 4th order scheme when calculating
    at the boundaries, but it does not seem as accurate when
    calculating values not on the boundary. 
    We can see that as Nx increases then the absolute difference between fourier and FD decreases.
    
    In fig9, I have plotted the time against Nx for Nx=256,512,1024. I have
    averaged the time over 10 iterations. We can see that the Fourier 
    method is faster than the FD method for all the values of Nx I have plotted. 
    
    In conclusion, DFT is a better method than FD as it is more accurate on the
    boundaries and it is faster.
    
    """
    def fourier_diff(g,Nx=256):
        c = np.fft.fft(g)/Nx
        n = np.fft.fftshift(np.arange(-Nx/2,Nx/2))
        k = 2j*np.pi*n/L
        dg = Nx*np.fft.ifft(k*c)
        return dg
    

    def FD_diff(g,Nx=256):
        alpha = 3/8
        a=25/16
        b=1/5
        c=-1/80
        N = g.size
        
    
        A = np.ones((3,N))
        A = np.array([[alpha,1,alpha]]).T*A
        A[0,0] = 1
        A[0,1] = 3
        A[-1,-1] = 1
        A[-1,-2] = 3
        
    
        hinv = 1/h
        zv = np.ones(N)
        zgv = np.zeros(N)
        zgv[0] = (-17/6)*hinv
        zgv[-1] = (17/6)*hinv
        agv = a*zv[1:]*hinv/2
        agv[0] = (3/2)*hinv
        bgv = b*zv[2:]*hinv/4
        bgv[0] = (3/2)*hinv
        cgv = c*zv[3:]*hinv/6
        cgv[0] = (-1/6)*hinv
        
        magv = -a*zv[1:]*hinv/2
        magv[-1] = (-3/2)*hinv
        mbgv = -b*zv[2:]*hinv/4
        mbgv[-1] = (-3/2)*hinv
        mcgv = -c*zv[3:]*hinv/6
        mcgv[-1] = (1/6)*hinv
        
        b = scipy.sparse.diags([[b*hinv/4, 0], [c*hinv/6, c*hinv/6, 0], mcgv,mbgv,magv,zgv,agv,bgv,cgv, [0, -c*hinv/6, -c*hinv/6], [0, -b*hinv/4]], [-N+2, -N+3, -3, -2, -1, 0, 1, 2, 3, N-3, N-2])
        
        B = b @ g
        dg = scipy.linalg.solve_banded((1,1),A, B)
        
        return dg
    

    L = 100
    Nt = 801
    T = 100
    
    Nx=256
    g = nwave(1-1j,1+2j, Nx=Nx, Nt=Nt, T=T)[-1]
    x = np.linspace(0,L,Nx+1)
    x= x[:-1]
    h = x[1] - x[0]
    dg = fourier_diff(g)
    dg_fdd = FD_diff(g)
    
    plt.figure()
    plt.plot(x,np.abs(dg),"rx--")
    plt.plot(x,np.abs(dg_fdd),"bx")
    plt.yscale("log")
    plt.title("Jimmy Yeung \n wavediff() \n Fourier vs Finite Difference")
    plt.xlabel('x')
    plt.ylabel('Amplitude')
    plt.legend(('Fourier','Finite Difference'))
    
    
    Nxarray = [256, 512, 1024]
    g0 = nwave(1 - 1j, 1 + 2j, Nx=Nxarray[0], Nt=Nt, T=T)[800]
    g1 = nwave(1 - 1j, 1 + 2j, Nx=Nxarray[1], Nt=Nt, T=T)[800]
    g2 = nwave(1 - 1j, 1 + 2j, Nx=Nxarray[2], Nt=Nt, T=T)[800]
    Glist = [g0, g1, g2]
    plt.figure()
    for i, j in enumerate(Nxarray):
        x = np.linspace(0, L, j+1)
        x = x[:-1]
        h = x[1] - x[0]
        dg1 = fourier_diff(Glist[i], Nx=j)
        dg2 = FD_diff(Glist[i], Nx=j)
        plt.plot(x, np.abs(dg1 - dg2))
    plt.title("Jimmy Yeung \n wavediff() \n Absolute difference between Amplitudes of Fourier and FD for different Nx")
    plt.xlabel("x")
    plt.ylabel("Absolute difference")
    plt.yscale("log")
    plt.legend(('Nx = 256', 'Nx = 512', 'Nx = 1024'))
    plt.grid()
    plt.show()
    
    
    fouriertimes=[]
    FDtimes=[]
    for Nx in [4,8,16,32,64,128,256,512,1024]:
        g = nwave(1-1j,1+2j, Nx=Nx, Nt=Nt, T=T)[-1]
        x = np.linspace(0,L,Nx+1)
        x= x[:-1]
        h = x[1] - x[0]  
        
        t1 = time.perf_counter()
        for i in range(10):
            dg = fourier_diff(g,Nx)
        t2 = time.perf_counter()
        fouriertimes.append((t2-t1)/10)
        
        t3 = time.perf_counter()
        for i in range(10):
            dg_fdd = FD_diff(g,Nx)
        t4 = time.perf_counter()
        FDtimes.append((t4-t3)/10)
    
    plt.figure()
    xaxis = [4,8,16,32,64,128,256,512,1024]
    plt.loglog(xaxis,fouriertimes)
    plt.loglog(xaxis,FDtimes)
    plt.title("Jimmy Yeung \n wavediff() \n Times for Fourier vs FD for different Nx")
    plt.xlabel("Nx")
    plt.ylabel("Time")
    plt.legend(('Fourier', 'FD'))
    plt.grid()
    plt.show()
    
    
    return None

if __name__=='__main__':
    analyze()
    wavediff()
