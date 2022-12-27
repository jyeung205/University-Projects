"""M3C 2018 Homework 1 Jimmy Yeung 01203261
"""
import numpy as np
import matplotlib.pyplot as plt
from random import choices

#creating a function to calculate the points of each M village - this will be called in simulate1    
def function(A,c,d):
    f = (c**A)*(d**(1-A))
    return f

def simulate1(N,Nt,b,e):
    """Simulate C vs. M competition on N x N grid over
    Nt generations. b and e are model parameters
    to be used in fitness calculations
    Output: S: Status of each gridpoint at tend of somulation, 0=M, 1=C
    fc: fraction of villages which are C at all Nt+1 times
    Do not modify input or return statement without instructor's permission.
    """

    #Set initial condition
    S  = np.ones((N,N),dtype=int) #Status of each gridpoint: 0=M, 1=C
    j = int((N-1)/2)
    S[j-1:j+2,j-1:j+2] = 0

    fc = np.zeros(Nt+1) #Fraction of points which are C
    fc[0] = S.sum()/(N*N)

    #creating a matrix with containings the number of neighbours for each village
    numberofneighbours = 8*np.ones((N,N))
    numberofneighbours[0,:] = 5
    numberofneighbours[N-1,:] = 5
    numberofneighbours[:,0] = 5
    numberofneighbours[:,N-1] = 5
    numberofneighbours[0,0] = 3
    numberofneighbours[0,N-1] = 3
    numberofneighbours[N-1,0] = 3
    numberofneighbours[N-1,N-1] = 3

    #initialising matrices used later
    S1 = np.zeros((N,N))
    S2 = np.zeros((N,N))
    S3 = np.zeros((N,N))
    S4 = np.zeros((N,N))
    S5 = np.zeros((N,N))
    S6 = np.zeros((N,N))
    S7 = np.zeros((N,N))     
    S8 = np.zeros((N,N))

    C1 = np.zeros((N,N))
    C2 = np.zeros((N,N))
    C3 = np.zeros((N,N))
    C4 = np.zeros((N,N))
    C5 = np.zeros((N,N))
    C6 = np.zeros((N,N))
    C7 = np.zeros((N,N))
    C8 = np.zeros((N,N))
     
    M1 = np.zeros((N,N))
    M2 = np.zeros((N,N))
    M3 = np.zeros((N,N))
    M4 = np.zeros((N,N))
    M5 = np.zeros((N,N))
    M6 = np.zeros((N,N))
    M7 = np.zeros((N,N))
    M8 = np.zeros((N,N))
    
    #iterating over Nt years
    for x in range(Nt):

        #creating matrices to calculate the points of each C village
        S1[1:N,1:N] = S[0:N-1,0:N-1]
        S2[1:N,0:N] = S[0:N-1,0:N]
        S3[1:N,0:N-1] = S[0:N-1,1:N]
        S4[0:N,1:N] = S[0:N,0:N-1]
        S5[0:N,0:N-1] = S[0:N,1:N]
        S6[0:N-1,1:N] = S[1:N,0:N-1]
        S7[0:N-1,0:N] = S[1:N,0:N]
        S8[0:N-1,0:N-1] = S[1:N,1:N]
        
        Cpoints = (S1+S2+S3+S4+S5+S6+S7+S8)*S


        #creating matrices to calculate the points of each M village
        M1 = function(S1,b,e)
        M1[0,:] = 0 
        M1[:,0] = 0

        M2 = function(S2,b,e)
        M2[0,:] = 0 

        M3 = function(S3,b,e)
        M3[0,:] = 0 
        M3[:,N-1] = 0

        M4 = function(S4,b,e) 
        M4[:,0] = 0

        M5 = function(S5,b,e)
        M5[:,N-1] = 0

        M6 = function(S6,b,e)
        M6[:,0] = 0 
        M6[N-1,:] = 0

        M7 = function(S7,b,e)
        M7[N-1,:] = 0

        M8 = function(S8,b,e)
        M6[:,N-1] = 0 
        M6[N-1,:] = 0

        Mpoints = (M1+M2+M3+M4+M5+M6+M7+M8)*(np.ones((N,N))-S)


        #calculating the fitness of each C and M village
        Cfitness = Cpoints / numberofneighbours
        Mfitness = Mpoints / numberofneighbours


        #creating matrices to calculate the total fitness of the community for each village
        C1[1:N,1:N] = Cfitness[0:N-1,0:N-1]
        C2[1:N,0:N] = Cfitness[0:N-1,0:N]
        C3[1:N,0:N-1] = Cfitness[0:N-1,1:N]
        C4[0:N,1:N] = Cfitness[0:N,0:N-1]
        C5[0:N,0:N-1] = Cfitness[0:N,1:N]
        C6[0:N-1,1:N] = Cfitness[1:N,0:N-1]
        C7[0:N-1,0:N] = Cfitness[1:N,0:N]
        C8[0:N-1,0:N-1] = Cfitness[1:N,1:N]
        
        Ccommfitness = Cfitness+C1+C2+C3+C4+C5+C6+C7+C8

        
        M1[1:N,1:N] = Mfitness[0:N-1,0:N-1]
        M2[1:N,0:N] = Mfitness[0:N-1,0:N]
        M3[1:N,0:N-1] = Mfitness[0:N-1,1:N]
        M4[0:N,1:N] = Mfitness[0:N,0:N-1]
        M5[0:N,0:N-1] = Mfitness[0:N,1:N]
        M6[0:N-1,1:N] = Mfitness[1:N,0:N-1]
        M7[0:N-1,0:N] = Mfitness[1:N,0:N]
        M8[0:N-1,0:N-1] = Mfitness[1:N,1:N]
        
        Mcommfitness = Mfitness+M1+M2+M3+M4+M5+M6+M7+M8


        totalfitness = Mcommfitness + Ccommfitness


        #calculating the probability of a village being C in the following year
        Cprobability = Ccommfitness / totalfitness

        #setting up S for the following year
        for i in range(N):
            for j in range(N):
                p = Cprobability[i][j]
                S[i][j] = choices([0,1],[1-p,p])[0]

        fc[x+1] = S.sum()/(N*N)                

    return S,fc

def plot_S(S):
    """Simple function to create plot from input S matrix
    """
    ind_s0 = np.where(S==0) #C locations
    ind_s1 = np.where(S==1) #M locations
    plt.plot(ind_s0[1],ind_s0[0],'rs')
    plt.hold(True)
    plt.plot(ind_s1[1],ind_s1[0],'bs')
    plt.hold(False)
    plt.show()
    plt.pause(0.05)
    return None


def simulate2(N,Nt,b,e):
    """simulate2 calculates an average of the fc arrays produced by simulate1 for each given b
    """
    fczeros = np.zeros(Nt+1)
   
    for a in range(100):
         _, fc = simulate1(N,Nt,b,e)
         fcsum = fczeros + fc
         fczeros = fcsum
         average = fczeros/100
    return average


def analyze():
    """analyze plots graphs of fc against Nt for different ranges of b"""
        
    #b from 1.05 to 1.25
    plt.figure(figsize=(12,4))
    for b in np.linspace(1.05,1.25,5):
        x = simulate2(21,300,b,0.01)
        plt.plot(range(300+1), x*100, label=b)
        plt.legend()
    plt.ylabel('% of villages which are C')
    plt.xlabel('Number of years past')
    plt.title('Evolution for the % of villages which are C for b from 1.05 to 1.25')
    plt.show()
        
    """
    For b = 1.05, the rate of which M gains territory is relatively constant over the 300 years.
    
    For b = 1.1, the rate at which M gains territory is relatively constant through years 50 and 200.
    The graph levels off at around year 200.
    
    For b = 1.15, the rate at which M gains territory is fastest between years 40 and 100 but begins to slow
    down after year 125. 
    M gains 100% territory around year 225.
    
    For b = 1.2, the rate at which M gains territory is quicker than for b = 1.15. M gains 100% of the territory 
    at around year 150.
    
    For b = 1.2, the rate at which M gains territory is quicker than for b = 1.2. M gains 100% of the territory 
    at around year 125.
    """
    
    #b from 1.3 to 1.5
    plt.figure(figsize=(12,4))
    for b in np.linspace(1.3,1.5,5):
        x = simulate2(21,300,b,0.01)
        plt.plot(range(300+1), x*100, label=b)
        plt.legend()
    plt.ylabel('% of villages which are C')
    plt.xlabel('Number of years past')
    plt.title('Evolution for the % of villages which are C for b from 1.30 to 1.50')
    plt.show()
    
    """
    M gains territory at a relatively quick rate for all values of b between 1.3 and 1.5.
    
    A larger b means M gains 100% of the territory in a shorter amount of time.
    
    M has gain 100% territory by 150 years
    """
    
    #b from 1.55 to 1.75
    plt.figure(figsize=(12,4))
    for b in np.linspace(1.55,1.75,5):
        x = simulate2(21,300,b,0.01)
        plt.plot(range(300+1), x*100, label=b)
        plt.legend()
    plt.ylabel('% of villages which are C')
    plt.xlabel('Number of years past')
    plt.title('Evolution for the % of villages which are C for b from 1.55 to 1.75')
    plt.show()
    
    """
    M gains territory at a relatively quick rate for all values of b between 1.55 and 1.75.
    M has gain 100% territory by 75 years.
    """
    
    #b from 1.8 to 2
    plt.figure(figsize=(12,4))
    for b in np.linspace(1.80,2,5):
        x = simulate2(21,300,b,0.01)
        plt.plot(range(300+1), x*100, label=b)
        plt.legend()
    plt.ylabel('% of villages which are C')
    plt.xlabel('Number of years past')
    plt.title('Evolution for the % of villages which are C for b from 1.80 to 2.0')
    plt.show()
    
    """
    There is no significant difference in the rate at which M gains territory for b in 1.8 to 2.
    Conclude that there is no significant difference in the rate at which M gains territory for
    b > 1.8.
    """
                


if __name__ == '__main__':
    #The code here should call analyze and generate the
    #figures that you are submitting with your code
    output = analyze()
