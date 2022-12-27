"""M345SC Homework 2, part 2
Jimmy Yeung 01203261
"""
import numpy as np
import networkx as nx
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def model1(G,x=0,params=(50,80,105,71,1,0),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf (see code below)
    display: A plot of S(t) for the infected node is generated when true

    x: node which is initially infected

    Output:
    S: Array containing S(t) for infected node
    """
    a,theta0,theta1,g,k,tau=params
    tarray = np.linspace(0,tf,Nt+1)
    S = np.zeros(Nt+1) 

    #Add code here
    S0 = 0.05
    I0 = 0.05
    V0 = 0.1
    
    y0 = [S0,I0,V0]
    
    def RHS(y,t):
        S = y[0]
        I = y[1]
        V = y[2]

        theta = theta0 + theta1*(1-np.sin(2*np.pi*t))
        dSdt = a*I - (g+k)*S
        dIdt = theta*S*V  - (k+a)*I
        dVdt = k*(1-V) - theta*S*V
        
        dydt = [dSdt,dIdt,dVdt]
        
        return dydt
    
    y = odeint(RHS,y0,tarray)
    S = y[:,0]
    
    if display:
        plt.figure()
        plt.plot(tarray,S)
        plt.title("Jimmy Yeung model1 \n S(t) against t")
        plt.xlabel('t')
        plt.ylabel('S(t)')
        plt.grid()

    return S

def modelN(G,x=0,params=(50,80,105,71,1,0.01),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf (see code below)
    display: A plot of S(t) for the infected node is generated when true

    x: node which is initially infected

    Output:
    Smean,Svar: Array containing mean and variance of S across network nodes at
                each time step.
    """
    a,theta0,theta1,g,k,tau=params
    tarray = np.linspace(0,tf,Nt+1)
    Smean = np.zeros(Nt+1)
    Svar = np.zeros(Nt+1)
    
    #Add code here
    N = nx.number_of_nodes(G)
    S0 = np.zeros(N)
    I0 = np.zeros(N)
    V0 = np.ones(N)
    S0[x] = 0.05
    I0[x] = 0.05
    V0[x] = 0.1 
    
    y0 = np.concatenate([S0,I0,V0])
    
    A = nx.adjacency_matrix(G)
    q = np.array([j for i,j in G.degree()])
    qA = np.diag(q)*A 
    sum_qA = np.diag(1/(q*A)) 
    F = tau*np.matmul(qA,sum_qA)                              
 
    def RHS(y,t):
        """Compute RHS of model at time t
        input: y should be a 3N x 1 array containing with
        y[:N],y[N:2*N],y[2*N:3*N] corresponding to
        S on nodes 0 to N-1, I on nodes 0 to N-1, and
        V on nodes 0 to N-1, respectively.
        output: dy: also a 3N x 1 array corresponding to dy/dt

        Discussion:
            
        Creating S from y would take - N operations.
            
        a*I is a scalar multiplied by a N-d array - N operations.
        
        (g+k)*S is a one addition then multiplied by a N-d array - N+1 operations.
        
        matmul(F,S) is taking the dot product of the rows in F with S, each
            dot product is N multiplications and N-1 additions. We do this N
            times - N(N+(N-1)) operations.
        
        np.sum(axis=0) makes N-1 additions, N times - N(N-1).
        
        sum(F,axis=0)*S makes N multiplications - N 
        
        adding these together - 3 operations.
        
        So the total number of operations is 
        
        N + N + (N+1) + N(N+(N-1)) + N(N-1) + N + 3
        
        = 3N^2 + 2N + 4
        
        Hence the time complexity is bounded by O(N^2).
        
        """
        S = y[:N]
        I = y[N:2*N]
        V = y[2*N:]
        
        theta = theta0 + theta1*(1-np.sin(2*np.pi*t))
        
        dSdt = a*I - (g+k)*S + np.matmul(F,S) - np.sum(F,axis = 0)*S           
        dIdt = theta*S*V - (k+a)*I + np.matmul(F,I) - np.sum(F,axis = 0)*I   
        dVdt = k*(1- V) - theta*S*V + np.matmul(F,V) - np.sum(F,axis = 0)*V                
        
        dydt = np.concatenate([dSdt,dIdt,dVdt])
        
        return dydt

    y = odeint(RHS,y0,tarray)
    S = y[:,:N]
    Smean = np.mean(S,axis = 1)
    Svar = np.var(S,axis = 1)
    
    if display:
        plt.figure()
        plt.plot(tarray,Smean)
        plt.title("Jimmy Yeung modelN \n Mean of S against t")
        plt.xlabel('t')
        plt.ylabel('<S>')
        plt.grid()
        
        plt.figure()
        plt.plot(tarray,Svar)
        plt.title("Jimmy Yeung modelN \n Variance of S against t")
        plt.xlabel('t')
        plt.ylabel('Var(S)')
        plt.grid()

    return Smean,Svar

def modeldiff(G,D,x=0,tf=6,Nt=40):
    
    tarray = np.linspace(0,tf,Nt+1)
    
    N = nx.number_of_nodes(G)
    S0 = np.zeros(N)
    I0 = np.zeros(N)
    V0 = np.ones(N)
    S0[x] = 0.05
    I0[x] = 0.05
    V0[x] = 0.1 
    
    y0 = np.concatenate([S0,I0,V0])
    
    Laplacian = nx.laplacian_matrix(G)
    
    def RHS_lap(y,t):
        S = y[:N]
        I = y[N:2*N]
        V = y[2*N:]
        
        dSdt = -D*Laplacian*S
        dIdt = -D*Laplacian*I
        dVdt = -D*Laplacian*V
        
        dydt = np.concatenate([dSdt,dIdt,dVdt])
        
        return dydt
    
    y = odeint(RHS_lap,y0,tarray)
    S_diff = y[:,:N]
    Smean_diff = np.mean(S_diff,axis = 1)
    Svar_diff = np.var(S_diff,axis = 1)

    return Smean_diff, Svar_diff



def diffusion(G,x=0,tf=30,Nt=1000):
    """Analyze similarities and differences
    between simplified infection model and linear diffusion on
    Barabasi-Albert networks.
    Modify input and output as needed.

    Discussion: 
        
    The mean is fixed at 0.0005 for both modeldiff and modelN for different tested values of tau and D. 
    This is our starting condition of S0 = 0.05 spread out over 100 nodes. 
    With the conditions set in the question the equation for S becomes:
        
        dS_i/dt = sum_j(Fij*Sj - Fji*Si)
        
    Taking the mean:
        
        d<S>/dt = sum_i(sum_j(Fij*Sj - Fji*Si)) / N
                = sum_i(sum_j(Fij*Sj)) - sum_i(sum_j(Fji*Si)) / N
                = sum_j(sum_i(Fij*Sj)) - sum_i(sum_j(Fji*Si)) / N
                = 0
       
    This shows that the mean does not change over time.     
        
    The variance of S will decrease as time goes on. This is because the spreaders will
    spread across the nodes as time progresses. Hence the variations of S across the nodes will be less
    and hence the variance will decrease. 
    Tau and D determine the rates at which S spreads. Larger tau and D means that S
    will spread more rapidly and hence the variance will decrease more rapidly.
    The variance tends to zero as t tends to infinity.
   
    """
    tarray = np.linspace(0,tf,Nt+1)
    
    Dlist = [0.01,0.1,0.5,0.8]
    
    plt.figure()
    for D in Dlist:
        Smean_diff, Svar_diff = modeldiff(G,D,x,tf,Nt)
        plt.plot(tarray,Smean_diff,label=D)
    plt.title("Jimmy Yeung \n modeldiff \n Mean for different values of D ")
    plt.xlabel('t')
    plt.ylabel('<S>')
    plt.legend()
    plt.grid()     
    
    plt.figure()
    for tau in Dlist:
        Smean_modelN, Svar_modelN = modelN(G,x,(0,80,0,0,0,tau),tf,Nt)
        plt.plot(tarray,Smean_modelN,label=tau)
    plt.title("Jimmy Yeung \n modelN \n Mean for different values of tau ")
    plt.xlabel('t')
    plt.ylabel('<S>')
    plt.legend()
    plt.grid()    
    
    plt.figure()
    for D in Dlist:
        Smean_diff, Svar_diff = modeldiff(G,D,x,tf,Nt)
        plt.plot(tarray,Svar_diff, label=D)
    plt.title("Jimmy Yeung \n modeldiff \n Variance for different values of D")
    plt.xlabel('t')
    plt.ylabel('Var(S)')
    plt.legend()
    plt.grid()
    
    plt.figure()
    for tau in Dlist:
        Smean_modelN, Svar_modelN = modelN(G,x,(0,80,0,0,0,tau),tf,Nt)
        plt.plot(tarray,Svar_modelN, label=tau)
    plt.title("Jimmy Yeung \n modelN \n Variance for different values of tau")
    plt.xlabel('t')
    plt.ylabel('Var(S)')
    plt.legend()
    plt.grid()
    
    return None 


if __name__=='__main__':

    G = nx.barabasi_albert_graph(100, 5)
    
    out = model1(G,display=True)
    out = modelN(G,display=True)
    
    out = diffusion(G) 
