"""M345SC Homework 3, part 2
Jimmy Yeung 01203261
"""
import numpy as np
import networkx as nx
from scipy.linalg import expm
import scipy.linalg as linalg


def growth1(G,params=(0.02,6,1,0.1,0.1),T=6):
    """
    Question 2.1
    Find maximum possible growth, G=e(t=T)/e(t=0) and corresponding initial
    condition.

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth
    y: 3N-element array containing computed initial condition where N is the
    number of nodes in G and y[:N], y[N:2*N], y[2*N:] correspond to S,I,V

    Discussion:
        
    The model can be expressed as dy/dt=My where M is a matrix and the solution
    is given by y = y0*exp(Mt). We can see that e(t)=y.T*y and so the
    problem is to maximise |y0*exp(Mt)|^2/|y0|^2 = |exp(Mt)|^2. The maximal growth
    corresponds to the largest eigenvalue and the initial conidtion is the corresponding 
    eigenvector. We find this using SVD.

    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------
    G=0
    y = np.zeros(3*N)

    #Add code here  
    zeros = np.zeros((N,N))
    M = np.block([[F-(tau+g+k)*np.eye(N), a*np.eye(N), zeros],
                 [theta*np.eye(N), F-(tau+k+a)*np.eye(N), zeros],
                 [-theta*np.eye(N),zeros,F+(k - tau)*np.eye(N)]])
  
    A = expm(M*T)
    u,s,vh = np.linalg.svd(A)
    G = s[0]**2
    y = vh[0]

    return G,y

def growth2(G,params=(0.02,6,1,0.1,0.1),T=6):
    """
    Question 2.2
    Find maximum possible growth, G=sum(Ii^2)/e(t=0) and corresponding initial
    condition.

    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth
    y: 3N-element array containing computed initial condition where N is the
    number of nodes in G and y[:N], y[N:2*N], y[2*N:] correspond to S,I,V

    Discussion: 
    
    This is similar to growth1 but we only want the maximal growth of I. The 
    maximal growth of I can be done by indexing rows N to 2N in the matrix
    A = exp(Mt) and solving using SVD to find the solution to find the largest 
    eigenvalue and eigenvector, similar to growth1.
    
    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------
    G=0
    y = np.zeros(3*N)

    #Add code here
    zeros = np.zeros((N,N))
    M = np.block([[F-(tau+g+k)*np.eye(N), a*np.eye(N), zeros],
                 [theta*np.eye(N), F-(tau+k+a)*np.eye(N), zeros],
                 [-theta*np.eye(N),zeros,F+(k - tau)*np.eye(N)]])
        
    A = expm(M*T)
    B = A[N:2*N,:]
    u,s,vh = np.linalg.svd(B)
    y = vh[0]
    G = s[0]**2
    
    return G,y


def growth3(G,params=(2,2.8,1,1.0,0.5),T=6):
    """
    Question 2.3
    Find maximum possible growth, G=sum(Si Vi)/e(t=0)
    Input:
    G: Networkx graph
    params: contains model parameters, see code below.
    T: time for energy growth

    Output:
    G: Maximum growth

    Discussion:
        
    First, we index rows 0 to N and 2N to 3N which correspond to S and V respectively.
    We want to construct a symmetric matrix.
    
    Let B1*y0 = S and B2*y0 = V.
    
    SV
    = SV / 2 + SV / 2
    = (B1*y0).T(B2*x) / 2 + (B2*x).T(B1*x) / 2
    = y0.T*B1.T*B2*y0 / 2 + y0.T*B2.T*B1*y0 / 2
    = y0.T(0.5B1.T*B2 + 0.5B2.T*B1)y0
    
    Then 0.5B1.T*B2 + 0.5*B2.T*B1 is our symmetric matrix which I have called M1.
    To find the maximal growth we find the largest eigenvalue of this symmetric matrix.
    
    """
    a,theta,g,k,tau=params
    N = G.number_of_nodes()

    #Construct flux matrix (use/modify as needed)
    Q = [t[1] for t in G.degree()]

    Pden = np.zeros(N)
    Pden_total = np.zeros(N)
    for j in G.nodes():
        for m in G.adj[j].keys():
            Pden[j] += Q[m]
    Pden = 1/Pden
    Q = tau*np.array(Q)
    F = nx.adjacency_matrix(G).toarray()
    F = F*np.outer(Q,Pden)
    #-------------------------------------

    #Add code here
    zeros = np.zeros((N,N))
    M = np.block([[F-(tau+g+k)*np.eye(N), a*np.eye(N), zeros],
                 [theta*np.eye(N), F-(tau+k+a)*np.eye(N), zeros],
                 [-theta*np.eye(N),zeros,F+(k - tau)*np.eye(N)]])
    A = expm(M*T)
    
    B1 = A[:N,:]
    B2 = A[2*N:,:]

    M1 = (np.matmul(B1.T, B2) + np.matmul(B2.T, B1))/2
    eigens = linalg.eigvalsh(M1)
    G = max(eigens)   
    
    return G


def Inew(D):
    """
    Question 2.4

    Input:
    D: N x M array, each column contains I for an N-node network

    Output:
    I: N-element array, approximation to D containing "large-variance"
    behavior

    Discussion: 
      
    I implemented the PCA algorithm on D. The 
    Variance of D can be found by summing the variance of the columns of the 
    readjusted matrix D - mean(Dcolumn).
    We can find the vector I that estimates variance in D by performing SVD 
    on D2 and finding the largest eigenvector which corresponds to the largest
    eigenvalue. I then multiply this by D2 to give I.
    The variance of I is 0.32. The actual variance of D is 0.40. 
        
    """
    D = np.loadtxt('q22test.txt')
    N,M = D.shape

    # Scale D to subtract mean of each corresponding column
    D2 = D - np.outer(np.ones(N), np.mean(D, axis=0))

    #Compute SVD
    U, S, VH = linalg.svd(D2)
    I = np.dot(D2, VH[0])
    return I


if __name__=='__main__':
    N,M = 100,5
    G = nx.barabasi_albert_graph(N,M,seed=1)
    G1, y1 = growth1(G,params=(0.02,6,1,0.1,0.1),T=6)
    G2, y2 = growth2(G,params=(0.02,6,1,0.1,0.1),T=6)
    G3 = growth3(G,params=(2,2.8,1,1.0,0.5),T=6)
    D = np.loadtxt("q22test.txt")
    I = Inew(D)
