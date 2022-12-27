"""M345SC Homework 1, part 2
Jimmy Yeung 01203261
"""
import numpy as np
import matplotlib.pyplot as plt
import time

def bsearch(L,x):
    #Set initial start and end indices for full list
    istart = 0
    iend = len(L)-1
    #Iterate and contract "active" portion of list
    while istart<=iend:

        imid = int(0.5*(istart+iend))
        if x==L[imid]:
            return imid
        elif x < L[imid]:
            iend = imid-1
        else:
            istart = imid+1

    return 

def isort(X):
    """Sort X using insertion sort algorithm and return sorted array
    """
    S = X.copy()
    for i,x in enumerate(X[1:],1):
        #place x appropriately in partially sorted array, S
        for j in range(i-1,-1,-1):
            if S[j+1]<S[j]:
                S[j],S[j+1] = S[j+1],S[j]
            else:
                break
    return S

def randL(N,M,P):
    L = []
    for i in range(N):
        L1 = list(np.random.randint(0,100,P))
        #sort L1
        L1 = isort(L1)
        L2 = list(np.random.randint(0,100,M-P))
        L3 = L1+L2
        L.append(L3)
    return L

def nsearch(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """
    
    Lout=[]
    M = len(L[0])
    for i in range(len(L)):
        L1 = L[i][:P]
        L2 = L[i][P:M]
        #binary search through L1
        position1 = bsearch(L1,target)
        if position1 != None:
            Lout.append([i,position1])
            #check if target occurs more than once
            position1_right = position1 + 1
            while position1_right <= P-1 and L1[position1_right] == target:
                Lout.append([i,position1_right])
                position1_right +=1
                
            position1_left = position1 - 1
            while position1_left >= 0 and L1[position1_left] == target:
                Lout.append([i,position1_left])
                position1_left -=1
            
        for j in range(M-P):
            if L2[j] == target:
                Lout.append([i,j+P])    

    return Lout


def nsearch_time():
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion:
    My algorithm loops through the list L and splits the sublists into
    the sorted part and the unsorted part. I then use binary search to search
    for the target within the sorted list and append the index of the sublist
    and the location within the sublist to Lout. Next, I need to check if 
    the target is found to the left or to the right of the location found in
    the binary search, and similarly append the index of the sublist and the
    location of the target to Lout.
    
    To search for the target in the unsorted part, I use linear search.
    If target is found, append the index of the sublist
    and the location of the target to Lout.
    
    The running time of the loop through the N sublists would be O(N). 
    
    Within the sublists, if the target only appears once in the sorted part, 
    then the running time would be O(log2(P)). However, the worst case would 
    be if the sorted part only contained the target and then 
    the running time would be O(log2(P)+(P-1)).
    
    The running time of the linear search through the unsorted part would be O(M-P). 
    So we have O(log2(P)+(P-1)) + O(M-P) = O(log2(P)+(P-1)+(M-P)).
    However M-P is large and so M >> P, hence the leading order running time
    to search the sublists is O(M).
    
    So overall the running time of my algorithm is O(M) * O(N) = O(MN).
    
    My algorithm can be considered efficient if we compare it to the running times 
    of other methods. Merge sorting the sublists would have a running time O(Mlog2(M)).
    With N sublists, the running time would be O(NMlog2(m)), which is greater than
    O(NM).  If we were to sort and then use binary search the running time would
    be O(Mlog2(M)+log2(M)) = O(Mlog2(M)). With N sublists, the running time is O(NMlog2(M)),
    which is greater than O(NM).
    
    Looking at fig1, with M fixed, we can see that the running time increases linearly as N 
    increases. As M is constant, O(NM) becomes O(N), which is what the figure shows.
    
    Looking at fig2, with N fixed, we can see that the running time increases linearly as M 
    increases. As N is constant, O(NM) becomes O(M).  
    
    Similarly for fig3, N fixed but now P increasing proportionally to M. We can see that 
    the running time also increases linearly as M increases. As N is constant, O(NM) becomes O(M).  
    """
    #Varying N
    x_axis = list(range(1000,6000,1000))
    times = list()
    M = 1000
    P = 500
    target = 40
    for N in x_axis:
        L = randL(N,M,P)
        t0 = time.time()
        nsearch(L,P,target)
        t1 = time.time()
        total = t1-t0
        times.append(total)

    plt.figure()
    plt.plot(x_axis,times,label="Time against N")
    plt.title("Jimmy Yeung\n nsearch_time")
    plt.xlabel('N')
    plt.ylabel('Time')
    plt.legend()
    plt.grid()
    plt.show()
    
    #Varying M
    x_axis = list(range(1000,6000,1000))
    times = list()
    N = 1000
    P = 500
    target = 40
    for M in x_axis:
        L = randL(N,M,P)
        t0 = time.time()
        nsearch(L,P,target)
        t1 = time.time()
        total = t1-t0
        times.append(total)

    plt.figure()
    plt.plot(x_axis,times,label="Time against M")
    plt.title("Jimmy Yeung\n nsearch_time")
    plt.xlabel('M')
    plt.ylabel('Time')
    plt.legend()
    plt.grid()
    plt.show()

    #Increasing M and P proportionally
    x_axis = list(range(1000,6000,1000))
    times = list()
    N = 1000
    target = 40
    for M in x_axis:
        P = int(M/2)
        L = randL(N,M,P)
        t0 = time.time()
        nsearch(L,P,target)
        t1 = time.time()
        total = t1-t0
        times.append(total)

    plt.figure()
    plt.plot(x_axis,times,label="Time against M with P increasing proportionally to M")
    plt.title("Jimmy Yeung\n nsearch_time")
    plt.xlabel('M')
    plt.ylabel('Time')
    plt.legend()
    plt.grid()
    plt.show()

    return None #Modify as needed


if __name__ == '__main__':

    nsearch_time() 

