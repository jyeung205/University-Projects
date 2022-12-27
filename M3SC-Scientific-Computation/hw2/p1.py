"""M345SC Homework 2, part 1
Jimmy Yeung 01203261
"""

def count_days(task_id):
    days = dict()
    if task_id in days:
    	return days[task_id]
    if not L[task_id]:
        days[task_id]=0
        return 0
    md = max([count_days(t) for t in L[task_id]])
    days[task_id] = md + 1
    return days[task_id]


def scheduler(L):
    """
    Question 1.1
    Schedule tasks using dependency list provided as input

    Input:
    L: Dependency list for tasks. L contains N sub-lists, and L[i] is a sub-list
    containing integers (the sub-list my also be empty). An integer, j, in this
    sub-list indicates that task j must be completed before task i can be started.

    Output:
    L: A list of integers corresponding to the schedule of tasks. L[i] indicates
    the day on which task i should be carried out. Days are numbered starting
    from 0.

    Discussion: 
    
    The number of days for a task to be completed will be the maximum of the
    number of days of the tasks it depends on plus 1. 
    If we compute the days recursively then we don't need to compute a new tree
    for each individual node. My function count_days counts the number of days needed
    for the task inputed. It does this by storing the number of days
    of each task in a dictionary, and for each task it finds the max
    days out of its dependent tasks. 
    
    The time complexity is O(N) as you only need to check each task once.
    
    """
    S=[]
    
    for task_id in range(len(L)):
        day = count_days(task_id)
        S.append(day)
    
    print(S)
    return S


def findPath(A,a0,amin,J1,J2):
    """
    Question 1.2 i)
    Search for feasible path for successful propagation of signal
    from node J1 to J2

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element tuples (the
    sub-list my also be empty) of the form (j,Lij). The integer, j, indicates that there is a link
    between nodes i and j and Lij is the loss parameter for the link.

    a0: Initial amplitude of signal at node J1

    amin: If a>=amin when the signal reaches a junction, it is boosted to a0.
    Otherwise, the signal is discarded and has not successfully
    reached the junction.

    J1: Signal starts at node J1 with amplitude, a0
    J2: Function should determine if the signal can successfully reach node J2 from node J1

    Output:
    L: A list of integers corresponding to a feasible path from J1 to J2.

    Discussion: 
    
    Here, I have used breadth first search adjusted so that 
    only the paths where Lij*a0 >= amin are considered. L4 is a list of paths
    where the index of the sublist corresponds to node number. 
    
    This is a BFS algorithm so the time complexity is O(N+M), where M is the 
    number of edges.
    
    """

    L2 = [0 for l in range(len(A))] 
    L3 = [-1000 for l in range(len(A))] 
    L4 = [[] for l in range(len(A))] #paths

    Q=[]
    Q.append(J1)
    L2[J1]=1
    L3[J1]=0
    L4[J1]=[J1]
    
    while len(Q)>0:
        x = Q.pop(0) 
        for tuples in A[x]:
            node = tuples[0]
            weight = tuples[1]
            if L2[node]==0 and weight*a0 >= amin:
                Q.append(node) 
                L2[node]=1
                L3[node]=1+L3[x]
                L4[node].extend(L4[x])
                L4[node].append(node)

    L = L4[J2]
    print(L)
    return L


def a0min(A,amin,J1,J2):
    """
    Question 1.2 ii)
    Find minimum initial amplitude needed for signal to be able to
    successfully propagate from node J1 to J2 in network (defined by adjacency list, A)

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element tuples (the
    sub-list my also be empty) of the form (j,Lij). The integer, j, indicates that there is a link
    between nodes i and j and Lij is the loss parameter for the link.

    amin: Threshold for signal boost
    If a>=amin when the signal reaches a junction, it is boosted to a0.
    Otherwise, the signal is discarded and has not successfully
    reached the junction.

    J1: Signal starts at node J1 with amplitude, a0
    J2: Function should determine min(a0) needed so the signal can successfully
    reach node J2 from node J1

    Output:
    (a0min,L) a two element tuple containing:
    a0min: minimum initial amplitude needed for signal to successfully reach J2 from J1
    L: A list of integers corresponding to a feasible path from J1 to J2 with
    a0=a0min
    If no feasible path exists for any a0, return output as shown below.

    Discussion: 
        
    The problem is to find the path whose minimum weight along that path is the maximum
    amongst all the paths. My algorithm keeps track of all possible paths to J2 and the 
    smallest weight in the paths in queue. 
    It looks at the neighbouring nodes of the current node and if the neighbours 
    are not in the path and the current node is not the target node J2, then 
    the neighbour nodes are appended to the path and the smallest weight is updated
    if necessary. If the current node is the target node and the weight of the path
    is less than the weight of the output path, then the output path is updated.
    Once a path reaches the target node, J2, this path is no longer updated.
    This continues until the queue becomes empty. 
        
    """
    path = [-1,[]]
    
    Q = []
    oldpath = [1,[J1]]
    Q.append(oldpath)
    
    while len(Q)>0:
        
        oldpath = Q.pop(0)     
        currentnode = oldpath[1][-1]
        
        if currentnode==J2 and oldpath[0]>path[0]: 
            path = oldpath
        
        else:
            for tuples in A[currentnode]:
                node=tuples[0]
                weight=tuples[1]
                if node not in oldpath[1] and currentnode != J2:
                    newpath = [oldpath[0],oldpath[1].copy()]
                    newpath[1].append(node)
                    if weight < newpath[0]:
                        newpath[0] = weight
                    Q.append(newpath)
    
    output = -1,[]    
    if path[1]: 
        a0min = amin/path[0]
        output = a0min, path[1]
    print(output)
    return output


if __name__=='__main__':

    L = [[1,2,3,4,5,7,8,9],[3,4,5,7,8,9],[5,6,9,10],[4,5,7,8,9],[5,8,9],[9],[10],[],[],[],[]]
    S = scheduler(L)
    
    A = [[(1,0.5),(2,0.3),(3,0.4)],
          [(0,0.5),(4,0.1)],
          [(0,0.3),(4,0.2),(5,0.3)],
          [(0,0.4),(5,0.4)],
          [(1,0.1),(2,0.2),(6,0.3)],
          [(2,0.3),(3,0.4),(6,0.15)],
          [(4,0.3),(5,0.15)]]
    
    a0 = 1
    amin = 0.2
    J1 = 0
    J2 = 6
    
    L = findPath(A,a0,amin,J1,J2)
    
    a0, path = a0min(A,amin,J1,J2)