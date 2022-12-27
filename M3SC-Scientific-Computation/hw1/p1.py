"""M345SC Homework 1, part 1
Jimmy Yeung 01203261
"""

def ksearch(S,k,f,x):
    """
    Search for frequently-occurring k-mers within a DNA sequence
    and find the number of point-x mutations for each frequently
    occurring k-mer.
    Input:
    S: A string consisting of A,C,G, and Ts
    k: The size of k-mer to search for (an integer)
    f: frequency parameter -- the search should identify k-mers which
    occur at least f times in S
    x: the location within each frequently-occurring k-mer where
    point-x mutations will differ.

    Output:
    L1: list containing the strings corresponding to the frequently-occurring k-mers
    L2: list containing the locations of the frequent k-mers
    L3: list containing the number of point-x mutations for each k-mer in L1.
    
    Discussion:
    My algorithm goes through the string S and splits it up into k-mers. If the k-mer is not in the dictionary,
    then add the k-mer to the dictionary as a key and set the value of the k-mer to be a list containing its position
    within S. If the k-mer is already in the dictionary then append the position to the list. 
    
    If the length of the list of positions of the k-mer is equal to the frequency parameter f,
    we append the k-mer into L1  (No need to check if the number of elements
    in the list is greater than f, as it will already be in L1 if it is). Similarly, append the locations 
    of the frequent k-mer into L2.
    
    I then created a function which returns the x-mutation of the inputed k-mer. Looping through L1, if 
    the x-mutations of the frequent k-mers are in the dictionary, we then append the lengths of the lists
    of positions to a list L, sum the list, and then append to L3.
    
    Let the length of S be n. The running time of the loop through S would
    be O(n-(k-1)) with k fixed. Let the length of L1 be m. The running time of the loop
    through L1 would be O(m). So O(n-(k-1)) + O(m) = O(n-(k-1)+m). As n >= m and k is constant, the leading
    order running time is O(n). The running time is proportional to the length of the string, n.
    """
    L1,L2,L3=[],[],[]    
    d = dict()
    
    for i in range(len(S)-(k-1)):
        key = S[i:i+k]
        if key in d:
            d[key].append(i)
        else:
            d[key] = [i]
        if len(d[key]) == f:
            L1.append(key)
            L2.append(d[key])
 
    for j in range(len(L1)):
        s1,s2,s3 = mutations(L1[j],x)
        L = []
        if s1 in d:
            L.append(len(d[s1]))
        if s2 in d:
            L.append(len(d[s2]))
        if s3 in d:
             L.append(len(d[s3]))
        L3.append(sum(L))     
     
    return L1,L2,L3

def mutations(s,x):
    
    s1 = list(s)
    s2 = s1.copy()
    s3 = s1.copy()
    if s[x] == "A":
        s1[x] = "C"
        s2[x] = "T"
        s3[x] = "G"
    if s[x] == "C":
        s1[x] = "A"
        s2[x] = "T"
        s3[x] = "G"
    if s[x] == "G":
        s1[x] = "C"
        s2[x] = "A"
        s3[x] = "T" 
    if s[x] == "T":
        s1[x] = "C"
        s2[x] = "A"
        s3[x] = "G"
    s1 = "".join(s1)
    s2 = "".join(s2)
    s3 = "".join(s3)
   
    return s1,s2,s3

if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    k=3
    x=2
    f=2
    L1,L2,L3=ksearch(S,k,f,x)

