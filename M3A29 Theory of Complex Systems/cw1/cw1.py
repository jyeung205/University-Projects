import numpy as np
import matplotlib.pyplot as plt

def SOBP(p,n,t):
    d = dict()
    p_list=[p]
    N = 2**(n+1)-1
    for i in range(t):
        L=[]
        for j in range(n):
            L.append([])
        
        initial = np.random.choice([0,1],None,p=[1-p,p])
        if initial == 1:
            #two active sites to start
            L[0].append(1)
            L[0].append(1)
            for b in range(n-1):
                for a in L[b]:
                    if a == 1:
                        site_status = np.random.choice([0,1],None,p=[1-p,p])
                        L[b+1].append(site_status)
                        L[b+1].append(site_status)              
        sigma = sum(L[n-1])
        p = p + (1-sigma)/N
        p_list.append(p)
        
        s=1
        for i1 in range(n):
            s = s + sum(L[i1])    

        key = s
        if key in d:
            d[key] += 1
        else:
            d[key] = 1
        
    return p_list, d

#part f
x_axis = range(3001)
y1,_ = SOBP(0.1,10,3000)
y2,_ = SOBP(0.6,10,3000)
plt.figure()
plt.plot(x_axis,y1,label="p(0)=0.1")
plt.plot(x_axis,y2,label="p(0)=0.6")
plt.title("Value of p(t) for a system with n = 10 generations")
plt.xlabel('t')
plt.ylabel('p(t)')
plt.legend()
plt.grid()
plt.show()   


#part g
time=[]
grads=[]
for n in range(8,16):
    p_list,_ = SOBP(0.1,n,30000)
    grad = (p_list[10]-p_list[0])/11
    grads.append(grad)
    t = 0
    while 0.5 - p_list[t] > 0.01:
        t = t+1
    time.append(t)

x_axis = range(8,16)
plt.figure()
plt.plot(x_axis,time,label="p(0)=0.1")
plt.title("Time t* at which p(t*) ~ 0.5 against n")
plt.xlabel('n')
plt.ylabel('t*')
plt.legend()
plt.grid()
plt.show()

print(grads)


#part h
plt.figure()
for n in [16,20,24,28]:
    _,d = SOBP(0.5,n,1000000)
    
    sizes=[]
    for key in d.keys():
        sizes.append(d[key])
    dist = np.asarray(sizes)/sum(sizes) 
     
    plt.plot([s for s in range(3,100)], [np.sqrt(2/np.pi)*s**-1.5 for s in range(3,100)])
    plt.loglog(d.keys(),dist,'.',label = n)   

plt.title("Avalance distribution D(s) for different system sizes")
plt.xlabel('s')
plt.ylabel('D(s)')
plt.legend()
plt.grid()
plt.show()    