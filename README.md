# M3SC-Scientific-Computation

M3SC Scientific Computation - Mathematics BSc at Imperial College London, 2019.

Grade: 70 (First Class)

### Overview
- Hw1: Searching Algorithms and Complexity Analysis
- Hw2: Graphs and Dijkstra's algorithm
- Hw3: Solving PDEs

|![](https://github.com/jyeung205/University-Projects/blob/main/M3SC-Scientific-Computation/hw3/fig1.png)|![](https://github.com/jyeung205/University-Projects/blob/main/M3SC-Scientific-Computation/hw3/fig7.png)|
|:-----------------------:|:-------------------:|
|![](https://github.com/jyeung205/University-Projects/blob/main/M3SC-Scientific-Computation/hw3/fig8.png)|![](https://github.com/jyeung205/University-Projects/blob/main/M3SC-Scientific-Computation/hw3/fig9.png)|

<details>
<summary> Code Snippets </summary>

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import hann
import scipy
import time


def nwave(alpha,beta,Nx=256,Nt=801,T=200,display=False):

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

```
</details>

---

# M3C-High-Performance-Computing

M3C High Performance Computing - Mathematics BSc at Imperial College London, 2019.

Grade: 73 (First Class)

### Overview:

- Hw1: Python Arrays and Vectorisation
- Hw2: Implementing Neural Networks using Fortran90
- Hw3: Fortran90 OpenMP
- Hw4: Parallel Computing of Fortran90 code using OpenMP + MPI

|![](https://github.com/jyeung205/University-Projects/blob/main/M3C-High-Performance-Computing/hw1/hw11.png)|![](https://github.com/jyeung205/University-Projects/blob/main/M3C-High-Performance-Computing/hw2/hw22.png)|
|:-----------------------:|:-------------------:|
|![](https://github.com/jyeung205/University-Projects/blob/main/M3C-High-Performance-Computing/hw4/part2/p31.png)|![](https://github.com/jyeung205/University-Projects/blob/main/M3C-High-Performance-Computing/hw3/hw322.png)|

<details><summary>Code Snippets</summary>
<p>

```python
import numpy as np
import matplotlib.pyplot as plt
from m1 import bmodel as bm #assumes p2.f90 has been compiled with: f2py3 -c p2.f90 -m m1
import time
from scipy import optimize

def simulate_jacobi(n,input_num=(10000,1e-8),input_mod=(1,1,1,2,1.5),display=False):
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

```

</p>
</details>

---

# MA407-Algorithms-And-Computation

## Overview
- MA407 Algorithms and Computation - Applicable Mathematics MSc at London School of Economics, 2020.
- Task: To generate Latex code (which draws Binary Search Trees) using Java.
- Grade: 100% (Highest Score in Cohort)

## Example
Example of a Binary Search Tree drawn using using Latex. The Latex code was generated using Java.

![](https://github.com/jyeung205/University-Projects/blob/main/MA407-Algorithms-And-Computation/png/java%20Tree%2020%203-1.png?raw=true)

<details>
<summary> Code Snippets </summary>

```java

public class Graph {
    Node nodeListHead;
    Node nodeListEnd;
    Edge edgeListHead;
    Edge edgeListEnd;
    final int maxCoordinate = 100;
    final double textWidth = 17.0;

    public static class Node {
        int x;
        int y;
        int label;
        String value;
        Node next;

        Node(int x, int y, int label, String value) {
            this.x = x;
            this.y = y;
            this.label = label;
            this.value = value;
        }
    }

    public static class Edge {
        int nodeLabel1;
        int nodeLabel2;
        Edge next;
    }

    void addNode(int x, int y, int label, String value) {

        // Check to see if coordinates exceed maxCoordinate value
        if (Math.abs(x) > maxCoordinate && Math.abs(y) > maxCoordinate) {
            System.out.println("% coordinate " + x + " and " + y + " in node (" + x + "," + y + ") is too large");
            return;
        }
        if (Math.abs(x) > maxCoordinate) {
            System.out.println("% coordinate " + x + " in node (" + x + "," + y + ") is too large");
            return;
        }
        if (Math.abs(y) > maxCoordinate) {
            System.out.println("% coordinate " + y + " in node (" + x + "," + y + ") is too large");
            return;
        }

        // If the nodeList is empty then add a new node to the list
        if (nodeListHead == null) {
            nodeListHead = new Node(x, y, label, value);
            nodeListEnd = nodeListHead;
        }
        else {
            // Check for repeated coordinates and update value if coordinates are repeated
            for (Node node = nodeListHead; node != null; node = node.next) {
                if (node.x == x && node.y == y) {
                    node.value = value;
                    return;
                }
            }
            // Else create a new node and add to the list
            Node node = new Node(x, y, label, value);
            nodeListEnd.next = node;
            nodeListEnd = node;
        }
    }

    void addEdge(int x1, int y1, int x2, int y2) {

        // find label of x1, y1 and x2, y2
        Edge edge = new Edge();
        boolean node1Present = false;
        boolean node2Present = false;
        for (Node node = nodeListHead; node != null; node = node.next) {
            if (node.x == x1 && node.y == y1) {
                edge.nodeLabel1 = node.label;
                node1Present = true;
            }
            if (node.x == x2 && node.y == y2) {
                edge.nodeLabel2 = node.label;
                node2Present = true;
            }
        }

        // Check to see if nodes connected by the edge are present before adding them to the edgeList
        if (!node1Present && !node2Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but both nodes are not present in the graph");
            return;
        }
        if (!node1Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but node with coordinates (" + x1 + "," + y1 + ") is not present in the graph");
            return;
        }
        if (!node2Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but node with coordinates (" + x2 + "," + y2 + ") is not present in the graph");
            return;
        }

        if (edgeListHead == null) {
            edgeListHead = edge;
        }
        else {
            edgeListEnd.next = edge;
        }
        edgeListEnd = edge;

    }

    void clear() {
        nodeListHead = null;
        nodeListEnd = null;
        edgeListHead = null;
        edgeListEnd = null;
    }

    void outHeader() {
        double oddSideMargin = -1* (2.54 - (21.0 - textWidth) / 2);

        System.out.print("\\documentclass[a4paper,11pt]{article}\n" +
                "\\usepackage{mathpazo}\n" +
                "\\usepackage{tikz}\n" +
                "\\usetikzlibrary{shapes}\n");

        System.out.print("\\oddsidemargin "); System.out.printf("%.2f", oddSideMargin); System.out.print("cm\n");

        System.out.print("\\textwidth "); System.out.printf("%.2f", textWidth); System.out.print("cm\n");

        System.out.println("\\textheight 24cm\n" +
                "\\topmargin -1.3cm\n" +
                "\\parindent 0pt\n" +
                "\\parskip 1ex\n" +
                "\\pagestyle{empty}\n" +
                "\\begin{document}\n" +
                "\\medskip\\hrule\\medskip\n");

    }

    void outFooter() {
        System.out.println("\\medskip\\hrule\\medskip\n" +
                "\\end{document}");
    }

    void outGraph() {

        // Loop through nodeList and output LaTex command
        for (Node node = nodeListHead; node != null; node = node.next) {
            int x = node.x;
            int y = node.y;
            int label = node.label;
            String value = node.value;

            System.out.print( "\\draw [thick] "+ "(");
            System.out.print(x + "," + y);
            System.out.print(") node[draw, rounded rectangle] ");
            System.out.println("(" + label + ")" + " {" + value + "};");
        }

        // Loop through edgeList and output LaTex command
        for (Edge edge = edgeListHead; edge != null; edge = edge.next) {
            int label1 = edge.nodeLabel1;
            int label2 = edge.nodeLabel2;
            System.out.println("\\draw [->, thick] (" + label1 + ") to (" + label2 + ");");
        }

    }

    // Find the smallest and largest x, y coordinates and return an array containing them
    int[] getGridCoordinates() {

        int xMax = -100;
        int yMax = -100;
        int xMin = 100;
        int yMin = 100;
        for (Node node = nodeListHead; node != null; node = node.next) {
            if (node.x < xMin && node.x >= -maxCoordinate) {
                xMin = node.x;
            }
            if (node.x > xMax && node.x <= maxCoordinate) {
                xMax = node.x;
            }
            if (node.y < yMin && node.y >= -maxCoordinate) {
                yMin = node.y;
            }
            if (node.y > yMax && node.y <= maxCoordinate) {
                yMax = node.y;
            }
        }

        int[] gridCoordinates = new int[4];
        gridCoordinates[0] = xMin;
        gridCoordinates[1] = yMin;
        gridCoordinates[2] = xMax;
        gridCoordinates[3] = yMax;
        return gridCoordinates;
    }

    // Outputs the LaTex commands to draw a grid between the smallest and largest coordinates
    void outGrid(int[] gridCoordinates) {

        int xMin = gridCoordinates[0];
        int yMin = gridCoordinates[1];
        int xMax = gridCoordinates[2];
        int yMax = gridCoordinates[3];
        System.out.println("\\draw [help lines, color=green] (" + xMin + "," + yMin + ") grid " + "(" + xMax + "," + yMax + ");\n");

    }

    // Method to output LaTex commands
    void latex(String[] args, int[] gridCoordinates, Tree tree) {

        if (args != null) {
            for (String arg : args) {
                System.out.print(arg + " ");
            }
            System.out.println("are inserted \n");
        }

        // printSorted and averageDepth are called when drawing a tree but not a graph
        if (tree != null) {
            System.out.print("In sorted order: ");
            tree.printSorted(tree.root);
            System.out.println("\n");

            tree.createDepthList(tree.root, 0);
            tree.averageDepth(tree.depthListHead, true);
        }

        double scale = 0.600;
        int xMin = gridCoordinates[0];
        int xMax = gridCoordinates[2];
        int width = 1;

        // Calculate width of grid
        if (xMin > 0 && xMax > 0 || xMax > 0 && xMin < 0) {
            width = xMax - xMin + 1;
        }
        if (xMin < 0 && xMax < 0) {
            width = Math.abs(xMin) - Math.abs(xMax) + 1;
        }

        // 28 is the max width for textWidth 17.0cm and scale 0.600cm
        if (width > 28) {
            scale = textWidth / width;
        }
        if (scale < 0.3) {
            scale = 0.3;
        }

        // Output LaTex commands to draw grid and graph
        System.out.print("\\begin{tikzpicture}");
        System.out.printf("[scale=%.3f]", scale);
        System.out.println();
        outGrid(gridCoordinates);
        outGraph();
        System.out.println("\n\\end{tikzpicture}\n");
    }

    public static void main(String[] args) {
        Graph graph = new Graph();

        // Take the command line inputs and add the nodes in the graph
        for (int i=0; i<args.length/2; i++) {
            int x = Integer.parseInt(args[2*i]);
            int y = Integer.parseInt(args[2*i+1]);
            graph.addNode(x, y, i, Integer.toString(i+1));
        }

        // Take the command line inputs and add the edges in the graph
        for (int i=0; i<args.length/2 - 1; i++) {
            int x1 = Integer.parseInt(args[2*i]);
            int y1 = Integer.parseInt(args[2*i+1]);
            int x2 = Integer.parseInt(args[2*i+2]);
            int y2 = Integer.parseInt(args[2*i+3]);
            graph.addEdge(x1, y1, x2, y2);
        }

        // Call methods to output LaTex commands
        graph.outHeader();

        // If no input arguments
        if (args.length == 0) {
            System.out.println("Please input some coordinates \n");
        }

        int[] gridCoordinates = new int[4];

        // If no nodes were added to the graph
        if (graph.nodeListHead == null) {
            System.out.println("The graph is empty");
        }
        else {
            gridCoordinates = graph.getGridCoordinates();
        }

        graph.latex(null, gridCoordinates, null);
        graph.outFooter();

    }
}
```   

</details>                             
