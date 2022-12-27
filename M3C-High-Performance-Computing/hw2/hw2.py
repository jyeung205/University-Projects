"""M3C 2018 Homework 2 Jimmy Yeung 01203261"""

import numpy as np
import matplotlib.pyplot as plt
from m1 import nmodel as nm #assumes that hw2_dev.f90 has been compiled with: f2py -c hw2_dev.f90 -m m1
# May also use scipy, scikit-learn, and time modules as needed
from scipy import optimize

def read_data(tsize=60000):
    """Read in image and label data from data.csv.
    The full image data is stored in a 784 x 70000 matrix, X
    and the corresponding labels are stored in a 70000 element array, y.
    The final 70000-tsize images and labels are stored in X_test and y_test, respectively.
    X,y,X_test, and y_test are all returned by the function.
    You are not required to use this function.
    """
    print("Reading data...") #may take 1-2 minutes
    Data=np.loadtxt('data.csv',delimiter=',')
    Data =Data.T
    X,y = Data[:-1,:]/255.,Data[-1,:].astype(int)%2 #rescale the image, convert the labels to 0s and 1s (For even and odd integers)
    Data = None

    # Extract testing data
    X_test = X[:,tsize:]
    y_test = y[tsize:]
    print("processed dataset")
    return X,y,X_test,y_test
#----------------------------

def snm_test(X,y,X_test,y_test,omethod,input=(None)):
    """Train single neuron model with input images and labels (i.e. use data in X and y), then compute and return testing error in test_error
    using X_test, y_test. The fitting parameters obtained via training should be returned in the 1-d array, fvec_f
    X: training image data, should be 784 x d with 1<=d<=60000
    y: training image labels, should contain d elements
    X_test,y_test: should be set as in read_data above
    omethod=1: use l-bfgs-b optimizer
    omethod=2: use stochastic gradient descent
    input: tuple, set if and as needed
    """
    """
    Fortran is a compiled programming language designed for numerical computing.
    An advantage of using the Python + Fortran approach is that Fortran will run faster,
    especially when doing computational heavy executions.
    There is also no need for me to vectorise my code in Fortran, unlike in Python.

    A disadvantage of using Python + Fortran are some of the differences in syntax.
    An example of this is that Fortran is case insensitive, whereas Python is case sensitive.
    This can cause problems when assigning variables.
    Another example is that the index in Fortran starts at 0, whereas in Python it starts
    at 1.
    """
    n = X.shape[0]
    fvec = np.random.randn(n+1) #initial fitting parameters

    #Add code to train SNM and evaluate testing test_error
    d = X.shape[1]

    nm.nm_x = X
    nm.nm_y = y

    if omethod == 1:
        def tomin(weights):
            c, cgrad = nm.snmodel(weights,d)
            return c
        result = optimize.minimize(tomin, fvec, method = "L-BFGS-B")
        fvec_f = result.x

    elif omethod == 2:
        fvec_f = nm.sgd(fvec,n,0,d,0.1)

    #compute test error
    nm.nm_x = X_test
    nm.nm_y = y_test
    a = nm.snforward(fvec_f,X_test.shape[1])
    n_correct = np.sum(np.round(a) == y_test)

    test_error = 1 - (n_correct/len(y_test)) #Modify to store testing error; see neural network notes for further details on definition of testing error
    output = (None) #output tuple, modify as needed
    return fvec_f,test_error,output
#--------------------------------------------

def nnm_test(X,y,X_test,y_test,m,omethod,input=(None)):
    """Train neural network model with input images and labels (i.e. use data in X and y), then compute and return testing error (in test_error)
    using X_test, y_test. The fitting parameters obtained via training should be returned in the 1-d array, fvec_f
    X: training image data, should be 784 x d with 1<=d<=60000
    y: training image labels, should contain d elements
    X_test,y_test: should be set as in read_data above
    m: number of neurons in inner layer
    omethod=1: use l-bfgs-b optimizer
    omethod=2: use stochastic gradient descent
    input: tuple, set if and as needed
    """
    n = X.shape[0]
    fvec = np.random.randn(m*(n+2)+1) #initial fitting parameters

    #Add code to train NNM and evaluate testing error, test_error
    d = X.shape[1]

    nm.nm_x = X
    nm.nm_y = y

    if omethod == 1:
        def tomin(weights):
            c, cgrad = nm.nnmodel(weights,n,m,d)
            return c
        result = optimize.minimize(tomin, fvec, method = "L-BFGS-B")
        fvec_f = result.x

    elif omethod == 2:
        fvec_f = nm.sgd(fvec,n,m,d,0.1)

    #compute test error
    nm.nm_x = X_test
    nm.nm_y = y_test
    a, a_outer = nm.nnforward(fvec_f,n,m,X_test.shape[1])
    n_correct = np.sum(np.round(a_outer) == y_test)

    test_error = 1 - (n_correct/len(y_test)) #Modify to store testing error; see neural network notes for further details on definition of testing error
    output = (None) #output tuple, modify as needed
    return fvec_f,test_error,output

#--------------------------------------------

def nm_analyze():
    """ Analyze performance of single neuron and neural network models
    on even/odd image classification problem
    Add input variables and modify return statement as needed.
    Should be called from
    name==main section below
    """
    """
    I've decided to investigate the effect of the training size,d, and the effect
    of the no. of inner neurons in NNM,m, on the testing error.
    Before investigating, you would expect the testing error to get smaller as the
    training size increases, and as the number of inner neurons increase.


    In hw21, I have plotted a graph of testing error against training size (where d goes from 1000 up to 30000,
    in increments of 1000) for both SNM and NNM.
    For NNM, the no. of inner neurons is set to 3.
    By analysing the graph, we can see that the testing error decreases as the training size increases, for both SNM and NNM.
    We can also see that the testing error for NNM is less than SNM for all of our plotted values of d between 1000 and 30000.

    The rate at which the Testing error for SNM decreases slows down as the training size increases, and it looks like it will stop.
    The graph suggests that the testing error for SNM may stay around the value 0.10, even
    if we increase the training size to a value larger than 30000.
    So it looks like there is a limit to how accurate SNM can become, no matter how large the training size.

    The Testing error for NNM decreases at a quicker rate than the testing error for SNM.
    The testing error begins to decrease at a slower rate as d gets to around 15000.
    The testing error doesn't change much for d from 20000 to 30000.
    This suggests that there may also be a limit to how accurate NNM for m = 3 can get,
    no matter how large the training size.

    We can conclude that the testing error for both SNM and NNM (m=3) decrease, as the training size increases.
    We can also conclude that NNM (m=3) is more accurate than SNM.
    But there is a limit to how accurate they both can become.
    This suggests that there is a limit to how accurate a NNM can become depending on the no. of neurons in
    the network, no matter how large our training size is.


    In hw22, I have plotted the testing error against the no. of inner neurons (from 1 to 14) to test whether increasing
    the number of inner neurons will make NNM more accurate. The training size is set to 5000.

    By analysing the graph, we can see that the testing error decreases as the number of inner neurons increases.
    However, the rate at which the Testing error decreases slows down at about 4/5 inner neurons.
    We can see that the testing error only decreases very slightly when the no. of inner neurons
    is increased from 4 to 14.
    This suggests there is limited benefit in increasing the no. of inner neurons
    to greater than 5, given the increased time it takes to compute - at least for a training size of 5000.
    Maybe the testing error would continue to decrease if I increased the training size.


    In hw23, I plotted the testing error against the no. of inner neurons(from 5 to 50 going up in 5s),
    but this time d is set to 10,000.

    We can see that the testing error does not decrease by much, even if the
    no. of inner neurons is increased a lot.
    Here the number of inner neurons increase tenfold, from 5 to 50.
    But the testing error only decrease by around 15%.


    In conclusion, increasing the training size and increasing the no. of inner neurons does decrease
    the testing error.
    But these both need to be increased proportionally.
    For each given training size, there is an optimal no. of inner neurons to make the
    computation most efficient (increasing the no. of inner neurons from this
    optimal number will just make the computation take longer, without actually decreasing the
    testing error by much).
    """
    X,y,X_test,y_test = read_data()

    nn_errors, sn_errors = list(), list()
    x_axis = list(range(1000, 31000, 1000))
    for d in x_axis:
        _, nn_error, _ = nnm_test(X[:,:d],y[:d],X_test,y_test,3,2)
        _, sn_error, _ = snm_test(X[:,:d],y[:d],X_test,y_test,2)
        nn_errors.append(nn_error)
        sn_errors.append(sn_error)
        print("calculating d=", d)
    plt.figure()
    plt.plot(x_axis,sn_errors, label="SNM")
    plt.plot(x_axis,nn_errors, label="NNM - no. of inner neurons is 3")
    plt.title("Testing error against training size")
    plt.xlabel('Training size, d')
    plt.ylabel('Testing Error')
    plt.legend()
    plt.show()

    nn_errors = list()
    x_axis = list(range(1,15,1))
    for m in x_axis:
        _, nn_error, _ = nnm_test(X[:,:5000],y[:5000],X_test,y_test,m,2)
        nn_errors.append(nn_error)
        print("calculating m=", m)
    plt.figure()
    plt.plot(x_axis,nn_errors, label="d = 5000")
    plt.title("Testing Error against No. of inner neurons")
    plt.xlabel('Number of Inner Neurons, m')
    plt.ylabel('Testing Error')
    plt.legend()
    plt.show()

    nn_errors = list()
    x_axis = list(range(5,55,5))
    for m in x_axis:
        _, nn_error, _ = nnm_test(X[:,:10000],y[:10000],X_test,y_test,m,2)
        nn_errors.append(nn_error)
        print("calculating m=", m)
    plt.figure()
    plt.plot(x_axis,nn_errors, label="d = 10000")
    plt.title("Testing Error against No. of inner neurons")
    plt.xlabel('Number of Inner Neurons, m')
    plt.ylabel('Testing Error')
    plt.legend()
    plt.show()

    return None
#--------------------------------------------

def display_image(X):
    """Displays image corresponding to input array of image data"""
    n2 = X.size
    n = np.sqrt(n2).astype(int) #Input array X is assumed to correspond to an n x n image matrix, M
    M = X.reshape(n,n)
    plt.figure()
    plt.imshow(M)
    #plt.show()
#--------------------------------------------
#--------------------------------------------


if __name__ == '__main__':
    #The code here should call analyze and generate the
    #figures that you are submitting with your code
    output = nm_analyze()
