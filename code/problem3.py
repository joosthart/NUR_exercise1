import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt

def f(x, a=2.2, b=0.5, c=3.1):
    """integrant"""
    return (x/b)**(a-1)*np.exp(-(x/b)**c)

def n(x, A, a=2.2, b=0.5, c=3.1, N=100):
    """Number density profile"""
    return A*N*(x/b)**(a-3)*np.exp(-(x/b)**c)

def romberg(f, start, end, m):
    """Romburg integration for a given function f from start to end. The oreder
    of the integration is given by m, i.e. the number of interval splits to use 
    to calculate the integral. 

    Args:
        f (func): Function taking as argument x coordinate and outputs 
            corresponding y
        start (int/float): Starting point for integration
        end (int/float): End point for integration
        m (int): Order of romburg integration

    Returns:
        float : Value of integral
    """
    # Initial stepsize
    h = end - start

    # Intialize array for estimates
    r = np.zeros(m)
    
    # 0th oreder estimate
    r[0] = 0.5*h*(f(start)+f(end))
    
    Np = 1 # Number of new points
    for i in range(1, m-1):
        delta = h # step size between points
        h *= 0.5 # decreasing stepsize
        x = start + h
        for _ in range(Np):
            r[i] += f(x) # new evaluation
            x += delta
        
        # Combine new point with previously calculated 
        r[i] = 0.5*(r[i-1]+ delta*r[i])
        
        # increment number of points for next run
        Np*=2
    
    # Upgrading previous results by combining them. 
    Np = 1 # Reset number of points
    for i in range(1, m-1):
        Np *= 4
        for j in range(0,m-i):
            # combining results j and j+1 and storing it in j
            r[j] = (Np*r[j+1] - r[j])/(Np-1)
    
    return r[0]

def bisect(x_range, x_i):
    """ Find the indices of elements in x_range that enclose x_i. If x_i outside 
        of x_range, the begin or end index is returned.

    Args:
        x_range (list): Sorted ascending list.
        x_i (float, int): Float or int to find encloding values for

    Raises:
        RuntimeError: No edges found.
    """
    if x_i <= x_range[0]:
        return 0
    if x_i >= x_range[-1]:
        return len(x_range)-1
    
    for i in range(len(x_range)-1):
        if x_i > x_range[i] and x_i < x_range[i+1]:
            return (i, i+1)
    raise RuntimeError("No edges found.")

def interpolate(y1, y2, x1, x2, x, log=False):
    """ Linear interpolator. 

    Args:
        y1 (int/float): y value of left interval boundry.
        y2 (int/float): y value of right interval boundry.
        x1 (int/flout): x value of left interval boundry.
        x2 (int/float): x value of right interval boundry.
        x (int/float): x value for which to calculate interpolation.
        log (bool, optional): If true, interpolate linear in log-space. Defaults
             to False.

    Returns:
        float: interpolated y(x)
    """
    if log:
        # transform values to log space
        y1 = np.log10(y1)
        y2 = np.log10(y2)
        x1 = np.log10(x1)
        x2 = np.log10(x2)
        x = np.log10(x)

    # Calculate y(x)
    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    
    if log:
        # transform result back
        y = 10**y
    return y

def poisson(k,l):
    """Poisson probability calculator. The caculation of e^-(l)l^k/k! is 
    performed in a loop. Thus, k! and l^k do not have to be calculated directly. 
    This prevents overflows up to higher k and l compared to direct caculations.

    Args:
        k (int): Integer for which to calculate probability
        l (float): Positive value for the mean of the distribution

    Returns:
        float: Poisson probability for P(k,l)

    Raises:
        ValueError: l should be larger than 0.
        ValueError: k should be an integer.
    """

    # Check l and k allowed values
    if l<0:
        raise ValueError("l should be larger than 0.")
    elif (
            not isinstance(k, int) and 
            not isinstance(k, np.int32) and 
            not isinstance(k, np.int64)
         ):
        raise ValueError("k should be an integer not {}.".format(type(k)))
    elif k < 0:
        ValueError("k should be larger than 0.")
    if k == 0:
        # l^k=1 and k! = 1
        return np.exp(-l)
    else:
        p = 1
        # Calculate l^k/k!*e^-l using a product
        e = np.exp(-l/k) # calculate ones
        for f in range(1, k+1):
            p*=l/f*e
        return p

if __name__ == '__main__':
    
    # 3a
    # parameters
    a=2.2 
    b=0.5
    c=3.1
    
    # Calculating integral using romberg integration
    integral = romberg(f, 0, 5, 10)

    # Using the value of the integral, calculate the estimate for A.
    A = 1/(4*np.pi*b**2*integral) # * r_vir^{-3}

    # Write to file
    with open('./output/3a_a.dat', 'w') as f:
        f.write('{:.8f}'.format(A))

    # 3b
    # Points used for interpolation
    points = np.array([1e-4, 1e-2, 1e-1, 1, 5])
    # Function evaluations at the above points
    evaluations = n(points, A)

    # Interpolating in logspace
    y_hat=[]
    for i in np.arange(1e-4, 5, 0.01):
        if i <= points[0]:
            idx1, idx2 = 0,1
        elif i >= points[-1]:
            idx1, idx2 = len(points)-2, len(points)-1
        else:
            idx1, idx2 = bisect(points, i)
        # liniarly interpolate in logspace
        y_inter = interpolate(
                      evaluations[idx1], 
                      evaluations[idx2], 
                      points[idx1], 
                      points[idx2], 
                      i, 
                      log=True
                  )
        y_hat.append(y_inter)
    
    # Plotting results
    plt.loglog(points, n(points, A), '.',ms=10, label='Points')
    plt.loglog(np.arange(1e-4, 5, 0.01), y_hat, label='Linear intepolation')
    plt.loglog(np.arange(1e-4, 4, 0.01), n(np.arange(1e-4, 4, 0.01), A), '--', 
               label=r'True $n(x)$')
    plt.axis(ymax =1e6, ymin=1e-5)
    plt.legend()
    plt.ylabel('Number Density')
    plt.xlabel(r'Relative Radius ($\frac{r}{r_\mathrm{vir}}$)')
    plt.tight_layout()
    plt.savefig('plots/3b_galaxy_satelite_numberdensity.png')

    # 3c
    # Combinations of (k,l)
    kl = [(0,1), (10,5), (21,3), (40,2.6), (200, 101)]

    doc = ''
    for k, l in kl:
        # Calculate Possion probabilities for (k,l) combinations.
        p = poisson(np.int64(k),np.int64(l))
        doc += 'P({},{}) = {:.7E}'.format(l,k,p)
        doc += '\n'

    print(doc)

    # Write to file
    with open('./output/3c_poisson.txt', 'w') as f:
        f.write(doc)
