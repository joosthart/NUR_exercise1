import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt

def f(x, a=2.2, b=0.5, c=3.1):
    return (x/b)**(a-1)*np.exp(-(x/b)**c)

def n(x, A, a=2.2, b=0.5, c=3.1, N=100):
    return A*N*(x/b)**(a-3)*np.exp(-(x/b)**c)

def romberg(f, start, end, m):
    h = end - start
    r = np.zeros(m)
    
    r[0] = 0.5*h*(f(start)+f(end))
    
    Np = 1
    for i in range(1, m-1):
        delta = h
        h *= 0.5
        x = start + h
        for _ in range(Np):
            r[i] += f(x)
            x += delta
            
        r[i] = 0.5*(r[i-1]+ delta*r[i])
        Np*=2
        
    Np=1
    for i in range(1, m-1):
        Np *= 4
        for j in range(0,m-i):
            r[j] = (Np*r[j+1] - r[j])/(Np-1)
    
    return r[0]

def bisect(x_range, x_i, order=1):
    if x_i <= x_range[0] and order==1:
        return 0
    if x_i >= x_range[-1] and order==1:
        return len(x_range)-1
    
    for i in range(len(x_range)-1):
        if x_i > x_range[i] and x_i < x_range[i+1]:
            return (i, i+1)
    raise RuntimeError

def interpolate(y1, y2, x1, x2, x, log=False):
    if log:
        y1 = np.log10(y1)
        y2 = np.log10(y2)
        x1 = np.log10(x1)
        x2 = np.log10(x2)
        x = np.log10(x)
    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    
    if log:
        y = 10**y
    return y

def poisson(k,l):
    p = 1
    for f in range(1, k+1):
        p*=l/f
    p*=np.exp(-l)
    return p

if __name__ == '__main__':
    
    # 3a

    # parameters
    a=2.2 
    b=0.5
    c=3.1
    
    # integrating
    integral = romberg(f, 0, 5, 10)

    A = 1/(4*np.pi*b**2*integral) # * r_vir^{-3}

    print('Solution for 3a:')
    print('Obtained value for A: {:.3f} r_vir^-3'.format(A))

    # 3b
    points = np.array([1e-4, 1e-2, 1e-1, 1, 5])
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

    print('Solution for 3c:')
    for k, l in kl:
        p = poisson(np.int64(k),np.int64(l))
        print('P({},{}) = {:.7E}'.format(l,k,p))
    