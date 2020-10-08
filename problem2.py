import numpy as np

def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")

def crout(M, separate=False):
    """ Generating LU-composition for input matrix M Using the improved Crout's 
    algorithm. L and U matrices are returned as 1 matrix, where the fact that 
    the diagnonal elements of the L matrix has been exploited. If the separate 
    flag is set to true, separate L and U matrices are returned. 

    Args:
        M (list, numpy.array): N x N array for which to calulate LU-decomposition
        separate (bool, optional): Return separate L and U matrix. Default set 
            to False.

    Returns:
        LU decomposition of M.
    """
    for k in range(M.shape[1]):
        for i in range(k+1,M.shape[0]):
            M[i,k] = M[i,k]/M[k,k]
            for j in range(k+1, M.shape[1]):
                M[i,j] = M[i,j]-M[i,k]*M[k,j]
    
    if separate:
        L = np.zeros(M.shape, dtype=np.dtype(M[0,0]))
        U = np.zeros(M.shape, dtype=np.dtype(M[0,0]))
        for i in range((M.shape[0])):
            L[i,:i+1] = M[i,:i+1]
            U[i,i:] = M[i,i:]
            L[i,i] = 1
        return L,U

    return M

def solve(LU, b):
    y = []
    for i in range(len(b)):
        y.append(b[i] - LU[i,:i].dot(y[:i]))
    for j in reversed(range(len(b))):
        b[j] = 1/LU[j,j] * (y[j] - LU[j,j+1:].dot(b[j+1:]))
    return b

if __name__ == '__main__':
    Wss = np.loadtxt('data/wss.dat', dtype=np.float32)
    Wgs = np.loadtxt('data/wgs.dat', dtype=np.float32)

    # 2a
    L, U = crout(np.copy(Wss), True)
    print('Solutions for 2a:')
    print('LU decomposition of Wgs:')
    print('L: {}'.format(L))
    print('U: {}'.format(U))
    x_hat = solve(LU, np.copy(Wgs))
    print('Solution for f: {}'.format(x_hat))
    print('Sum of f: {}'.format(sum(x_hat)))

    # 2b
    print('Solutions for 2b:')
    b_res = Wss.dot(x_hat) - Wgs
    x_delta = solve(LU, np.copy(b_res))
    x_hat_second = x_hat - x_delta
    print('Improved f: {}'.format(x_hat_second))
    print('Sum of improved f {}'.format(sum(x_hat_second)))
