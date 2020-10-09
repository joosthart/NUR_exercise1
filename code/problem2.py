import numpy as np

def write_mat_to_file(M, filename):
    """ Writing Matrix to a file with 5 decimal numbers.

    Args:
        M (list): n x m list/array
        filename (str): filename where to store matrix.
    """
    file=''
    for i in range(M.shape[0]):
        line = ''
        for j in range(M.shape[1]):
            line += '{: 3.5f} '.format(M[i,j])
        file += line 
        file += '\n'
    with open(filename, 'w') as f:
        f.write(file)

def write_vec_to_file(f, filename):
    """ Writing list as column vector to a file with 9 decimal numbers.

    Args:
        f (list): list/array
        filename (str): filename where to store vector.
    """
    file=''
    for i in range(f.shape[0]):
        file += '{:.9f}\n'.format(f[i])
    with open(filename, 'w') as f:
        f.write(file.strip())

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
            # Calculate L part of LU matrix
            M[i,k] = M[i,k]/M[k,k]
            for j in range(k+1, M.shape[1]):
                # Calculate U part of LU matrix
                M[i,j] = M[i,j]-M[i,k]*M[k,j]
    
    if separate:
        # Splitting LU matrix intp two separate matrices, L and U.
        L = np.zeros(M.shape, dtype=np.dtype(M[0,0]))
        U = np.zeros(M.shape, dtype=np.dtype(M[0,0]))
        for i in range((M.shape[0])):
            L[i,:i+1] = M[i,:i+1]
            U[i,i:] = M[i,i:]
            L[i,i] = 1
        return L,U

    return M

def solve(LU, b):
    """ Solve LUx=b for x using LU decomposed matrix. L should have 1s on the 
    diagonal and L and U should be stored in the same array. The algorithm 
    solves the system by first solving Ly = b for y (forward pass) and than 
    Ux = y for x (backward pass).

    Args:
        LU (numpy.array): Numpy array (n x n) containing LU decomposition. 
        b (np.array): Numpy array (n) containg column vector.

    Returns:
        np.array: solution vector
    """
    # list since it will be appended and no numpy attributes required from y.
    y = []
    # Forward pass through matrix to solve Ly = b for y.
    for i in range(len(b)):
        # Dot product used instead of sum
        y.append(b[i] - LU[i,:i].dot(y[:i]))
    # Backward pass through matrix to solve Ux = y for x.
    for j in reversed(range(len(b))):
        # Dot product used instead of sum
        b[j] = 1/LU[j,j] * (y[j] - LU[j,j+1:].dot(b[j+1:]))
    return b

if __name__ == '__main__':
    Wss = np.loadtxt('data/wss.dat', dtype=np.float32)
    Wgs = np.loadtxt('data/wgs.dat', dtype=np.float32)

    # 2a
    
    # Calculating separate L and U matrix.
    L, U = crout(np.copy(Wss), separate=True)

    # Writing matrices to file
    write_mat_to_file(L, 'output/2a_L_matrix.txt')
    write_mat_to_file(U, 'output/2a_U_matrix.txt')


    # The LU is calculated again, since the solve function expects the LU-matries
    # in 1 matrix. Calculating it again is a waste of computing power, but it is
    # in order to print the L and U part separate, as asked in the exercise. 
    LU = crout(np.copy(Wss), separate=False)
    # Solving LU*x=Wgs for x.
    x_hat = solve(LU, np.copy(Wgs))    
    
    # Write sum of f to file
    with open('output/2a_f_vector_sum.dat', 'w') as f:
        f.write('{:.9f}'.format(sum(x_hat)))

    # Write vector to file
    write_vec_to_file(x_hat, 'output/2a_f_vector.txt')

    # 2b

    # Calculate residual vector, b_res.
    b_res = Wss.dot(x_hat) - Wgs
    # Solving LU*x_delta = b_res for x_delta.
    x_delta = solve(LU, np.copy(b_res))
    # Calculate first iterative improvement of x
    x_hat_second = x_hat - x_delta

    # Write sum of single iterative improvement f to file
    with open('output/2b_f_vector_sum.dat', 'w') as f:
        f.write('{:.9f}'.format(sum(x_hat_second)))

    # write vector of single iterative improvement to file
    write_vec_to_file(x_hat_second, 'output/2b_f_vector.txt')
