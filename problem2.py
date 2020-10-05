import numpy as np

def crout(M):
    LU = np.copy(M)
    for k in range(M.shape[1]):
        for i in range(k+1,M.shape[0]):
            LU[i,k] = LU[i,k]/LU[k,k]
            for j in range(k+1, M.shape[1]):
                LU[i,j] = LU[i,j]-LU[i,k]*LU[k,j]
    return LU

def solve(LU, b):
    x = np.copy(b)
    y = np.zeros(b.shape, dtype=np.float32)
    for i in range(len(b)):
        y[i] = b[i] - LU[i,:i].dot(y[:i])
    for i in reversed(range(len(b))):
        x[i] = 1/LU[i,i] * (y[i] - LU[i,i+1:].dot(x[i+1:]))
    return x

if __name__ == '__main__':
    Wss = np.loadtxt('data/wss.dat', dtype=np.float32)
    Wgs = np.loadtxt('data/wgs.dat', dtype=np.float32)

    # 2a
    LU = crout(Wss)
    print('Solutions for 2a:')
    print('LU decomposition of Wgs:')
    print(LU)
    x_hat = solve(LU, Wgs)
    print('Solution for f: {}'.format(x_hat))
    print('Sum of f: {}'.format(sum(x_hat)))

    # 2b
    print('Solutions for 2a:')
    b_res = Wss.dot(x_hat) - Wgs
    x_delta = solve(LU, b_res)
    x_hat_second = x_hat - x_delta
    print('Improved f: {}'.format(x_hat_second))
    print('Sum of improved f {}'.format(sum(x_hat_second)))
