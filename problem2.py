import numpy as np

def crout(LU):
    for k in range(LU.shape[1]):
        for i in range(k+1,LU.shape[0]):
            LU[i,k] = LU[i,k]/LU[k,k]
            for j in range(k+1, LU.shape[1]):
                LU[i,j] = LU[i,j]-LU[i,k]*LU[k,j]
    return LU

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
    LU = crout(np.copy(Wss))
    print('Solutions for 2a:')
    print('LU decomposition of Wgs:')
    print(LU)
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
