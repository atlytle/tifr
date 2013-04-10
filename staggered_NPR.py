import numpy as np
from numpy import dot, trace

# Define gamma matrices.
gamma = [None, None, None, None]

gamma[0] = np.array([[0, 0, 0, 1j],
                     [0, 0, 1j, 0],
                     [0, -1j, 0, 0],
                     [-1j, 0, 0, 0]])

gamma[1] = np.array([[0, 0, 0, -1],
                     [0, 0, 1, 0],
                     [0, 1, 0, 0],
                     [-1, 0, 0, 0]])

gamma[2] = np.array([[0, 0, 1j, 0],
                     [0, 0, 0, -1j],
                     [-1j, 0, 0, 0],
                     [0, 1j, 0, 0]])

gamma[3] = np.array([[0, 0, 1, 0],
                     [0, 0, 0, 1],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0]])

id4 = np.identity(4, dtype=complex)


# QDP enumeration of gamma matrices.
def hc(M):
    "Hermitian conjugate."
    return np.conjugate(np.transpose(M))

def dtob4(d):
    "Convert a decimal integer 0..15 to binary."
    result = [None, None, None, None]
    if d in range(16):
        for i in range(4):
            r = d % 2
            result[i] = r
            d = (d-r)/2
        return result
    else:
        raise ValueError(str(d) + ' not in 0..15')

def one_or_M(M, num):
    if M.shape == (4, 4): 
        if num is 0:
            return id4
        if num is 1:
            return M
        else:
            raise ValueError("num must be 0 or 1")
    else:
        raise TypeError("M must be a 4x4 array")

def Gamma(n):
    '''QDP enumeration of gamma matrices.'''
    binary = dtob4(n)
    gammas = [gamma[0], gamma[1], gamma[2], gamma[3]]
    gammas = map(one_or_M, gammas, binary)  # gamma[0]^n[0]..gamma[3]^n[3]
    return reduce(np.dot, gammas)

G = [Gamma(n) for n in range(16)]

# Hypercubic matrices.
def hypercubic(S, F, A, B):
    '''Hypercubic matrix elements.'''
    return (1./4)*trace(reduce(dot, [hc(G[A]), G[S], G[B], hc(G[F])]))

def hypercubicM(S, F):
    '''Hypercubic staggered matrices.'''
    return np.array([[hypercubic(S,F,A,B) for B in range(16)] 
                      for A in range(16)])

# Helper functions.
def dotp(A,B):
    return np.dot(dtob4(A), dtob4(B))
    
def sfun(A, B, C, D):
    return (-1)**(dotp(A,C) + dotp(D,B))
    
dotpM = np.array([[np.dot(dtob4(A), dtob4(B)) for B in range(16)]  
                                              for A in range(16)])
    
def sfunM(A,B):
    '''Fast implementation making use of cached dotpM;
    S(A,B)_{CD} == (-1)**(dotp(A,C) + dotp(D,B))'''
    
    tmp = dotpM[A,:] + dotpM[:,B].reshape((16,1))
    return np.cos(np.pi*tmp)
    
# Polespace matrices.
def polespace(S, F, A, B):
    '''Naive implementation of polespace matrix elements.'''
    result = 0.
    H = hypercubicM(S,F)
    for C in range(16):
        for D in range(16):
            sign = sfun(A, B, C, D)  # Can speed up.
            result += sign*H[C,D]
    return (1./16)*result
 
def polespaceM(S, F):
    '''Fast implementation of polespace matrices;
    
    P(S,F)_{A,B} == H(S,F)_{CD} S(A,B)_{CD}.'''
    
    H = hypercubicM(S,F)
    tmp = np.array([[np.sum(H*sfunM(A,B)) for B in range(16)]
                                          for A in range(16)])
    return (1./16)*tmp
    
# Some stats stuff.
def JKsample(list):
    '''Yield the jackknife samples of the elements in list.'''
    sample = []
    for dummy in range(len(list)):
        x, xs = list[0], list[1:]
        sample.append(xs[:])
        xs.append(x)
        list = xs
    return sample

def JKsigma(JKvals, ave):
    N = len(JKvals)
    diffs = [(JKval - ave)*(JKval - ave) for JKval in JKvals]
    return np.sqrt(sum(diffs)*(1-1./N))

def bootstrap_sample(data, N):
    '''Yield N boostrap samples of the elements in data.'''
    ri = random.randint
    L = len(data)
    sample = []
    for foo in range(N):
        sample.append([data[ri(0, L-1)] for bar in range(L)])
    return sample
                      
def test():
    print 'calculating polespace naive'
    ps = np.array([[polespace(6,11,A,B) for B in range(16)] for A in range(16)])
    print 'calculating polespace fast'
    ps3 = polespaceM(6,11)
    print 'Same result? ',
    print ps.all() == ps3.all()
    
    import cProfile
    cProfile.run('polespaceM(1,1)')
    
if __name__ == "__main__":
    test()
    
        
    
    

