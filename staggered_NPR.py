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


# QDP enumeration.

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
    binary = dtob4(n)
    gammas = [gamma[0], gamma[1], gamma[2], gamma[3]]
    gammas = map(one_or_M, gammas, binary)  # gamma[0]^n[0]..gamma[3]^n[3]
    return reduce(np.dot, gammas)

G = [Gamma(n) for n in range(16)]

# Staggered matrices.

def hypercubic(S, F, A, B):
    '''Hypercubic matrix elements.'''
    return (1./4)*trace(reduce(dot, [hc(G[A]), G[S], G[B], hc(G[F])]))
    
hypercubic_cache = np.array([[[[hypercubic(S,F,A,B) for B in range(16)]
                                                    for A in range(16)] 
                                                    for F in range(16)]
                                                    for S in range(16)])

def hypercubicM(S, F):
    '''Hypercubic staggered matrices.'''
    return np.array([[hypercubic(S,F,A,B) for B in range(16)] 
                      for A in range(16)])
           
def dotp(A,B):
    return np.dot(dtob4(A), dtob4(B))
    
def polespace(S, F, A, B):
    result = 0.
    for C in range(16):
        for D in range(16):
            sign = (-1)**(dotp(A,C) + dotp(D,B))  # Can speed up.
            result += sign*hypercubic_cache[S,F,C,D]
    return (1./16)*result
    
polespace_cache = np.array([[[[polespace(S,F,A,B) for B in range(16)]
                                                    for A in range(16)] 
                                                    for F in range(16)]
                                                    for S in range(16)])
    
def polespaceM(S,F):
    '''Polespace matrices.'''
    return np.array([[polespace_cache[S,F,A,B] for B in range(16)] 
                      for A in range(16)])
    
        
    
    

