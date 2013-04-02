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
    
#hypercubic_cache = np.array([[[[hypercubic(S,F,A,B) for B in range(16)]
#                                                    for A in range(16)] 
#                                                    for F in range(16)]
#                                                    for S in range(16)])

def hypercubicM(S, F):
    '''Hypercubic staggered matrices.'''
    return np.array([[hypercubic(S,F,A,B) for B in range(16)] 
                      for A in range(16)])
           
def dotp(A,B):
    return np.dot(dtob4(A), dtob4(B))
    
dotpM = np.array([[np.dot(dtob4(A), dtob4(B)) for B in range(16)]
                                              for A in range(16)])
    
def sfun(A, B, C, D):
    return (-1)**(dotp(A,C) + dotp(D,B))
    
def sfunM(A,B):
    tmp = np.array([[dotp(A,C) + dotp(D,B) for D in range(16)]
                                           for C in range(16)]) # fromfunction?
    return np.cos(np.pi*tmp)
    
def sfunM2(A,B):
    tmp = dotpM[A,:] + dotpM[:,B].reshape((16,1))
    return np.cos(np.pi*tmp)
    
def polespace(S, F, A, B):
    result = 0.
    H = hypercubicM(S,F)
    #Sfun = sfunM(A,B)
    for C in range(16):
        for D in range(16):
            sign = (-1)**(dotp(A,C) + dotp(D,B))  # Can speed up.
            result += sign*H[C,D]
    return (1./16)*result
    
def polespace2(S, F, A, B):
    H = hypercubicM(S,F)
    Sfun = sfunM(A,B)
    return (1./16)*np.sum(H*Sfun)
    
def polespace3(S, F):
    H = hypercubicM(S,F)
    tmp = np.array([[np.sum(H*sfunM2(A,B)) for B in range(16)]
                                          for A in range(16)])
    return (1./16)*tmp
    
#polespace_cache = np.array([[[[polespace(S,F,A,B) for B in range(16)]
#                                                    for A in range(16)] 
#                                                    for F in range(16)]
#                                                    for S in range(16)])
    
def polespaceM(S,F):
    '''Polespace matrices.'''
    return np.array([[polespace_cache[S,F,A,B] for B in range(16)] 
                      for A in range(16)])
                      
def main():
#    print 'calculating hypercubicM 1 1'
#    hc = hypercubicM(1,1)
#    print 'calculating SfunM 1 1'
#    sf = sfunM(1,1)
    print 'calculating polespace 1 1'
    ps = np.array([[polespace(6,11,A,B) for B in range(16)] for A in range(16)])
#    print 'calculating polespace2 1 1'
#    ps2 = np.array([[polespace2(1,1,A,B) for B in range(16)] for A in range(16)])
    print 'calculating polespace3 1 1'
    ps3 = polespace3(6,11)
     
    print ps.all() == ps3.all()
    #import cProfile
    #cProfile.run('polespace3(1,1)')
    
if __name__ == "__main__":
    main()
    
        
    
    

