import numpy as np
from numpy import dot, trace
from numpy.linalg import inv

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
    tmp = (1./16)*np.array([[np.sum(H*sfunM(A,B)) for B in range(16)]
                                                  for A in range(16)])
    return np.kron(tmp, np.identity(3)) # Color dof.
    
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
    
def JKcompute(f, samples):
    '''Perform computation 'f' on the samples and JKsamples.'''
    return f(samples), map(f, JKsample(samples))

def JKsigma(ave, JKvals):
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
    
# Measurements.
def ps_trace(S, F, M):
    '''Trace M with polespace matrix.'''
    return (1./48)*trace(dot(polespaceM(S,F), M))
    
def Zq(prop, ap):
    Sinv = (1./2)*inv(prop) # 1/2 from action def.
    tr = ps_trace
    traces = tr(1,0,Sinv), tr(2,0,Sinv), tr(4,0,Sinv), tr(8,0,Sinv)
    Zinv = np.dot(ap, traces).imag/np.dot(ap, ap)
    return 1./Zinv
    
def bilinear_Lambda(prop, bilinear, S, F):
    '''Trace of amputated bilinear correlation function.'''
    iprop = inv(prop)
    amputated = reduce(dot, [iprop, bilinear, iprop])
    return ps_trace(S, F, amputated)
    
def bilinear_Z(prop, bilinear, ap, S, F):
    '''Bilinear Z-factor in RI' scheme.'''
    # Should add Cosine factors to make general.
    return Zq(prop, ap)/bilinear_Lambda(prop, bilinear, S, F)
    
def Lambda_V(prop, bilinear, ap, muhat):
    '''Traced vector bilinear divided by tree-level value.'''
    S = 2**muhat
    # 1/2 comes from coding of vector operator.
    tr = (1./2)*bilinear_Lambda(prop, bilinear, S, 0)
    return tr/np.cos(ap[muhat])

def Lambda_Vave(prop, bl0, bl1, bl2, bl3, ap):
    LV = lambda bl, muhat: Lambda_V(prop, bl, ap, muhat)
    return (1/4.)*(LV(bl0,0) + LV(bl1,1) + LV(bl2,2) + LV(bl3,3))
                      
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
    
        
    
    

