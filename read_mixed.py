import itertools

nx, ny, nz, nt = 24, 24, 24, 64
V = nx*ny*nz
nc = 3
    
def reorder_HISQ(c1, c2, z, y, x):
    return c1*nz*ny*nx*nc + z*ny*nx*nc + y*nx*nc + x*nc + c2
                                              
new_indices = [reorder_HISQ2(c1,c2,z,y,x) 
               for c1, c2, z, y, x in itertools.product(range(nc), range(nc),
                                              range(nz), range(ny), range(nx))]
                                              
print new_indices[0:10]
    

