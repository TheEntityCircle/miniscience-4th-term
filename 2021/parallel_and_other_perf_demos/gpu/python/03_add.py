import numpy as np
from numba import cuda

@cuda.jit
def add(a, b):
    pos = cuda.grid(1)
    if pos < a.size:
         b[pos] = a[pos] + b[pos]

def main():
    vec_size = int(1e7)

    a = b = np.full(vec_size, 1.0, dtype=np.float32)

    d_A = cuda.to_device(a)
    d_B = cuda.to_device(b) 

    # Configure the blocks
    threadsperblock = 32
    blockspergrid = (a.shape[0] + (threadsperblock - 1)) // threadsperblock
    
    # Start the kernel 
    add[blockspergrid, threadsperblock](d_A, d_B)

    # Get the results
    b = d_B.copy_to_host()

#    is_ok = True
#    for i in range(b.size):
#        is_ok = is_ok and (b[i] == 2)
#    print("We are ok: " + str(is_ok))


if __name__ == '__main__':
    main()
