import numpy as np
from numba import guvectorize, cuda

@guvectorize(['void(float32[:], float32[:])'], '(n)->(n)', target='cuda', nopython=True)
def add(a, b):
    for i in range(a.shape[0]):
         b[i] = a[i] + b[i]

def main():
    vec_size = int(1e7)

    a = b = np.full(vec_size, 1.0, dtype=np.float32)

    d_A = cuda.to_device(a)
    d_B = cuda.to_device(b) 
    add(d_A, d_B)
    b = d_B.copy_to_host()

#    is_ok = True
#    for i in range(b.size):
#        is_ok = is_ok and (b[i] == 2)
#    print("We are ok: " + str(is_ok))


if __name__ == '__main__':
    main()
