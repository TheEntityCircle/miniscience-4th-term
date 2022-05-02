from cython.parallel import prange

cpdef long long test(int x):
    cdef long long y = 1
    cdef long long i
    for i in prange(1, x+1, nogil=True):
        y += i
    return y

