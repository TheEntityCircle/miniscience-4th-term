cpdef long long test(int x):
    cdef long long y = 1
    cdef long long i
    for i in range(1, x+1):
        y += i
    return y

