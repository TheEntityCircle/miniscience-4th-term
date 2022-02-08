import plane_python
import using_cython
import using_cython_and_openmp
import time

number = 100000000

start = time.time()
print(plane_python.test(number))
end =  time.time()

py_time = end - start
print("Python time = {}".format(py_time))

start = time.time()
print(using_cython.test(number))
end =  time.time()

cy_time = end - start
print("Cython time = {}".format(cy_time))

start = time.time()
print(using_cython_and_openmp.test(number))
end =  time.time()

parcy_time = end - start
print("Cython && OpenMP time = {}".format(parcy_time))

