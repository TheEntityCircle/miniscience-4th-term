from mpi4py import MPI
import os, sys
from math import cos, pi

host = os.uname()[1]

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

def f(x):
    return cos(x)

def integrate(a, n, h):
    print("Starting from %f and doing %d steps with step %f" % (a, n, h))
    res = 0
    h2 = h/2

    for i in range(0, n):
        x = a + i * h
        res += f(x + h2) * h

    return res

a = 0.0
b = pi / 2
n = 100000
if len(sys.argv) > 1:
    n = int(sys.argv[1])

print("We are going to do %d steps" % n)

h = (b-a)/n
local_n = n/p

local_a = a + my_rank * local_n * h
local_b = local_a + local_n * h

integral = integrate(local_a, int(local_n), h)


if my_rank == 0:
    total = integral
    for s in range(1, p):
        integral = comm.recv(source=s)
        total = total + integral
else:
    comm.send(integral, dest = 0)

if my_rank == 0:
    print("res = %f" % total)

MPI.Finalize
