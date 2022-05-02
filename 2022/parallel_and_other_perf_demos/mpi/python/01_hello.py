from mpi4py import MPI
import os

host = os.uname()[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
p = comm.Get_size()

print("Hello world, I'm process number: %d on host %s" % (rank, host))
