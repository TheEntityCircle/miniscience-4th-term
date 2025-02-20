from dolfin import *

# Create mesh and define function space
mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundaries

def boundary_D1(x, on_boundary):
    return on_boundary and (x[0] < -1.0 + DOLFIN_EPS)

def boundary_D2(x, on_boundary):
    return on_boundary and (x[0] > 1.0 - DOLFIN_EPS)

def boundary_Z(x, on_boundary):
    return on_boundary and (x[1] < -1.0 + DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)

# Define boundary condition
u_l = Constant(-1.0)
u_r = Constant(1.0)
u_z = Expression('x[0]', degree=2)
bc1 = DirichletBC(V, u_l, boundary_D1)
bc2 = DirichletBC(V, u_r, boundary_D2)
bc3 = DirichletBC(V, u_z, boundary_Z)

# Define variational problem
u = Function(V)
v = TestFunction(V)
la = Expression("1/(100 * exp(-25 * (x[0]*x[0] + x[1]*x[1])) + 1)", degree=2)

a = inner(grad(u), grad(v))*la*dx

# Compute solution
solve(a == 0, u, [bc1, bc2, bc3])

# Save solution in VTK format
W = VectorFunctionSpace(mesh, "Lagrange", 1)

file = File("test_phi.pvd")
file << u

E = project(-grad(u), W)
file = File("test_E.pvd")
file << E

j = project(-la*grad(u), W)
file = File("test_j.pvd")
file << j
