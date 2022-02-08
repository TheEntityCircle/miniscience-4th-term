import numpy as np

def add(a, b):
    for i in range(a.size):
         b[i] = a[i] + b[i]

def main():
    vec_size = int(1e7)

    a = b = np.full(vec_size, 1.0, dtype=np.float32)

    add(a, b)

#    is_ok = True
#    for i in range(b.size):
#        is_ok = is_ok and (b[i] == 2)
#    print("We are ok: " + str(is_ok))

if __name__ == '__main__':
    main()
