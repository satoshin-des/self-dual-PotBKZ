import ctypes, random, os, platform, sys
import numpy as np

pf = platform.system()

#N = int(input("lattice dimension = "))
print("lattice dimension = ")
N = 90

if pf == 'Linux':
    SDPB = ctypes.cdll.LoadLibrary("./SelfDualPotBKZ.so")
elif pf == 'Windows':
    os.add_dll_directory(os.getcwd())
    SDPB = ctypes.cdll.LoadLibrary('SelfDualPotBKZ.dll')
else:
    print(f"Platform {pf} is not supported.")
    sys.exit(0)

def PotLLL(b, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    for i in range(N):
        for j in range(N):
            pp[i][j] = ctypes.c_long(b[i, j])

    SDPB.PotLLL.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.PotLLL.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.PotLLL(pp, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]
    

def DualPotLLL(b, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotLLL.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotLLL.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.DualPotLLL(pp, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def PotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    SDPB.PotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.PotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    bb = SDPB.PotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def DualPotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.DualPotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def SelfDualPotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.SelfDualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.SelfDualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.SelfDualPotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

if __name__ == '__main__':
    b = np.eye(N, dtype=int)
    
    is_svp_challenge = True
    
    if is_svp_challenge:
        with open(f'svp_challenge_list/SVP-{N}-{0}.svp') as f:
            basis = np.array(f.read().split()[1:], dtype=int)
        
        for i in range(N):
            for j in range(N):
                b[i, j] = basis[j + i * N]
        c = b.copy()
    else:
        for i in range(N):
            b[i, 0] = random.randint(100, 1000)
        c = b.copy()

    print(b)
    print("PotLLL-reduce:")
    PotLLL(c, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Dual-PotLLL-reduced:")
    DualPotLLL(c, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("PotBKZ-reduced:")
    PotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Dual-PotBKZ-reduced:")
    #DualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Self-Dual-PotBKZ-reduced:")
    SelfDualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)
    
    
