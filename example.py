import ctypes, random, os, platform, sys
import numpy as np

pf = platform.system()

if pf == 'Linux':
    SDPB = ctypes.cdll.LoadLibrary("./libsdpb.so")
elif pf == 'Windows':
    os.add_dll_directory(os.getcwd())
    SDPB = ctypes.cdll.LoadLibrary('SDPB.dll')
else:
    print(f"Platform {pf} is not supported.")
    sys.exit(0)

def PotLLL(b, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.PotLLL.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.PotLLL(pp, d, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]
    

def DualPotLLL(b, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotLLL.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotLLL(pp, d, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]

def BKZ(b, beta, d, lp):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.BKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
    SDPB.BKZ(pp, beta, d, lp, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]

def PotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.PotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.PotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]

def DualPotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]

def SelfDualPotBKZ(b, beta, d):
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.SelfDualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.SelfDualPotBKZ(pp, beta, d, n, m)

    for i in range(N):
        for j in range(N): b[i, j] = pp[i][j]

if __name__ == '__main__':
    N = 80

    b = np.eye(N).astype(np.int64)
    for i in range(N): b[i, 0] = random.randint(100, 1000)
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

    #c = b.copy()
    #print("BKZ-reduced:")
    #BKZ(c, 40, 0.99, 10)
    #print(np.linalg.norm(c[0]))
    #print(c)

    c = b.copy()
    print("PotBKZ-reduced:")
    PotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Dual-PotBKZ-reduced:")
    DualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Self-Dual-PotBKZ-reduced:")
    SelfDualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)
