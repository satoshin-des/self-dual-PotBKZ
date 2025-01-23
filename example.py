import ctypes, random, platform, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

is_graphs_output_mode = True # output graphs of potential or not

N = int(input("lattice dimension = "))

if platform.system() == 'Linux':
    SDPB = ctypes.cdll.LoadLibrary("./SelfDualPotBKZ.so")
else:
    print(f"Platform {platform.system()} is not supported.")
    sys.exit(0)

def PotLLL(b: np.ndarray, d: float) -> None:
    """PotLLL reduction in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        d (float): reduction parameter
    """
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
    
def BKZ(b: np.ndarray, block_size: int, d: float, max_loop: int) -> None:
    """BKZ reduction in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
        max_loop (int): limit number of tours
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    for i in range(N):
        for j in range(N):
            pp[i][j] = ctypes.c_long(b[i, j])

    SDPB.BKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
    SDPB.BKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.BKZ(pp, block_size, d, max_loop, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def DualPotLLL(b: np.ndarray, d: float) -> None:
    """dual version of PotLLL reduction in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotLLL.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotLLL.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.DualPotLLL(pp, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def PotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """PotBKZ reduction in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    SDPB.PotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.PotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    bb = SDPB.PotBKZ(pp, block_size, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def DualPotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """DualPotBKZ reduction in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.DualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.DualPotBKZ(pp, block_size, d, n, m)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def SelfDualPotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """SelfDualPotBKZ in SelfDualPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.SelfDualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
    SDPB.SelfDualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.SelfDualPotBKZ(pp, block_size, d, n, m)

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

    print(np.linalg.norm(b[0]))
    print(b)    
    
    print("PotLLL-reduce:")
    PotLLL(c, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("BKZ-reduced:")
    BKZ(c, 40, 0.99, 10)
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
    DualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    c = b.copy()
    print("Self-Dual-PotBKZ-reduced:")
    SelfDualPotBKZ(c, 40, 0.99)
    print(np.linalg.norm(c[0]))
    print(c)

    if is_graphs_output_mode:
        bkz_potential = pd.read_csv(".data/potential_of_BKZ.csv")['Potential']
        bkz_x = np.arange(len(bkz_potential)) / N
        dual_pot_bkz_potential = pd.read_csv(".data/potential_of_DualPotBKZ.csv")['Potential']
        dual_pot_bkz_x = np.arange(len(dual_pot_bkz_potential)) / N
        pot_bkz_potential = pd.read_csv(".data/potential_of_PotBKZ.csv")['Potential']
        pot_bkz_x = np.arange(len(pot_bkz_potential)) / N
        self_dual_pot_bkz_potential = pd.read_csv(".data/potential_of_SelfDualPotBKZ.csv")['Potential']
        self_dual_pot_bkz_x = np.arange(len(self_dual_pot_bkz_potential)) / N

        fig, ax = plt.subplots()
        ax.set_xlabel("number of tours")
        ax.set_ylabel("log value of potential")
        ax.plot(bkz_x, bkz_potential, marker = "", label="Potential of BKZ")
        ax.plot(pot_bkz_x, pot_bkz_potential, marker="", label="Potential of PotBKZ")
        ax.plot(dual_pot_bkz_x, dual_pot_bkz_potential, marker="", label="Potential of dual PotBKZ")
        ax.plot(self_dual_pot_bkz_x, self_dual_pot_bkz_potential, marker = "", label="Potential of self dual PotBKZ")
        plt.tick_params()
        plt.legend()
        plt.show()
        plt.savefig(f'potential_graph/SVP_{N}_{0}.png')
    