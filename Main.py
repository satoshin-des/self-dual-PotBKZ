import ctypes, random, platform, sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

is_graphs_output_mode = True # output graphs of potential or not

N, seed, mode = list(map(int, input().split()))

if platform.system() == 'Linux':
    SDPB = ctypes.cdll.LoadLibrary("./build/libSDPotBKZ.so")
else:
    print(f"Platform {platform.system()} is not supported.")
    sys.exit(0)

def interpolation(curves: list, upper_bounds: float, interval: int) -> tuple[np.ndarray, np.ndarray]:
    x_common = np.linspace(0, upper_bounds, interval)
    ys_interp = []
    for x, y in curves:
        y_interp = np.interp(x_common, x, y)
        ys_interp.append(y_interp)

    return x_common, np.mean(ys_interp, axis=0)


def meanization(curves: list) -> tuple[np.ndarray, np.ndarray]:
    max = 0
    mean_curv = []
    for curv in curves:
        if len(curv) > max:
            max = len(curv)
    for curv in curves:
        mean_curv.append(np.concatenate([curv, np.zeros(max - len(curv))], 0))
        print(len(mean_curv[-1]))

    return np.arange(max) + 1, np.mean(mean_curv, axis=0)

def make_graph(alg: str, x_axis: str) -> tuple[np.ndarray, np.ndarray]:
    """Make data to draw graphs

    Args:
        alg (str): algorithm to draw

    Returns:
        tuple[np.ndarray, np.ndarray]: x-axis and y-axis
    """
    max = 0
    csv_file = glob.glob(f'./.data/potential_of_{alg}{N}_*.csv')
    curves_potential = []
    for s in csv_file:
        potential = np.array(pd.read_csv(s)['Potential'])
        if x_axis == "tour":
            num = np.arange(len(potential))
            curves_potential.append((num, potential))
            if num[-1] > max:
                max = num[-1]
        elif x_axis == "time":
            time = np.array(pd.read_csv(s)['Time'])
            idx = len(time)
            for i in range(len(time)):
                if time[i] > 60:
                    idx = i
                    break
            curves_potential.append((time[:idx], potential[:idx]))
            if time[idx - 1] > max:
                max = time[idx - 1]
    return interpolation(curves_potential, max, int(max) * 10)

def BKZ(b: np.ndarray, block_size: int, d: float, max_loop: int) -> None:
    """BKZ reduction in libSDPotBKZ.so

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

    SDPB.BKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int
    SDPB.BKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.BKZ(pp, block_size, d, max_loop, n, m, seed)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def PotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """PotBKZ reduction in libSDPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    SDPB.PotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.PotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
    bb = SDPB.PotBKZ(pp, block_size, d, n, m, seed)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def DualPotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """DualPotBKZ reduction in libSDPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.DualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
    SDPB.DualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.DualPotBKZ(pp, block_size, d, n, m, seed)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

def SelfDualPotBKZ(b: np.ndarray, block_size: int, d: float) -> None:
    """libSDPotBKZ in libSDPotBKZ.so

    Args:
        b (np.ndarray): lattice basis
        block_size (int): block size
        d (float): reduction parameter
    """
    n, m = b.shape

    ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
    pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

    SDPB.SelfDualPotBKZ.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int
    SDPB.SelfDualPotBKZ.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
    bb = SDPB.SelfDualPotBKZ(pp, block_size, d, n, m, seed)

    for i in range(N):
        for j in range(N):
            b[i, j] = bb[i][j]

if __name__ == '__main__':
    if mode == 0:
        b = np.eye(N, dtype=int)
    
        ptrs = [array.ctypes.data_as(ctypes.POINTER(ctypes.c_long)) for array in b]
        pp = (ctypes.POINTER(ctypes.c_long) * N)(*ptrs)

        SDPB.generator.argtypes = ctypes.POINTER(ctypes.POINTER(ctypes.c_long)), ctypes.c_int, ctypes.c_int
        SDPB.generator.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long))
        bb = SDPB.generator(pp, N, seed)

        for i in range(N):
            for j in range(N):
                b[i, j] = bb[i][j]

        print(np.linalg.norm(b[0]))
        print(b)

        c = b.copy()
        print("BKZ-reduced:")
        BKZ(c, 40, 0.99, 20)
        print(np.linalg.norm(c[0]))
        print(c)

        c = b.copy()
        print("PotBKZ-reduced:")
        PotBKZ(c, 40, 0.99)
        print(np.linalg.norm(c[0]))
        print(c)

        c = b.copy()
        print("Self-Dual-PotBKZ-reduced:")
        SelfDualPotBKZ(c, 40, 0.99)
        print(np.linalg.norm(c[0]))
        print(c)
        
    elif mode == 1:
        # Tour
        fig, ax = plt.subplots()
        ax.set_xlabel("tours")
        ax.set_ylabel("logarithm values of potential")
        
        # Potential graph of BKZ
        x_BKZ, y_BKZ = make_graph("BKZ", "tour")
        ax.plot(x_BKZ / N, y_BKZ, marker="", label="BKZ", lw=1.7)

        # Potential graph of PotBKZ
        x_PotBKZ, y_PotBKZ = make_graph("PotBKZ", "tour")
        ax.plot(x_PotBKZ / N, y_PotBKZ, marker="", label="PotBKZ", lw=1.7)

        # Potential graph of SelfDualPotBKZ
        x_SDPotBKZ, y_SDPotBKZ = make_graph("SelfDualPotBKZ", "tour")
        ax.plot(x_SDPotBKZ / N, y_SDPotBKZ, marker="", label="self-dual PotBKZ", lw=1.7)
        
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tick_params()
        plt.legend()
        fig.set_size_inches(4 * 1.7, 3 * 1.7)
        plt.savefig(f'graph/{N}_compare_potential_tour.png')
        plt.savefig(f'graph/{N}_compare_potential_tour.pdf')
        
        # Time
        fig, ax = plt.subplots()
        ax.set_xlabel("Run-times[secs]")
        ax.set_ylabel("logarithm values of potential")
        
        # Potential graph of BKZ
        x_BKZ, y_BKZ = make_graph("BKZ", "time")
        ax.plot(x_BKZ, y_BKZ, marker="", label="BKZ", lw=1.7)

        # Potential graph of PotBKZ
        x_PotBKZ, y_PotBKZ = make_graph("PotBKZ", "time")
        ax.plot(x_PotBKZ, y_PotBKZ, marker="", label="PotBKZ", lw=1.7)

        # Potential graph of SelfDualPotBKZ
        x_SDPotBKZ, y_SDPotBKZ = make_graph("SelfDualPotBKZ", "time")
        ax.plot(x_SDPotBKZ, y_SDPotBKZ, marker="", label="self-dual PotBKZ", lw=1.7)
        
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tick_params()
        plt.legend()
        fig.set_size_inches(4 * 1.7, 3 * 1.7)
        plt.savefig(f'graph/{N}_compare_potential_time.png')
        plt.savefig(f'graph/{N}_compare_potential_time.pdf')
