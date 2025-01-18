import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
N = 90

if 1:
    potential1 = pd.read_csv(".data/dual_potential_graph.csv")['Potential']
    X1 = np.arange(len(potential1)) / N

    potentialBKZ = pd.read_csv(".data/potential_graph_DualENUM_and_PrimalLLL.csv")['Potential']
    XBKZ = np.arange(len(potentialBKZ)) / N

    potential2 = pd.read_csv(".data/potential_graph.csv")['Potential']
    X2 = np.arange(len(potential2)) / N
    
    #potential3 = pd.read_csv("potential_graph_primalENUM_and_DualLLL.csv")['Potential']
    #X3 = np.arange(len(potential3)) / N

    #potential4 = pd.read_csv("potential_graph_DualENUM_and_PrimalLLL.csv")['Potential']
    #X4 = np.arange(len(potential4)) / N

    #T = range(int(np.ceil(np.max(list(X1) + list(X2) + list(X3) + list(X4)))) + 1)

    fig, ax = plt.subplots()
    ax.set_xlabel("ツアー数")
    ax.set_ylabel("ポテンシャル量の対数値")
    
    ax.plot(XBKZ, potentialBKZ, marker="", label="BKZのポテンシャル量", color="green", lw=1.7)
    #ax.plot(X1, potential1, marker = "", label="PotBKZのポテンシャル量", color="blue", lw=1.7)
    ax.plot(X2, potential2, marker = "", label="自己双対型PotBKZのポテンシャル量", color="red", lw=1.7)
    #ax.plot(X3, potential3, marker = "", label="Potential of the combination of primal PotENUM and dual PotLLL")
    #ax.plot(X4, potential4, marker = "", label="Potential of the combination of dual PotENUM and primal PotLLL")
    
    #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #for i in tour_data:
    #    plt.vlines(i / N, np.min(np.array(potential1)), np.max(np.array(potential1)), color='g', linestyles='dotted')
    #ax.set_xticks(T)
    plt.tick_params()
    #ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
    plt.legend()
    fig.set_size_inches(4 * 1.7, 3 * 1.7)
    plt.show()
    plt.savefig(f'../POTENTIAL_GRAPH/{N}_{S}_compare_potential.eps')