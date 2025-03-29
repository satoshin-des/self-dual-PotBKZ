# self-dual-PotBKZ

## How to Install This
First, you need to clone this repository to your computer:
```shell
$ git clone https://github.com/satoshin-des/self-dual-PotBKZ.git
```

Using ``cd`` command to change the directory ``self-dual-PotBKZ``:
```shell
$ cd self-dual-PotBKZ
```

Next, use ``make`` command to compile the source codes:

```shell
$ make
g++ -shared -fPIC -O3 -fopenmp -mtune=native -march=native -mfpmath=both -o libSDPotBKZ.so src/SelfDualPotBKZ.cpp -lntl
```

You can execute the algorithms, BKZ, PotBKZ, and self-dual PotBKZ, with python( the below example is 90 dimensional svp-challenge lattice):

```shell
$ python example.py
lattice dimension = 90
5885.448411123829
[[-374  501 -527 ...  -63    0    0]
 [-176 -163  964 ...    0    0    0]
 [-606 -757 -977 ...   63    0    0]
 ...
 [ -15  524 -520 ...   41   33   -4]
 [ 177 -326 -795 ...   27  123  -93]
 [ -88  -57  582 ...  257  -78  145]]
PotLLL-reduce:
3785.9412303943654
[[-431 -258  328 ... -176 -377  247]
 [ 399  353 -136 ... -215  264 -147]
 [-151 -172 -262 ... -134 -317  240]
 ...
 [ -65  852  579 ... -622 -327   12]
 [ 163 -266  728 ...  -89 -253  -62]
 [-719 -125 -616 ... -260  489 -680]]
Dual-PotLLL-reduced:
4237.737014020573
[[   669    822    172 ...   -383    765   -112]
 [  -629    190    -99 ...   -132   -832    148]
 [   151    531    -33 ...    163    848    293]
 ...
 [-17324  74122 -46921 ...  22175  67305  95860]
 [-10527   1273   2064 ... -10667  -4711  12943]
 [ -3418  -1106   2538 ...  -7146  -7995   6704]]
PotBKZ-reduced:
3695.132879884024
[[-646   39 -268 ... -497  444   19]
 [-616 -307 -220 ... -504  482  145]
 [-349  346  417 ...   -7  444  150]
 ...
 [ 381 -491 -523 ...   86 -108 -198]
 [ -87 -124  159 ... -375  486  342]
 [ 158  388 -271 ...  -92 -107  773]]
Dual-PotBKZ-reduced:
4164.9043206297065
[[   399   -180    138 ...   -258   -231    294]
 [  -218      6    149 ...    -67   -677   -292]
 [  -523    121   -231 ...    135   -122    662]
 ...
 [-26751 -19563  -3702 ...  -2910 -13387  20780]
 [ -9134  -5972    843 ...  22412  14581 -17352]
 [  5008   3722   1887 ...   3958   4143  -8773]]
Self-Dual-PotBKZ-reduced:
3276.4476189922525
[[  505   341   448 ...  -142  -185    72]
 [  185   540   223 ...  -542   -90   471]
 [  325  -366   268 ...    39  -407  -562]
 ...
 [  161   101  -491 ...   467  -429  -259]
 [ -629 -1408   125 ...  -352  -172    71]
 [ -591  -220   899 ...   -23  -646   247]]
```

The sources in this repository are research results that was presented in SCIS2025 at Kita-kyushu city in Japan. The research results are new variants of BKZ with provably termination **PotBKZ** and its dual version and self-dual version.<br>Specifically, BKZ algorithm[SE94] outputs much better lattice basis than, for example, LLL algorithm[LLL82]. But as trade-off, BKZ algorithm is exponentially on lattice dimension and is not guaranteed terminating.<br>PotBKZ, its dual version, and its self-dual version are new variants of BKZ algorithm that terminates polynomaial tour on lattice dimension by monotonically decreasing potential of lattice basis.

- [LLL82] Arjen Klaas Lenstra and Hendrik Willem Lenstra and László Lovász, "Factoring polynomials with rational coefficients", 1982
- [SE94] Claus-Peter Schnorr and Martin Euchner, "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", 1994
