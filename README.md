# self-dual-PotBKZ

## Setup
This libray need NTL library and Eigen library.
If you have not installed these libraries, please install with the below command
```shell
# install NTL library
$ sudo apt-get install -y libntl-dev
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following additional packages will be installed:
  libgf2x3 libntl44
The following NEW packages will be installed:
  libgf2x3 libntl-dev libntl44
0 upgraded, 3 newly installed, 0 to remove and 29 not upgraded.
Need to get 2,236 kB of archives.
After this operation, 11.9 MB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu jammy/universe amd64 libgf2x3 amd64 1.3.0-2 [27.9 kB]
Get:2 http://archive.ubuntu.com/ubuntu jammy/universe amd64 libntl44 amd64 11.5.1-1 [838 kB]
Get:3 http://archive.ubuntu.com/ubuntu jammy/universe amd64 libntl-dev amd64 11.5.1-1 [1,370 kB]
Fetched 2,236 kB in 1s (2,323 kB/s)
Selecting previously unselected package libgf2x3:amd64.
(Reading database ... 126210 files and directories currently installed.)
Preparing to unpack .../libgf2x3_1.3.0-2_amd64.deb ...
Unpacking libgf2x3:amd64 (1.3.0-2) ...
Selecting previously unselected package libntl44:amd64.
Preparing to unpack .../libntl44_11.5.1-1_amd64.deb ...
Unpacking libntl44:amd64 (11.5.1-1) ...
Selecting previously unselected package libntl-dev.
Preparing to unpack .../libntl-dev_11.5.1-1_amd64.deb ...
Unpacking libntl-dev (11.5.1-1) ...
Setting up libgf2x3:amd64 (1.3.0-2) ...
Setting up libntl44:amd64 (11.5.1-1) ...
Setting up libntl-dev (11.5.1-1) ...
Processing triggers for libc-bin (2.35-0ubuntu3.8) ...


# install Eigen library
$ sudo apt install libeigen3-dev
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Suggested packages:
  libeigen3-doc libmpfrc++-dev
The following NEW packages will be installed:
  libeigen3-dev
0 upgraded, 1 newly installed, 0 to remove and 29 not upgraded.
Need to get 1,056 kB of archives.
After this operation, 9,081 kB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu jammy/universe amd64 libeigen3-dev all 3.4.0-2ubuntu2 [1,056 kB]
Fetched 1,056 kB in 0s (9,706 kB/s)
Selecting previously unselected package libeigen3-dev.
(Reading database ... 126493 files and directories currently installed.)
Preparing to unpack .../libeigen3-dev_3.4.0-2ubuntu2_all.deb ...
Unpacking libeigen3-dev (3.4.0-2ubuntu2) ...
Setting up libeigen3-dev (3.4.0-2ubuntu2) ...
```


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

If you do not have Ubuntu environment etc., you can execute this on Google Colab. The link is [here](https://colab.research.google.com/drive/14V0Hhf8cbboGDOs2phoJpCYRIRrr8Gis?usp=sharing).

The sources in this repository are research results that was presented in SCIS2025 at Kita-kyushu city in Japan. The research results are new variants of BKZ with provably termination **PotBKZ** and its dual version and self-dual version.<br>Specifically, BKZ algorithm[SE94] outputs much better lattice basis than, for example, LLL algorithm[LLL82]. But as trade-off, BKZ algorithm is exponentially on lattice dimension and is not guaranteed terminating.<br>PotBKZ, its dual version, and its self-dual version are new variants of BKZ algorithm that terminates polynomaial tour on lattice dimension by monotonically decreasing potential of lattice basis.

- [LLL82] Arjen Klaas Lenstra and Hendrik Willem Lenstra and László Lovász, "Factoring polynomials with rational coefficients", 1982
- [SE94] Claus-Peter Schnorr and Martin Euchner, "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", 1994
