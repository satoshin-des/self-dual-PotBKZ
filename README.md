# self-dual-PotBKZ

## Setup
This libray need NTL library and Eigen library.
If you have not installed these libraries, please install with the below command
```shell
# install NTL library
$ sudo apt-get install -y libntl-dev

# install Eigen library
$ sudo apt install libeigen3-dev
```


## How to Install This
First, you need to clone this repository to your computer:

```bash
$ git clone https://github.com/satoshin-des/self-dual-PotBKZ.git
```

Using ``cd`` command to change the directory ``self-dual-PotBKZ``:

```shell
$ cd self-dual-PotBKZ
```

Next, running `Make.sh` to make

```bash
self-dual-PotBKZ$ ./Make.sh
```

Below is result of running above command:

```shell
-- The C compiler identification is GNU 11.4.0
-- The CXX compiler identification is GNU 11.4.0
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /usr/bin/cc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found OpenMP_C: -fopenmp (found version "4.5") 
-- Found OpenMP_CXX: -fopenmp (found version "4.5") 
-- Found OpenMP: TRUE (found version "4.5")  
-- Configuring done
-- Generating done
-- Build files have been written to: hoge/self-dual-PotBKZ/build
[  6%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/BKZ.cpp.o
[ 13%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/DualPotBKZ.cpp.o
[ 20%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/DualPotENUM.cpp.o
[ 26%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/DualPotLLL.cpp.o
[ 33%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/ENUM.cpp.o
[ 40%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/GSO.cpp.o
[ 46%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/Insert.cpp.o
[ 53%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/Lattice.cpp.o
[ 60%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/PotBKZ.cpp.o
[ 66%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/PotENUM.cpp.o
[ 73%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/PotLLL.cpp.o
[ 80%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/SDPotBKZ.cpp.o
[ 86%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/SelfDualPotBKZ.cpp.o
[ 93%] Building CXX object CMakeFiles/SDPotBKZ.dir/src/generator.cpp.o
[100%] Linking CXX shared library libSDPotBKZ.so
[100%] Built target SDPotBKZ
```

You can execute the algorithms, BKZ, PotBKZ, and self-dual PotBKZ, with running `Main.sh`:

```shell
self-dual-PotBKZ$ ./Main.sh
```

The sources in this repository are research results that was presented in SCIS2025 at Kita-kyushu city in Japan. The research results are new variants of BKZ with provably termination <b>PotBKZ</b> and its dual version and self-dual version.<br>Specifically, BKZ algorithm[SE94] outputs much better lattice basis than, for example, LLL algorithm[LLL82]. But as trade-off, BKZ algorithm is exponentially on lattice dimension and is not guaranteed terminating.<br>PotBKZ, its dual version, and its self-dual version are new variants of BKZ algorithm that terminates polynomaial tour on lattice dimension by monotonically decreasing potential of lattice basis.

- [LLL82] Arjen Klaas Lenstra and Hendrik Willem Lenstra and László Lovász, "Factoring polynomials with rational coefficients", 1982
- [SE94] Claus-Peter Schnorr and Martin Euchner, "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", 1994
